import logging
import astropy.units as u
from astropy.coordinates import SkyCoord
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from approved_data_sources import (
    OFFICIAL_SPHEREX_COLLECTIONS,
    discover_spherex_collections,
    get_approved_sources_manifest,
)
from config import (
    ALL_OFFICIAL_SPHEREX_COLLECTIONS,
    BANDS,
    COLLECTIONS,
    CORS_ALLOWED_ORIGINS,
    DEFAULT_REFERENCE_EPOCH_JYEAR,
    DEFAULT_REFERENCE_MJD,
    TARGETS,
)
from spherex_core import (
    audit_all_benchmark_targets,
    build_cube,
    compute_spectrum,
    get_cache_status,
    propagate_coords,
    resolve_object_name,
    validate_query_parameters,
)
from serializers import package_cube_json

logger = logging.getLogger(__name__)

app = FastAPI(
    title="Project Astraea SPHEREx API",
    description=(
        "Scientific backend for NASA SPHEREx multi-epoch infrared image alignment, "
        "background normalization, coverage-aware difference imaging, and provenance validation."
    ),
    version="2.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=CORS_ALLOWED_ORIGINS,
    allow_methods=["GET", "OPTIONS"],
    allow_headers=["*"],
)


def _error_response(
    status_code: int, error_code: str, detail: str, extra: dict | None = None
) -> JSONResponse:
    content = {
        "status": "error",
        "error_code": error_code,
        "detail": detail,
        "message": detail,
    }
    if extra:
        content.update(extra)
    return JSONResponse(status_code=status_code, content=content)


@app.exception_handler(HTTPException)
async def structured_http_exception_handler(_request: Request, exc: HTTPException):
    if isinstance(exc.detail, dict):
        detail_msg = str(
            exc.detail.get("detail") or exc.detail.get("message") or exc.detail
        )
        error_code = str(exc.detail.get("error_code") or "http_error")
    else:
        detail_msg = str(exc.detail)
        error_code = (
            "object_not_found"
            if exc.status_code == 404
            else "invalid_parameters"
            if exc.status_code == 400
            else "http_error"
        )
    return _error_response(exc.status_code, error_code, detail_msg)


@app.exception_handler(Exception)
async def unhandled_exception_handler(request: Request, exc: Exception):
    logger.exception("Unhandled exception on %s: %s", request.url.path, exc)
    return _error_response(
        500,
        "internal_server_error",
        "An unexpected internal server error occurred.",
    )


def _find_nearby_benchmark_target(
    ra_deg: float,
    dec_deg: float,
    max_sep_deg: float = 0.1,
    to_mjd: float = DEFAULT_REFERENCE_MJD,
):
    query_coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    best_match = None
    best_sep = max_sep_deg
    for key, t in TARGETS.items():
        ra_now, dec_now = propagate_coords(
            t["ra"],
            t["dec"],
            t["pm_ra_arcsec_yr"],
            t["pm_dec_arcsec_yr"],
            t["epoch"],
            to_mjd=float(to_mjd),
        )
        prop_coord = SkyCoord(ra=ra_now * u.deg, dec=dec_now * u.deg, frame="icrs")
        cat_coord = SkyCoord(ra=t["ra"] * u.deg, dec=t["dec"] * u.deg, frame="icrs")
        sep_prop_deg = float(query_coord.separation(prop_coord).deg)
        sep_cat_deg = float(query_coord.separation(cat_coord).deg)
        sep_deg = min(sep_prop_deg, sep_cat_deg)
        matched_epoch = (
            "propagated_reference_epoch"
            if sep_prop_deg <= sep_cat_deg
            else "catalog_epoch"
        )
        if sep_deg <= best_sep:
            best_sep = sep_deg
            best_match = {
                "target_key": key,
                "name": t["name"],
                "category": t.get("category"),
                "test_case": t.get("test_case"),
                "data_status": t.get("data_status"),
                "separation_arcsec": round(sep_deg * 3600.0, 3),
                "matched_against": matched_epoch,
                "ra_2026": round(ra_now, 6),
                "dec_2026": round(dec_now, 6),
                "propagated_epoch_jyear": DEFAULT_REFERENCE_EPOCH_JYEAR,
                "propagated_to_mjd": float(to_mjd),
                "catalog_ra": t["ra"],
                "catalog_dec": t["dec"],
                "catalog_epoch_jyear": t["epoch"],
                "pm_convention": t.get("pm_convention", "pm_ra_cosdec"),
            }
    return best_match


@app.get("/api/targets")
def get_targets():
    out = {}
    for key, t in TARGETS.items():
        ra_now, dec_now = propagate_coords(
            t["ra"],
            t["dec"],
            t["pm_ra_arcsec_yr"],
            t["pm_dec_arcsec_yr"],
            t["epoch"],
            to_mjd=DEFAULT_REFERENCE_MJD,
        )
        out[key] = {
            **t,
            "ra_2026": round(ra_now, 6),
            "dec_2026": round(dec_now, 6),
            "propagated_epoch_jyear": DEFAULT_REFERENCE_EPOCH_JYEAR,
            "propagated_to_mjd": DEFAULT_REFERENCE_MJD,
            "catalog_epoch_jyear": t["epoch"],
            "pm_convention": t.get("pm_convention", "pm_ra_cosdec"),
        }
    return {
        "status": "ok",
        "reference_epoch": {
            "jyear": DEFAULT_REFERENCE_EPOCH_JYEAR,
            "mjd": DEFAULT_REFERENCE_MJD,
            "iso_utc": "2026-01-01T00:00:00.000",
            "pm_convention": "pm_ra_cosdec",
        },
        "official_spherex_collections": list(OFFICIAL_SPHEREX_COLLECTIONS),
        "default_queried_collections": list(COLLECTIONS),
        "targets": out,
        "bands_um": BANDS,
    }


@app.get("/api/sources")
def get_approved_sources(
    discover_live: bool = Query(
        False,
        description="Whether to probe live IRSA SIA across all 4 official SPHEREx collections",
    ),
    ra: float = Query(269.446294, description="RA (deg) for optional live collection discovery"),
    dec: float = Query(4.767881, description="Dec (deg) for optional live collection discovery"),
    fov: float = Query(0.1, description="FOV (deg) for optional live collection discovery"),
):
    """Return the approved scientific data source ecosystem, dataset separation rules,
    Gaia DR3 astrometric reference table, NEOSSat/CADC adapter specification, and
    SPHEREx collection discovery summary.
    """
    coord = SkyCoord(ra=float(ra) * u.deg, dec=float(dec) * u.deg, frame="icrs")
    discovery = discover_spherex_collections(
        coord=coord,
        fov_deg=float(fov),
        candidate_collections=OFFICIAL_SPHEREX_COLLECTIONS,
        allow_network=bool(discover_live),
    )
    return {
        "status": "ok",
        **get_approved_sources_manifest(),
        "spherex_collection_discovery": discovery,
    }


@app.get("/api/targets/audit")
def get_targets_audit(
    query_live_irsa: bool = Query(
        False,
        description="Whether to query live IRSA SIA services in addition to local cache inspection",
    ),
):
    return {
        "status": "ok",
        **audit_all_benchmark_targets(query_live_irsa=query_live_irsa),
    }


@app.get("/api/search/coordinates")
def search_coordinates(
    ra: str = Query(..., description="Right Ascension in degrees or sexagesimal"),
    dec: str = Query(..., description="Declination in degrees or sexagesimal"),
    fov: float = Query(0.1, description="Field of view in degrees (0.03 to 0.30)"),
    detectors: str = Query("5,6", description="Comma-separated detector numbers, e.g. 5,6"),
    max_frames: int = Query(30, description="Maximum frames to load (1 to 100)"),
):
    """Validate and normalize sky coordinates and check nearby NASA benchmark targets.

    Explicitly does not query IRSA for observation availability; call ``/api/frames``
    to query and retrieve SPHEREx observations.
    """
    try:
        validated = validate_query_parameters(
            ra=ra, dec=dec, fov_deg=fov, detectors=detectors, max_frames=max_frames
        )
    except ValueError as exc:
        return _error_response(400, "invalid_parameters", str(exc))

    coord = SkyCoord(ra=validated["ra"] * u.deg, dec=validated["dec"] * u.deg, frame="icrs")
    nearby = _find_nearby_benchmark_target(
        validated["ra"],
        validated["dec"],
        max_sep_deg=0.1,
        to_mjd=DEFAULT_REFERENCE_MJD,
    )
    return {
        "status": "ok",
        "operation": "coordinate_validation_and_benchmark_lookup",
        "irsa_queried": False,
        "data_availability_checked": False,
        "valid": True,
        "ra": round(validated["ra"], 6),
        "dec": round(validated["dec"], 6),
        "coordinate_frame": "ICRS",
        "hms_dms": coord.to_string("hmsdms", precision=3),
        "fov_deg": validated["fov_deg"],
        "detectors": list(validated["detectors"]),
        "max_frames": validated["max_frames"],
        "nearby_benchmark_target": nearby,
        "note": (
            "Coordinates validated and checked against benchmark catalog positions. "
            "IRSA archive availability is not queried by this endpoint; use /api/frames."
        ),
    }


@app.get("/api/search/object")
@app.get("/api/resolve")
def search_object(
    name: str = Query(..., description="Astronomical object name or benchmark alias"),
    to_mjd: float = Query(
        DEFAULT_REFERENCE_MJD,
        description="Target epoch MJD for proper-motion propagation (default 61041.0 = J2026.0)",
    ),
    allow_external: bool = Query(True, description="Allow fallback to Astropy Sesame resolver"),
):
    try:
        resolved = resolve_object_name(name=name, to_mjd=to_mjd, allow_external=allow_external)
        return {"status": "ok", **resolved}
    except ValueError as exc:
        return _error_response(400, "invalid_object_name", str(exc))
    except LookupError as exc:
        return _error_response(404, "object_not_found", str(exc))


@app.get("/api/frames")
def get_frames(
    ra: str = Query(..., description="Right Ascension in degrees or sexagesimal"),
    dec: str = Query(..., description="Declination in degrees or sexagesimal"),
    fov: float = Query(0.1, description="Field of view in degrees"),
    detectors: str = Query("5,6", description="Comma-separated detector numbers, e.g. 5,6"),
    max_frames: int = Query(30, description="Maximum frames to load"),
):
    try:
        validated = validate_query_parameters(
            ra=ra, dec=dec, fov_deg=fov, detectors=detectors, max_frames=max_frames
        )
    except ValueError as exc:
        return _error_response(400, "invalid_parameters", str(exc))

    cube = build_cube(
        validated["ra"],
        validated["dec"],
        fov_deg=validated["fov_deg"],
        detectors=validated["detectors"],
        max_frames=validated["max_frames"],
    )
    payload = package_cube_json(cube)
    if payload.get("error_code") in {
        "upstream_service_unavailable",
        "partial_upstream_failure",
    }:
        return JSONResponse(status_code=503, content=payload)
    return payload


@app.get("/api/candidates")
def get_candidates(
    ra: str = Query(..., description="Right Ascension in degrees or sexagesimal"),
    dec: str = Query(..., description="Declination in degrees or sexagesimal"),
    fov: float = Query(0.1, description="Field of view in degrees"),
    detectors: str = Query("5,6", description="Comma-separated detector numbers, e.g. 5,6"),
    max_frames: int = Query(30, description="Maximum frames to load"),
):
    try:
        validated = validate_query_parameters(
            ra=ra, dec=dec, fov_deg=fov, detectors=detectors, max_frames=max_frames
        )
    except ValueError as exc:
        return _error_response(400, "invalid_parameters", str(exc))

    cube = build_cube(
        validated["ra"],
        validated["dec"],
        fov_deg=validated["fov_deg"],
        detectors=validated["detectors"],
        max_frames=validated["max_frames"],
    )
    payload = package_cube_json(cube)
    body = {
        "status": payload.get("status", "ok"),
        "error_code": payload.get("error_code"),
        "message": payload.get("message"),
        "center": payload.get("center", {"ra": validated["ra"], "dec": validated["dec"]}),
        "candidates_by_detector": payload.get("experimental_analytics", {}).get(
            "candidates_by_detector", {}
        ),
        "sky_health": payload.get("sky_health"),
        "validation_summary": payload.get("validation_summary", {}),
        "disclaimer": (
            "Heuristic candidate detections derived from multi-epoch difference imaging; "
            "require human inspection of quality flags and coverage masks. Spatial proximity "
            "to a catalog or benchmark coordinate is never treated as independent confirmation."
        ),
    }
    if payload.get("error_code") in {
        "upstream_service_unavailable",
        "partial_upstream_failure",
    }:
        return JSONResponse(status_code=503, content=body)
    return body


@app.get("/api/cache/status")
def cache_status():
    return {"status": "ok", **get_cache_status()}


@app.get("/api/spectrum")
def get_spectrum(
    ra: str = Query(..., description="Right Ascension in degrees or sexagesimal"),
    dec: str = Query(..., description="Declination in degrees or sexagesimal"),
    fov: float = Query(0.1, description="Field of view in degrees"),
):
    try:
        validated = validate_query_parameters(
            ra=ra, dec=dec, fov_deg=fov, detectors=(1, 2, 3, 4, 5, 6), max_frames=30
        )
    except ValueError as exc:
        return _error_response(400, "invalid_parameters", str(exc))

    spec = compute_spectrum(validated["ra"], validated["dec"], fov_deg=validated["fov_deg"])
    return {
        "status": "ok" if spec else "empty",
        "center": {"ra": validated["ra"], "dec": validated["dec"]},
        "fov_deg": validated["fov_deg"],
        "spectrum": spec,
        "spectral_metadata": {
            "unit": "MJy / sr",
            "measurement_type": "aperture_surface_brightness_excess",
            "aperture_radius_px": 2.5,
            "calibration_grade": False,
            "wavelength_source": "WCS-WAVE approximate lookup table",
            "band_count_queried": 6,
            "scientific_note": (
                "SPHEREx uses Linear Variable Filters (LVFs) across 6 detectors rather than "
                "102 simultaneous static-band images per exposure. Reported wavelengths are "
                "approximate per-detector averages from WCS-WAVE, not science-grade calibrated "
                "point-source spectra."
            ),
        },
    }
