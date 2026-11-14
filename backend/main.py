import astropy.units as u
from astropy.coordinates import SkyCoord
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from config import TARGETS, BANDS
from spherex_core import (
    build_cube,
    compute_spectrum,
    get_cache_status,
    propagate_coords,
    resolve_object_name,
    validate_query_parameters,
)
from serializers import package_cube_json

app = FastAPI(title="Project Astraea SPHEREx API")

# Enable CORS so React+Vite (localhost:5173) can call this backend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


def _error_response(status_code: int, error_code: str, detail: str) -> JSONResponse:
    return JSONResponse(
        status_code=status_code,
        content={
            "status": "error",
            "error_code": error_code,
            "detail": detail,
            "message": detail,
        },
    )


@app.exception_handler(HTTPException)
async def structured_http_exception_handler(_request: Request, exc: HTTPException):
    if isinstance(exc.detail, dict):
        detail_msg = str(exc.detail.get("detail") or exc.detail.get("message") or exc.detail)
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


def _find_nearby_benchmark_target(ra_deg: float, dec_deg: float, max_sep_deg: float = 0.1, to_mjd: float = 61000.0):
    query_coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    best_match = None
    best_sep = max_sep_deg
    for key, t in TARGETS.items():
        ra_now, dec_now = propagate_coords(
            t["ra"], t["dec"], t["pm_ra_arcsec_yr"], t["pm_dec_arcsec_yr"], t["epoch"], to_mjd=to_mjd
        )
        prop_coord = SkyCoord(ra=ra_now * u.deg, dec=dec_now * u.deg, frame="icrs")
        cat_coord = SkyCoord(ra=t["ra"] * u.deg, dec=t["dec"] * u.deg, frame="icrs")
        sep_deg = min(
            float(query_coord.separation(prop_coord).deg),
            float(query_coord.separation(cat_coord).deg),
        )
        if sep_deg <= best_sep:
            best_sep = sep_deg
            best_match = {
                "target_key": key,
                "name": t["name"],
                "category": t.get("category"),
                "test_case": t.get("test_case"),
                "separation_arcsec": round(sep_deg * 3600.0, 3),
                "ra_2026": round(ra_now, 6),
                "dec_2026": round(dec_now, 6),
            }
    return best_match


@app.get("/api/targets")
def get_targets():
    out = {}
    for key, t in TARGETS.items():
        ra_now, dec_now = propagate_coords(
            t["ra"], t["dec"], t["pm_ra_arcsec_yr"], t["pm_dec_arcsec_yr"], t["epoch"]
        )
        out[key] = {**t, "ra_2026": round(ra_now, 6), "dec_2026": round(dec_now, 6)}
    return {"status": "ok", "targets": out, "bands_um": BANDS}


@app.get("/api/search/coordinates")
def search_coordinates(
    ra: str = Query(..., description="Right Ascension in degrees or sexagesimal"),
    dec: str = Query(..., description="Declination in degrees or sexagesimal"),
    fov: float = Query(0.1, description="Field of view in degrees (0.03 to 0.30)"),
    detectors: str = Query("5,6", description="Comma-separated detector numbers, e.g. 5,6"),
    max_frames: int = Query(30, description="Maximum frames to load (1 to 100)"),
):
    try:
        validated = validate_query_parameters(
            ra=ra, dec=dec, fov_deg=fov, detectors=detectors, max_frames=max_frames
        )
    except ValueError as exc:
        return _error_response(400, "invalid_parameters", str(exc))

    coord = SkyCoord(ra=validated["ra"] * u.deg, dec=validated["dec"] * u.deg, frame="icrs")
    nearby = _find_nearby_benchmark_target(validated["ra"], validated["dec"], max_sep_deg=0.1)
    return {
        "status": "ok",
        "valid": True,
        "ra": round(validated["ra"], 6),
        "dec": round(validated["dec"], 6),
        "hms_dms": coord.to_string("hmsdms", precision=3),
        "fov_deg": validated["fov_deg"],
        "detectors": list(validated["detectors"]),
        "max_frames": validated["max_frames"],
        "nearby_benchmark_target": nearby,
    }


@app.get("/api/search/object")
@app.get("/api/resolve")
def search_object(
    name: str = Query(..., description="Astronomical object name or benchmark alias"),
    to_mjd: float = Query(61000.0, description="Target epoch MJD for proper-motion propagation"),
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
    return package_cube_json(cube)


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
    return {
        "status": payload.get("status", "ok"),
        "center": payload.get("center", {"ra": validated["ra"], "dec": validated["dec"]}),
        "candidates_by_detector": payload.get("experimental_analytics", {}).get(
            "candidates_by_detector", {}
        ),
        "sky_health": payload.get("sky_health"),
        "disclaimer": (
            "Heuristic candidate detections derived from multi-epoch difference imaging; "
            "require human inspection of quality flags and coverage masks."
        ),
    }


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
