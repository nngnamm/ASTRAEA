import hashlib
import json
import logging
import os
import re
import threading
import time
import warnings
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import astropy.units as u
import numpy as np
import requests
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.stats import mad_std, sigma_clip
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.ipac.irsa import Irsa
from approved_data_sources import (
    OFFICIAL_SPHEREX_COLLECTIONS,
    discover_spherex_collections,
    get_astrometric_reference,
    infer_spherex_collection_from_url,
    propagate_astrometry_with_parallax,
    validate_dataset_separation,
)
from config import (
    ALL_OFFICIAL_SPHEREX_COLLECTIONS,
    BANDS,
    CACHE_DIR,
    COLLECTIONS,
    DEFAULT_REFERENCE_EPOCH_JYEAR,
    DEFAULT_REFERENCE_MJD,
    HTTP_TIMEOUT_SECONDS,
    MIN_COMMON_COVERAGE,
    MIN_EPOCH_COVERAGE,
    TARGETS,
)
from image_validation import build_operation_flag_mask, validate_image_payload
from observation_metadata import extract_observation_metadata
from reproject import reproject_interp

try:
    from erfa import ErfaWarning
except ImportError:
    ErfaWarning = Warning

warnings.filterwarnings("ignore", category=ErfaWarning)
logging.getLogger("astropy").setLevel(logging.ERROR)
logger = logging.getLogger(__name__)


class SiaObservationList(list):
    """List of SIA observation records carrying structured query diagnostics.

    Behaves as a standard Python ``list`` for backward compatibility with existing
    callers, while exposing ``.diagnostics`` for upstream error and completeness
    reporting.
    """

    def __init__(self, iterable=(), *, diagnostics=None):
        super().__init__(iterable)
        self.diagnostics = diagnostics or {}


def _compute_file_sha256(path: str) -> str | None:
    """Compute SHA-256 hex digest of a file on disk."""
    hasher = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            while True:
                chunk = fh.read(65536)
                if not chunk:
                    break
                hasher.update(chunk)
        return hasher.hexdigest()
    except OSError:
        return None


def propagate_coords(
    ra0,
    dec0,
    pm_ra_arcsec_yr=0.0,
    pm_dec_arcsec_yr=0.0,
    epoch0=2000.0,
    to_mjd=DEFAULT_REFERENCE_MJD,
    *,
    to_epoch_jyear=None,
):
    """Propagate ICRS sky coordinates from ``epoch0`` (Julian year) to target epoch.

    Proper motion convention:
    - ``pm_ra_arcsec_yr`` is ``pm_ra_cosdec`` (mu_alpha* = mu_alpha * cos(delta)) in arcsec/yr.
    - ``pm_dec_arcsec_yr`` is ``mu_delta`` in arcsec/yr.
    - Uses ``astropy.coordinates.SkyCoord.apply_space_motion`` for exact spherical
      space-motion propagation.
    - Defaults to ``to_mjd=DEFAULT_REFERENCE_MJD`` (61041.0 = J2026.0, 2026-01-01T00:00:00 UTC)
      unless an explicit ``to_mjd`` or ``to_epoch_jyear`` is supplied.
    """
    warnings.filterwarnings("ignore", category=ErfaWarning)
    if to_epoch_jyear is not None:
        target_time = Time(float(to_epoch_jyear), format="jyear")
    else:
        target_time = Time(float(to_mjd), format="mjd")

    if pm_ra_arcsec_yr == 0.0 and pm_dec_arcsec_yr == 0.0:
        return float(ra0), float(dec0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c0 = SkyCoord(
            ra=float(ra0) * u.deg,
            dec=float(dec0) * u.deg,
            pm_ra_cosdec=float(pm_ra_arcsec_yr) * u.arcsec / u.yr,
            pm_dec=float(pm_dec_arcsec_yr) * u.arcsec / u.yr,
            obstime=Time(float(epoch0), format="jyear"),
            frame="icrs",
        )
        c1 = c0.apply_space_motion(new_obstime=target_time)
    return float(c1.ra.deg), float(c1.dec.deg)


def propagate_coords_with_metadata(
    ra0,
    dec0,
    pm_ra_arcsec_yr=0.0,
    pm_dec_arcsec_yr=0.0,
    epoch0=2000.0,
    to_mjd=DEFAULT_REFERENCE_MJD,
    *,
    to_epoch_jyear=None,
) -> dict:
    """Propagate coordinates and return explicit epoch and proper-motion metadata."""
    if to_epoch_jyear is not None:
        target_time = Time(float(to_epoch_jyear), format="jyear")
    else:
        target_time = Time(float(to_mjd), format="mjd")

    ra1, dec1 = propagate_coords(
        ra0,
        dec0,
        pm_ra_arcsec_yr=pm_ra_arcsec_yr,
        pm_dec_arcsec_yr=pm_dec_arcsec_yr,
        epoch0=epoch0,
        to_mjd=float(target_time.mjd),
    )
    return {
        "ra": round(ra1, 6),
        "dec": round(dec1, 6),
        "catalog_ra": float(ra0),
        "catalog_dec": float(dec0),
        "catalog_epoch_jyear": float(epoch0),
        "propagated_to_mjd": round(float(target_time.mjd), 6),
        "propagated_epoch_jyear": round(float(target_time.jyear), 6),
        "propagated_iso_utc": str(target_time.utc.isot),
        "pm_ra_arcsec_yr": float(pm_ra_arcsec_yr),
        "pm_dec_arcsec_yr": float(pm_dec_arcsec_yr),
        "pm_convention": "pm_ra_cosdec",
        "proper_motion_applied": bool(
            pm_ra_arcsec_yr != 0.0 or pm_dec_arcsec_yr != 0.0
        ),
    }


def _parse_coordinate_deg(value, axis_name, default_sexagesimal_unit):
    if value is None or isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{axis_name} must be a valid coordinate; received {value!r}")

    if isinstance(value, (int, float, np.number)):
        deg = float(value)
        if not np.isfinite(deg):
            raise ValueError(f"{axis_name} must be finite; received {value!r}")
        return deg

    if isinstance(value, str):
        text = value.strip()
        if not text:
            raise ValueError(f"{axis_name} must not be empty")
        try:
            deg = float(text)
            if not np.isfinite(deg):
                raise ValueError(f"{axis_name} must be finite; received {value!r}")
            return deg
        except ValueError:
            pass
        try:
            try:
                angle = Angle(text)
            except u.UnitsError:
                angle = Angle(text, unit=default_sexagesimal_unit)
            deg = float(angle.to(u.deg).value)
            if not np.isfinite(deg):
                raise ValueError(f"{axis_name} must be finite; received {value!r}")
            return deg
        except Exception as exc:
            raise ValueError(
                f"Invalid {axis_name} coordinate {value!r}: {exc}"
            ) from exc

    raise ValueError(f"Unsupported {axis_name} coordinate type: {type(value).__name__}")


def validate_query_parameters(
    ra, dec, fov_deg=0.1, detectors=(5, 6), max_frames=30
):
    """Validate and normalize sky coordinates, FOV, detector IDs, and max_frames."""
    ra_deg = _parse_coordinate_deg(ra, "RA", u.hourangle)
    if ra_deg < 0.0 or ra_deg >= 360.0:
        raise ValueError(
            f"RA must satisfy 0.0 <= ra < 360.0 deg; received {ra_deg}"
        )

    dec_deg = _parse_coordinate_deg(dec, "Dec", u.deg)
    if dec_deg < -90.0 or dec_deg > 90.0:
        raise ValueError(
            f"Dec must satisfy -90.0 <= dec <= 90.0 deg; received {dec_deg}"
        )

    if fov_deg is None or isinstance(fov_deg, (bool, np.bool_)):
        raise ValueError(f"fov_deg must be a finite number in [0.03, 0.30]; received {fov_deg!r}")
    try:
        fov_val = float(fov_deg)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid fov_deg {fov_deg!r}") from exc
    if not np.isfinite(fov_val) or fov_val < 0.03 or fov_val > 0.30:
        raise ValueError(
            f"fov_deg must satisfy 0.03 <= fov_deg <= 0.30 deg; received {fov_val}"
        )

    if detectors is None or isinstance(detectors, (bool, np.bool_)):
        raise ValueError("detectors must be a non-empty collection of detector IDs in 1..6")
    if isinstance(detectors, str):
        tokens = [tok.strip() for tok in detectors.split(",") if tok.strip()]
        if not tokens:
            raise ValueError("detectors must not be empty")
        raw_items = tokens
    elif isinstance(detectors, (int, np.integer)):
        raw_items = [detectors]
    else:
        try:
            raw_items = list(detectors)
        except TypeError as exc:
            raise ValueError(f"Invalid detectors specification: {detectors!r}") from exc

    if not raw_items:
        raise ValueError("At least one detector in {1, 2, 3, 4, 5, 6} is required")

    parsed_detectors = set()
    for item in raw_items:
        if isinstance(item, (bool, np.bool_)):
            raise ValueError(f"Invalid detector ID: {item!r}")
        try:
            det_int = int(str(item).strip()) if isinstance(item, str) else int(item)
            if isinstance(item, (float, np.floating)) and not float(item).is_integer():
                raise ValueError(f"Detector ID must be an integer: {item!r}")
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid detector ID {item!r}; expected integer in 1..6") from exc
        if det_int not in {1, 2, 3, 4, 5, 6}:
            raise ValueError(
                f"Detector ID {det_int} is out of valid SPHEREx range {{1, 2, 3, 4, 5, 6}}"
            )
        parsed_detectors.add(det_int)

    det_tuple = tuple(sorted(parsed_detectors))

    if max_frames is None or isinstance(max_frames, (bool, np.bool_)):
        raise ValueError(f"max_frames must be an integer in [1, 100]; received {max_frames!r}")
    try:
        if isinstance(max_frames, (float, np.floating)) and not float(max_frames).is_integer():
            raise ValueError(f"max_frames must be an integer; received {max_frames!r}")
        max_frames_int = int(max_frames)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid max_frames {max_frames!r}") from exc
    if max_frames_int < 1 or max_frames_int > 100:
        raise ValueError(
            f"max_frames must satisfy 1 <= max_frames <= 100; received {max_frames_int}"
        )

    return {
        "ra": float(ra_deg),
        "dec": float(dec_deg),
        "fov_deg": float(fov_val),
        "detectors": det_tuple,
        "max_frames": int(max_frames_int),
    }


def _normalize_target_token(text: str) -> str:
    cleaned = re.sub(r"['`’]", "", text.strip().lower())
    cleaned = re.sub(r"[\s_]+", " ", cleaned)
    return cleaned


def resolve_object_name(
    name: str, to_mjd: float = DEFAULT_REFERENCE_MJD, allow_external: bool = True
) -> dict:
    """Resolve an astronomical object name against NASA benchmark targets or Sesame.

    Returns explicit epoch and proper-motion convention metadata and never applies
    proper motion twice.
    """
    if name is None or not isinstance(name, str) or not name.strip():
        raise ValueError("Object name must not be empty")

    raw_lower = name.strip().lower()
    norm_query = _normalize_target_token(name)

    for key, target in TARGETS.items():
        candidates = {
            key.lower(),
            _normalize_target_token(key),
            target["name"].lower(),
            _normalize_target_token(target["name"]),
        }
        for alias in target.get("aliases", []):
            candidates.add(alias.strip().lower())
            candidates.add(_normalize_target_token(alias))

        if raw_lower in candidates or norm_query in candidates:
            prop_meta = propagate_coords_with_metadata(
                target["ra"],
                target["dec"],
                target["pm_ra_arcsec_yr"],
                target["pm_dec_arcsec_yr"],
                target["epoch"],
                to_mjd=float(to_mjd),
            )
            return {
                "status": "ok",
                "query": name,
                "resolved_name": target["name"],
                "target_key": key,
                "source": "nasa_benchmark_catalog",
                "benchmark_target": True,
                "category": target.get("category"),
                "test_case": target.get("test_case"),
                "ra": prop_meta["ra"],
                "dec": prop_meta["dec"],
                "catalog_ra": target["ra"],
                "catalog_dec": target["dec"],
                "catalog_epoch": target["epoch"],
                "catalog_epoch_jyear": target["epoch"],
                "pm_ra_arcsec_yr": target["pm_ra_arcsec_yr"],
                "pm_dec_arcsec_yr": target["pm_dec_arcsec_yr"],
                "pm_convention": target.get("pm_convention", "pm_ra_cosdec"),
                "coordinate_source": target.get("coordinate_source"),
                "propagated_to_mjd": prop_meta["propagated_to_mjd"],
                "propagated_epoch_jyear": prop_meta["propagated_epoch_jyear"],
                "propagated_iso_utc": prop_meta["propagated_iso_utc"],
                "proper_motion_applied": True,
                "astrometric_reference": get_astrometric_reference(key),
                "parallax_mas": target.get("parallax_mas"),
                "annual_parallax_amplitude_arcsec": (
                    round(float(target["parallax_mas"]) / 1000.0, 4)
                    if target.get("parallax_mas") is not None
                    else None
                ),
                "radial_velocity_km_s": target.get("radial_velocity_km_s"),
                "gaia_dr3_source_id": target.get("gaia_dr3_source_id"),
            }

    if not allow_external:
        raise LookupError(f"Could not resolve astronomical object '{name}'")

    try:
        coord = SkyCoord.from_name(name.strip())
    except Exception as exc:
        logger.warning("Sesame name resolution failed for %r: %s", name, exc)
        raise LookupError(f"Could not resolve astronomical object '{name}'") from exc

    return {
        "status": "ok",
        "query": name,
        "resolved_name": name.strip(),
        "target_key": None,
        "source": "astropy_sesame",
        "benchmark_target": False,
        "category": None,
        "ra": round(float(coord.ra.deg), 6),
        "dec": round(float(coord.dec.deg), 6),
        "catalog_epoch_jyear": 2000.0,
        "propagated_to_mjd": None,
        "propagated_epoch_jyear": None,
        "pm_convention": None,
        "proper_motion_applied": False,
        "coordinate_epoch_note": (
            "ICRS catalog coordinates returned by CDS Sesame; space motion not applied "
            "because basic name resolution does not supply proper-motion vectors."
        ),
    }


def get_cache_status() -> dict:
    """Return a summary of cached FITS cutouts and SIA metadata queries.

    Does not leak absolute host filesystem paths in API responses.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    fits_count = 0
    sia_count = 0
    corrupt_quarantined_count = 0
    total_bytes = 0
    try:
        for entry in os.scandir(CACHE_DIR):
            if not entry.is_file():
                continue
            try:
                size = entry.stat().st_size
            except OSError:
                size = 0
            if entry.name.endswith(".fits"):
                fits_count += 1
                total_bytes += size
            elif entry.name.startswith("sia_") and entry.name.endswith(".json"):
                sia_count += 1
                total_bytes += size
            elif ".corrupt" in entry.name:
                corrupt_quarantined_count += 1
    except OSError as exc:
        logger.warning("Could not scan cache directory: %s", exc)
    safe_dir_label = os.path.basename(os.path.normpath(str(CACHE_DIR))) or "cache"
    return {
        "cache_dir": safe_dir_label,
        "fits_files_count": fits_count,
        "sia_queries_count": sia_count,
        "corrupt_quarantined_count": corrupt_quarantined_count,
        "total_size_bytes": int(total_bytes),
        "total_size_mb": round(float(total_bytes) / (1024.0 * 1024.0), 3),
    }


def find_images(
    coord,
    fov_deg,
    use_cache=True,
    return_diagnostics=False,
    collections=None,
):
    """Query IRSA SIA across configured SPHEREx collections with completeness tracking.

    Distinguishes:
    - ``complete_with_observations``: all configured collections succeeded and >=1 row returned.
    - ``complete_empty``: all configured collections succeeded and 0 rows returned.
    - ``partial_failure``: at least one collection succeeded and at least one failed.
    - ``upstream_failure``: all configured collections failed.

    Incomplete multi-collection queries (``partial_failure`` or ``upstream_failure``)
    are NEVER written to the SIA cache file as complete results.
    """
    os.makedirs(CACHE_DIR, exist_ok=True)
    collections_list = list(collections) if collections is not None else list(COLLECTIONS)
    collections_str = ",".join(collections_list)
    sia_key = hashlib.md5(
        f"sia_{collections_str}_{coord.ra.deg:.5f}_{coord.dec.deg:.5f}_{fov_deg:.4f}".encode()
    ).hexdigest()
    sia_cache_path = f"{CACHE_DIR}/sia_{sia_key}.json"

    if use_cache and os.path.exists(sia_cache_path):
        try:
            with open(sia_cache_path, "r", encoding="utf-8") as fh:
                cached = json.load(fh)
            if isinstance(cached, dict):
                is_complete = bool(cached.get("complete") is True)
                cached_cols = cached.get("collections_queried")
                observations = cached.get("observations")
                if (
                    is_complete
                    and isinstance(observations, list)
                    and (cached_cols is None or list(cached_cols) == collections_list)
                ):
                    diag = {
                        "query_status": cached.get(
                            "query_status",
                            "complete_with_observations"
                            if observations
                            else "complete_empty",
                        ),
                        "complete": True,
                        "from_cache": True,
                        "collections_queried": collections_list,
                        "collection_status": cached.get("collection_status", {}),
                        "observations": observations,
                        "total_observations": len(observations),
                    }
                    if return_diagnostics:
                        return diag
                    return SiaObservationList(observations, diagnostics=diag)
                else:
                    logger.warning(
                        "Ignoring incomplete or mismatched SIA cache file %s; re-querying.",
                        os.path.basename(sia_cache_path),
                    )
            elif isinstance(cached, list):
                diag = {
                    "query_status": "complete_with_observations"
                    if cached
                    else "complete_empty",
                    "complete": True,
                    "from_cache": True,
                    "collections_queried": collections_list,
                    "collection_status": {
                        col: {
                            "status": "ok",
                            "count": sum(
                                1 for item in cached if item.get("collection") == col
                            ),
                            "error": None,
                        }
                        for col in collections_list
                    },
                    "observations": cached,
                    "total_observations": len(cached),
                }
                if return_diagnostics:
                    return diag
                return SiaObservationList(cached, diagnostics=diag)
        except Exception as exc:
            logger.warning(
                "Corrupted SIA cache file %s (%s); evicting and re-querying.",
                os.path.basename(sia_cache_path),
                exc,
            )
            try:
                os.remove(sia_cache_path)
            except OSError:
                pass

    found = []
    collection_status = {}
    n_succeeded = 0
    n_failed = 0

    for col in collections_list:
        try:
            r = Irsa.query_sia(pos=(coord, (fov_deg / 2) * u.deg), collection=col)
            t = r.to_table() if hasattr(r, "to_table") else r
            n_succeeded += 1
        except Exception as exc:
            n_failed += 1
            err_msg = f"{type(exc).__name__}: {exc}"
            logger.warning(
                "IRSA SIA query failed for collection=%s at ra=%.5f dec=%.5f fov=%.4f: %s",
                col,
                coord.ra.deg,
                coord.dec.deg,
                fov_deg,
                err_msg,
            )
            collection_status[col] = {
                "status": "error",
                "count": 0,
                "error": err_msg,
            }
            continue

        col_count = 0
        for row in t:
            url = str(row["access_url"])
            m = re.search(r"_\d+D(\d)_spx", url)
            verified_col = infer_spherex_collection_from_url(url, fallback_collection=col)
            found.append(
                {
                    "url": url,
                    "mjd": float(row["t_min"]),
                    "detector": int(m.group(1)) if m else None,
                    "collection": verified_col,
                    "queried_collection": col,
                    "dataset_identity": "SPHEREx",
                }
            )
            col_count += 1

        collection_status[col] = {
            "status": "ok",
            "count": col_count,
            "error": None,
        }

    found.sort(key=lambda d: d["mjd"])
    complete = bool(len(collections_list) > 0 and n_failed == 0)

    if complete:
        query_status = (
            "complete_with_observations" if len(found) > 0 else "complete_empty"
        )
    elif n_succeeded > 0:
        query_status = "partial_failure"
    else:
        query_status = "upstream_failure"

    diagnostics = {
        "query_status": query_status,
        "complete": complete,
        "from_cache": False,
        "collections_queried": collections_list,
        "collection_status": collection_status,
        "observations": found,
        "total_observations": len(found),
    }

    # Only cache when ALL configured collections succeeded. Never cache partial
    # or total upstream failures as complete results.
    if complete:
        cache_payload = {
            "schema_version": "2.0",
            "complete": True,
            "query_status": query_status,
            "queried_at_utc": datetime.now(timezone.utc).isoformat(),
            "collections_queried": collections_list,
            "query_parameters": {
                "ra_deg": round(float(coord.ra.deg), 6),
                "dec_deg": round(float(coord.dec.deg), 6),
                "fov_deg": round(float(fov_deg), 4),
            },
            "collection_status": collection_status,
            "observations": found,
        }
        tmp_path = f"{sia_cache_path}.tmp.{os.getpid()}.{threading.get_ident()}"
        try:
            with open(tmp_path, "w", encoding="utf-8") as fh:
                json.dump(cache_payload, fh)
            os.replace(tmp_path, sia_cache_path)
        except Exception as exc:
            logger.warning("Failed to write SIA cache %s: %s", os.path.basename(sia_cache_path), exc)
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except OSError:
                pass

    if return_diagnostics:
        return diagnostics
    return SiaObservationList(found, diagnostics=diagnostics)


def _validate_and_parse_fits_file(
    file_path: str,
    meta: dict,
    coord: SkyCoord,
    *,
    from_cache: bool,
    cache_filename: str | None = None,
) -> tuple[dict | None, dict]:
    """Open, structurally validate, and extract a frame from a FITS file on disk.

    Used both to validate newly downloaded ``.tmp`` files BEFORE promoting them
    to ``cache/{key}.fits`` and to verify already-cached FITS files.
    """
    sha256_hex = _compute_file_sha256(file_path)
    reported_filename = cache_filename or os.path.basename(file_path)

    try:
        with fits.open(file_path, memmap=False) as hdul:
            try:
                image_hdu = hdul["IMAGE"]
            except (KeyError, IndexError):
                report = {
                    "valid": False,
                    "status": "invalid",
                    "errors": ["Required IMAGE extension is missing."],
                    "warnings": [],
                    "detector": meta.get("detector"),
                    "collection": meta.get("collection"),
                    "from_cache": bool(from_cache),
                    "sha256": sha256_hex,
                    "cache_filename": reported_filename,
                }
                return None, report

            try:
                flags_hdu = hdul["FLAGS"]
                flags_data = flags_hdu.data
                flags_header = flags_hdu.header
            except (KeyError, IndexError):
                flags_data = None
                flags_header = None

            def _extension_data(name):
                try:
                    return hdul[name].data
                except (KeyError, IndexError):
                    return None

            variance_data = _extension_data("VARIANCE")
            zodi_data = _extension_data("ZODI")
            epsf_data = _extension_data("EPSF")
            wcs_wave_data = _extension_data("WCS-WAVE")

            report = validate_image_payload(
                image=image_hdu.data,
                flags=flags_data,
                header=image_hdu.header,
                expected_detector=meta.get("detector"),
                variance=variance_data,
                zodi=zodi_data,
                epsf=epsf_data,
                wcs_wave=wcs_wave_data,
                flags_header=flags_header,
                require_full_product=True,
            )
            report["collection"] = meta.get("collection")
            report["source_mjd"] = meta.get("mjd")
            report["from_cache"] = bool(from_cache)
            report["sha256"] = sha256_hex
            report["cache_filename"] = reported_filename
            if not report["valid"]:
                return None, report

            header = image_hdu.header.copy()
            image = np.asarray(image_hdu.data, dtype=np.float32).copy()
            flags = np.asarray(flags_data).copy()
            variance = np.asarray(variance_data, dtype=np.float32).copy()
            zodi = np.asarray(zodi_data, dtype=np.float32).copy()
            detector = int(header["DETECTOR"])
            observation_mjd = float(header["MJD-AVG"])

            try:
                target_x, target_y = WCS(header, relax=True).celestial.world_to_pixel(
                    coord
                )
            except Exception:
                target_x, target_y = np.nan, np.nan

            try:
                spectral_wcs = WCS(header, fobj=hdul, key="W")
                spectral_wcs.sip = None
                wave_value = float(
                    spectral_wcs.pixel_to_world(target_x, target_y)[0].to(u.um).value
                )
            except Exception:
                wave_value = np.nan

            if not np.isfinite(wave_value) and detector in BANDS:
                wave_value = float(np.mean(BANDS[detector]))

            observation_metadata = extract_observation_metadata(
                header,
                image.shape,
                collection=meta.get("collection"),
                wcs_valid=report.get("wcs_valid"),
                wcs_wave=wcs_wave_data,
                sha256=sha256_hex,
                cache_filename=reported_filename,
                access_url=meta.get("url"),
                from_cache=bool(from_cache),
            )

            frame = {
                "time": observation_mjd,
                "detector": detector,
                "wave_um": wave_value,
                "image": image,
                "flags": flags,
                "variance": variance,
                "zodi": zodi,
                "header": header,
                "collection": meta.get("collection"),
                "from_cache": bool(from_cache),
                "sha256": sha256_hex,
                "cache_filename": reported_filename,
                "metadata": observation_metadata,
                "validation": report,
            }
            return frame, report

    except Exception as exc:
        report = {
            "valid": False,
            "status": "invalid",
            "errors": [f"FITS could not be read: {type(exc).__name__}: {exc}"],
            "warnings": [],
            "detector": meta.get("detector"),
            "collection": meta.get("collection"),
            "from_cache": bool(from_cache),
            "sha256": sha256_hex,
            "cache_filename": reported_filename,
        }
        return None, report


def _load_frame_with_validation(meta, coord, fov_deg, max_retries=4):
    """Load one FITS cutout and return (frame, validation_report).

    Cache integrity guarantees:
    1. If a cached FITS file exists at ``cache/{key}.fits``, it is validated on read.
       If the cached file is corrupted or fails structural validation, it is evicted
       from the cache directory and a fresh download is attempted.
    2. Downloads are written to a process/thread-isolated temporary file and must pass
       full structural FITS validation BEFORE atomic ``os.replace(tmp_path, path)``.
       Incomplete or invalid downloads are deleted and never promoted into the cache.
    """
    key = meta.get("cache_key") or hashlib.md5(
        f"{meta['url']}{coord.ra.deg:.5f}{coord.dec.deg:.5f}{fov_deg}".encode()
    ).hexdigest()
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = f"{CACHE_DIR}/{key}.fits"
    cache_filename = f"{key}.fits"
    corrupted_cache_evicted = False
    cached_failure_report = None

    if os.path.exists(path):
        frame, report = _validate_and_parse_fits_file(
            path,
            meta,
            coord,
            from_cache=True,
            cache_filename=cache_filename,
        )
        if frame is not None:
            return frame, report

        # Existing cached file is corrupt or invalid: evict it so it cannot poison
        # future requests, and attempt recovery via re-download.
        corrupted_cache_evicted = True
        cached_failure_report = report
        logger.warning(
            "Evicting corrupted/invalid cached FITS file %s: %s",
            cache_filename,
            "; ".join(report.get("errors", [])),
        )
        try:
            os.remove(path)
        except OSError as exc:
            logger.warning("Failed to remove corrupted cache file %s: %s", cache_filename, exc)

    url = f"{meta['url']}?center={coord.ra.deg},{coord.dec.deg}&size={fov_deg}"
    last_error = "download failed"
    last_validation_report = None
    tmp_path = f"{path}.tmp.{os.getpid()}.{threading.get_ident()}"

    for attempt in range(max(1, int(max_retries))):
        try:
            response = requests.get(url, timeout=HTTP_TIMEOUT_SECONDS)
            if response.status_code == 503:
                last_error = "HTTP 503 Service Unavailable"
                if attempt + 1 < max_retries:
                    time.sleep(min(2**attempt, 4))
                continue
            response.raise_for_status()
            content = response.content
            if not content:
                raise ValueError("Empty response body received for FITS download")

            with open(tmp_path, "wb") as output:
                output.write(content)

            # Validate temporary file BEFORE atomic promotion to final cache path
            frame, report = _validate_and_parse_fits_file(
                tmp_path,
                meta,
                coord,
                from_cache=False,
                cache_filename=cache_filename,
            )
            if frame is None:
                last_validation_report = report
                last_error = "; ".join(report.get("errors", ["Invalid FITS payload"]))
                try:
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
                except OSError:
                    pass
                # Structural rejection of a complete FITS response does not need
                # blind network retries unless the file could not be read at all.
                if not any(
                    "FITS could not be read" in err for err in report.get("errors", [])
                ):
                    break
                if attempt + 1 < max_retries:
                    time.sleep(min(2**attempt, 4))
                continue

            os.replace(tmp_path, path)
            if corrupted_cache_evicted:
                report["recovered_from_corrupted_cache"] = True
                report["corrupted_cache_evicted"] = True
            return frame, report

        except Exception as exc:
            last_error = f"{type(exc).__name__}: {exc}"
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except OSError:
                pass
            if attempt + 1 < max_retries:
                time.sleep(min(2**attempt, 4))

    if last_validation_report is not None:
        if corrupted_cache_evicted:
            last_validation_report["corrupted_cache_evicted"] = True
        return None, last_validation_report

    if cached_failure_report is not None:
        cached_failure_report["corrupted_cache_evicted"] = True
        cached_failure_report["errors"] = list(cached_failure_report.get("errors", [])) + [
            f"Re-download recovery failed: {last_error}"
        ]
        return None, cached_failure_report

    report = {
        "valid": False,
        "status": "invalid",
        "errors": [f"FITS download failed: {last_error}"],
        "warnings": [],
        "detector": meta.get("detector"),
        "collection": meta.get("collection"),
        "from_cache": False,
        "corrupted_cache_evicted": corrupted_cache_evicted,
    }
    return None, report


def load_frame(meta, coord, fov_deg):
    """Backward-compatible loader; returns a frame or None."""
    frame, report = _load_frame_with_validation(meta, coord, fov_deg)
    if frame is None:
        logger.warning(
            "Rejected FITS observation (detector=%s, collection=%s): %s",
            meta.get("detector"),
            meta.get("collection"),
            "; ".join(report.get("errors", [])),
        )
    return frame


def compute_reference_pixel_scale(frames):
    """Compute a deterministic median pixel scale (arcsec/px) across loaded frames."""
    if not frames:
        raise ValueError(
            "Cannot compute reference pixel scale from an empty frame list"
        )

    frame_scales = []
    for frame in frames:
        header = frame.get("header")
        if header is None:
            continue
        celestial_wcs = WCS(header, relax=True).celestial
        if not celestial_wcs.has_celestial:
            continue
        scales = proj_plane_pixel_scales(celestial_wcs) * 3600.0
        mean_scale = float(np.mean(scales))
        if np.isfinite(mean_scale) and mean_scale > 0:
            frame_scales.append(mean_scale)

    if not frame_scales:
        raise ValueError("No valid celestial WCS pixel scales found in loaded frames")

    return float(np.median(frame_scales))


def make_grid(coord, fov_deg, pix_arcsec):
    """Construct a North-up ICRS TAN reference WCS grid centered on ``coord``."""
    if not np.isfinite(fov_deg) or fov_deg <= 0:
        raise ValueError(
            f"fov_deg must be a positive finite number; received {fov_deg}"
        )
    if not np.isfinite(pix_arcsec) or pix_arcsec <= 0:
        raise ValueError(
            f"pix_arcsec must be a positive finite number; received {pix_arcsec}"
        )

    n = int(round(fov_deg * 3600.0 / pix_arcsec))
    if n < 2:
        raise ValueError(f"Computed grid size n={n} is too small for fov_deg={fov_deg}")

    w = WCS(naxis=2)
    w.wcs.radesys = "ICRS"
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.cunit = ["deg", "deg"]
    w.wcs.crval = [float(coord.ra.deg), float(coord.dec.deg)]
    w.wcs.crpix = [(n + 1) / 2.0, (n + 1) / 2.0]
    w.wcs.cdelt = [-pix_arcsec / 3600.0, pix_arcsec / 3600.0]
    return w, n


def align_frame_to_grid(frame, coord, grid, shape_out):
    """Reproject one validated frame onto the common target WCS grid.

    Preserves uncovered pixels as ``np.nan`` (never zero), stores the boolean
    ``valid_mask`` and ``footprint``, reprojects ``variance`` when present, and
    records native-vs-reprojected target coverage diagnostics without rejecting
    partial-footprint edge exposures.
    """
    image = np.asarray(frame["image"], dtype=float)
    ny, nx = image.shape
    celestial_wcs = WCS(frame["header"], relax=True).celestial
    if not celestial_wcs.has_celestial or celestial_wcs.pixel_n_dim != 2:
        raise ValueError("Frame header does not contain a valid 2-D celestial WCS")

    native_scales = proj_plane_pixel_scales(celestial_wcs) * 3600.0
    cx_native, cy_native = (nx - 1) / 2.0, (ny - 1) / 2.0
    center_sky = celestial_wcs.pixel_to_world(cx_native, cy_native)
    rx, ry = celestial_wcs.world_to_pixel(center_sky)
    roundtrip_error_px = float(np.hypot(rx - cx_native, ry - cy_native))
    if not np.isfinite(roundtrip_error_px) or roundtrip_error_px > 1e-2:
        raise ValueError(
            f"Celestial WCS round-trip error too large: {roundtrip_error_px} px"
        )

    tx, ty = celestial_wcs.world_to_pixel(coord)
    tx_f, ty_f = float(tx), float(ty)
    target_inside_native = bool(-0.5 <= tx_f < nx - 0.5 and -0.5 <= ty_f < ny - 0.5)

    aligned, footprint = reproject_interp(
        (image, celestial_wcs),
        grid,
        shape_out=shape_out,
    )
    aligned = np.asarray(aligned, dtype=float)
    footprint = np.asarray(footprint, dtype=float)

    if aligned.shape != tuple(shape_out) or footprint.shape != tuple(shape_out):
        raise ValueError("Reprojected image or footprint has an unexpected shape")

    valid_mask = np.isfinite(aligned) & np.isfinite(footprint) & (footprint > 0)
    if not np.any(valid_mask):
        raise ValueError("Reprojected image contains no valid covered pixels")

    aligned = aligned.copy()
    aligned[~valid_mask] = np.nan

    # Propagate VARIANCE onto the aligned grid when available
    variance = frame.get("variance")
    if variance is not None:
        var_arr = np.asarray(variance, dtype=float)
        if var_arr.shape == image.shape and np.isfinite(var_arr).any():
            try:
                aligned_var, var_fp = reproject_interp(
                    (var_arr, celestial_wcs),
                    grid,
                    shape_out=shape_out,
                )
                aligned_var = np.asarray(aligned_var, dtype=float)
                var_valid = (
                    valid_mask
                    & np.isfinite(aligned_var)
                    & (np.asarray(var_fp) > 0)
                    & (aligned_var > 0)
                )
                aligned_var_clean = np.full_like(aligned, np.nan, dtype=float)
                aligned_var_clean[var_valid] = aligned_var[var_valid]
                frame["aligned_variance"] = aligned_var_clean
            except Exception:
                pass

    finite_fraction = float(np.isfinite(aligned).mean())
    covered_fraction = float((np.isfinite(footprint) & (footprint > 0)).mean())

    gx, gy = grid.world_to_pixel(coord)
    ix, iy = int(round(float(gx))), int(round(float(gy)))
    out_ny, out_nx = shape_out
    target_covered_reprojected = bool(
        0 <= ix < out_nx and 0 <= iy < out_ny and valid_mask[iy, ix]
    )

    frame["aligned"] = aligned
    frame["footprint"] = footprint
    frame["valid_mask"] = valid_mask

    validation = frame.setdefault("validation", {})
    validation["aligned_shape"] = list(aligned.shape)
    validation["aligned_finite_pixel_fraction"] = round(finite_fraction, 6)
    validation["reprojection_covered_fraction"] = round(covered_fraction, 6)
    validation["alignment"] = {
        "status": "aligned",
        "native_shape_yx": [int(ny), int(nx)],
        "aligned_shape_yx": list(aligned.shape),
        "native_pixel_scale_arcsec": [round(float(s), 6) for s in native_scales],
        "wcs_roundtrip_error_px": round(roundtrip_error_px, 9),
        "target_native_pixel_xy": [round(tx_f, 4), round(ty_f, 4)],
        "target_inside_native_footprint": target_inside_native,
        "target_covered_after_reprojection": target_covered_reprojected,
        "aligned_finite_pixel_fraction": round(finite_fraction, 6),
        "reprojection_covered_fraction": round(covered_fraction, 6),
        "partial_coverage": covered_fraction < 0.5,
    }
    return frame


def evaluate_cross_frame_registration(
    frames: list[dict],
    *,
    pix_arcsec: float | None = None,
    snr_threshold: float = 5.0,
    window_radius_px: int = 3,
) -> dict:
    """Measure empirical cross-frame source centroid registration residuals.

    Does not rely solely on single-frame WCS self-round-trip checks. Within each
    detector and visit epoch (where proper motion over a few days is sub-pixel),
    measures the flux-weighted centroid of the brightest high-SNR compact source
    in the common grid and reports the empirical intra-epoch RMS registration
    scatter in pixels and arcseconds. If no suitable high-SNR compact reference
    source exists across >=2 frames, explicitly marks independent astrometric
    validation as ``incomplete_no_reference_sources``.
    """
    if not frames or len(frames) < 2:
        return {
            "status": "incomplete_insufficient_frames",
            "independent_astrometric_validation": False,
            "frames_evaluated": len(frames) if frames else 0,
            "note": "At least two aligned frames are required for cross-frame registration checks.",
        }

    by_det = group_epochs_by_detector(frames, gap_days=30)
    intra_epoch_offsets_px = []
    centroids_measured = 0
    epochs_validated = 0

    for _det, epochs in by_det.items():
        for epoch in epochs:
            if len(epoch) < 2:
                continue
            try:
                epoch_coadd = np.asarray(coadd(epoch), dtype=float)
            except Exception:
                continue
            if not np.isfinite(epoch_coadd).any():
                continue

            ny, nx = epoch_coadd.shape
            cy_grid, cx_grid = ny // 2, nx // 2
            r_search = max(window_radius_px + 2, min(ny, nx) // 5)
            y0 = max(window_radius_px + 2, cy_grid - r_search)
            y1 = min(ny - window_radius_px - 2, cy_grid + r_search + 1)
            x0 = max(window_radius_px + 2, cx_grid - r_search)
            x1 = min(nx - window_radius_px - 2, cx_grid + r_search + 1)
            if y1 <= y0 or x1 <= x0:
                continue

            coadd_sub = epoch_coadd[y0:y1, x0:x1]
            if not np.isfinite(coadd_sub).any():
                continue
            ref_rel_y, ref_rel_x = np.unravel_index(
                np.nanargmax(coadd_sub), coadd_sub.shape
            )
            ref_py, ref_px = y0 + int(ref_rel_y), x0 + int(ref_rel_x)

            epoch_centroids = []
            for frame in epoch:
                arr = np.asarray(
                    frame.get("normalized", frame.get("aligned")), dtype=float
                )
                valid_mask = np.isfinite(arr)
                if int(np.count_nonzero(valid_mask)) < 25:
                    continue

                finite_vals = arr[valid_mask]
                sky_med = float(np.median(finite_vals))
                sky_sigma = float(mad_std(finite_vals))
                if not np.isfinite(sky_sigma) or sky_sigma <= 0:
                    continue

                # Search within +/-2 px of the epoch reference source position
                ly0 = max(window_radius_px, ref_py - 2)
                ly1 = min(ny - window_radius_px, ref_py + 3)
                lx0 = max(window_radius_px, ref_px - 2)
                lx1 = min(nx - window_radius_px, ref_px + 3)
                local_box = arr[ly0:ly1, lx0:lx1]
                if not np.isfinite(local_box).all():
                    continue

                l_ry, l_rx = np.unravel_index(np.argmax(local_box), local_box.shape)
                py, px = ly0 + int(l_ry), lx0 + int(l_rx)
                peak_snr = (float(arr[py, px]) - sky_med) / sky_sigma
                if not np.isfinite(peak_snr) or peak_snr < snr_threshold:
                    continue

                wy0, wy1 = py - window_radius_px, py + window_radius_px + 1
                wx0, wx1 = px - window_radius_px, px + window_radius_px + 1
                if wy0 < 0 or wy1 > ny or wx0 < 0 or wx1 > nx:
                    continue

                stamp = arr[wy0:wy1, wx0:wx1] - sky_med
                if not np.isfinite(stamp).all():
                    continue
                weights = np.clip(stamp, 0.0, None)
                total_w = float(np.sum(weights))
                if total_w <= 0:
                    continue

                yy, xx = np.mgrid[wy0:wy1, wx0:wx1]
                cx_meas = float(np.sum(xx * weights) / total_w)
                cy_meas = float(np.sum(yy * weights) / total_w)
                epoch_centroids.append((cx_meas, cy_meas))
                centroids_measured += 1

            if len(epoch_centroids) >= 2:
                epochs_validated += 1
                pts = np.asarray(epoch_centroids, dtype=float)
                mean_pt = np.mean(pts, axis=0)
                offsets = np.hypot(pts[:, 0] - mean_pt[0], pts[:, 1] - mean_pt[1])
                intra_epoch_offsets_px.extend(float(v) for v in offsets)

    if not intra_epoch_offsets_px:
        return {
            "status": "incomplete_no_reference_sources",
            "independent_astrometric_validation": False,
            "frames_evaluated": len(frames),
            "centroids_measured": centroids_measured,
            "snr_threshold": float(snr_threshold),
            "note": (
                "No multi-frame visit epoch contained a compact reference source above "
                f"SNR>={snr_threshold:.1f} for independent empirical centroid registration verification."
            ),
        }

    offsets_arr = np.asarray(intra_epoch_offsets_px, dtype=float)
    rms_px = float(np.sqrt(np.mean(offsets_arr**2)))
    max_px = float(np.max(offsets_arr))
    rms_arcsec = (
        round(rms_px * float(pix_arcsec), 4)
        if pix_arcsec is not None and np.isfinite(pix_arcsec)
        else None
    )
    return {
        "status": "empirically_validated",
        "independent_astrometric_validation": True,
        "frames_evaluated": len(frames),
        "centroids_measured": centroids_measured,
        "epochs_validated": epochs_validated,
        "snr_threshold": float(snr_threshold),
        "intra_epoch_centroid_rms_px": round(rms_px, 4),
        "intra_epoch_centroid_max_offset_px": round(max_px, 4),
        "intra_epoch_centroid_rms_arcsec": rms_arcsec,
        "note": (
            "Empirical flux-weighted source centroid dispersion measured across aligned "
            "frames within each visit epoch on the common TAN grid."
        ),
    }


def normalize_frame_background(frame, grid=None, shape_out=None):
    """Scientifically normalize an aligned frame by removing ZODI and residual sky pedestal.

    Preserves physical surface-brightness units (MJy/sr) via additive background
    subtraction only (never dividing by flux dispersion, which would distort
    variable-source amplitudes). Excludes pixels flagged as sources or detector
    defects when estimating the sky pedestal (when sufficient unflagged sky pixels
    exist). Uncovered pixels remain ``np.nan``.
    """
    if "aligned" not in frame:
        raise ValueError("Frame must be aligned before background normalization")

    aligned = np.asarray(frame["aligned"], dtype=float)
    valid_mask = frame.get("valid_mask")
    if valid_mask is None:
        valid_mask = np.isfinite(aligned)
    else:
        valid_mask = np.asarray(valid_mask, dtype=bool) & np.isfinite(aligned)

    if not np.any(valid_mask):
        raise ValueError("Cannot normalize a frame with no valid finite pixels")

    raw_valid_pixels = aligned[valid_mask]
    raw_median = float(np.median(raw_valid_pixels))

    zodi_subtracted = False
    zodi_median_removed = 0.0
    working = aligned.copy()

    zodi = frame.get("zodi")
    header = frame.get("header")
    if (
        zodi is not None
        and header is not None
        and grid is not None
        and shape_out is not None
    ):
        zodi_arr = np.asarray(zodi, dtype=float)
        if (
            zodi_arr.shape == np.asarray(frame.get("image", aligned)).shape
            and np.isfinite(zodi_arr).any()
        ):
            try:
                celestial_wcs = WCS(header, relax=True).celestial
                aligned_zodi, zodi_fp = reproject_interp(
                    (zodi_arr, celestial_wcs),
                    grid,
                    shape_out=shape_out,
                )
                aligned_zodi = np.asarray(aligned_zodi, dtype=float)
                zodi_valid = (
                    valid_mask & np.isfinite(aligned_zodi) & (np.asarray(zodi_fp) > 0)
                )
                if np.any(zodi_valid):
                    zodi_median_removed = float(np.median(aligned_zodi[zodi_valid]))
                    working = np.where(zodi_valid, aligned - aligned_zodi, np.nan)
                    valid_mask = zodi_valid
                    zodi_subtracted = True
            except Exception:
                zodi_subtracted = False

    # Construct operation-specific quality-flag exclusion mask for sky background estimation
    flags = frame.get("flags")
    flag_defs = (frame.get("validation") or {}).get("flag_definitions")
    sky_mask = valid_mask.copy()
    sky_mask_policy = "all_valid_pixels_sigma_clipped"

    if flags is not None:
        flags_arr = np.asarray(flags)
        if flags_arr.ndim == 2 and np.issubdtype(flags_arr.dtype, np.integer):
            native_exclude_all = build_operation_flag_mask(
                flags_arr, flag_defs, operation="background_estimation"
            )
            native_exclude_defects = build_operation_flag_mask(
                flags_arr, flag_defs, operation="detector_defects_only"
            )

            def _map_flag_mask(native_mask):
                if native_mask is None:
                    return None
                if native_mask.shape == aligned.shape:
                    return np.asarray(native_mask, dtype=bool)
                if header is not None and grid is not None and shape_out is not None:
                    try:
                        celestial_wcs = WCS(header, relax=True).celestial
                        reproj_flag, _ = reproject_interp(
                            (native_mask.astype(float), celestial_wcs),
                            grid,
                            shape_out=shape_out,
                            order="nearest-neighbor",
                        )
                        return np.isfinite(reproj_flag) & (reproj_flag > 0.5)
                    except Exception:
                        return None
                return None

            aligned_exclude_all = _map_flag_mask(native_exclude_all)
            min_sky_px = min(25, max(4, int(0.10 * np.count_nonzero(valid_mask))))
            if (
                aligned_exclude_all is not None
                and int(np.count_nonzero(valid_mask & ~aligned_exclude_all)) >= min_sky_px
            ):
                sky_mask = valid_mask & ~aligned_exclude_all
                sky_mask_policy = "unflagged_sky_pixels"
            else:
                aligned_exclude_defects = _map_flag_mask(native_exclude_defects)
                if (
                    aligned_exclude_defects is not None
                    and int(np.count_nonzero(valid_mask & ~aligned_exclude_defects)) >= min_sky_px
                ):
                    sky_mask = valid_mask & ~aligned_exclude_defects
                    sky_mask_policy = "detector_defect_free_pixels_fallback"

    # Estimate residual sky pedestal using sigma-clipped median of unflagged sky pixels
    valid_working = working[sky_mask]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clipped_sky = sigma_clip(
            valid_working,
            sigma=3.0,
            cenfunc="median",
            stdfunc=mad_std,
            masked=True,
        )
    sky_pixels = (
        valid_working[~clipped_sky.mask] if np.any(~clipped_sky.mask) else valid_working
    )
    residual_pedestal = float(np.median(sky_pixels))
    sky_mad_sigma = float(mad_std(sky_pixels)) if sky_pixels.size > 1 else 0.0

    normalized = np.full_like(aligned, np.nan, dtype=float)
    normalized[valid_mask] = working[valid_mask] - residual_pedestal

    frame["normalized"] = normalized
    frame["valid_mask"] = valid_mask

    method = (
        "zodi_model_plus_sigma_clipped_sky_median"
        if zodi_subtracted
        else "sigma_clipped_sky_median"
    )
    valid_count = int(np.count_nonzero(valid_mask))
    sky_candidate_count = int(valid_working.size)
    validation = frame.setdefault("validation", {})
    validation["normalization"] = {
        "status": "normalized",
        "method": method,
        "sky_estimation_mask_policy": sky_mask_policy,
        "flag_aware_sky_estimation": sky_mask_policy != "all_valid_pixels_sigma_clipped",
        "valid_pixel_count": valid_count,
        "sky_estimation_pixel_count": sky_candidate_count,
        "flagged_pixels_excluded_from_sky": int(valid_count - sky_candidate_count),
        "units_preserved": frame.get("metadata", {})
        .get("product", {})
        .get("image_unit")
        or "MJy / sr",
        "zodi_model_subtracted": zodi_subtracted,
        "raw_aligned_median_mjy_sr": round(raw_median, 6),
        "zodi_median_removed_mjy_sr": round(zodi_median_removed, 6),
        "residual_pedestal_removed_mjy_sr": round(residual_pedestal, 6),
        "total_background_removed_mjy_sr": round(
            zodi_median_removed + residual_pedestal, 6
        ),
        "normalized_sky_median_mjy_sr": round(float(np.nanmedian(normalized)), 6),
        "sky_noise_mad_mjy_sr": round(sky_mad_sigma, 6),
        "limitations_note": (
            "Additive ZODI + sigma-clipped sky pedestal subtraction preserves MJy/sr "
            "point-source flux amplitudes; in crowded fields, extended nebulosity, or "
            "cutouts dominated by bright PSF wings, pedestal estimation can absorb "
            "diffuse extended emission."
        ),
    }
    return frame


def build_cube(ra, dec, fov_deg=0.1, detectors=(5, 6), max_frames=40):
    coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    raw_sia = find_images(coord, fov_deg)
    sia_diag = getattr(raw_sia, "diagnostics", None) or {
        "query_status": "complete_with_observations" if raw_sia else "complete_empty",
        "complete": True,
        "from_cache": False,
        "collections_queried": list(COLLECTIONS),
        "collection_status": {},
    }
    imgs = [m for m in raw_sia if m["detector"] in detectors]

    if not imgs:
        query_status = sia_diag.get("query_status", "complete_empty")
        if query_status == "upstream_failure":
            return {
                "status": "error",
                "error_code": "upstream_service_unavailable",
                "frames": [],
                "message": "Upstream IRSA SIA archive service failed across all queried collections.",
                "validation_summary": {
                    "query_records": 0,
                    "valid_fits": 0,
                    "aligned_frames": 0,
                    "rejected_frames": 0,
                    "rejections": [],
                    "sia_query_diagnostics": sia_diag,
                    "quality_flags_policy": "header_decoded_no_mask_applied_operation_specific_policy_required",
                },
            }
        if query_status == "partial_failure":
            return {
                "status": "error",
                "error_code": "partial_upstream_failure",
                "frames": [],
                "message": "IRSA SIA query partially failed and returned no observations for the requested detectors.",
                "validation_summary": {
                    "query_records": 0,
                    "valid_fits": 0,
                    "aligned_frames": 0,
                    "rejected_frames": 0,
                    "rejections": [],
                    "sia_query_diagnostics": sia_diag,
                    "quality_flags_policy": "header_decoded_no_mask_applied_operation_specific_policy_required",
                },
            }
        return {
            "status": "empty",
            "error_code": "no_observations_found",
            "frames": [],
            "message": "no SPHEREx images cover this position",
            "validation_summary": {
                "query_records": 0,
                "valid_fits": 0,
                "aligned_frames": 0,
                "rejected_frames": 0,
                "rejections": [],
                "sia_query_diagnostics": sia_diag,
                "quality_flags_policy": "header_decoded_no_mask_applied_operation_specific_policy_required",
            },
        }

    if len(imgs) > max_frames:
        imgs = [imgs[i] for i in np.linspace(0, len(imgs) - 1, max_frames).astype(int)]

    with ThreadPoolExecutor(max_workers=3) as executor:
        load_results = list(
            executor.map(
                lambda item: _load_frame_with_validation(item, coord, fov_deg), imgs
            )
        )

    loaded = [frame for frame, report in load_results if frame is not None]
    rejected_reports = [report for frame, report in load_results if frame is None]

    if not loaded:
        return {
            "status": "empty",
            "error_code": "validation_rejected",
            "frames": [],
            "message": "no FITS observations passed structural validation",
            "validation_summary": {
                "query_records": len(imgs),
                "valid_fits": 0,
                "aligned_frames": 0,
                "rejected_frames": len(rejected_reports),
                "rejections": [
                    {
                        "detector": report.get("detector"),
                        "collection": report.get("collection"),
                        "errors": report.get("errors", []),
                    }
                    for report in rejected_reports
                ],
                "sia_query_diagnostics": sia_diag,
                "quality_flags_policy": "header_decoded_no_mask_applied_operation_specific_policy_required",
            },
        }

    try:
        pix_arcsec = compute_reference_pixel_scale(loaded)
        grid, n = make_grid(coord, fov_deg, pix_arcsec)
    except Exception as exc:
        return {
            "status": "empty",
            "error_code": "alignment_failed",
            "frames": [],
            "message": f"Could not construct the common image grid: {exc}",
            "validation_summary": {
                "query_records": len(imgs),
                "valid_fits": len(loaded),
                "aligned_frames": 0,
                "rejected_frames": len(rejected_reports) + len(loaded),
                "rejections": [
                    {
                        "detector": report.get("detector"),
                        "collection": report.get("collection"),
                        "errors": report.get("errors", []),
                    }
                    for report in rejected_reports
                ]
                + [
                    {"detector": frame.get("detector"), "errors": [str(exc)]}
                    for frame in loaded
                ],
                "sia_query_diagnostics": sia_diag,
                "quality_flags_policy": "header_decoded_no_mask_applied_operation_specific_policy_required",
            },
        }

    aligned_frames = []
    alignment_rejections = []
    for frame in loaded:
        try:
            align_frame_to_grid(frame, coord, grid, shape_out=(n, n))
            normalize_frame_background(frame, grid=grid, shape_out=(n, n))
            aligned_frames.append(frame)
        except Exception as exc:
            alignment_rejections.append(
                {
                    "detector": frame.get("detector"),
                    "collection": frame.get("collection"),
                    "errors": [
                        f"Reprojection validation failed: {type(exc).__name__}: {exc}"
                    ],
                }
            )

    all_rejections = [
        {
            "detector": report.get("detector"),
            "collection": report.get("collection"),
            "errors": report.get("errors", []),
        }
        for report in rejected_reports
    ] + alignment_rejections

    coverages = [
        frame["validation"]["reprojection_covered_fraction"] for frame in aligned_frames
    ]
    target_in_native_count = sum(
        1
        for frame in aligned_frames
        if frame["validation"]["alignment"]["target_inside_native_footprint"]
    )
    partial_coverage_count = sum(
        1
        for frame in aligned_frames
        if frame["validation"]["alignment"]["partial_coverage"]
    )

    cross_frame_reg = evaluate_cross_frame_registration(
        aligned_frames, pix_arcsec=pix_arcsec
    )

    validation_summary = {
        "query_records": len(imgs),
        "valid_fits": len(loaded),
        "aligned_frames": len(aligned_frames),
        "rejected_frames": len(all_rejections),
        "rejections": all_rejections,
        "sia_query_diagnostics": sia_diag,
        "quality_flags_policy": "header_decoded_no_mask_applied_operation_specific_policy_required",
        "alignment_summary": {
            "reference_frame": "ICRS",
            "projection": "TAN",
            "grid_shape_yx": [n, n],
            "reference_pixel_scale_arcsec": round(pix_arcsec, 6),
            "target_inside_native_footprint_count": target_in_native_count,
            "partial_coverage_frame_count": partial_coverage_count,
            "min_covered_fraction": round(float(min(coverages)), 6)
            if coverages
            else 0.0,
            "median_covered_fraction": round(float(np.median(coverages)), 6)
            if coverages
            else 0.0,
            "max_covered_fraction": round(float(max(coverages)), 6)
            if coverages
            else 0.0,
            "cross_frame_registration": cross_frame_reg,
        },
        "normalization_summary": {
            "method": "zodi_model_plus_sigma_clipped_sky_median",
            "units_preserved": "MJy / sr",
            "zodi_subtracted_frames": sum(
                1
                for f in aligned_frames
                if f["validation"]["normalization"]["zodi_model_subtracted"]
            ),
            "flag_aware_sky_frames": sum(
                1
                for f in aligned_frames
                if f["validation"]["normalization"].get("flag_aware_sky_estimation")
            ),
        },
    }

    aligned_frames.sort(key=lambda frame: frame["time"])
    if not aligned_frames:
        return {
            "status": "empty",
            "error_code": "alignment_failed",
            "frames": [],
            "message": "no observations passed alignment validation",
            "validation_summary": validation_summary,
        }

    return {
        "status": "ok",
        "frames": aligned_frames,
        "grid": grid,
        "n": n,
        "pix_arcsec": pix_arcsec,
        "center": (ra, dec),
        "validation_summary": validation_summary,
    }


def group_epochs(frames, gap_days=30):
    """Group one detector's observations into epochs separated by large gaps.

    Frames are sorted by observation time (MJD-AVG) inside this function. Empty
    input returns an empty list instead of raising IndexError.
    """
    if not frames:
        return []

    ordered = sorted(frames, key=lambda frame: float(frame["time"]))
    epochs = []
    current = [ordered[0]]

    for previous, following in zip(ordered[:-1], ordered[1:]):
        gap = float(following["time"]) - float(previous["time"])
        if gap > gap_days:
            epochs.append(current)
            current = []
        current.append(following)

    epochs.append(current)
    return epochs


def group_epochs_by_detector(frames, gap_days=30):
    """Group observations independently for each detector.

    SPHEREx detectors represent different wavelength channels. Keeping their
    timelines separate prevents a temporal epoch or coadd from mixing bands.
    Returns ``{detector_number: [epoch, ...]}`` in detector-number order.
    """
    by_detector = {}
    detector_ids = sorted({int(frame["detector"]) for frame in frames})

    for detector in detector_ids:
        detector_frames = [
            frame for frame in frames if int(frame["detector"]) == detector
        ]
        by_detector[detector] = group_epochs(detector_frames, gap_days=gap_days)

    return by_detector


def coadd(frames, use_normalized=True):
    """Robustly combine aligned frames that share a detector and grid.

    Uses background-normalized arrays (``frame['normalized']``) when available
    so partial-footprint frames with different zodiacal backgrounds do not
    imprint step-edge seams into the coadded epoch image.
    """
    if not frames:
        raise ValueError("Cannot coadd an empty observation list")

    arrays = [
        np.asarray(
            frame["normalized"]
            if (use_normalized and "normalized" in frame)
            else frame["aligned"],
            dtype=float,
        )
        for frame in frames
    ]
    reference_shape = arrays[0].shape
    if any(array.shape != reference_shape for array in arrays):
        raise ValueError("All aligned images in a coadd must have the same shape")

    if len(arrays) == 1:
        return arrays[0].copy()

    stack = np.stack(arrays, axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clipped = sigma_clip(
            stack,
            sigma=2.5,
            cenfunc="median",
            stdfunc=mad_std,
            axis=0,
            masked=True,
        )
        with np.errstate(all="ignore"):
            return np.nanmedian(clipped.filled(np.nan), axis=0)


def _coadd_variance_map(frames, target_shape: tuple[int, ...]) -> np.ndarray | None:
    """Propagate per-frame variance maps into an epoch median-coadd variance map.

    For N independent measurements at a pixel with positive finite variances v_i,
    the asymptotic variance of the sample median is (pi / 2) * sum(v_i) / N^2
    for N > 1, and v_1 for N = 1. Returns None if no frame provides a valid
    variance array matching ``target_shape``.
    """
    var_arrays = []
    for frame in frames:
        v = frame.get("aligned_variance")
        if v is None:
            v = frame.get("variance")
        if v is None:
            continue
        v_arr = np.asarray(v, dtype=float)
        if v_arr.shape != tuple(target_shape):
            continue
        clean = np.where(np.isfinite(v_arr) & (v_arr > 0), v_arr, np.nan)
        if np.isfinite(clean).any():
            var_arrays.append(clean)

    if not var_arrays:
        return None

    stack = np.stack(var_arrays, axis=0)
    valid_counts = np.sum(np.isfinite(stack), axis=0)
    with np.errstate(all="ignore"):
        sum_var = np.nansum(stack, axis=0)
        coadd_var = np.full(target_shape, np.nan, dtype=float)
        single = valid_counts == 1
        multi = valid_counts > 1
        coadd_var[single] = sum_var[single]
        coadd_var[multi] = (np.pi / 2.0) * sum_var[multi] / (valid_counts[multi] ** 2)
    return coadd_var


def _sigma_clipped_sky_stats(values):
    """Compute sigma-clipped (sigma=3.0, median, mad_std) background pedestal and MAD."""
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        raise ValueError("Cannot compute background statistics on an empty pixel set")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clipped = sigma_clip(
            arr,
            sigma=3.0,
            cenfunc="median",
            stdfunc=mad_std,
            masked=True,
        )
    unmasked = arr[~clipped.mask] if np.any(~clipped.mask) else arr
    pedestal = float(np.median(unmasked))
    noise_mad = float(mad_std(unmasked)) if unmasked.size > 1 else 0.0
    return pedestal, noise_mad


def _extract_epoch_wavelength_stats(frames):
    """Extract [min_wave, max_wave] and mean wavelength (um) across epoch frames."""
    waves = []
    for frame in frames:
        w = frame.get("wave_um")
        if w is not None:
            try:
                wf = float(w)
                if np.isfinite(wf):
                    waves.append(wf)
            except (TypeError, ValueError):
                continue
    if not waves:
        return None, None
    return [round(float(min(waves)), 4), round(float(max(waves)), 4)], float(
        np.mean(waves)
    )


def _extract_epoch_observation_ids(frames):
    """Return ordered list of OBSID strings from an epoch's frames."""
    obs_ids = []
    for frame in frames:
        obs_id = (
            frame.get("observation_id")
            or frame.get("obsid")
            or frame.get("OBSID")
            or (frame.get("metadata") or {}).get("observation", {}).get("observation_id")
            or (frame.get("header") or {}).get("OBSID")
        )
        if obs_id is not None:
            text = str(obs_id).strip()
            if text:
                obs_ids.append(text)
    return obs_ids


def compute_difference_with_diagnostics(
    frames_a, frames_b, detector=None, min_common_coverage=MIN_COMMON_COVERAGE
):
    """Compute a coverage-aware, background-subtracted difference image and diagnostics.

    Scientific integrity guarantees:
    - Never mixes detectors (raises ValueError if ``detector=None`` and multiple
      detectors are present across ``frames_a`` and ``frames_b``).
    - Computes the joint validity mask ``common_mask = np.isfinite(img_a) & np.isfinite(img_b)``.
    - Computes sigma-clipped (sigma=3.0, median, mad_std) background pedestals for
      both epochs strictly over ``common_mask`` so non-overlapping edge regions
      cannot bias the background subtraction.
    - Propagates per-frame ``VARIANCE`` into a difference variance/uncertainty map
      when available, or explicitly documents its absence.
    - Preserves all pixels outside ``common_mask`` as ``np.nan`` (never 0.0).
    """
    if detector is not None:
        resolved_detector = int(detector)
        frames_a = [f for f in frames_a if int(f["detector"]) == resolved_detector]
        frames_b = [f for f in frames_b if int(f["detector"]) == resolved_detector]
    else:
        detectors_present = {
            int(f["detector"]) for f in list(frames_a) + list(frames_b) if "detector" in f
        }
        if len(detectors_present) > 1:
            raise ValueError(
                f"Cannot mix multiple detectors in difference image: {sorted(detectors_present)}"
            )
        resolved_detector = int(next(iter(detectors_present))) if detectors_present else None

    if not frames_a or not frames_b:
        raise ValueError("Both comparison epochs must contain observations")

    detectors_a = {int(f["detector"]) for f in frames_a if "detector" in f}
    detectors_b = {int(f["detector"]) for f in frames_b if "detector" in f}
    if len(detectors_a | detectors_b) > 1:
        raise ValueError(
            f"Cannot mix multiple detectors in difference image: {sorted(detectors_a | detectors_b)}"
        )

    img_a = np.asarray(coadd(frames_a), dtype=float)
    img_b = np.asarray(coadd(frames_b), dtype=float)
    if img_a.shape != img_b.shape:
        raise ValueError("Comparison epochs must use the same aligned image grid")

    finite_a = np.isfinite(img_a)
    finite_b = np.isfinite(img_b)
    if not finite_a.any() or not finite_b.any():
        raise ValueError("Comparison epoch contains no finite image pixels")

    common_mask = finite_a & finite_b
    total_pixels = int(img_a.size)
    common_valid_pixels = int(np.count_nonzero(common_mask))
    common_valid_pixel_fraction = (
        float(common_valid_pixels / total_pixels) if total_pixels > 0 else 0.0
    )

    if common_valid_pixels == 0 or common_valid_pixel_fraction < float(
        min_common_coverage
    ):
        raise ValueError(
            f"Insufficient common valid pixel coverage ({common_valid_pixel_fraction:.4f} < {min_common_coverage:.4f})"
        )

    pedestal_a, _ = _sigma_clipped_sky_stats(img_a[common_mask])
    pedestal_b, _ = _sigma_clipped_sky_stats(img_b[common_mask])

    img_a_sub = np.full_like(img_a, np.nan, dtype=float)
    img_b_sub = np.full_like(img_b, np.nan, dtype=float)
    diff = np.full_like(img_a, np.nan, dtype=float)

    img_a_sub[common_mask] = img_a[common_mask] - pedestal_a
    img_b_sub[common_mask] = img_b[common_mask] - pedestal_b
    diff[common_mask] = img_b_sub[common_mask] - img_a_sub[common_mask]

    if not np.isfinite(diff).any():
        raise ValueError("Difference image contains no finite pixels")

    _, diff_noise_mad = _sigma_clipped_sky_stats(diff[common_mask])

    # Propagate per-epoch variance maps when available
    var_a = _coadd_variance_map(frames_a, img_a.shape)
    var_b = _coadd_variance_map(frames_b, img_b.shape)
    variance_propagated = False
    median_prop_sigma = None
    empirical_to_prop_ratio = None
    if var_a is not None and var_b is not None:
        joint_var_mask = common_mask & np.isfinite(var_a) & np.isfinite(var_b) & ((var_a + var_b) > 0)
        if np.any(joint_var_mask):
            variance_propagated = True
            diff_sigma_map = np.sqrt(var_a[joint_var_mask] + var_b[joint_var_mask])
            median_prop_sigma = round(float(np.median(diff_sigma_map)), 6)
            if median_prop_sigma > 0:
                empirical_to_prop_ratio = round(float(diff_noise_mad / median_prop_sigma), 4)

    ny, nx = common_mask.shape
    target_covered_in_difference = bool(common_mask[ny // 2, nx // 2])

    early_wave_range, early_wave_mean = _extract_epoch_wavelength_stats(frames_a)
    late_wave_range, late_wave_mean = _extract_epoch_wavelength_stats(frames_b)
    if early_wave_mean is not None and late_wave_mean is not None:
        mean_wavelength_shift_um = round(float(abs(late_wave_mean - early_wave_mean)), 4)
    else:
        mean_wavelength_shift_um = 0.0

    diagnostics = {
        "detector": resolved_detector,
        "early_coverage_fraction": round(float(finite_a.mean()), 6),
        "late_coverage_fraction": round(float(finite_b.mean()), 6),
        "common_valid_pixel_fraction": round(common_valid_pixel_fraction, 6),
        "common_valid_pixels": common_valid_pixels,
        "target_covered_in_difference": target_covered_in_difference,
        "early_sky_median_removed_mjy_sr": round(pedestal_a, 6),
        "late_sky_median_removed_mjy_sr": round(pedestal_b, 6),
        "difference_noise_mad_mjy_sr": round(diff_noise_mad, 6),
        "variance_propagated": variance_propagated,
        "variance_propagation_method": (
            "quadrature_sum_of_epoch_coadd_median_variances"
            if variance_propagated
            else "unavailable_using_empirical_mad_noise_only"
        ),
        "median_propagated_sigma_mjy_sr": median_prop_sigma,
        "empirical_to_propagated_noise_ratio": empirical_to_prop_ratio,
        "early_wavelength_range_um": early_wave_range,
        "late_wavelength_range_um": late_wave_range,
        "mean_wavelength_shift_um": mean_wavelength_shift_um,
        "early_observation_ids": _extract_epoch_observation_ids(frames_a),
        "late_observation_ids": _extract_epoch_observation_ids(frames_b),
        "uncovered_pixel_policy": "preserved_as_nan",
    }

    return diff, img_a_sub, img_b_sub, diagnostics


def difference(frames_a, frames_b, detector=None):
    """Backward-compatible wrapper returning ``(diff, img_a_sub, img_b_sub)``."""
    diff, img_a_sub, img_b_sub, _ = compute_difference_with_diagnostics(
        frames_a,
        frames_b,
        detector=detector,
        min_common_coverage=0.0,
    )
    return diff, img_a_sub, img_b_sub


def compute_spectrum(ra, dec, fov_deg=0.1, max_per_det=8, r_ap=2.5):
    coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    all_imgs = find_images(coord, fov_deg)
    if not all_imgs:
        return []

    picked = []
    for d in range(1, 7):
        d_imgs = [m for m in all_imgs if m["detector"] == d]
        if len(d_imgs) > max_per_det:
            d_imgs = [
                d_imgs[i]
                for i in np.linspace(0, len(d_imgs) - 1, max_per_det).astype(int)
            ]
        picked.extend(d_imgs)

    with ThreadPoolExecutor(max_workers=3) as ex:
        loaded = [
            f for f in ex.map(lambda m: load_frame(m, coord, fov_deg), picked) if f
        ]
    if not loaded:
        return []

    pix_arcsec = compute_reference_pixel_scale(loaded)
    grid, n = make_grid(coord, fov_deg, pix_arcsec)
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]

    def aper_mean(img, cx, cy):
        m = np.hypot(xx - cx, yy - cy) <= r_ap
        v = img[m]
        return float(np.nanmean(v)) if np.isfinite(v).sum() > 0.5 * m.sum() else np.nan

    spectrum = []
    for d in range(1, 7):
        g = [f for f in loaded if f["detector"] == d]
        if len(g) < 2:
            continue
        for f in g:
            align_frame_to_grid(f, coord, grid, shape_out=(n, n))
            normalize_frame_background(f, grid=grid, shape_out=(n, n))
        img = coadd(g)

        vals = [
            aper_mean(img, x, y)
            for x in range(8, n - 8, 3)
            for y in range(8, n - 8, 3)
            if np.hypot(x - c, y - c) > 8
        ]
        vals = np.array([v for v in vals if np.isfinite(v)])
        center_flux = aper_mean(img, c, c)
        if len(vals) == 0 or not np.isfinite(center_flux):
            continue

        med = float(np.median(vals))
        sig = float(1.4826 * np.median(np.abs(vals - med)))
        if not np.isfinite(sig) or sig <= 0:
            continue

        excess = float(center_flux - med)
        finite_waves = [
            float(f["wave_um"])
            for f in g
            if f.get("wave_um") is not None and np.isfinite(float(f["wave_um"]))
        ]
        wave = float(np.mean(finite_waves)) if finite_waves else float(np.mean(BANDS[d]))
        sampled_range = (
            [round(float(min(finite_waves)), 4), round(float(max(finite_waves)), 4)]
            if finite_waves
            else list(BANDS[d])
        )

        spectrum.append(
            {
                "detector": d,
                "n_frames": len(g),
                "wave_um": round(wave, 3),
                "flux_mjy_sr": round(excess, 4),
                "err_mjy_sr": round(sig, 4),
                "snr": round(excess / sig, 1),
                "band_range_um": list(BANDS[d]),
                "sampled_wavelength_range_um": sampled_range,
                "calibration_grade": False,
                "wavelength_source": "WCS-WAVE_approximate",
                "aperture_radius_px": float(r_ap),
                "observation_ids": _extract_epoch_observation_ids(g),
            }
        )
    return spectrum


def audit_benchmark_target(
    target_key: str,
    *,
    fov_deg: float = 0.1,
    to_mjd: float = DEFAULT_REFERENCE_MJD,
    query_live_irsa: bool = False,
) -> dict:
    """Run a reproducible scientific provenance and readiness audit for one benchmark target.

    Never reports a failed or skipped IRSA query as zero observations. Evaluates both
    live/cached SIA query status and any locally cached FITS products matching the
    target's sky position.
    """
    if target_key not in TARGETS:
        raise KeyError(f"Unknown benchmark target key: {target_key!r}")

    target = TARGETS[target_key]
    prop_meta = propagate_coords_with_metadata(
        target["ra"],
        target["dec"],
        target["pm_ra_arcsec_yr"],
        target["pm_dec_arcsec_yr"],
        target["epoch"],
        to_mjd=float(to_mjd),
    )
    coord = SkyCoord(prop_meta["ra"] * u.deg, prop_meta["dec"] * u.deg, frame="icrs")

    sia_diag = None
    if query_live_irsa:
        sia_diag = find_images(coord, fov_deg=fov_deg, use_cache=True, return_diagnostics=True)
    else:
        # Check if a complete SIA JSON cache file already exists on disk without hitting network
        collections_str = ",".join(COLLECTIONS)
        sia_key = hashlib.md5(
            f"sia_{collections_str}_{coord.ra.deg:.5f}_{coord.dec.deg:.5f}_{fov_deg:.4f}".encode()
        ).hexdigest()
        sia_path = os.path.join(CACHE_DIR, f"sia_{sia_key}.json")
        if os.path.exists(sia_path):
            sia_diag = find_images(coord, fov_deg=fov_deg, use_cache=True, return_diagnostics=True)
        else:
            sia_diag = {
                "query_status": "not_queried_offline",
                "complete": False,
                "from_cache": False,
                "collections_queried": list(COLLECTIONS),
                "collection_status": {},
                "observations": None,
                "total_observations": None,
            }

    # Scan CACHE_DIR for locally cached FITS cutouts covering or near this target
    local_fits_frames = []
    local_fits_rejections = []
    if os.path.isdir(CACHE_DIR):
        for entry in sorted(os.scandir(CACHE_DIR), key=lambda e: e.name):
            if not entry.is_file() or not entry.name.endswith(".fits"):
                continue
            try:
                with fits.open(entry.path, memmap=False) as hdul:
                    if "IMAGE" not in hdul:
                        continue
                    hdr = hdul["IMAGE"].header
                    ny, nx = hdul["IMAGE"].data.shape
                    w = WCS(hdr, relax=True).celestial
                    center_sky = w.pixel_to_world((nx - 1) / 2.0, (ny - 1) / 2.0)
                    if float(center_sky.separation(coord).deg) > max(fov_deg * 2.0, 0.25):
                        continue
                    obsid = str(hdr.get("OBSID", ""))
                    # Infer collection from campaign prefix (2025W.. / 2026W01..19 -> dr1; 2026W33.. -> qr3)
                    col = "spherex_qr3" if obsid.startswith("2026W3") else "spherex_dr1"
                    meta = {
                        "url": f"local_cache://{entry.name}",
                        "mjd": float(hdr.get("MJD-BEG", hdr.get("MJD-AVG", 0.0))),
                        "detector": int(hdr["DETECTOR"]) if "DETECTOR" in hdr else None,
                        "collection": col,
                    }
            except Exception:
                continue

            frame, report = _validate_and_parse_fits_file(
                entry.path,
                meta,
                coord,
                from_cache=True,
                cache_filename=entry.name,
            )
            if frame is not None:
                local_fits_frames.append(frame)
            else:
                local_fits_rejections.append(report)

    local_det_counts = dict(
        sorted(Counter(int(f["detector"]) for f in local_fits_frames).items())
    )
    local_col_counts = dict(
        sorted(Counter(str(f.get("collection")) for f in local_fits_frames).items())
    )

    alignment_and_comparison = {
        "alignable": False,
        "aligned_frame_count": 0,
        "by_detector": {},
    }
    if local_fits_frames:
        try:
            pix_arcsec = compute_reference_pixel_scale(local_fits_frames)
            grid, n = make_grid(coord, fov_deg, pix_arcsec)
            aligned_ok = []
            for f in local_fits_frames:
                try:
                    align_frame_to_grid(f, coord, grid, shape_out=(n, n))
                    normalize_frame_background(f, grid=grid, shape_out=(n, n))
                    aligned_ok.append(f)
                except Exception:
                    pass
            grouped = group_epochs_by_detector(aligned_ok, gap_days=30)
            det_summaries = {}
            for det, epochs in grouped.items():
                covs = [
                    round(float(np.isfinite(coadd(ep)).mean()), 4) for ep in epochs
                ]
                eligible_epochs = [c for c in covs if c >= MIN_EPOCH_COVERAGE]
                det_summaries[str(det)] = {
                    "n_frames": sum(len(ep) for ep in epochs),
                    "n_epochs_total": len(epochs),
                    "n_epochs_meeting_coverage_threshold": len(eligible_epochs),
                    "epoch_coverages": covs,
                    "comparison_eligible": len(eligible_epochs) >= 2,
                }
            alignment_and_comparison = {
                "alignable": len(aligned_ok) > 0,
                "aligned_frame_count": len(aligned_ok),
                "grid_shape_yx": [n, n],
                "reference_pixel_scale_arcsec": round(pix_arcsec, 4),
                "by_detector": det_summaries,
            }
        except Exception as exc:
            alignment_and_comparison["error"] = str(exc)

    mjds = sorted(float(f["time"]) for f in local_fits_frames)
    any_comparison_eligible = any(
        d.get("comparison_eligible", False)
        for d in alignment_and_comparison.get("by_detector", {}).values()
    )

    if sia_diag["query_status"] == "complete_with_observations" and local_fits_frames:
        verification_status = "verified_live_and_cached"
    elif local_fits_frames:
        verification_status = "verified_local_cache_subset"
    elif sia_diag["query_status"] == "complete_with_observations":
        verification_status = "verified_sia_metadata_only"
    elif sia_diag["query_status"] == "complete_empty":
        verification_status = "verified_empty"
    else:
        verification_status = "unverified_offline"

    return {
        "target_key": target_key,
        "name": target["name"],
        "category": target.get("category"),
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "verification_status": verification_status,
        "catalog_coordinates": {
            "ra_deg": target["ra"],
            "dec_deg": target["dec"],
            "epoch_jyear": target["epoch"],
            "pm_ra_arcsec_yr": target["pm_ra_arcsec_yr"],
            "pm_dec_arcsec_yr": target["pm_dec_arcsec_yr"],
            "pm_convention": target.get("pm_convention", "pm_ra_cosdec"),
            "coordinate_source": target.get("coordinate_source"),
            "proper_motion_relevant": bool(
                abs(target["pm_ra_arcsec_yr"]) > 0.05
                or abs(target["pm_dec_arcsec_yr"]) > 0.05
            ),
        },
        "propagated_coordinates": prop_meta,
        "query_parameters": {
            "ra_deg": prop_meta["ra"],
            "dec_deg": prop_meta["dec"],
            "fov_deg": float(fov_deg),
            "collections": list(COLLECTIONS),
        },
        "sia_query": {
            "query_status": sia_diag["query_status"],
            "complete": sia_diag["complete"],
            "from_cache": sia_diag["from_cache"],
            "collection_status": sia_diag["collection_status"],
            "total_observations": sia_diag["total_observations"],
        },
        "local_cache": {
            "valid_fits_count": len(local_fits_frames),
            "rejected_fits_count": len(local_fits_rejections),
            "by_detector": local_det_counts,
            "by_collection": local_col_counts,
            "mjd_min": round(mjds[0], 6) if mjds else None,
            "mjd_max": round(mjds[-1], 6) if mjds else None,
            "baseline_days": round(mjds[-1] - mjds[0], 2) if len(mjds) > 1 else 0.0,
        },
        "alignment_and_comparison": alignment_and_comparison,
        "comparison_eligible": any_comparison_eligible,
    }


def audit_all_benchmark_targets(
    *,
    fov_deg: float = 0.1,
    to_mjd: float = DEFAULT_REFERENCE_MJD,
    query_live_irsa: bool = False,
) -> dict:
    """Audit all five NASA benchmark targets independently."""
    results = {}
    for key in TARGETS:
        results[key] = audit_benchmark_target(
            key,
            fov_deg=fov_deg,
            to_mjd=to_mjd,
            query_live_irsa=query_live_irsa,
        )
    return {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "query_live_irsa": bool(query_live_irsa),
        "reference_mjd": float(to_mjd),
        "targets": results,
    }
