import hashlib
import json
import logging
import os
import re
import time
import warnings
from concurrent.futures import ThreadPoolExecutor

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
from config import BANDS, CACHE_DIR, COLLECTIONS, TARGETS
from image_validation import validate_image_payload
from observation_metadata import extract_observation_metadata
from reproject import reproject_interp

try:
    from erfa import ErfaWarning
except ImportError:
    ErfaWarning = Warning

warnings.filterwarnings("ignore", category=ErfaWarning)
logging.getLogger("astropy").setLevel(logging.ERROR)


def propagate_coords(
    ra0, dec0, pm_ra_arcsec_yr=0.0, pm_dec_arcsec_yr=0.0, epoch0=2000.0, to_mjd=61000.0
):
    warnings.filterwarnings("ignore", category=ErfaWarning)
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
        c1 = c0.apply_space_motion(new_obstime=Time(float(to_mjd), format="mjd"))
    return float(c1.ra.deg), float(c1.dec.deg)


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
    name: str, to_mjd: float = 61000.0, allow_external: bool = True
) -> dict:
    """Resolve an astronomical object name against NASA benchmark targets or Sesame."""
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
            ra_now, dec_now = propagate_coords(
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
                "ra": round(ra_now, 6),
                "dec": round(dec_now, 6),
                "catalog_ra": target["ra"],
                "catalog_dec": target["dec"],
                "catalog_epoch": target["epoch"],
                "pm_ra_arcsec_yr": target["pm_ra_arcsec_yr"],
                "pm_dec_arcsec_yr": target["pm_dec_arcsec_yr"],
                "propagated_to_mjd": float(to_mjd),
                "proper_motion_applied": True,
            }

    if not allow_external:
        raise LookupError(f"Could not resolve astronomical object '{name}'")

    try:
        coord = SkyCoord.from_name(name.strip())
    except Exception as exc:
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
        "proper_motion_applied": False,
    }


def get_cache_status() -> dict:
    """Return a summary of cached FITS cutouts and SIA metadata queries."""
    os.makedirs(CACHE_DIR, exist_ok=True)
    fits_count = 0
    sia_count = 0
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
    except OSError:
        pass
    return {
        "cache_dir": CACHE_DIR,
        "fits_files_count": fits_count,
        "sia_queries_count": sia_count,
        "total_size_bytes": int(total_bytes),
        "total_size_mb": round(float(total_bytes) / (1024.0 * 1024.0), 3),
    }


def find_images(coord, fov_deg, use_cache=True):
    os.makedirs(CACHE_DIR, exist_ok=True)
    collections_str = ",".join(COLLECTIONS)
    sia_key = hashlib.md5(
        f"sia_{collections_str}_{coord.ra.deg:.5f}_{coord.dec.deg:.5f}_{fov_deg:.4f}".encode()
    ).hexdigest()
    sia_cache_path = f"{CACHE_DIR}/sia_{sia_key}.json"

    if use_cache and os.path.exists(sia_cache_path):
        try:
            with open(sia_cache_path, "r", encoding="utf-8") as fh:
                cached = json.load(fh)
            if isinstance(cached, list):
                return cached
        except Exception:
            pass

    found = []
    any_query_succeeded = False
    for col in COLLECTIONS:
        try:
            r = Irsa.query_sia(pos=(coord, (fov_deg / 2) * u.deg), collection=col)
            t = r.to_table() if hasattr(r, "to_table") else r
            any_query_succeeded = True
        except Exception as e:
            print(col, "query failed:", str(e)[:60])
            continue
        for row in t:
            url = str(row["access_url"])
            m = re.search(r"_\d+D(\d)_spx", url)
            found.append(
                {
                    "url": url,
                    "mjd": float(row["t_min"]),
                    "detector": int(m.group(1)) if m else None,
                    "collection": col,
                }
            )
    found.sort(key=lambda d: d["mjd"])

    if any_query_succeeded or found:
        tmp_path = f"{sia_cache_path}.tmp.{os.getpid()}"
        try:
            with open(tmp_path, "w", encoding="utf-8") as fh:
                json.dump(found, fh)
            os.replace(tmp_path, sia_cache_path)
        except Exception:
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except OSError:
                pass

    return found


def _load_frame_with_validation(meta, coord, fov_deg):
    """Load one FITS cutout and return (frame, validation_report).

    Validation is structural: image geometry, numerical content, flags-array
    shape, essential metadata, and celestial WCS are checked. The numerical
    meaning of individual FLAGS bits is deliberately not interpreted here.
    """
    key = hashlib.md5(
        f"{meta['url']}{coord.ra.deg:.5f}{coord.dec.deg:.5f}{fov_deg}".encode()
    ).hexdigest()
    os.makedirs(CACHE_DIR, exist_ok=True)
    path = f"{CACHE_DIR}/{key}.fits"
    existed_before_download = os.path.exists(path)

    if not existed_before_download:
        url = f"{meta['url']}?center={coord.ra.deg},{coord.dec.deg}&size={fov_deg}"
        downloaded = False
        last_error = "download failed"
        tmp_path = f"{path}.tmp.{os.getpid()}"
        for attempt in range(4):
            try:
                response = requests.get(url, timeout=120)
                if response.status_code == 503:
                    time.sleep(2**attempt)
                    continue
                response.raise_for_status()
                with open(tmp_path, "wb") as output:
                    output.write(response.content)
                os.replace(tmp_path, path)
                downloaded = True
                break
            except Exception as exc:
                last_error = str(exc)
                try:
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
                except OSError:
                    pass
                time.sleep(2**attempt)
        if not downloaded:
            report = {
                "valid": False,
                "status": "invalid",
                "errors": [f"FITS download failed: {last_error}"],
                "warnings": [],
                "detector": meta.get("detector"),
                "collection": meta.get("collection"),
                "from_cache": False,
            }
            return None, report

    try:
        with fits.open(path, memmap=False) as hdul:
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
                    "from_cache": bool(existed_before_download),
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
            report["from_cache"] = bool(existed_before_download)
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
                "from_cache": bool(existed_before_download),
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
            "from_cache": bool(existed_before_download),
        }
        return None, report


def load_frame(meta, coord, fov_deg):
    """Backward-compatible loader; returns a frame or None."""
    frame, report = _load_frame_with_validation(meta, coord, fov_deg)
    if frame is None:
        logging.warning(
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
    ``valid_mask`` and ``footprint``, and records native-vs-reprojected target
    coverage diagnostics without rejecting partial-footprint edge exposures.
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


def normalize_frame_background(frame, grid=None, shape_out=None):
    """Scientifically normalize an aligned frame by removing ZODI and residual sky pedestal.

    Preserves physical surface-brightness units (MJy/sr) via additive background
    subtraction only (never dividing by flux dispersion, which would distort
    variable-source amplitudes). Uncovered pixels remain ``np.nan``.
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

    # Estimate residual sky pedestal using sigma-clipped median of valid pixels
    valid_working = working[valid_mask]
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
    validation = frame.setdefault("validation", {})
    validation["normalization"] = {
        "status": "normalized",
        "method": method,
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
    }
    return frame


def build_cube(ra, dec, fov_deg=0.1, detectors=(5, 6), max_frames=40):
    coord = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    imgs = [m for m in find_images(coord, fov_deg) if m["detector"] in detectors]
    if not imgs:
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

    validation_summary = {
        "query_records": len(imgs),
        "valid_fits": len(loaded),
        "aligned_frames": len(aligned_frames),
        "rejected_frames": len(all_rejections),
        "rejections": all_rejections,
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
        },
        "normalization_summary": {
            "method": "zodi_model_plus_sigma_clipped_sky_median",
            "units_preserved": "MJy / sr",
            "zodi_subtracted_frames": sum(
                1
                for f in aligned_frames
                if f["validation"]["normalization"]["zodi_model_subtracted"]
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

    Frames are sorted by observation time inside this function. Empty input
    returns an empty list instead of raising IndexError.
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
    frames_a, frames_b, detector=None, min_common_coverage=0.10
):
    """Compute a coverage-aware, background-subtracted difference image and diagnostics.

    Scientific integrity guarantees:
    - Never mixes detectors (raises ValueError if ``detector=None`` and multiple
      detectors are present across ``frames_a`` and ``frames_b``).
    - Computes the joint validity mask ``common_mask = np.isfinite(img_a) & np.isfinite(img_b)``.
    - Computes sigma-clipped (sigma=3.0, median, mad_std) background pedestals for
      both epochs strictly over ``common_mask`` so non-overlapping edge regions
      cannot bias the background subtraction.
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
