"""Extract a stable, JSON-safe metadata record from a SPHEREx IMAGE HDU.

The IMAGE header is the canonical source for observation timing, identifiers,
image units, and celestial-WCS reference metadata. WCS-WAVE is reported only
as an approximate wavelength/bandpass lookup; it is not presented as a
calibration-grade spectral solution.
"""
from __future__ import annotations

from typing import Any, Mapping

import numpy as np


def _text(header: Mapping[str, Any], key: str) -> str | None:
    value = header.get(key)
    if value is None:
        return None
    value = str(value).strip()
    return value or None


def _number(header: Mapping[str, Any], key: str) -> float | None:
    value = header.get(key)
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if np.isfinite(result) else None


def _integer(header: Mapping[str, Any], key: str) -> int | None:
    value = _number(header, key)
    if value is None or not value.is_integer():
        return None
    return int(value)


def _boolean(header: Mapping[str, Any], key: str) -> bool | None:
    value = header.get(key)
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().lower()
    if text in {"t", "true", "1", "yes"}:
        return True
    if text in {"f", "false", "0", "no"}:
        return False
    return None


def _range(values: np.ndarray) -> list[float] | None:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    if finite.size == 0:
        return None
    return [float(np.min(finite)), float(np.max(finite))]


def _approximate_spectral_metadata(wcs_wave: Any) -> dict[str, Any]:
    if wcs_wave is None:
        return {
            "available": False,
            "source_extension": "WCS-WAVE",
            "calibration_grade": False,
            "warning": "WCS-WAVE extension is unavailable.",
        }

    try:
        if len(wcs_wave) != 1:
            raise ValueError("WCS-WAVE must contain exactly one lookup-table row.")
        values = np.asarray(wcs_wave[0]["VALUES"], dtype=float)
        if values.ndim != 3 or values.shape[-1] != 2:
            raise ValueError(f"Unexpected WCS-WAVE VALUES shape: {values.shape}")
        wavelength = values[..., 0]
        bandwidth = values[..., 1]
        wave_range = _range(wavelength)
        band_range = _range(bandwidth)
        return {
            "available": wave_range is not None and band_range is not None,
            "source_extension": "WCS-WAVE",
            "calibration_grade": False,
            "wavelength_range_um_approx": wave_range,
            "bandwidth_range_um_approx": band_range,
            "lookup_shape": list(values.shape),
            "warning": (
                "Approximate per-pixel wavelength/bandpass metadata for context; "
                "use the Spectral WCS calibration product CWAVE/CBAND for science-grade values."
            ),
        }
    except (KeyError, TypeError, ValueError, IndexError) as exc:
        return {
            "available": False,
            "source_extension": "WCS-WAVE",
            "calibration_grade": False,
            "warning": f"WCS-WAVE could not be summarized: {exc}",
        }


def extract_observation_metadata(
    image_header: Mapping[str, Any],
    image_shape: tuple[int, ...] | list[int],
    *,
    collection: str | None = None,
    wcs_valid: bool | None = None,
    wcs_wave: Any = None,
) -> dict[str, Any]:
    """Return compact metadata for one SPHEREx spectral-image observation.

    Missing optional FITS keywords are returned as ``None``. Values are taken
    from this observation's own IMAGE header, never from a cross-file default.
    """
    try:
        shape = [int(value) for value in image_shape]
    except (TypeError, ValueError):
        shape = []

    # These counters are copied as header provenance only. Their scope is not
    # assumed to be the requested cutout or its pixel count.
    l2_counters: dict[str, int | float] = {}
    for key in image_header.keys():
        normalized = str(key).strip()
        if not normalized.startswith("L2 N_"):
            continue
        value = image_header.get(key)
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(number):
            l2_counters[normalized] = int(number) if number.is_integer() else number

    return {
        "schema_version": "1.0",
        "observation": {
            "observation_id": _text(image_header, "OBSID"),
            "exposure_id": _text(image_header, "EXPIDN"),
            "detector": _integer(image_header, "DETECTOR"),
            "collection": collection,
        },
        "time": {
            "time_system": _text(image_header, "TIMESYS"),
            "observation_mjd": _number(image_header, "MJD-OBS"),
            "observation_utc": _text(image_header, "DATE-OBS"),
            "start_mjd": _number(image_header, "MJD-BEG"),
            "start_utc": _text(image_header, "DATE-BEG"),
            "midpoint_mjd": _number(image_header, "MJD-AVG"),
            "midpoint_utc": _text(image_header, "DATE-AVG"),
            "end_mjd": _number(image_header, "MJD-END"),
            "end_utc": _text(image_header, "DATE-END"),
        },
        "exposure": {
            "integration_seconds": _number(image_header, "XPOSURE"),
            "elapsed_seconds": _number(image_header, "TELAPSE"),
            "sample_interval_seconds": _number(image_header, "TSAMP"),
        },
        "product": {
            "created_utc": _text(image_header, "DATE"),
            "image_unit": _text(image_header, "BUNIT"),
            "image_shape_yx": shape,
        },
        "spatial_reference": {
            "coordinate_system": _text(image_header, "RADESYS"),
            "axis_types": [
                _text(image_header, "CTYPE1"),
                _text(image_header, "CTYPE2"),
            ],
            "reference_ra_deg": _number(image_header, "CRVAL1"),
            "reference_dec_deg": _number(image_header, "CRVAL2"),
            "wcs_valid": wcs_valid,
        },
        "spectral_approximation": _approximate_spectral_metadata(wcs_wave),
        "pipeline_quality": {
            "level1_qa": _text(image_header, "L1DQAFLG"),
            "level2_qa": _text(image_header, "L2DQAFLG"),
            "astrometry_solved": _boolean(image_header, "ASTRONETSOLVED"),
            "matched_stars": _integer(image_header, "MATCH_NSTARS"),
            "psf_fwhm": _number(image_header, "PSF_FWHM"),
            "header_reported_l2_counters": l2_counters,
        },
        "provenance_note": (
            "Header-reported L2 counters are retained as provenance and are not "
            "assumed to describe only this spatial cutout."
        ),
    }
