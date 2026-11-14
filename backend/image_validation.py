"""Structural and metadata validation for SPHEREx spectral-image FITS products.

The validator separates four distinct readiness dimensions:
1. Structural validity (required HDUs, dimensions, dtypes, header keywords).
2. Numerical validity (finite pixel fractions, variance positivity, mask consistency).
3. Alignment validity (celestial WCS invertibility and reprojection readiness).
4. Scientific-operation suitability (operation-specific readiness for visual display,
   background estimation, epoch coaddition, difference imaging, and photometry).

It reads flag definitions from each FITS FLAGS header and reports both per-bit counts
and exact bitwise-union counts (preventing double-counting when a pixel carries
multiple flags). It also provides operation-specific flag masks on demand without
applying a destructive universal mask at load time.
"""
from __future__ import annotations

import re
from typing import Any

import numpy as np
from astropy.wcs import WCS


_REQUIRED_EPSF_COLUMNS = {
    "BINX", "BINY", "XCENTER", "YCENTER", "XWIDTH", "YWIDTH",
    "NSTAR", "NEFF_MEAN", "CWAVE", "CBAND", "EPSF",
}
_REQUIRED_WAVE_COLUMNS = {"X", "Y", "VALUES"}

# SPHEREx Level-2 astrophysical / source-related flag names that do not indicate
# detector hardware defects.
ASTROPHYSICAL_FLAG_NAMES = {"SOURCE", "TRANSIENT"}

# Severe detector/pipeline defect flag names that should be excluded from sky
# background estimation and flagged during photometric quality checks.
SEVERE_DEFECT_FLAG_NAMES = {
    "OVERFLOW",
    "NONFUNC",
    "MISSING_DATA",
    "HOT",
    "COLD",
    "OUTLIER",
    "PERSIST",
    "ghost",
    "GHOST",
}

# Fallback SPHEREx Level-2 bitmask for severe detector defects when FLAGS header
# cards (MSKNnnnn) are not supplied:
# Bit 1 (OVERFLOW), Bit 2 (SNOWBALL), Bit 4 (PERSIST), Bit 6 (NONFUNC),
# Bit 7 (OUTLIER), Bit 9 (MISSING_DATA), Bit 10 (HOT), Bit 11 (COLD), Bit 19 (PHAN_MISS).
# Deliberately excludes Bit 0 (TRANSIENT) and Bit 21 (SOURCE).
DEFAULT_NON_ASTROPHYSICAL_DEFECT_BITMASK = int(
    (1 << 1)
    | (1 << 2)
    | (1 << 4)
    | (1 << 6)
    | (1 << 7)
    | (1 << 9)
    | (1 << 10)
    | (1 << 11)
    | (1 << 19)
)


def _table_column_names(table: Any) -> set[str]:
    """Return binary-table column names for an Astropy FITS_rec/table array."""
    names = getattr(getattr(table, "dtype", None), "names", None)
    if names:
        return {str(name).upper() for name in names}
    columns = getattr(table, "columns", None)
    if columns is not None:
        try:
            return {str(name).upper() for name in columns.names}
        except Exception:
            pass
    return set()


def build_operation_flag_mask(
    flags: np.ndarray | None,
    flag_definitions: list[dict[str, Any]] | None = None,
    operation: str = "background_estimation",
) -> np.ndarray | None:
    """Construct a boolean mask (True = flagged/exclude) for a specific scientific operation.

    Operations:
    - ``"background_estimation"``: excludes both known astrophysical sources (``SOURCE``,
      ``TRANSIENT``) and detector/pipeline defect bits so the sky pedestal is estimated
      from unflagged background pixels where possible.
    - ``"detector_defects_only"``: excludes only non-astrophysical defect bits (preserving
      ``SOURCE`` and ``TRANSIENT`` pixels so stars and transients are never masked out).
    """
    if flags is None:
        return None
    flag_array = np.asarray(flags)
    if flag_array.ndim != 2 or not np.issubdtype(flag_array.dtype, np.integer):
        return None

    flags_u32 = flag_array.astype(np.uint32, copy=False)
    bitmask = np.uint32(0)

    if flag_definitions:
        for item in flag_definitions:
            bit = item.get("bit")
            if bit is None or not (0 <= int(bit) < 32):
                continue
            name = str(item.get("name", "")).strip().upper()
            if operation == "background_estimation":
                bitmask |= np.uint32(1 << int(bit))
            elif operation == "detector_defects_only":
                if name not in ASTROPHYSICAL_FLAG_NAMES:
                    bitmask |= np.uint32(1 << int(bit))
    else:
        if operation == "background_estimation":
            bitmask = np.uint32(DEFAULT_NON_ASTROPHYSICAL_DEFECT_BITMASK | (1 << 21) | (1 << 0))
        elif operation == "detector_defects_only":
            bitmask = np.uint32(DEFAULT_NON_ASTROPHYSICAL_DEFECT_BITMASK)

    return (flags_u32 & bitmask) != 0


def _decode_flag_header(flags: np.ndarray | None, flags_header: Any) -> dict[str, Any]:
    """Decode per-bit names/descriptions from MSKNnnnn/MSKDnnnn FITS cards.

    Per-bit counts are non-exclusive: a pixel can carry several bits. To prevent
    downstream consumers from summing non-exclusive fractions, this function also
    computes exact pixelwise bitwise-union counts for all flags and for
    non-astrophysical defect flags.
    """
    if flags is None or flags_header is None:
        return {
            "decoded": False,
            "definitions": [],
            "active_bits_without_definitions": [],
            "active_bits_outside_32bit_range": [],
            "any_flagged_pixels": None,
            "any_flagged_fraction": None,
            "non_astrophysical_flagged_pixels": None,
            "non_astrophysical_flagged_fraction": None,
            "non_astrophysical_bitmask": None,
        }

    flag_array = np.asarray(flags)
    if not np.issubdtype(flag_array.dtype, np.integer):
        return {
            "decoded": False,
            "definitions": [],
            "active_bits_without_definitions": [],
            "active_bits_outside_32bit_range": [],
            "any_flagged_pixels": None,
            "any_flagged_fraction": None,
            "non_astrophysical_flagged_pixels": None,
            "non_astrophysical_flagged_fraction": None,
            "non_astrophysical_bitmask": None,
        }

    # SPHEREx image FLAGS in this dataset are 32-bit integer bitmaps. Casting
    # signed int32 to uint32 preserves the underlying bit pattern.
    flags_u32 = flag_array.astype(np.uint32, copy=False)
    definitions_by_bit: dict[int, tuple[str, str]] = {}

    for key in flags_header.keys():
        match = re.fullmatch(r"MSKN(\d{4})", str(key).upper())
        if not match:
            continue
        bit = int(match.group(1))
        name = str(flags_header[key]).strip()
        description = str(flags_header.get(f"MSKD{bit:04d}", "")).strip()
        definitions_by_bit[bit] = (name, description)

    active_bits: list[int] = []
    for bit in range(32):
        mask = np.uint32(1 << bit)
        if np.any((flags_u32 & mask) != 0):
            active_bits.append(bit)

    details = []
    non_astro_mask = np.uint32(0)
    for bit, (name, description) in sorted(definitions_by_bit.items()):
        if bit >= 32:
            details.append({
                "bit": bit,
                "name": name,
                "description": description,
                "flagged_pixels": None,
                "flagged_fraction": None,
                "definition_only": True,
            })
            continue
        mask = np.uint32(1 << bit)
        if name.upper() not in ASTROPHYSICAL_FLAG_NAMES:
            non_astro_mask |= mask
        count = int(np.count_nonzero((flags_u32 & mask) != 0))
        details.append({
            "bit": bit,
            "name": name,
            "description": description,
            "flagged_pixels": count,
            "flagged_fraction": round(count / flag_array.size, 8) if flag_array.size else 0.0,
            "definition_only": count == 0,
        })

    any_flagged_pixels = int(np.count_nonzero(flags_u32 != 0))
    non_astro_pixels = int(np.count_nonzero((flags_u32 & non_astro_mask) != 0))

    return {
        "decoded": bool(definitions_by_bit),
        "definitions": details,
        "active_bits_without_definitions": [
            bit for bit in active_bits if bit not in definitions_by_bit
        ],
        "active_bits_outside_32bit_range": [
            bit for bit in definitions_by_bit if bit >= 32
        ],
        "any_flagged_pixels": any_flagged_pixels,
        "any_flagged_fraction": round(any_flagged_pixels / flag_array.size, 8) if flag_array.size else 0.0,
        "non_astrophysical_flagged_pixels": non_astro_pixels,
        "non_astrophysical_flagged_fraction": round(non_astro_pixels / flag_array.size, 8) if flag_array.size else 0.0,
        "non_astrophysical_bitmask": int(non_astro_mask),
    }


def validate_image_payload(
    image,
    flags,
    header,
    expected_detector=None,
    *,
    variance=None,
    zodi=None,
    epsf=None,
    wcs_wave=None,
    flags_header=None,
    require_full_product=False,
):
    """Validate image geometry, numerical arrays, FITS product extensions and WCS.

    ``valid`` means the core image product passed structural checks. EPSF and
    WCS-WAVE are reported as capability extensions: malformed/missing versions
    produce warnings and make those capabilities unavailable, but do not by
    themselves invalidate a usable 2-D image. Set ``require_full_product`` to
    require VARIANCE and ZODI, which are needed by the scientific processing
    path. This function never applies a quality-flag mask.
    """
    errors: list[str] = []
    structural_errors: list[str] = []
    numerical_errors: list[str] = []
    warnings: list[str] = []
    extension_reports: dict[str, Any] = {}

    # IMAGE ---------------------------------------------------------------
    image_array = None
    if image is None:
        msg = "IMAGE extension has no data."
        errors.append(msg)
        structural_errors.append(msg)
    else:
        try:
            image_array = np.asarray(image)
            if image_array.ndim != 2:
                msg = f"IMAGE data must be 2-D; received {image_array.ndim}-D."
                errors.append(msg)
                structural_errors.append(msg)
            elif image_array.size == 0 or min(image_array.shape) < 2:
                msg = "IMAGE data is empty or too small."
                errors.append(msg)
                structural_errors.append(msg)
            elif not np.issubdtype(image_array.dtype, np.number):
                msg = "IMAGE data must have a numeric data type."
                errors.append(msg)
                structural_errors.append(msg)
        except Exception as exc:
            msg = f"IMAGE data could not be interpreted: {type(exc).__name__}: {exc}"
            errors.append(msg)
            structural_errors.append(msg)

    image_valid_shape = (
        image_array is not None and image_array.ndim == 2 and image_array.size > 0
    )
    image_finite = None
    image_finite_fraction = None
    if image_valid_shape:
        image_finite = np.isfinite(image_array)
        image_finite_fraction = float(image_finite.mean())
        if not image_finite.any():
            msg = "IMAGE data contains no finite pixels."
            errors.append(msg)
            numerical_errors.append(msg)
        elif image_finite_fraction < 1.0:
            warnings.append(
                "IMAGE contains non-finite pixels; preserve this mask in all downstream processing."
            )

    # FLAGS ---------------------------------------------------------------
    flags_array = None
    flags_shape = None
    if flags is None:
        msg = "Required FLAGS extension is missing or has no data."
        errors.append(msg)
        structural_errors.append(msg)
    else:
        try:
            flags_array = np.asarray(flags)
            flags_shape = list(flags_array.shape)
            if flags_array.ndim != 2:
                msg = f"FLAGS data must be 2-D; received {flags_array.ndim}-D."
                errors.append(msg)
                structural_errors.append(msg)
            if not np.issubdtype(flags_array.dtype, np.integer):
                msg = "FLAGS data must have an integer data type."
                errors.append(msg)
                structural_errors.append(msg)
            if image_valid_shape and flags_array.shape != image_array.shape:
                msg = f"FLAGS shape {flags_array.shape} does not match IMAGE shape {image_array.shape}."
                errors.append(msg)
                structural_errors.append(msg)
        except Exception as exc:
            msg = f"FLAGS data could not be interpreted: {type(exc).__name__}: {exc}"
            errors.append(msg)
            structural_errors.append(msg)

    # DETECTOR / TIME ------------------------------------------------------
    detector = None
    mjd_avg = None
    try:
        detector = int(header["DETECTOR"])
        if detector not in range(1, 7):
            msg = f"DETECTOR must be between 1 and 6; received {detector}."
            errors.append(msg)
            structural_errors.append(msg)
        if expected_detector is not None and detector != int(expected_detector):
            msg = f"FITS DETECTOR={detector} does not match queried detector={expected_detector}."
            errors.append(msg)
            structural_errors.append(msg)
    except Exception:
        msg = "IMAGE header is missing a valid DETECTOR value."
        errors.append(msg)
        structural_errors.append(msg)

    try:
        mjd_avg = float(header["MJD-AVG"])
        if not np.isfinite(mjd_avg) or mjd_avg <= 0:
            msg = f"MJD-AVG must be a finite positive number; received {mjd_avg}."
            errors.append(msg)
            structural_errors.append(msg)
    except Exception:
        msg = "IMAGE header is missing a valid MJD-AVG value."
        errors.append(msg)
        structural_errors.append(msg)

    # IMAGE / VARIANCE / ZODI arrays --------------------------------------
    variance_array = None
    variance_finite_fraction = None
    variance_zero_count = None
    variance_negative_count = None
    mask_consistency = None
    if variance is None:
        extension_reports["VARIANCE"] = {"present": False, "valid": False}
        if require_full_product:
            msg = "Required VARIANCE extension is missing or has no data."
            errors.append(msg)
            structural_errors.append(msg)
    else:
        try:
            variance_array = np.asarray(variance)
            variance_errors = []
            variance_warnings = []
            if variance_array.ndim != 2 or not np.issubdtype(variance_array.dtype, np.number):
                variance_errors.append("VARIANCE must be a numeric 2-D array.")
            if image_valid_shape and variance_array.shape != image_array.shape:
                variance_errors.append(
                    f"VARIANCE shape {variance_array.shape} does not match IMAGE shape {image_array.shape}."
                )
            if variance_array.ndim == 2 and np.issubdtype(variance_array.dtype, np.number):
                variance_finite = np.isfinite(variance_array)
                variance_finite_fraction = float(variance_finite.mean())
                finite_variance_values = variance_array[variance_finite]
                variance_zero_count = int(np.count_nonzero(finite_variance_values == 0))
                variance_negative_count = int(np.count_nonzero(finite_variance_values < 0))
                if not variance_finite.any():
                    variance_errors.append("VARIANCE contains no finite pixels.")
                elif variance_finite_fraction < 1.0:
                    variance_warnings.append("VARIANCE contains non-finite pixels.")
                if variance_zero_count:
                    variance_warnings.append(
                        f"VARIANCE contains {variance_zero_count} zero-valued pixels; inverse-variance weighting must handle them explicitly."
                    )
                if variance_negative_count:
                    variance_warnings.append(
                        f"VARIANCE contains {variance_negative_count} negative pixels; those pixels are unsuitable for variance-weighted operations."
                    )
                if image_finite is not None and variance_array.shape == image_array.shape:
                    image_bad = ~image_finite
                    variance_bad = ~variance_finite
                    mask_consistency = {
                        "image_only_nonfinite": int(np.count_nonzero(image_bad & ~variance_bad)),
                        "variance_only_nonfinite": int(np.count_nonzero(variance_bad & ~image_bad)),
                        "both_nonfinite": int(np.count_nonzero(image_bad & variance_bad)),
                        "mask_disagreements": int(np.count_nonzero(image_bad != variance_bad)),
                    }
                    if mask_consistency["mask_disagreements"]:
                        variance_warnings.append(
                            "IMAGE and VARIANCE non-finite masks differ; construct a joint validity mask before uncertainty-aware processing."
                        )
            errors.extend(variance_errors)
            structural_errors.extend(variance_errors)
            warnings.extend(variance_warnings)
            extension_reports["VARIANCE"] = {
                "present": True,
                "valid": not variance_errors,
                "shape": list(variance_array.shape),
                "dtype": str(variance_array.dtype),
                "finite_fraction": round(variance_finite_fraction, 6) if variance_finite_fraction is not None else None,
                "zero_count": variance_zero_count,
                "negative_count": variance_negative_count,
                "image_mask_consistency": mask_consistency,
                "errors": variance_errors,
                "warnings": variance_warnings,
            }
        except Exception as exc:
            msg = f"VARIANCE data could not be interpreted: {type(exc).__name__}: {exc}"
            errors.append(msg)
            structural_errors.append(msg)
            extension_reports["VARIANCE"] = {"present": True, "valid": False}

    if zodi is None:
        extension_reports["ZODI"] = {"present": False, "valid": False}
        if require_full_product:
            msg = "Required ZODI extension is missing or has no data."
            errors.append(msg)
            structural_errors.append(msg)
    else:
        try:
            zodi_array = np.asarray(zodi)
            zodi_errors = []
            zodi_warnings = []
            if zodi_array.ndim != 2 or not np.issubdtype(zodi_array.dtype, np.number):
                zodi_errors.append("ZODI must be a numeric 2-D array.")
            if image_valid_shape and zodi_array.shape != image_array.shape:
                zodi_errors.append(
                    f"ZODI shape {zodi_array.shape} does not match IMAGE shape {image_array.shape}."
                )
            finite_fraction = (
                float(np.isfinite(zodi_array).mean())
                if zodi_array.ndim == 2 and np.issubdtype(zodi_array.dtype, np.number)
                else None
            )
            if finite_fraction is not None and finite_fraction < 1.0:
                zodi_warnings.append("ZODI contains non-finite pixels.")
            errors.extend(zodi_errors)
            structural_errors.extend(zodi_errors)
            warnings.extend(zodi_warnings)
            extension_reports["ZODI"] = {
                "present": True,
                "valid": not zodi_errors,
                "shape": list(zodi_array.shape),
                "dtype": str(zodi_array.dtype),
                "finite_fraction": round(finite_fraction, 6) if finite_fraction is not None else None,
                "errors": zodi_errors,
                "warnings": zodi_warnings,
            }
        except Exception as exc:
            msg = f"ZODI data could not be interpreted: {type(exc).__name__}: {exc}"
            errors.append(msg)
            structural_errors.append(msg)
            extension_reports["ZODI"] = {"present": True, "valid": False}

    # EPSF capability extension ------------------------------------------
    epsf_report: dict[str, Any] = {"present": epsf is not None, "valid": False}
    if epsf is None:
        epsf_report["warnings"] = ["EPSF extension is unavailable; PSF-assisted processing is unavailable."]
        warnings.extend(epsf_report["warnings"])
    else:
        try:
            names = _table_column_names(epsf)
            missing_columns = sorted(_REQUIRED_EPSF_COLUMNS - names)
            epsf_report["columns"] = sorted(names)
            epsf_report["missing_columns"] = missing_columns
            epsf_report["rows"] = int(len(epsf))
            epsf_warnings = []
            if missing_columns:
                epsf_warnings.append(f"EPSF is missing expected columns: {missing_columns}.")
            if len(epsf) != 441:
                epsf_warnings.append(f"Expected 441 EPSF rows for this documented product structure; found {len(epsf)}.")
            if "EPSF" in names:
                psf_values = np.asarray(epsf["EPSF"])
                epsf_report["array_shape"] = list(psf_values.shape)
                epsf_report["finite_fraction"] = float(np.isfinite(psf_values).mean()) if psf_values.size else 0.0
                if psf_values.shape != (441, 33, 33):
                    epsf_warnings.append(
                        f"EPSF array shape is {psf_values.shape}; documented DR1 structure is (441, 33, 33)."
                    )
                if not np.isfinite(psf_values).all():
                    epsf_warnings.append("EPSF contains non-finite values.")
            else:
                epsf_report["array_shape"] = None
            epsf_report["valid"] = not missing_columns and not epsf_warnings
            epsf_report["warnings"] = epsf_warnings
            warnings.extend(epsf_warnings)
        except Exception as exc:
            epsf_report["warnings"] = [f"EPSF table could not be inspected: {type(exc).__name__}: {exc}"]
            warnings.extend(epsf_report["warnings"])
    extension_reports["EPSF"] = epsf_report

    # WCS-WAVE capability extension --------------------------------------
    wave_report: dict[str, Any] = {"present": wcs_wave is not None, "valid": False}
    if wcs_wave is None:
        wave_report["warnings"] = ["WCS-WAVE extension is unavailable; lookup-table spectral metadata is unavailable."]
        warnings.extend(wave_report["warnings"])
    else:
        try:
            names = _table_column_names(wcs_wave)
            missing_columns = sorted(_REQUIRED_WAVE_COLUMNS - names)
            wave_warnings = []
            wave_report["columns"] = sorted(names)
            wave_report["missing_columns"] = missing_columns
            wave_report["rows"] = int(len(wcs_wave))
            if missing_columns:
                wave_warnings.append(f"WCS-WAVE is missing expected columns: {missing_columns}.")
            if len(wcs_wave) != 1:
                wave_warnings.append(f"Expected one WCS-WAVE row; found {len(wcs_wave)}.")
            if not missing_columns and len(wcs_wave) == 1:
                row = wcs_wave[0]
                x = np.asarray(row["X"])
                y = np.asarray(row["Y"])
                values = np.asarray(row["VALUES"])
                expected_shape = (len(y), len(x), 2)
                wave_report.update({
                    "x_shape": list(x.shape),
                    "y_shape": list(y.shape),
                    "values_shape": list(values.shape),
                })
                if x.ndim != 1 or y.ndim != 1 or values.shape != expected_shape:
                    wave_warnings.append(
                        f"WCS-WAVE shape mismatch: expected VALUES {expected_shape}, got {values.shape}."
                    )
                if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isfinite(values).all():
                    wave_warnings.append("WCS-WAVE contains non-finite lookup values.")
                if values.ndim == 3 and values.shape[-1] == 2:
                    wavelength = values[..., 0]
                    bandwidth = values[..., 1]
                    wave_report["wavelength_min_um"] = float(np.nanmin(wavelength)) if np.isfinite(wavelength).any() else None
                    wave_report["wavelength_max_um"] = float(np.nanmax(wavelength)) if np.isfinite(wavelength).any() else None
                    wave_report["bandwidth_min_um"] = float(np.nanmin(bandwidth)) if np.isfinite(bandwidth).any() else None
                    wave_report["bandwidth_max_um"] = float(np.nanmax(bandwidth)) if np.isfinite(bandwidth).any() else None
                    if np.any(np.isfinite(wavelength) & (wavelength <= 0)):
                        wave_warnings.append("WCS-WAVE contains non-positive wavelength values.")
                    if np.any(np.isfinite(bandwidth) & (bandwidth <= 0)):
                        wave_warnings.append("WCS-WAVE contains non-positive bandwidth values.")
            wave_report["valid"] = not wave_warnings
            wave_report["warnings"] = wave_warnings
            warnings.extend(wave_warnings)
        except Exception as exc:
            wave_report["warnings"] = [f"WCS-WAVE table could not be inspected: {type(exc).__name__}: {exc}"]
            warnings.extend(wave_report["warnings"])
    extension_reports["WCS-WAVE"] = wave_report

    # Celestial WCS --------------------------------------------------------
    wcs_valid = False
    if image_valid_shape:
        try:
            celestial_wcs = WCS(header, relax=True).celestial
            if not celestial_wcs.has_celestial or celestial_wcs.pixel_n_dim != 2:
                raise ValueError("WCS does not describe a two-dimensional celestial image.")
            center_x = (image_array.shape[1] - 1) / 2.0
            center_y = (image_array.shape[0] - 1) / 2.0
            world = celestial_wcs.pixel_to_world_values(center_x, center_y)
            if len(world) < 2 or not all(np.isfinite(float(value)) for value in world[:2]):
                raise ValueError("WCS maps the image center to non-finite coordinates.")
            wcs_valid = True
        except Exception as exc:
            errors.append(f"Invalid celestial WCS: {exc}")

    # Header-driven flag decoding -----------------------------------------
    flag_report = _decode_flag_header(flags_array, flags_header)
    if flags_array is not None and flags_header is None:
        warnings.append("FLAGS header was not supplied; flag bits could not be decoded.")
    if flag_report["active_bits_without_definitions"]:
        warnings.append(
            "Active FLAGS bits have no definition in this FITS header: "
            + ", ".join(map(str, flag_report["active_bits_without_definitions"]))
        )

    valid = not errors
    structurally_valid = not structural_errors
    numerically_valid = not numerical_errors and bool(
        image_finite_fraction is not None and image_finite_fraction > 0.0
    )
    variance_ok = extension_reports.get("VARIANCE", {}).get("valid", False)
    zodi_ok = extension_reports.get("ZODI", {}).get("valid", False)

    return {
        "valid": valid,
        "status": "valid" if valid else "invalid",
        "errors": errors,
        "warnings": warnings,
        "detector": detector,
        "mjd_avg": mjd_avg,
        "obsid": header.get("OBSID"),
        "product_date": header.get("DATE"),
        "image_bunit": header.get("BUNIT"),
        "image_shape": list(image_array.shape) if image_array is not None else None,
        "image_dtype": str(image_array.dtype) if image_array is not None else None,
        "finite_pixel_fraction": round(image_finite_fraction, 6) if image_finite_fraction is not None else None,
        "flags_shape": flags_shape,
        "quality_flags_decoded": bool(flag_report["decoded"]),
        "flag_definitions": flag_report["definitions"],
        "any_flagged_pixels": flag_report["any_flagged_pixels"],
        "any_flagged_fraction": flag_report["any_flagged_fraction"],
        "non_astrophysical_flagged_pixels": flag_report["non_astrophysical_flagged_pixels"],
        "non_astrophysical_flagged_fraction": flag_report["non_astrophysical_flagged_fraction"],
        "non_astrophysical_bitmask": flag_report["non_astrophysical_bitmask"],
        "active_bits_without_definitions": flag_report["active_bits_without_definitions"],
        "image_variance_mask_consistency": mask_consistency,
        "extensions": extension_reports,
        "wcs_valid": wcs_valid,
        "structural_validity": {
            "valid": structurally_valid,
            "errors": structural_errors,
        },
        "numerical_validity": {
            "valid": numerically_valid,
            "errors": numerical_errors,
            "image_finite_fraction": round(image_finite_fraction, 6) if image_finite_fraction is not None else None,
            "variance_finite_fraction": round(variance_finite_fraction, 6) if variance_finite_fraction is not None else None,
            "variance_zero_count": variance_zero_count,
            "variance_negative_count": variance_negative_count,
        },
        "alignment_validity": {
            "wcs_valid": wcs_valid,
            "reprojection_eligible": bool(structurally_valid and numerically_valid and wcs_valid),
        },
        "scientific_operation_suitability": {
            "visual_display": bool(valid and numerically_valid),
            "background_estimation": bool(valid and numerically_valid),
            "epoch_coadd": bool(valid and numerically_valid and wcs_valid),
            "difference_imaging": bool(valid and numerically_valid and wcs_valid),
            "uncertainty_weighted_photometry": bool(
                valid
                and numerically_valid
                and wcs_valid
                and variance_ok
                and (variance_negative_count or 0) == 0
            ),
            "zodi_background_subtraction": bool(valid and zodi_ok),
            "psf_fitting": bool(valid and wcs_valid and epsf_report.get("valid", False)),
            "wavelength_lookup": bool(valid and wave_report.get("valid", False)),
        },
        "flag_mask_policy": "none_applied; operation_specific_policy_required",
    }
