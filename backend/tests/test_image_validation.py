import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from image_validation import validate_image_payload


def make_valid_header(detector=5, mjd=61000.0):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [269.0, 4.0]
    wcs.wcs.crpix = [2.0, 2.0]
    wcs.wcs.cdelt = [-0.001, 0.001]
    header = wcs.to_header()
    header["DETECTOR"] = detector
    header["MJD-AVG"] = mjd
    header["OBSID"] = "TEST_OBSERVATION"
    return header


def make_flag_header():
    header = fits.Header()
    header["MSKN0000"] = "TRANSIENT"
    header["MSKD0000"] = "Transient detected during SUR"
    header["MSKN0021"] = "SOURCE"
    header["MSKD0021"] = "Pixel mapped to a known source"
    return header


def make_epsf():
    dtype = [
        ("BINX", "i4"), ("BINY", "i4"),
        ("XCENTER", "f4"), ("YCENTER", "f4"),
        ("XWIDTH", "f4"), ("YWIDTH", "f4"),
        ("NSTAR", "i4"), ("NEFF_MEAN", "f4"),
        ("CWAVE", "f4"), ("CBAND", "f4"),
        ("EPSF", "f4", (33, 33)),
    ]
    table = np.zeros(441, dtype=dtype)
    table["EPSF"] = 1.0
    table["CWAVE"] = 4.0
    table["CBAND"] = 0.04
    return table


def make_wcs_wave(n=9):
    dtype = [("X", "i4", (n,)), ("Y", "i4", (n,)), ("VALUES", "f4", (n, n, 2))]
    table = np.zeros(1, dtype=dtype)
    table["X"][0] = np.arange(n)
    table["Y"][0] = np.arange(n)
    table["VALUES"][0, ..., 0] = 4.5
    table["VALUES"][0, ..., 1] = 0.04
    return table


def test_valid_image_flags_headers_and_wcs_pass():
    image = np.ones((4, 4), dtype=np.float32)
    flags = np.zeros((4, 4), dtype=np.int32)
    result = validate_image_payload(
        image, flags, make_valid_header(), expected_detector=5,
        flags_header=make_flag_header(),
    )
    assert result["valid"] is True
    assert result["wcs_valid"] is True
    assert result["finite_pixel_fraction"] == 1.0
    assert result["quality_flags_decoded"] is True
    assert {item["name"] for item in result["flag_definitions"]} == {"TRANSIENT", "SOURCE"}


def test_non_finite_image_pixels_are_reported():
    image = np.ones((4, 4), dtype=np.float32)
    image[0, 0] = np.nan
    flags = np.zeros((4, 4), dtype=np.int32)
    result = validate_image_payload(image, flags, make_valid_header(), expected_detector=5)
    assert result["valid"] is True
    assert result["finite_pixel_fraction"] < 1.0
    assert any("non-finite" in warning for warning in result["warnings"])


def test_flags_shape_mismatch_fails_validation():
    image = np.ones((4, 4), dtype=np.float32)
    flags = np.zeros((3, 4), dtype=np.int32)
    result = validate_image_payload(image, flags, make_valid_header(), expected_detector=5)
    assert result["valid"] is False
    assert any("does not match IMAGE shape" in error for error in result["errors"])


def test_detector_mismatch_fails_validation():
    image = np.ones((4, 4), dtype=np.float32)
    flags = np.zeros((4, 4), dtype=np.int32)
    result = validate_image_payload(image, flags, make_valid_header(detector=6), expected_detector=5)
    assert result["valid"] is False
    assert any("does not match queried detector" in error for error in result["errors"])


def test_missing_wcs_fails_validation():
    image = np.ones((4, 4), dtype=np.float32)
    flags = np.zeros((4, 4), dtype=np.int32)
    header = fits.Header()
    header["DETECTOR"] = 5
    header["MJD-AVG"] = 61000.0
    result = validate_image_payload(image, flags, header, expected_detector=5)
    assert result["valid"] is False
    assert result["wcs_valid"] is False
    assert any("Invalid celestial WCS" in error for error in result["errors"])


def test_missing_flags_fails_validation():
    image = np.ones((4, 4), dtype=np.float32)
    result = validate_image_payload(image, None, make_valid_header(), expected_detector=5)
    assert result["valid"] is False
    assert any("FLAGS extension" in error for error in result["errors"])


def test_full_product_extensions_and_dynamic_wcs_wave_grid_validate():
    image = np.ones((4, 4), dtype=np.float32)
    flags = np.zeros((4, 4), dtype=np.int32)
    flags[0, 0] = (1 << 0) | (1 << 21)
    result = validate_image_payload(
        image, flags, make_valid_header(), expected_detector=5,
        variance=np.ones((4, 4), dtype=np.float32),
        zodi=np.ones((4, 4), dtype=np.float32),
        epsf=make_epsf(),
        wcs_wave=make_wcs_wave(11),
        flags_header=make_flag_header(),
        require_full_product=True,
    )
    assert result["valid"] is True
    assert result["extensions"]["EPSF"]["valid"] is True
    assert result["extensions"]["WCS-WAVE"]["valid"] is True
    assert result["extensions"]["WCS-WAVE"]["values_shape"] == [11, 11, 2]
    counts = {item["name"]: item["flagged_pixels"] for item in result["flag_definitions"]}
    assert counts["TRANSIENT"] == 1
    assert counts["SOURCE"] == 1


def test_image_variance_nonfinite_masks_are_compared():
    image = np.ones((4, 4), dtype=np.float32)
    variance = np.ones((4, 4), dtype=np.float32)
    image[0, 0] = np.nan
    variance[0, 0] = np.nan
    flags = np.zeros((4, 4), dtype=np.int32)
    result = validate_image_payload(
        image, flags, make_valid_header(), variance=variance,
        zodi=np.ones((4, 4)), flags_header=make_flag_header(),
    )
    assert result["valid"] is True
    assert result["image_variance_mask_consistency"]["both_nonfinite"] == 1
    assert result["image_variance_mask_consistency"]["mask_disagreements"] == 0


def test_source_flag_is_reported_but_not_automatically_masked():
    image = np.ones((4, 4), dtype=np.float32)
    flags = np.zeros((4, 4), dtype=np.int32)
    flags[1, 1] = 1 << 21
    result = validate_image_payload(
        image, flags, make_valid_header(), flags_header=make_flag_header()
    )
    assert result["valid"] is True
    assert result["flag_mask_policy"] == "none_applied; operation_specific_policy_required"
    source = next(item for item in result["flag_definitions"] if item["name"] == "SOURCE")
    assert source["flagged_pixels"] == 1


def test_variance_shape_mismatch_fails_full_product_validation():
    image = np.ones((4, 4), dtype=np.float32)
    flags = np.zeros((4, 4), dtype=np.int32)
    result = validate_image_payload(
        image, flags, make_valid_header(), variance=np.ones((3, 4)),
        zodi=np.ones((4, 4)), require_full_product=True,
    )
    assert result["valid"] is False
    assert any("VARIANCE shape" in error for error in result["errors"])


def test_negative_and_zero_variance_reported_and_suitability_separated():
    image = np.ones((4, 4), dtype=np.float32)
    flags = np.zeros((4, 4), dtype=np.int32)
    variance = np.ones((4, 4), dtype=np.float32)
    variance[0, 0] = 0.0
    variance[0, 1] = -1.5

    result = validate_image_payload(
        image,
        flags,
        make_valid_header(),
        variance=variance,
        zodi=np.ones((4, 4), dtype=np.float32),
        require_full_product=True,
    )
    assert result["valid"] is True
    assert result["structural_validity"]["valid"] is True
    assert result["numerical_validity"]["valid"] is True
    assert result["numerical_validity"]["variance_zero_count"] == 1
    assert result["numerical_validity"]["variance_negative_count"] == 1
    assert result["alignment_validity"]["wcs_valid"] is True
    # Negative variance makes uncertainty_weighted_photometry unsuitable while visual_display remains suitable
    assert result["scientific_operation_suitability"]["visual_display"] is True
    assert result["scientific_operation_suitability"]["uncertainty_weighted_photometry"] is False
