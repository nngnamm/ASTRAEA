import json

import numpy as np

from observation_metadata import extract_observation_metadata


def _sample_wave_table():
    dtype = [
        ("X", np.int32, (3,)),
        ("Y", np.int32, (2,)),
        ("VALUES", np.float32, (2, 3, 2)),
    ]
    table = np.zeros(1, dtype=dtype)
    table[0]["X"] = [0, 1, 2]
    table[0]["Y"] = [0, 1]
    table[0]["VALUES"][..., 0] = [[4.4, 4.5, 4.6], [4.41, 4.51, 4.61]]
    table[0]["VALUES"][..., 1] = 0.035
    return table


def test_extracts_observation_and_product_timestamps_separately():
    header = {
        "OBSID": "2025W33_2C_0321_1",
        "EXPIDN": 202533203038646,
        "DETECTOR": 6,
        "MJD-AVG": 60902.865380738,
        "DATE-AVG": "2025-08-15T20:46:08.000",
        "MJD-BEG": 60902.8647,
        "DATE-BEG": "2025-08-15T20:45:10.000",
        "MJD-END": 60902.8660,
        "DATE-END": "2025-08-15T20:47:05.000",
        "MJD-OBS": 60902.8647,
        "DATE-OBS": "2025-08-15T20:45:10.000",
        "TIMESYS": "UTC",
        "DATE": "2026-08-16T13:49:47.923",
        "XPOSURE": 113.5826,
        "TELAPSE": 118.1873,
        "TSAMP": 1.5349,
        "BUNIT": "MJy / sr",
        "RADESYS": "ICRS",
        "CTYPE1": "RA---TAN-SIP",
        "CTYPE2": "DEC--TAN-SIP",
        "CRVAL1": 269.9,
        "CRVAL2": 6.18,
        "L1DQAFLG": "Pass",
        "L2DQAFLG": "Pass",
        "ASTRONETSOLVED": True,
        "MATCH_NSTARS": 1889,
        "PSF_FWHM": 9.04,
        "L2 N_SOURCE": 3357522,
    }

    result = extract_observation_metadata(
        header,
        (58, 58),
        collection="spherex_dr1",
        wcs_valid=True,
        wcs_wave=_sample_wave_table(),
    )

    assert result["observation"]["observation_id"] == "2025W33_2C_0321_1"
    assert result["observation"]["detector"] == 6
    assert result["time"]["midpoint_mjd"] == 60902.865380738
    assert result["product"]["created_utc"].startswith("2026-")
    assert result["time"]["midpoint_utc"].startswith("2025-")
    assert result["exposure"]["integration_seconds"] == 113.5826
    assert result["product"]["image_unit"] == "MJy / sr"
    assert result["pipeline_quality"]["header_reported_l2_counters"]["L2 N_SOURCE"] == 3357522
    assert result["spectral_approximation"]["calibration_grade"] is False
    assert np.allclose(
        result["spectral_approximation"]["wavelength_range_um_approx"],
        [4.4, 4.61],
        rtol=1e-6,
    )
    json.dumps(result)  # Ensure the result is JSON-serializable.


def test_missing_optional_keywords_return_none():
    result = extract_observation_metadata({}, (5, 7))

    assert result["observation"]["observation_id"] is None
    assert result["time"]["midpoint_mjd"] is None
    assert result["product"]["image_unit"] is None
    assert result["pipeline_quality"]["astrometry_solved"] is None
    assert result["spectral_approximation"]["available"] is False
    json.dumps(result)


def test_approximate_wavelength_metadata_is_not_overstated():
    result = extract_observation_metadata(
        {"DETECTOR": 5},
        (58, 59),
        wcs_wave=_sample_wave_table(),
    )

    spectral = result["spectral_approximation"]
    assert spectral["available"] is True
    assert spectral["calibration_grade"] is False
    assert "approx" in spectral["warning"].lower()


def test_bad_wcs_wave_table_is_reported_as_unavailable():
    bad_table = np.zeros(1, dtype=[("VALUES", np.float32, (2, 2))])
    result = extract_observation_metadata({}, (10, 10), wcs_wave=bad_table)

    assert result["spectral_approximation"]["available"] is False
    assert "could not be summarized" in result["spectral_approximation"]["warning"]
