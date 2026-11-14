import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time
from fastapi.testclient import TestClient

from config import DEFAULT_REFERENCE_EPOCH_JYEAR, DEFAULT_REFERENCE_MJD, TARGETS
from main import app
import spherex_core
from spherex_core import (
    propagate_coords,
    propagate_coords_with_metadata,
    resolve_object_name,
    validate_query_parameters,
)

client = TestClient(app)


def test_propagate_coords_exact_spherical_space_motion():
    for key in ("barnard", "wise0855", "eps_ind_b"):
        target = TARGETS[key]
        to_mjd = 61000.0
        ra1, dec1 = propagate_coords(
            target["ra"],
            target["dec"],
            target["pm_ra_arcsec_yr"],
            target["pm_dec_arcsec_yr"],
            epoch0=target["epoch"],
            to_mjd=to_mjd,
        )
        c0 = SkyCoord(
            ra=float(target["ra"]) * u.deg,
            dec=float(target["dec"]) * u.deg,
            pm_ra_cosdec=float(target["pm_ra_arcsec_yr"]) * u.arcsec / u.yr,
            pm_dec=float(target["pm_dec_arcsec_yr"]) * u.arcsec / u.yr,
            obstime=Time(float(target["epoch"]), format="jyear"),
            frame="icrs",
        )
        c1 = c0.apply_space_motion(new_obstime=Time(to_mjd, format="mjd"))
        assert np.isclose(ra1, float(c1.ra.deg), atol=1e-10)
        assert np.isclose(dec1, float(c1.dec.deg), atol=1e-10)

    # Verify Barnard's Star moves northward by ~269 arcsec (~0.0748 deg) from J2000 to MJD 61000
    ra_b, dec_b = propagate_coords(
        TARGETS["barnard"]["ra"],
        TARGETS["barnard"]["dec"],
        TARGETS["barnard"]["pm_ra_arcsec_yr"],
        TARGETS["barnard"]["pm_dec_arcsec_yr"],
        epoch0=2000.0,
        to_mjd=61000.0,
    )
    assert np.isclose(ra_b, 269.446294, atol=1e-5)
    assert np.isclose(dec_b, 4.767881, atol=1e-5)


def test_propagate_coords_default_epoch_is_exact_j2026_not_mjd_61000():
    # J2026.0 is MJD 61041.0 (2026-01-01T00:00:00.000 UTC), 41 days after MJD 61000.0
    assert np.isclose(float(Time(2026.0, format="jyear").mjd), DEFAULT_REFERENCE_MJD)
    meta = propagate_coords_with_metadata(
        TARGETS["barnard"]["ra"],
        TARGETS["barnard"]["dec"],
        TARGETS["barnard"]["pm_ra_arcsec_yr"],
        TARGETS["barnard"]["pm_dec_arcsec_yr"],
        epoch0=2000.0,
    )
    assert meta["propagated_to_mjd"] == 61041.0
    assert meta["propagated_epoch_jyear"] == DEFAULT_REFERENCE_EPOCH_JYEAR == 2026.0
    assert meta["pm_convention"] == "pm_ra_cosdec"
    assert meta["proper_motion_applied"] is True
    # At J2026.0, Barnard's Star is ~1.16 arcsec further North than at MJD 61000.0
    assert meta["dec"] > 4.7681


def test_validate_query_parameters_accepts_decimal_and_sexagesimal():
    res_dec = validate_query_parameters(
        ra=269.446294,
        dec=4.767881,
        fov_deg=0.1,
        detectors="5,6",
        max_frames=30,
    )
    assert np.isclose(res_dec["ra"], 269.446294)
    assert np.isclose(res_dec["dec"], 4.767881)
    assert res_dec["fov_deg"] == 0.1
    assert res_dec["detectors"] == (5, 6)
    assert res_dec["max_frames"] == 30

    # WISE 0855 sexagesimal J2010.5 position: 08h55m10.831s, -7d14m42.53s
    res_sex = validate_query_parameters(
        ra="08h55m10.831s",
        dec="-7d14m42.53s",
        fov_deg=0.15,
        detectors=[6, 5, 1],
        max_frames=25,
    )
    assert np.isclose(res_sex["ra"], 133.795129, atol=1e-5)
    assert np.isclose(res_sex["dec"], -7.245147, atol=1e-5)
    assert res_sex["detectors"] == (1, 5, 6)


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"ra": 365.0, "dec": 0.0}, "RA must satisfy"),
        ({"ra": -1.0, "dec": 0.0}, "RA must satisfy"),
        ({"ra": 180.0, "dec": -95.0}, "Dec must satisfy"),
        ({"ra": 180.0, "dec": 91.0}, "Dec must satisfy"),
        ({"ra": 180.0, "dec": 0.0, "fov_deg": 0.01}, "fov_deg must satisfy"),
        ({"ra": 180.0, "dec": 0.0, "fov_deg": 0.50}, "fov_deg must satisfy"),
        ({"ra": 180.0, "dec": 0.0, "detectors": [0]}, "out of valid SPHEREx range"),
        ({"ra": 180.0, "dec": 0.0, "detectors": [7]}, "out of valid SPHEREx range"),
        ({"ra": 180.0, "dec": 0.0, "detectors": "abc"}, "Invalid detector ID"),
    ],
)
def test_validate_query_parameters_rejects_invalid_inputs(kwargs, match):
    with pytest.raises(ValueError, match=match):
        validate_query_parameters(**kwargs)


@pytest.mark.parametrize(
    "query_alias,expected_key,expected_name",
    [
        ("Barnard's Star", "barnard", "Barnard's Star"),
        ("GJ 699", "barnard", "Barnard's Star"),
        ("WISE 0855", "wise0855", "WISE 0855-0714"),
        ("eps Ind B", "eps_ind_b", "Epsilon Indi B"),
        ("RW Aur", "rw_aur", "RW Aurigae"),
        ("R Dra", "r_dra", "R Draconis"),
    ],
)
def test_resolve_object_name_benchmark_catalog_offline(
    query_alias, expected_key, expected_name
):
    resolved = resolve_object_name(query_alias, to_mjd=61000.0, allow_external=False)
    assert resolved["target_key"] == expected_key
    assert resolved["resolved_name"] == expected_name
    assert resolved["source"] == "nasa_benchmark_catalog"
    assert resolved["benchmark_target"] is True
    assert resolved["pm_convention"] == "pm_ra_cosdec"
    assert resolved["proper_motion_applied"] is True
    assert 0.0 <= resolved["ra"] < 360.0
    assert -90.0 <= resolved["dec"] <= 90.0


def test_resolve_object_name_sesame_fallback_does_not_apply_proper_motion_twice(monkeypatch):
    monkeypatch.setattr(
        spherex_core.SkyCoord,
        "from_name",
        lambda name: SkyCoord(10.684708 * u.deg, 41.26875 * u.deg, frame="icrs"),
    )
    resolved = resolve_object_name("M31", allow_external=True)
    assert resolved["source"] == "astropy_sesame"
    assert resolved["benchmark_target"] is False
    assert resolved["proper_motion_applied"] is False
    assert resolved["propagated_to_mjd"] is None


def test_resolve_object_name_unknown_and_empty_handling():
    with pytest.raises(ValueError, match="Object name must not be empty"):
        resolve_object_name("   ", allow_external=False)

    with pytest.raises(LookupError, match="Could not resolve astronomical object"):
        resolve_object_name("NonExistentObjectXYZ-9999", allow_external=False)


def test_fastapi_search_and_validation_endpoints():
    resp_obj = client.get(
        "/api/search/object",
        params={"name": "Barnard's Star", "to_mjd": 61000.0, "allow_external": False},
    )
    assert resp_obj.status_code == 200
    body_obj = resp_obj.json()
    assert body_obj["target_key"] == "barnard"
    assert body_obj["source"] == "nasa_benchmark_catalog"
    assert np.isclose(body_obj["ra"], 269.446294, atol=1e-5)

    resp_not_found = client.get(
        "/api/search/object",
        params={"name": "UnknownTarget9999", "allow_external": False},
    )
    assert resp_not_found.status_code == 404

    resp_coord = client.get(
        "/api/search/coordinates",
        params={"ra": "08h55m10.831s", "dec": "-7d14m42.53s", "fov": 0.1, "detectors": "5,6"},
    )
    assert resp_coord.status_code == 200
    body_coord = resp_coord.json()
    assert body_coord["valid"] is True
    assert body_coord["operation"] == "coordinate_validation_and_benchmark_lookup"
    assert body_coord["irsa_queried"] is False
    assert body_coord["data_availability_checked"] is False
    assert np.isclose(body_coord["ra"], 133.795129, atol=1e-5)
    assert np.isclose(body_coord["dec"], -7.245147, atol=1e-5)
    assert body_coord["nearby_benchmark_target"]["target_key"] == "wise0855"
    assert body_coord["nearby_benchmark_target"]["propagated_epoch_jyear"] == 2026.0
    assert body_coord["nearby_benchmark_target"]["propagated_to_mjd"] == 61041.0

    for bad_params in (
        {"ra": 365.0, "dec": 4.7},
        {"ra": 269.4, "dec": -95.0},
        {"ra": 269.4, "dec": 4.7, "fov": 0.9},
        {"ra": 269.4, "dec": 4.7, "detectors": "7"},
    ):
        resp_bad = client.get("/api/frames", params=bad_params)
        assert resp_bad.status_code == 400
