"""Tests for Approved Scientific Data Sources, Dataset Separation, Gaia DR3
astrometric reference metadata, dynamic SPHEREx collection discovery, and the
independent NEOSSat/CADC instrument adapter.
"""
from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from fastapi.testclient import TestClient

from approved_data_sources import (
    DATASET_IDENTITIES,
    GAIA_DR3_BENCHMARK_ASTROMETRY,
    NEOSSatCADCAdapter,
    OFFICIAL_SPHEREX_COLLECTIONS,
    discover_spherex_collections,
    infer_spherex_collection_from_url,
    mark_synthetic_fixture,
    propagate_astrometry_with_parallax,
    validate_dataset_separation,
)
from main import app
from spherex_core import resolve_object_name

client = TestClient(app)


def test_official_spherex_collections_and_url_provenance():
    assert OFFICIAL_SPHEREX_COLLECTIONS == (
        "spherex_dr1",
        "spherex_dr1_deep",
        "spherex_qr3",
        "spherex_qr3_deep",
    )
    assert (
        infer_spherex_collection_from_url(
            "https://irsa.ipac.caltech.edu/ibe/data/spherex/dr1_deep/level2/cutout.fits",
            fallback_collection="spherex_dr1",
        )
        == "spherex_dr1_deep"
    )
    assert (
        infer_spherex_collection_from_url(
            "https://irsa.ipac.caltech.edu/ibe/data/spherex/qr3_deep/level2/cutout.fits",
            fallback_collection="spherex_qr3",
        )
        == "spherex_qr3_deep"
    )
    assert (
        infer_spherex_collection_from_url(
            "https://irsa.ipac.caltech.edu/ibe/data/spherex/qr3/level2/cutout.fits",
            fallback_collection="spherex_dr1",
        )
        == "spherex_qr3"
    )


def test_discover_spherex_collections_dynamic_probing():
    coord = SkyCoord(ra=269.446294 * u.deg, dec=4.767881 * u.deg, frame="icrs")

    # 1. Offline mode does not invent results
    offline = discover_spherex_collections(coord=coord, allow_network=False)
    assert offline["network_queried"] is False
    assert offline["official_collections_investigated"] == list(
        OFFICIAL_SPHEREX_COLLECTIONS
    )
    for col in OFFICIAL_SPHEREX_COLLECTIONS:
        assert offline["collection_status"][col]["status"] == "unverified_offline"

    # 2. Live/mocked mode records per-collection availability without assuming all 4 have data
    def fake_query_sia(pos, collection):
        if collection in {"spherex_dr1", "spherex_qr3"}:
            return Table(
                {
                    "access_url": [
                        f"https://irsa.ipac.caltech.edu/ibe/data/spherex/{collection.replace('spherex_', '')}/level2/2025_1D5_spx.fits"
                    ],
                    "t_min": [60769.6],
                }
            )
        return Table({"access_url": [], "t_min": []})

    probed = discover_spherex_collections(
        coord=coord,
        fov_deg=0.1,
        allow_network=True,
        query_fn=fake_query_sia,
    )
    assert probed["network_queried"] is True
    assert probed["complete"] is True
    assert probed["collections_with_observations"] == ["spherex_dr1", "spherex_qr3"]
    assert probed["collection_status"]["spherex_dr1_deep"]["status"] == "ok_empty"
    assert probed["collection_status"]["spherex_qr3_deep"]["status"] == "ok_empty"


def test_gaia_dr3_astrometric_reference_and_non_gaia_fallback():
    # Barnard's Star has authoritative Gaia DR3 source_id, parallax, and radial velocity
    barnard = GAIA_DR3_BENCHMARK_ASTROMETRY["barnard"]
    assert barnard["dataset_identity"] == "Gaia"
    assert barnard["is_spherex_measurement"] is False
    assert barnard["gaia_dr3_present"] is True
    assert barnard["source_id"] == "4472832130942575872"
    assert barnard["parallax_mas"] == pytest.approx(546.9759)
    assert barnard["radial_velocity_km_s"] == pytest.approx(-110.51)
    assert barnard["astrometric_quality"]["ruwe"] == pytest.approx(1.08)

    # WISE 0855-0714 is absent from Gaia DR3 and explicitly documents its independent IR reference
    wise = GAIA_DR3_BENCHMARK_ASTROMETRY["wise0855"]
    assert wise["gaia_dr3_present"] is False
    assert wise["source_id"] is None
    assert "Luhman" in wise["independent_reference_source"]
    assert wise["reference_confidence"] == "high_peer_reviewed_ir_astrometry"
    assert wise["parallax_mas"] == pytest.approx(439.0)

    # Object resolution attaches astrometric_reference and parallax metadata
    res = resolve_object_name("Barnard's Star", to_mjd=61041.0, allow_external=False)
    assert res["gaia_dr3_source_id"] == "4472832130942575872"
    assert res["parallax_mas"] == pytest.approx(546.9759)
    assert res["annual_parallax_amplitude_arcsec"] == pytest.approx(0.547, abs=1e-3)
    assert res["astrometric_reference"]["is_spherex_measurement"] is False


def test_propagate_astrometry_with_parallax_accounts_for_3d_space_motion():
    barnard = GAIA_DR3_BENCHMARK_ASTROMETRY["barnard"]
    prop = propagate_astrometry_with_parallax(
        barnard["ra_deg"],
        barnard["dec_deg"],
        pm_ra_cosdec_mas_yr=barnard["pm_ra_cosdec_mas_yr"],
        pm_dec_mas_yr=barnard["pm_dec_mas_yr"],
        reference_epoch_jyear=barnard["reference_epoch_jyear"],
        to_mjd=61041.0,
        parallax_mas=barnard["parallax_mas"],
        radial_velocity_km_s=barnard["radial_velocity_km_s"],
    )
    assert prop["dataset_identity"] == "Gaia"
    assert prop["is_spherex_measurement"] is False
    assert prop["parallax_accounted_in_space_motion"] is True
    assert prop["annual_parallax_amplitude_arcsec"] == pytest.approx(0.547, abs=1e-3)
    assert prop["max_peak_to_peak_parallax_shift_arcsec"] == pytest.approx(1.094, abs=1e-3)


def test_neossat_cadc_adapter_and_cross_mission_separation():
    # Construct a synthetic NEOSSat-format FITS HDUList (optical broadband CCD in ADU)
    hdr = fits.Header()
    hdr["INSTRUME"] = "NEOSSat"
    hdr["TELESCOP"] = "NEOSSat / CSA"
    hdr["OBSID"] = "NEOSSAT_2025_00142"
    hdr["BUNIT"] = "ADU"
    hdr["DATE-OBS"] = "2025-04-04T14:30:00"
    hdr["EXPTIME"] = 100.0
    hdr["CTYPE1"] = "RA---TAN"
    hdr["CTYPE2"] = "DEC--TAN"
    hdr["CRVAL1"] = 269.4463
    hdr["CRVAL2"] = 4.7679
    hdr["CRPIX1"] = 16.0
    hdr["CRPIX2"] = 16.0
    hdr["CDELT1"] = -3.0 / 3600.0
    hdr["CDELT2"] = 3.0 / 3600.0

    img = np.ones((32, 32), dtype=np.float32) * 120.0
    hdul = fits.HDUList([fits.PrimaryHDU(data=img, header=hdr)])

    meta = NEOSSatCADCAdapter.inspect_fits_product(
        hdul,
        download_url="https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/NEOSSAT/NEOSSAT_2025_00142.fits",
    )
    assert meta["valid"] is True
    assert meta["dataset_identity"] == "NEOSSat"
    assert meta["is_spherex_measurement"] is False
    assert meta["archive_collection"] == "NEOSSAT"
    assert meta["observation_id"] == "NEOSSAT_2025_00142"
    assert meta["image_unit"] == "ADU"
    assert meta["pixel_scale_arcsec"] == pytest.approx(3.0, abs=0.05)
    assert "Open Government Licence" in meta["usage_and_citation"]["licence"]

    # Cross-mission comparison explicitly prohibits direct pixel subtraction or flux equivalence
    cmp_summary = NEOSSatCADCAdapter.build_cross_mission_comparison_summary(
        {"mjd": 60769.61, "product": {"image_unit": "MJy / sr"}},
        meta,
    )
    assert cmp_summary["direct_pixel_subtraction_allowed"] is False
    assert cmp_summary["flux_equivalence_claimed"] is False
    assert cmp_summary["timing_comparison"]["absolute_time_offset_days"] == pytest.approx(
        0.0058, abs=0.01
    )


def test_dataset_separation_and_synthetic_fixture_labeling():
    # Pure SPHEREx collection passes separation check
    res = validate_dataset_separation(
        [{"dataset_identity": "SPHEREx"}, {"dataset_identity": "SPHEREx"}]
    )
    assert res["valid_separation"] is True
    assert res["dataset_identity"] == "SPHEREx"

    # Mixing SPHEREx and NEOSSat raises ValueError
    with pytest.raises(ValueError, match="Dataset separation violation"):
        validate_dataset_separation(
            [{"dataset_identity": "SPHEREx"}, {"dataset_identity": "NEOSSat"}]
        )

    # Mixing Synthetic or WISE into an observational collection raises ValueError
    syn = mark_synthetic_fixture({"data": [1, 2, 3]}, purpose="unit_test_dipole")
    assert syn["dataset_identity"] == "Synthetic"
    assert syn["is_synthetic_test_fixture"] is True
    assert syn["allowed_in_production_responses"] is False
    with pytest.raises(ValueError, match="Dataset separation violation"):
        validate_dataset_separation([syn])


def test_sources_endpoint_and_wiseview_policy():
    resp = client.get("/api/sources")
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "ok"
    assert body["official_spherex_collections"] == list(OFFICIAL_SPHEREX_COLLECTIONS)
    assert set(body["approved_datasets"].keys()) == {
        "SPHEREx",
        "Gaia",
        "NEOSSat",
        "WISE",
        "Synthetic",
    }
    assert body["approved_datasets"]["WISE"]["is_substitute_for_spherex"] is False
    assert body["dataset_separation_rules"]["silent_cross_mission_merging_prohibited"] is True
