import hashlib
import json
import os
from pathlib import Path

import astropy.units as u
import pytest
from astropy.coordinates import SkyCoord

import spherex_core

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
PROVENANCE_JSON = REPO_ROOT / "audit" / "data_provenance.json"
SAMPLE_RESPONSE_JSON = REPO_ROOT / "backend" / "sample_barnard_response.json"


def test_all_30_cached_fits_files_match_provenance_ledger_and_validate():
    assert PROVENANCE_JSON.exists(), "audit/data_provenance.json must exist"
    ledger = json.loads(PROVENANCE_JSON.read_text(encoding="utf-8"))
    products = ledger["verified_cached_fits_products"]
    assert len(products) == 30

    col_counts = {"spherex_dr1": 0, "spherex_qr3": 0}
    det_counts = {5: 0, 6: 0}

    for entry in products:
        rel_path = entry["filename"]
        fits_path = REPO_ROOT / rel_path
        assert fits_path.exists(), f"Missing cached FITS file {rel_path}"

        # 1. Verify cryptographic SHA-256 checksum matches file on disk
        actual_sha256 = hashlib.sha256(fits_path.read_bytes()).hexdigest()
        assert actual_sha256 == entry["sha256"]

        # 2. Verify MD5 cache filename matches upstream IRSA IBE access_url + query params
        qp = entry["cutout_query_parameters"]
        preimage = (
            f"{entry['access_url']}"
            f"{qp['center_ra_deg']:.5f}"
            f"{qp['center_dec_deg']:.5f}"
            f"{qp['size_deg']}"
        )
        expected_md5 = hashlib.md5(preimage.encode("utf-8")).hexdigest()
        assert expected_md5 == entry["cache_key_md5"]

        # 3. Verify full SPHEREx Level-2 FITS structural & quality validation
        coord = SkyCoord(
            qp["center_ra_deg"] * u.deg, qp["center_dec_deg"] * u.deg, frame="icrs"
        )
        meta = {
            "url": entry["access_url"],
            "mjd": entry["timing"]["mjd_beg"],
            "detector": entry["detector"],
            "collection": entry["collection"],
            "cache_key": entry["cache_key_md5"],
        }
        frame, report = spherex_core._load_frame_with_validation(
            meta, coord, fov_deg=qp["size_deg"], max_retries=1
        )
        assert frame is not None
        assert report["valid"] is True
        assert report["from_cache"] is True
        assert report["sha256"] == entry["sha256"]
        assert frame["detector"] == entry["detector"]
        assert frame["metadata"]["observation"]["observation_id"] == entry["observation_id"]

        col_counts[entry["collection"]] += 1
        det_counts[entry["detector"]] += 1

    assert col_counts == {"spherex_dr1": 19, "spherex_qr3": 11}
    assert det_counts == {5: 16, 6: 14}


def test_sample_barnard_response_has_true_collections_and_unconfirmed_candidate_labels():
    assert SAMPLE_RESPONSE_JSON.exists()
    payload = json.loads(SAMPLE_RESPONSE_JSON.read_text(encoding="utf-8"))
    assert payload["status"] == "ok"
    assert payload["n_frames_total"] == 30

    cols = [f["collection"] for f in payload["frames"]]
    assert cols.count("spherex_dr1") == 19
    assert cols.count("spherex_qr3") == 11
    assert "spherex_qr2" not in cols

    for diff in payload["differences"]:
        assert diff["comparison_eligible"] is True
        assert diff["diagnostics"]["variance_propagated"] is True
        for cand in diff["candidates"]:
            assert "Confirmed Benchmark Object" not in cand["classification"]
            assert cand["independent_confirmation"] is False
            assert cand["is_calibrated_probability"] is False


@pytest.mark.live_irsa
@pytest.mark.skipif(
    os.environ.get("ASTRAEA_RUN_LIVE_IRSA") != "1",
    reason="Live IRSA network integration test disabled by default (set ASTRAEA_RUN_LIVE_IRSA=1 to run)",
)
def test_live_irsa_benchmark_audit():
    report = spherex_core.audit_all_benchmark_targets(query_live_irsa=True)
    assert report["query_live_irsa"] is True
    for key, item in report["targets"].items():
        assert item["sia_query"]["query_status"] in {
            "complete_with_observations",
            "complete_empty",
        }, f"Live IRSA query failed for {key}: {item['sia_query']}"
