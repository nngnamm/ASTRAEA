import base64
import io

import matplotlib.image as mpimg
import numpy as np
from fastapi.testclient import TestClient

import main as main_module
from config import DEFAULT_REFERENCE_EPOCH_JYEAR, DEFAULT_REFERENCE_MJD
from main import app
from serializers import array_to_b64png, package_cube_json

client = TestClient(app, raise_server_exceptions=False)


def test_targets_endpoint_reports_honest_provenance_and_explicit_j2026_epoch():
    resp = client.get("/api/targets")
    assert resp.status_code == 200
    payload = resp.json()
    assert payload["status"] == "ok"

    ref_epoch = payload["reference_epoch"]
    assert ref_epoch["jyear"] == DEFAULT_REFERENCE_EPOCH_JYEAR == 2026.0
    assert ref_epoch["mjd"] == DEFAULT_REFERENCE_MJD == 61041.0
    assert ref_epoch["pm_convention"] == "pm_ra_cosdec"

    targets = payload["targets"]
    assert set(targets.keys()) == {
        "barnard",
        "eps_ind_b",
        "wise0855",
        "rw_aur",
        "r_dra",
    }

    # Only barnard has locally cached, SHA-256 verified FITS cutouts in this repository
    assert targets["barnard"]["data_status"] == "verified_local_cache_subset"
    barnard_audit = targets["barnard"]["provenance_audit"]
    assert barnard_audit["cached_fits_count"] == 30
    assert barnard_audit["collections_verified_in_cache"] == {
        "spherex_dr1": 19,
        "spherex_qr3": 11,
    }
    assert barnard_audit["by_detector_in_cache"] == {"5": 16, "6": 14}

    # Targets without local FITS cache must never be labeled verified_available
    for key in ("eps_ind_b", "wise0855", "rw_aur", "r_dra"):
        t = targets[key]
        assert t["data_status"] == "unverified_offline"
        assert "irsa_audit" not in t
        assert t["provenance_audit"]["cached_fits_count"] == 0
        assert t["propagated_epoch_jyear"] == 2026.0
        assert t["propagated_to_mjd"] == 61041.0
        assert t["pm_convention"] == "pm_ra_cosdec"


def test_targets_audit_endpoint_runs_reproducible_local_audit():
    resp = client.get("/api/targets/audit", params={"query_live_irsa": False})
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "ok"
    assert body["query_live_irsa"] is False
    assert body["reference_mjd"] == 61041.0

    barnard = body["targets"]["barnard"]
    assert barnard["verification_status"] == "verified_local_cache_subset"
    assert barnard["local_cache"]["valid_fits_count"] == 30
    assert barnard["comparison_eligible"] is True

    for key in ("eps_ind_b", "wise0855", "rw_aur", "r_dra"):
        target_audit = body["targets"][key]
        assert target_audit["verification_status"] == "unverified_offline"
        assert target_audit["sia_query"]["query_status"] == "not_queried_offline"
        # Offline skipped query must never be reported as 0 observations
        assert target_audit["sia_query"]["total_observations"] is None


def test_package_cube_json_empty_state_contract():
    empty_cube = {
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
        },
    }
    out = package_cube_json(empty_cube)
    assert out["status"] == "empty"
    assert out["error_code"] == "no_observations_found"
    assert out["error"] == "no SPHEREx images cover this position"
    assert out["frames"] == []
    assert out["epochs"] == []
    assert out["differences"] == []
    assert out["validation_summary"]["query_records"] == 0


def test_package_cube_json_includes_epoch_grouping_policy_and_epoch_bounds():
    img1 = np.ones((8, 8), dtype=float)
    img2 = np.ones((8, 8), dtype=float) * 2.0
    img2[4, 4] = 10.0

    cube = {
        "status": "ok",
        "center": (269.446294, 4.767881),
        "n": 8,
        "pix_arcsec": 6.15,
        "frames": [
            {
                "time": 60900.0,
                "detector": 5,
                "wave_um": 4.10,
                "aligned": img1,
                "normalized": img1 - 1.0,
                "collection": "spherex_dr1",
                "metadata": {"observation": {"observation_id": "OBS-1"}},
                "validation": {},
            },
            {
                "time": 60902.5,
                "detector": 5,
                "wave_um": 4.12,
                "aligned": img1,
                "normalized": img1 - 1.0,
                "collection": "spherex_dr1",
                "metadata": {"observation": {"observation_id": "OBS-2"}},
                "validation": {},
            },
            {
                "time": 61265.0,
                "detector": 5,
                "wave_um": 4.11,
                "aligned": img2,
                "normalized": img2 - 2.0,
                "collection": "spherex_qr3",
                "metadata": {"observation": {"observation_id": "OBS-3"}},
                "validation": {},
            },
        ],
        "validation_summary": {"aligned_frames": 3},
    }

    out = package_cube_json(cube)
    assert out["status"] == "ok"
    assert out["wcs_summary"]["reference_frame"] == "ICRS"
    assert out["wcs_summary"]["projection"] == "TAN"
    assert "icrs_tan_wcs_reprojection" in out["processing_history"]

    policy = out["epoch_grouping_policy"]
    assert policy["method"] == "time_gap_heuristic"
    assert policy["gap_days"] == 30
    assert policy["official_mission_epoch"] is False

    ep0 = out["epochs"][0]
    assert ep0["grouping_heuristic_gap_days"] == 30
    assert ep0["mjd_start"] == 60900.0
    assert ep0["mjd_end"] == 60902.5
    assert ep0["invalid_pixel_rendering"] == "transparent_rgba_alpha_0"


def test_array_to_b64png_renders_nan_pixels_as_transparent_alpha_zero():
    arr = np.ones((6, 6), dtype=float) * 5.0
    arr[2, 2] = 20.0
    arr[0, :] = np.nan  # Top row in array (becomes bottom row after flipud)

    uri = array_to_b64png(arr, vmin=0.0, vmax=20.0, cmap_name="gray")
    assert uri.startswith("data:image/png;base64,")
    raw_bytes = base64.b64decode(uri.split(",", 1)[1])
    rgba = mpimg.imread(io.BytesIO(raw_bytes))

    # After flipud, row 0 of input is row -1 of output PNG
    assert np.allclose(rgba[-1, :, 3], 0.0)
    # Finite rows must have full opacity (alpha = 1.0 in float RGBA)
    assert np.allclose(rgba[:-1, :, 3], 1.0)


def test_structured_error_responses_on_invalid_queries():
    resp_bad_frames = client.get(
        "/api/frames",
        params={"ra": 400.0, "dec": 4.7, "fov": 0.1, "detectors": "5,6"},
    )
    assert resp_bad_frames.status_code == 400
    err_frames = resp_bad_frames.json()
    assert err_frames["status"] == "error"
    assert err_frames["error_code"] == "invalid_parameters"
    assert "RA must satisfy" in err_frames["detail"]
    assert err_frames["message"] == err_frames["detail"]

    resp_unknown_obj = client.get(
        "/api/search/object",
        params={"name": "CompletelyUnknownObject-12345", "allow_external": False},
    )
    assert resp_unknown_obj.status_code == 404
    err_obj = resp_unknown_obj.json()
    assert err_obj["status"] == "error"
    assert err_obj["error_code"] == "object_not_found"
    assert "Could not resolve" in err_obj["detail"]


def test_frames_and_candidates_distinguish_upstream_failure_from_empty(monkeypatch):
    upstream_fail_cube = {
        "status": "error",
        "error_code": "upstream_service_unavailable",
        "frames": [],
        "message": "Upstream IRSA SIA archive service failed across all queried collections.",
        "validation_summary": {
            "query_records": 0,
            "sia_query_diagnostics": {"query_status": "upstream_failure", "complete": False},
        },
    }
    monkeypatch.setattr(
        main_module, "build_cube", lambda *args, **kwargs: upstream_fail_cube
    )

    resp_frames = client.get("/api/frames", params={"ra": 10.0, "dec": 20.0})
    assert resp_frames.status_code == 503
    body_frames = resp_frames.json()
    assert body_frames["status"] == "error"
    assert body_frames["error_code"] == "upstream_service_unavailable"

    resp_cands = client.get("/api/candidates", params={"ra": 10.0, "dec": 20.0})
    assert resp_cands.status_code == 503
    body_cands = resp_cands.json()
    assert body_cands["status"] == "error"
    assert body_cands["error_code"] == "upstream_service_unavailable"


def test_unhandled_exception_returns_sanitized_500(monkeypatch):
    def boom(*_args, **_kwargs):
        raise RuntimeError("/secret/internal/path/database_ failure")

    monkeypatch.setattr(main_module, "get_cache_status", boom)
    resp = client.get("/api/cache/status")
    assert resp.status_code == 500
    body = resp.json()
    assert body["status"] == "error"
    assert body["error_code"] == "internal_server_error"
    assert "/secret/internal/path" not in body["detail"]
