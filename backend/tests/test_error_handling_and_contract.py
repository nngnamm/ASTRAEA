import numpy as np
from fastapi.testclient import TestClient

from main import app
from serializers import package_cube_json

client = TestClient(app)


def test_targets_endpoint_includes_verified_availability_and_irsa_audit():
    resp = client.get("/api/targets")
    assert resp.status_code == 200
    payload = resp.json()
    assert payload["status"] == "ok"
    targets = payload["targets"]
    assert set(targets.keys()) == {
        "barnard",
        "eps_ind_b",
        "wise0855",
        "rw_aur",
        "r_dra",
    }
    expected_totals = {
        "barnard": 386,
        "eps_ind_b": 486,
        "wise0855": 263,
        "rw_aur": 342,
        "r_dra": 578,
    }
    for key, t in targets.items():
        assert t["data_status"] == "verified_available"
        assert "recommended_detectors" in t
        audit = t["irsa_audit"]
        assert audit["total_observations"] == expected_totals[key]
        assert audit["baseline_days"] > 350.0
        assert "by_detector" in audit
        assert "visit_groups_30d" in audit


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
                "collection": "spherex_qr2",
                "metadata": {"observation": {"observation_id": "OBS-1"}},
                "validation": {},
            },
            {
                "time": 60902.5,
                "detector": 5,
                "wave_um": 4.12,
                "aligned": img1,
                "normalized": img1 - 1.0,
                "collection": "spherex_qr2",
                "metadata": {"observation": {"observation_id": "OBS-2"}},
                "validation": {},
            },
            {
                "time": 61265.0,
                "detector": 5,
                "wave_um": 4.11,
                "aligned": img2,
                "normalized": img2 - 2.0,
                "collection": "spherex_qr2",
                "metadata": {"observation": {"observation_id": "OBS-3"}},
                "validation": {},
            },
        ],
        "validation_summary": {"aligned_frames": 3},
    }

    out = package_cube_json(cube)
    assert out["status"] == "ok"
    policy = out["epoch_grouping_policy"]
    assert policy["method"] == "time_gap_heuristic"
    assert policy["gap_days"] == 30
    assert policy["official_mission_epoch"] is False

    ep0 = out["epochs"][0]
    assert ep0["grouping_heuristic_gap_days"] == 30
    assert ep0["mjd_start"] == 60900.0
    assert ep0["mjd_end"] == 60902.5


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
