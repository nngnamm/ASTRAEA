import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from fastapi.testclient import TestClient

import main as main_module
from analytics import compute_sky_health, detect_anomalies
from main import app
from spherex_core import make_grid

client = TestClient(app)


def test_compute_sky_health_does_not_penalize_source_flags():
    aligned = np.ones((20, 20), dtype=float)
    flags = np.zeros((20, 20), dtype=np.int32)
    cube = {
        "frames": [
            {
                "time": 60900.0,
                "aligned": aligned,
                "flags": flags,
                "validation": {
                    "flag_definitions": [
                        {"bit": 10, "name": "SOURCE", "flagged_fraction": 0.45},
                        {"bit": 11, "name": "TRANSIENT", "flagged_fraction": 0.05},
                        {"bit": 1, "name": "HOT", "flagged_fraction": 0.0},
                    ]
                },
            },
            {
                "time": 61080.0,
                "aligned": aligned,
                "flags": flags,
                "validation": {
                    "flag_definitions": [
                        {"bit": 10, "name": "SOURCE", "flagged_fraction": 0.40},
                        {"bit": 2, "name": "NONFUNC", "flagged_fraction": 0.0},
                    ]
                },
            },
        ]
    }
    health = compute_sky_health(cube)
    assert health["experimental_heuristic"] is True
    assert health["clean_detector_pct"] == 100.0
    assert health["valid_coverage_pct"] == 100.0
    assert health["score_pct"] == 100.0
    assert "not an official NASA SPHEREx quality metric" in health["provenance_note"]


def test_detect_anomalies_identifies_dipole_brightening_and_fading():
    coord = SkyCoord(269.446294 * u.deg, 4.767881 * u.deg, frame="icrs")
    grid, n = make_grid(coord, fov_deg=0.1, pix_arcsec=6.0)
    rng = np.random.default_rng(42)

    img_early = rng.normal(0.0, 0.02, size=(n, n))
    img_late = rng.normal(0.0, 0.02, size=(n, n))

    # 1. Moving dipole near (30, 30) -> (33, 30)
    img_early[30, 30] += 8.0
    img_late[33, 30] += 8.0

    # 2. Brightening variable at (15, 45)
    img_late[15, 45] += 7.0

    # 3. Fading variable at (45, 15)
    img_early[45, 15] += 7.0

    diff = img_late - img_early
    candidates = detect_anomalies(diff, img_early, img_late, grid, z_thresh=4.5, border=5)

    types = {c["type"] for c in candidates}
    assert "Moving Object Candidate" in types
    assert "Variable Source Candidate (Brightening)" in types
    assert "Variable Source Candidate (Fading)" in types

    for c in candidates:
        assert c["experimental_heuristic"] is True
        assert "score_disclaimer" in c
        assert c["anomaly_score"] == c["confidence"]


def test_detect_anomalies_ignores_partial_coverage_nan_boundary():
    coord = SkyCoord(269.446294 * u.deg, 4.767881 * u.deg, frame="icrs")
    grid, n = make_grid(coord, fov_deg=0.1, pix_arcsec=6.0)
    rng = np.random.default_rng(123)

    img_early = np.full((n, n), np.nan, dtype=float)
    img_late = np.full((n, n), np.nan, dtype=float)

    # Half-covered image with a non-zero background offset that would create a step if NaNs were zero-filled
    img_early[:, : n // 2] = 5.0 + rng.normal(0.0, 0.01, size=(n, n // 2))
    img_late[:, : n // 2] = 10.0 + rng.normal(0.0, 0.01, size=(n, n // 2))
    diff = img_late - img_early

    candidates = detect_anomalies(diff, img_early, img_late, grid, z_thresh=6.0, border=5)
    assert candidates == []


def test_candidates_endpoint_returns_structured_payload(monkeypatch):
    coord = SkyCoord(269.446294 * u.deg, 4.767881 * u.deg, frame="icrs")
    grid, n = make_grid(coord, fov_deg=0.1, pix_arcsec=6.0)
    rng = np.random.default_rng(7)
    early_img = rng.normal(0.0, 0.02, size=(n, n))
    late_img = rng.normal(0.0, 0.02, size=(n, n))
    early_img[n // 2, n // 2] += 10.0
    late_img[n // 2 + 2, n // 2] += 10.0

    fake_cube = {
        "status": "ok",
        "center": (269.446294, 4.767881),
        "grid": grid,
        "n": n,
        "pix_arcsec": 6.0,
        "frames": [
            {
                "time": 60900.0,
                "detector": 5,
                "wave_um": 4.10,
                "aligned": early_img,
                "normalized": early_img,
                "flags": np.zeros((n, n), dtype=np.int32),
                "collection": "spherex_qr2",
                "metadata": {"observation": {"observation_id": "OBS-E"}},
                "validation": {},
            },
            {
                "time": 61265.0,
                "detector": 5,
                "wave_um": 4.11,
                "aligned": late_img,
                "normalized": late_img,
                "flags": np.zeros((n, n), dtype=np.int32),
                "collection": "spherex_qr2",
                "metadata": {"observation": {"observation_id": "OBS-L"}},
                "validation": {},
            },
        ],
        "validation_summary": {"aligned_frames": 2},
    }
    monkeypatch.setattr(main_module, "build_cube", lambda *args, **kwargs: fake_cube)

    resp = client.get(
        "/api/candidates",
        params={"ra": 269.446294, "dec": 4.767881, "fov": 0.1, "detectors": "5"},
    )
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "ok"
    assert "5" in body["candidates_by_detector"]
    assert len(body["candidates_by_detector"]["5"]) >= 1
    assert body["sky_health"]["experimental_heuristic"] is True
    assert "Heuristic candidate detections" in body["disclaimer"]
