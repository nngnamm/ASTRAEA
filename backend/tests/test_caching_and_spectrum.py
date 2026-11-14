import astropy.units as u
from astropy.coordinates import SkyCoord
from fastapi.testclient import TestClient

import main as main_module
import spherex_core
from main import app

client = TestClient(app)


def test_find_images_writes_and_reads_local_sia_json_cache(monkeypatch, tmp_path):
    monkeypatch.setattr(spherex_core, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(spherex_core, "COLLECTIONS", ["spherex_dr1"])

    call_count = {"n": 0}

    def fake_query_sia(pos, collection):
        call_count["n"] += 1
        return [
            {
                "access_url": "https://irsa.ipac.caltech.edu/SIA/level2_2025W34_1B_0557_2D5_spx.fits",
                "t_min": 60913.1,
            },
            {
                "access_url": "https://irsa.ipac.caltech.edu/SIA/level2_2025W33_2C_0321_1D6_spx.fits",
                "t_min": 60906.2,
            },
        ]

    monkeypatch.setattr(spherex_core.Irsa, "query_sia", fake_query_sia)

    coord = SkyCoord(269.446294 * u.deg, 4.767881 * u.deg, frame="icrs")
    first = spherex_core.find_images(coord, fov_deg=0.1, use_cache=True)
    assert call_count["n"] == 1
    assert len(first) == 2
    assert [item["detector"] for item in first] == [6, 5]

    # Second call with use_cache=True must read from the JSON cache without calling Irsa.query_sia
    def fail_on_network(*_args, **_kwargs):
        raise AssertionError("Expected SIA cache hit without network query")

    monkeypatch.setattr(spherex_core.Irsa, "query_sia", fail_on_network)
    second = spherex_core.find_images(coord, fov_deg=0.1, use_cache=True)
    assert second == first


def test_get_cache_status_and_endpoint():
    status = spherex_core.get_cache_status()
    assert "cache_dir" in status
    assert status["fits_files_count"] >= 0
    assert status["sia_queries_count"] >= 0
    assert status["total_size_bytes"] >= 0
    assert status["total_size_mb"] >= 0.0

    resp = client.get("/api/cache/status")
    assert resp.status_code == 200
    body = resp.json()
    assert body["status"] == "ok"
    assert body["fits_files_count"] == status["fits_files_count"]


def test_spectrum_endpoint_returns_spectral_metadata_and_preserves_contract(monkeypatch):
    sample_spectrum = [
        {
            "detector": 5,
            "n_frames": 6,
            "wave_um": 4.142,
            "flux_mjy_sr": 18.45,
            "err_mjy_sr": 0.12,
            "snr": 153.8,
            "band_range_um": [3.83, 4.41],
            "sampled_wavelength_range_um": [3.8545, 4.4443],
            "calibration_grade": False,
            "wavelength_source": "WCS-WAVE_approximate",
            "aperture_radius_px": 2.5,
            "observation_ids": ["2025W34_1B_0557_2"],
        }
    ]
    monkeypatch.setattr(
        main_module, "compute_spectrum", lambda ra, dec, fov_deg=0.1: sample_spectrum
    )

    resp = client.get(
        "/api/spectrum", params={"ra": 269.446294, "dec": 4.767881, "fov": 0.1}
    )
    assert resp.status_code == 200
    payload = resp.json()
    assert payload["status"] == "ok"
    assert payload["center"] == {"ra": 269.446294, "dec": 4.767881}
    assert payload["fov_deg"] == 0.1
    assert payload["spectrum"] == sample_spectrum
    meta = payload["spectral_metadata"]
    assert meta["unit"] == "MJy / sr"
    assert meta["calibration_grade"] is False
    assert "Linear Variable Filters" in meta["scientific_note"]
