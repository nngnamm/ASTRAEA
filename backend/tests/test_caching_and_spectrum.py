import io
import os
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from fastapi.testclient import TestClient

import main as main_module
import spherex_core
from main import app

client = TestClient(app)


def _build_valid_fits_bytes(detector=5, mjd=61000.0, shape=(8, 8)) -> bytes:
    wcs = WCS(naxis=2)
    wcs.wcs.radesys = "ICRS"
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = ["deg", "deg"]
    wcs.wcs.crval = [269.446294, 4.767881]
    wcs.wcs.crpix = [4.5, 4.5]
    wcs.wcs.cdelt = [-6.15 / 3600.0, 6.15 / 3600.0]
    img_hdr = wcs.to_header()
    img_hdr["DETECTOR"] = int(detector)
    img_hdr["MJD-AVG"] = float(mjd)
    img_hdr["MJD-BEG"] = float(mjd - 0.001)
    img_hdr["MJD-END"] = float(mjd + 0.001)
    img_hdr["OBSID"] = f"TEST_OBS_D{detector}"
    img_hdr["BUNIT"] = "MJy / sr"

    flags_hdr = fits.Header()
    flags_hdr["MSKN0001"] = "OVERFLOW"
    flags_hdr["MSKD0001"] = "Overflow pixel"
    flags_hdr["MSKN0021"] = "SOURCE"
    flags_hdr["MSKD0021"] = "Source pixel"

    hdul = fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(data=np.ones(shape, dtype=np.float32) * 2.0, header=img_hdr, name="IMAGE"),
            fits.ImageHDU(data=np.zeros(shape, dtype=np.int32), header=flags_hdr, name="FLAGS"),
            fits.ImageHDU(data=np.ones(shape, dtype=np.float32) * 0.04, name="VARIANCE"),
            fits.ImageHDU(data=np.ones(shape, dtype=np.float32) * 0.5, name="ZODI"),
        ]
    )
    buf = io.BytesIO()
    hdul.writeto(buf)
    return buf.getvalue()


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
    assert first.diagnostics["query_status"] == "complete_with_observations"
    assert first.diagnostics["complete"] is True

    # Second call with use_cache=True must read from the JSON cache without calling Irsa.query_sia
    def fail_on_network(*_args, **_kwargs):
        raise AssertionError("Expected SIA cache hit without network query")

    monkeypatch.setattr(spherex_core.Irsa, "query_sia", fail_on_network)
    second = spherex_core.find_images(coord, fov_deg=0.1, use_cache=True)
    assert second == first
    assert second.diagnostics["from_cache"] is True


def test_find_images_does_not_cache_partial_multi_collection_failure(monkeypatch, tmp_path):
    monkeypatch.setattr(spherex_core, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(spherex_core, "COLLECTIONS", ["spherex_dr1", "spherex_qr3"])

    def partial_query_sia(pos, collection):
        if collection == "spherex_dr1":
            return [
                {
                    "access_url": "https://irsa.ipac.caltech.edu/SIA/level2_2025W34_1B_0557_2D5_spx.fits",
                    "t_min": 60913.1,
                }
            ]
        raise ConnectionError("IRSA spherex_qr3 service timeout")

    monkeypatch.setattr(spherex_core.Irsa, "query_sia", partial_query_sia)

    coord = SkyCoord(269.446294 * u.deg, 4.767881 * u.deg, frame="icrs")
    res = spherex_core.find_images(coord, fov_deg=0.1, use_cache=True)
    assert len(res) == 1
    assert res.diagnostics["query_status"] == "partial_failure"
    assert res.diagnostics["complete"] is False
    assert res.diagnostics["collection_status"]["spherex_dr1"]["status"] == "ok"
    assert res.diagnostics["collection_status"]["spherex_qr3"]["status"] == "error"

    # Partial failure must NEVER be cached to disk as a complete SIA result
    sia_files = list(tmp_path.glob("sia_*.json"))
    assert sia_files == []


def test_total_upstream_failure_is_not_mistaken_for_empty_observations(monkeypatch, tmp_path):
    monkeypatch.setattr(spherex_core, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(spherex_core, "COLLECTIONS", ["spherex_dr1", "spherex_qr3"])

    def failing_query_sia(pos, collection):
        raise ConnectionError(f"{collection} unreachable")

    monkeypatch.setattr(spherex_core.Irsa, "query_sia", failing_query_sia)

    cube = spherex_core.build_cube(150.0, 20.0, fov_deg=0.1, detectors=(5, 6))
    assert cube["status"] == "error"
    assert cube["error_code"] == "upstream_service_unavailable"
    assert (
        cube["validation_summary"]["sia_query_diagnostics"]["query_status"]
        == "upstream_failure"
    )
    assert list(tmp_path.glob("sia_*.json")) == []


def test_complete_empty_query_is_cached_and_reported_as_no_observations(monkeypatch, tmp_path):
    monkeypatch.setattr(spherex_core, "CACHE_DIR", str(tmp_path))
    monkeypatch.setattr(spherex_core, "COLLECTIONS", ["spherex_dr1", "spherex_qr3"])
    monkeypatch.setattr(spherex_core.Irsa, "query_sia", lambda pos, collection: [])

    cube = spherex_core.build_cube(150.0, 20.0, fov_deg=0.1, detectors=(5, 6))
    assert cube["status"] == "empty"
    assert cube["error_code"] == "no_observations_found"
    assert (
        cube["validation_summary"]["sia_query_diagnostics"]["query_status"]
        == "complete_empty"
    )
    assert len(list(tmp_path.glob("sia_*.json"))) == 1


def test_corrupted_cached_fits_is_evicted_and_recovered_via_redownload(monkeypatch, tmp_path):
    monkeypatch.setattr(spherex_core, "CACHE_DIR", str(tmp_path))
    coord = SkyCoord(269.446294 * u.deg, 4.767881 * u.deg, frame="icrs")
    meta = {
        "url": "https://irsa.ipac.caltech.edu/SIA/test_corrupt_recovery_D5_spx.fits",
        "mjd": 61000.0,
        "detector": 5,
        "collection": "spherex_dr1",
    }

    # Pre-populate cache with a corrupted/truncated FITS file
    key = spherex_core.hashlib.md5(
        f"{meta['url']}{coord.ra.deg:.5f}{coord.dec.deg:.5f}{0.1}".encode()
    ).hexdigest()
    corrupt_path = tmp_path / f"{key}.fits"
    corrupt_path.write_bytes(b"NOT_A_VALID_FITS_FILE_TRUNCATED")

    valid_bytes = _build_valid_fits_bytes(detector=5, mjd=61000.0)

    class FakeResponse:
        status_code = 200
        content = valid_bytes

        def raise_for_status(self):
            pass

    monkeypatch.setattr(spherex_core.requests, "get", lambda *args, **kwargs: FakeResponse())

    frame, report = spherex_core._load_frame_with_validation(meta, coord, fov_deg=0.1, max_retries=1)
    assert frame is not None
    assert report["valid"] is True
    assert report["corrupted_cache_evicted"] is True
    assert report["recovered_from_corrupted_cache"] is True
    assert report["sha256"] is not None
    assert corrupt_path.exists()
    assert corrupt_path.read_bytes() == valid_bytes


def test_invalid_fits_download_is_never_promoted_to_cache(monkeypatch, tmp_path):
    monkeypatch.setattr(spherex_core, "CACHE_DIR", str(tmp_path))
    coord = SkyCoord(269.446294 * u.deg, 4.767881 * u.deg, frame="icrs")
    meta = {
        "url": "https://irsa.ipac.caltech.edu/SIA/test_bad_download_D5_spx.fits",
        "mjd": 61000.0,
        "detector": 5,
        "collection": "spherex_dr1",
    }

    class BadHtmlResponse:
        status_code = 200
        content = b"<html>502 Bad Gateway</html>"

        def raise_for_status(self):
            pass

    monkeypatch.setattr(spherex_core.requests, "get", lambda *args, **kwargs: BadHtmlResponse())

    frame, report = spherex_core._load_frame_with_validation(meta, coord, fov_deg=0.1, max_retries=1)
    assert frame is None
    assert report["valid"] is False
    # Neither a .fits nor a leftover .tmp file may exist in CACHE_DIR
    assert list(tmp_path.glob("*.fits")) == []
    assert list(tmp_path.glob("*.tmp*")) == []


def test_get_cache_status_and_endpoint():
    status = spherex_core.get_cache_status()
    assert "cache_dir" in status
    # Must not leak absolute filesystem paths
    assert not os.path.isabs(status["cache_dir"])
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
