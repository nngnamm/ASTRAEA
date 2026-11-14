import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from serializers import compute_detector_display_limits
from spherex_core import (
    align_frame_to_grid,
    coadd,
    make_grid,
    normalize_frame_background,
)
from tests.test_image_alignment import _make_wcs_header


def test_normalize_frame_removes_zodi_and_preserves_source_flux_amplitude():
    coord = SkyCoord(ra=269.446294 * u.deg, dec=4.767881 * u.deg, frame="icrs")
    grid, n = make_grid(coord, fov_deg=0.05, pix_arcsec=6.15)
    hdr = _make_wcs_header(coord.ra.deg, coord.dec.deg, nx=n, ny=n, pix_arcsec=6.15)

    # Sky has 0.65 MJy/sr ZODI + 0.12 MJy/sr residual pedestal; star has +50.0 MJy/sr excess
    zodi = np.full((n, n), 0.65, dtype=np.float32)
    image = zodi + np.float32(0.12)
    c = n // 2
    image[c, c] += 50.0

    frame = {"image": image, "zodi": zodi, "header": hdr, "validation": {}}
    align_frame_to_grid(frame, coord, grid, shape_out=(n, n))
    normalize_frame_background(frame, grid=grid, shape_out=(n, n))

    norm_meta = frame["validation"]["normalization"]
    assert norm_meta["zodi_model_subtracted"] is True
    assert np.isclose(norm_meta["total_background_removed_mjy_sr"], 0.77, atol=1e-2)
    assert np.isclose(np.nanmedian(frame["normalized"]), 0.0, atol=1e-3)
    # Physical amplitude of the star above sky must be preserved in MJy/sr
    assert np.isclose(np.nanmax(frame["normalized"]), 50.0, rtol=1e-2)


def test_normalize_frame_excludes_flagged_source_and_defect_pixels_from_sky_pedestal():
    # Create a 10x10 cutout where 60% of pixels are contaminated by a +20.0 MJy/sr source halo
    # that is marked in FLAGS (bit 21 = SOURCE), while true unflagged sky is at 1.5 MJy/sr.
    aligned = np.full((10, 10), 1.5, dtype=float)
    flags = np.zeros((10, 10), dtype=np.int32)
    aligned[:6, :] = 21.5
    flags[:6, :] = 1 << 21

    frame = {
        "aligned": aligned,
        "valid_mask": np.isfinite(aligned),
        "flags": flags,
        "validation": {
            "flag_definitions": [
                {"bit": 21, "name": "SOURCE", "flagged_fraction": 0.60},
            ]
        },
    }
    normalize_frame_background(frame)
    norm_meta = frame["validation"]["normalization"]
    assert norm_meta["flag_aware_sky_estimation"] is True
    assert norm_meta["sky_estimation_mask_policy"] == "unflagged_sky_pixels"
    assert norm_meta["flagged_pixels_excluded_from_sky"] == 60
    # True sky (1.5 MJy/sr) is removed, not the biased median (21.5 MJy/sr)
    assert np.isclose(norm_meta["residual_pedestal_removed_mjy_sr"], 1.5, atol=1e-6)
    assert np.allclose(frame["normalized"][6:, :], 0.0, atol=1e-6)
    assert np.allclose(frame["normalized"][:6, :], 20.0, atol=1e-6)


def test_normalize_frame_preserves_nan_outside_valid_mask():
    aligned = np.full((10, 10), 2.5, dtype=float)
    aligned[:4, :] = np.nan
    valid_mask = np.isfinite(aligned)

    frame = {"aligned": aligned, "valid_mask": valid_mask, "validation": {}}
    normalize_frame_background(frame)

    assert np.all(np.isnan(frame["normalized"][:4, :]))
    assert np.allclose(frame["normalized"][4:, :], 0.0)


def test_coadd_uses_normalized_frames_to_prevent_partial_footprint_seams():
    # Frame 1: full coverage with background 0.20 MJy/sr
    f1_aligned = np.full((10, 10), 0.20, dtype=float)
    f1 = {"detector": 6, "time": 100.0, "aligned": f1_aligned, "validation": {}}
    normalize_frame_background(f1)

    # Frame 2: partial coverage (top half NaN) with high seasonal background 0.80 MJy/sr
    f2_aligned = np.full((10, 10), 0.80, dtype=float)
    f2_aligned[:5, :] = np.nan
    f2 = {"detector": 6, "time": 101.0, "aligned": f2_aligned, "validation": {}}
    normalize_frame_background(f2)

    combined = coadd([f1, f2], use_normalized=True)
    # Both top half and bottom half should sit smoothly at 0.0 MJy/sr with no step seam
    assert np.allclose(combined, 0.0, atol=1e-6)


def test_detector_display_limits_ignore_partial_coverage_frames():
    full_norm = np.zeros((10, 10), dtype=float)
    full_norm[5, 5] = 100.0
    partial_norm = np.random.default_rng(0).normal(0.0, 0.05, size=(10, 10))

    frames = [
        {
            "detector": 5,
            "aligned": full_norm,
            "normalized": full_norm,
            "validation": {"alignment": {"partial_coverage": False}},
        },
        {
            "detector": 5,
            "aligned": partial_norm,
            "normalized": partial_norm,
            "validation": {"alignment": {"partial_coverage": True}},
        },
    ]
    limits = compute_detector_display_limits(frames, low_pct=5.0, high_pct=99.0)
    vmin, vmax = limits[5]
    assert (
        vmax >= 1.0
    )  # Driven by full frame with the star, not the empty partial strip
