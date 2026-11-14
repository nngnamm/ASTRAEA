import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from spherex_core import (
    align_frame_to_grid,
    compute_reference_pixel_scale,
    evaluate_cross_frame_registration,
    make_grid,
)


def _make_wcs_header(ra_deg, dec_deg, nx=20, ny=20, pix_arcsec=6.15, angle_deg=0.0):
    wcs = WCS(naxis=2)
    wcs.wcs.radesys = "ICRS"
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.cunit = ["deg", "deg"]
    wcs.wcs.crval = [ra_deg, dec_deg]
    wcs.wcs.crpix = [(nx + 1) / 2.0, (ny + 1) / 2.0]
    scale_deg = pix_arcsec / 3600.0
    theta = np.radians(angle_deg)
    cos_t, sin_t = np.cos(theta), np.sin(theta)
    # Apply rotation to CD matrix (with negative RA axis convention)
    wcs.wcs.cd = np.array(
        [
            [-scale_deg * cos_t, scale_deg * sin_t],
            [-scale_deg * sin_t, scale_deg * cos_t],
        ]
    )
    return wcs.to_header()


def test_make_grid_constructs_icrs_north_up_tan_wcs():
    coord = SkyCoord(ra=269.446294 * u.deg, dec=4.767881 * u.deg, frame="icrs")
    grid, n = make_grid(coord, fov_deg=0.1, pix_arcsec=6.15)

    assert n == 59
    assert grid.wcs.radesys == "ICRS"
    assert list(grid.wcs.ctype) == ["RA---TAN", "DEC--TAN"]
    cx, cy = (n - 1) / 2.0, (n - 1) / 2.0
    sky = grid.pixel_to_world(cx, cy)
    assert sky.separation(coord).arcsec < 1e-6


def test_compute_reference_pixel_scale_uses_median_across_frames():
    frames = [
        {
            "header": _make_wcs_header(269.44, 4.76, pix_arcsec=6.50)
        },  # outlier first frame
        {"header": _make_wcs_header(269.44, 4.76, pix_arcsec=6.15)},
        {"header": _make_wcs_header(269.44, 4.76, pix_arcsec=6.15)},
    ]
    scale = compute_reference_pixel_scale(frames)
    assert np.isclose(scale, 6.15, rtol=1e-4)


def test_align_frame_to_grid_preserves_star_position_under_rotation():
    coord = SkyCoord(ra=269.446294 * u.deg, dec=4.767881 * u.deg, frame="icrs")
    grid, n = make_grid(coord, fov_deg=0.05, pix_arcsec=6.15)

    # Create a frame rotated by 35 degrees with a bright source at the target coord
    hdr = _make_wcs_header(
        coord.ra.deg, coord.dec.deg, nx=31, ny=31, pix_arcsec=6.15, angle_deg=35.0
    )
    image = np.ones((31, 31), dtype=np.float32) * 0.5
    image[15, 15] = 100.0  # center pixel corresponds to crval (coord)

    frame = {"image": image, "header": hdr, "validation": {}}
    align_frame_to_grid(frame, coord, grid, shape_out=(n, n))

    assert frame["aligned"].shape == (n, n)
    assert frame["valid_mask"].shape == (n, n)
    assert frame["validation"]["alignment"]["target_inside_native_footprint"] is True
    assert frame["validation"]["alignment"]["target_covered_after_reprojection"] is True

    # Peak after reprojection must sit at the center of the output grid
    peak_y, peak_x = np.unravel_index(
        np.nanargmax(frame["aligned"]), frame["aligned"].shape
    )
    center_idx = (n - 1) / 2.0
    assert np.hypot(peak_x - center_idx, peak_y - center_idx) <= 1.0


def test_partial_footprint_edge_cutout_aligns_and_preserves_nans():
    coord = SkyCoord(ra=269.446294 * u.deg, dec=4.767881 * u.deg, frame="icrs")
    grid, n = make_grid(coord, fov_deg=0.1, pix_arcsec=6.15)

    # Offset frame so target is outside the small (10, 58) strip, but strip overlaps grid edge
    offset_ra = coord.ra.deg
    offset_dec = coord.dec.deg - 0.035
    hdr = _make_wcs_header(offset_ra, offset_dec, nx=58, ny=10, pix_arcsec=6.15)
    image = np.full((10, 58), 5.0, dtype=np.float32)

    frame = {"image": image, "header": hdr, "validation": {}}
    align_frame_to_grid(frame, coord, grid, shape_out=(n, n))

    alignment_meta = frame["validation"]["alignment"]
    assert alignment_meta["target_inside_native_footprint"] is False
    assert alignment_meta["partial_coverage"] is True
    assert 0.05 < alignment_meta["reprojection_covered_fraction"] < 0.50

    # Uncovered pixels must be NaN, never 0.0
    uncovered = ~frame["valid_mask"]
    assert np.any(uncovered)
    assert np.all(np.isnan(frame["aligned"][uncovered]))


def test_evaluate_cross_frame_registration_empirically_validates_and_handles_empty_sky():
    rng = np.random.default_rng(99)
    img1 = rng.normal(0.0, 0.05, size=(25, 25))
    img2 = rng.normal(0.0, 0.05, size=(25, 25))
    img1[12, 12] += 15.0
    img2[12, 12] += 15.0

    frames_with_star = [
        {"detector": 5, "time": 60900.0, "aligned": img1, "normalized": img1},
        {"detector": 5, "time": 60901.0, "aligned": img2, "normalized": img2},
    ]
    reg_ok = evaluate_cross_frame_registration(frames_with_star, pix_arcsec=6.15)
    assert reg_ok["status"] == "empirically_validated"
    assert reg_ok["independent_astrometric_validation"] is True
    assert reg_ok["intra_epoch_centroid_rms_px"] < 0.25

    # Pure noise frames with no reference source must explicitly report incomplete_no_reference_sources
    noise1 = rng.normal(0.0, 0.05, size=(25, 25))
    noise2 = rng.normal(0.0, 0.05, size=(25, 25))
    frames_no_star = [
        {"detector": 5, "time": 60900.0, "aligned": noise1, "normalized": noise1},
        {"detector": 5, "time": 60901.0, "aligned": noise2, "normalized": noise2},
    ]
    reg_incomplete = evaluate_cross_frame_registration(
        frames_no_star, pix_arcsec=6.15, snr_threshold=10.0
    )
    assert reg_incomplete["status"] == "incomplete_no_reference_sources"
    assert reg_incomplete["independent_astrometric_validation"] is False
