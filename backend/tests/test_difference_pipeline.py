import numpy as np
import pytest

from serializers import _difference_for_detector
from spherex_core import (
    coadd,
    compute_difference_with_diagnostics,
    difference,
    group_epochs_by_detector,
)


def make_frame(
    detector: int,
    mjd: float,
    image: np.ndarray,
    obsid: str | None = None,
    variance: np.ndarray | None = None,
) -> dict:
    frame = {
        "detector": detector,
        "time": mjd,
        "aligned": np.asarray(image, dtype=float),
    }
    if variance is not None:
        frame["aligned_variance"] = np.asarray(variance, dtype=float)
    if obsid is not None:
        frame["metadata"] = {"observation": {"observation_id": obsid}}
        frame["header"] = {"OBSID": obsid}
    return frame


def test_epochs_are_grouped_independently_per_detector():
    image = np.array([[0.0, 1.0], [2.0, 3.0]])
    frames = [
        make_frame(5, 100.0, image),
        make_frame(6, 101.0, image),
        make_frame(5, 101.0, image),
        make_frame(6, 200.0, image),
        make_frame(5, 200.0, image),
    ]

    grouped = group_epochs_by_detector(frames, gap_days=30)

    assert set(grouped) == {5, 6}
    assert [len(epoch) for epoch in grouped[5]] == [2, 1]
    assert [len(epoch) for epoch in grouped[6]] == [1, 1]
    assert all(
        all(int(frame["detector"]) == detector for frame in epoch)
        for detector, epochs in grouped.items()
        for epoch in epochs
    )


def test_single_frame_epoch_can_be_coadded():
    image = np.array([[1.0, 2.0], [3.0, 4.0]])
    result = coadd([make_frame(5, 100.0, image)])
    np.testing.assert_array_equal(result, image)


def test_difference_filters_to_requested_detector():
    pattern = np.array([[0.0, 0.0], [0.0, 2.0]])
    early = [
        make_frame(5, 100.0, pattern),
        make_frame(6, 100.0, pattern * 100.0),
    ]
    late = [
        make_frame(5, 200.0, pattern * 2.0),
        make_frame(6, 200.0, pattern * 500.0),
    ]

    result, _, _ = difference(early, late, detector=5)
    expected, _, _ = difference(
        [early[0]], [late[0]], detector=5
    )
    np.testing.assert_allclose(result, expected)


def test_empty_epochs_are_rejected():
    try:
        difference([], [], detector=5)
    except ValueError as exc:
        assert "Both comparison epochs" in str(exc)
    else:
        raise AssertionError("Expected an empty-epoch ValueError")


def test_difference_preserves_nan_outside_common_mask_and_records_diagnostics():
    img_a = np.full((4, 4), np.nan, dtype=float)
    img_b = np.full((4, 4), np.nan, dtype=float)
    # Rows 0..2 valid in A (12 pixels), rows 1..3 valid in B (12 pixels) -> rows 1..2 (8 pixels) in common_mask
    img_a[0:3, :] = 10.0
    img_b[1:4, :] = 15.0
    img_b[2, 2] = 25.0  # +10.0 source increment at center (2, 2)

    early = [
        make_frame(5, 100.0, img_a, obsid="OBS-EARLY-1"),
        make_frame(5, 100.5, img_a, obsid="OBS-EARLY-2"),
    ]
    late = [
        make_frame(5, 250.0, img_b, obsid="OBS-LATE-1"),
    ]

    diff, img_a_sub, img_b_sub, diag = compute_difference_with_diagnostics(
        early, late, detector=5, min_common_coverage=0.10
    )

    common_mask = np.isfinite(img_a) & np.isfinite(img_b)
    assert np.isnan(diff[~common_mask]).all()
    assert np.isnan(img_a_sub[~common_mask]).all()
    assert np.isnan(img_b_sub[~common_mask]).all()
    assert np.isfinite(diff[common_mask]).all()
    assert np.isclose(diff[2, 2], 10.0)

    assert diag["detector"] == 5
    assert diag["common_valid_pixels"] == 8
    assert np.isclose(diag["common_valid_pixel_fraction"], 8 / 16)
    assert np.isclose(diag["early_coverage_fraction"], 12 / 16)
    assert np.isclose(diag["late_coverage_fraction"], 12 / 16)
    assert diag["target_covered_in_difference"] is True
    assert diag["early_observation_ids"] == ["OBS-EARLY-1", "OBS-EARLY-2"]
    assert diag["late_observation_ids"] == ["OBS-LATE-1"]
    assert diag["uncovered_pixel_policy"] == "preserved_as_nan"
    assert diag["variance_propagated"] is False


def test_difference_propagates_variance_when_available():
    img_a = np.ones((6, 6), dtype=float) * 5.0
    img_b = np.ones((6, 6), dtype=float) * 5.0
    img_b[3, 3] = 12.0
    var_a = np.full((6, 6), 0.09, dtype=float)  # sigma_a = 0.3
    var_b = np.full((6, 6), 0.16, dtype=float)  # sigma_b = 0.4 -> quadrature diff sigma = 0.5

    _, _, _, diag = compute_difference_with_diagnostics(
        [make_frame(5, 100.0, img_a, variance=var_a)],
        [make_frame(5, 200.0, img_b, variance=var_b)],
        detector=5,
    )
    assert diag["variance_propagated"] is True
    assert np.isclose(diag["median_propagated_sigma_mjy_sr"], 0.5, atol=1e-5)


def test_difference_for_detector_does_not_silently_fallback_to_low_coverage_epochs():
    # Epoch 1 has 100% coverage; Epoch 2 has only 20% coverage (< 50% threshold)
    full_img = np.ones((10, 10), dtype=float)
    low_cov_img = np.full((10, 10), np.nan, dtype=float)
    low_cov_img[:2, :] = 1.0  # 20/100 = 0.20 coverage

    epochs = [
        [make_frame(5, 100.0, full_img)],
        [make_frame(5, 200.0, low_cov_img)],
    ]
    res = _difference_for_detector(5, epochs, min_epoch_coverage=0.50)
    assert res["status"] == "ineligible_low_coverage"
    assert res["comparison_eligible"] is False
    assert res["n_epochs"] == 2
    assert res["n_eligible_epochs"] == 1
    assert res["png_b64"] is None


def test_background_pedestal_computed_strictly_over_common_mask():
    # Rows 0..1 are covered only in epoch A and contain a bright edge artifact (+5000 MJy/sr).
    # Rows 2..5 are covered in both epochs with steady sky (10.0 in A, 20.0 in B).
    img_a = np.full((6, 6), 10.0, dtype=float)
    img_a[0:2, :] = 5000.0

    img_b = np.full((6, 6), np.nan, dtype=float)
    img_b[2:6, :] = 20.0

    diff, _, _, diag = compute_difference_with_diagnostics(
        [make_frame(5, 100.0, img_a)],
        [make_frame(5, 200.0, img_b)],
        detector=5,
    )

    common_mask = np.isfinite(img_a) & np.isfinite(img_b)
    assert np.isclose(diag["early_sky_median_removed_mjy_sr"], 10.0)
    assert np.isclose(diag["late_sky_median_removed_mjy_sr"], 20.0)
    np.testing.assert_allclose(diff[common_mask], 0.0, atol=1e-12)
    assert np.isnan(diff[~common_mask]).all()


def test_disjoint_or_low_overlap_footprints_raise_value_error():
    img_a = np.full((10, 10), np.nan, dtype=float)
    img_b = np.full((10, 10), np.nan, dtype=float)
    img_a[0:4, :] = 1.0
    img_b[6:10, :] = 1.0  # 0% overlap

    with pytest.raises(ValueError, match="Insufficient common valid pixel coverage"):
        compute_difference_with_diagnostics(
            [make_frame(5, 100.0, img_a)],
            [make_frame(5, 200.0, img_b)],
            detector=5,
            min_common_coverage=0.10,
        )


def test_mixing_detectors_when_detector_is_none_raises_value_error():
    image = np.ones((4, 4), dtype=float)
    with pytest.raises(ValueError, match="Cannot mix multiple detectors"):
        compute_difference_with_diagnostics(
            [make_frame(5, 100.0, image)],
            [make_frame(6, 200.0, image)],
            detector=None,
        )
