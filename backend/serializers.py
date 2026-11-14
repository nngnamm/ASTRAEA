import base64
import io

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from analytics import compute_sky_health, detect_anomalies
from config import MIN_COMMON_COVERAGE, MIN_EPOCH_COVERAGE
from spherex_core import (
    coadd,
    compute_difference_with_diagnostics,
    difference,
    group_epochs_by_detector,
)


def array_to_b64png(arr, vmin=None, vmax=None, cmap_name="gray"):
    """Encode a 2-D array as a PNG data URI with transparent RGBA for NaN pixels.

    Uncovered or invalid (non-finite) pixels are rendered with RGBA ``[0, 0, 0, 0]``
    (alpha = 0) so they are never mistaken for zero-flux black sky or negative
    difference residuals.
    """
    image = np.asarray(arr, dtype=float)
    if image.ndim != 2:
        raise ValueError("PNG image input must be a 2-D array")

    finite = image[np.isfinite(image)]
    if finite.size == 0:
        raise ValueError("Cannot render an image with no finite pixels")

    if vmin is None:
        vmin = float(np.percentile(finite, 5))
    if vmax is None:
        vmax = float(np.percentile(finite, 99))

    if not np.isfinite(vmin) or not np.isfinite(vmax):
        raise ValueError("Image display limits must be finite")

    if vmax <= vmin:
        center = float(np.median(finite))
        scale = max(abs(center) * 1e-6, 1e-12)
        vmin, vmax = center - scale, center + scale

    flipped = np.flipud(image)
    invalid_mask = ~np.isfinite(flipped)
    with np.errstate(all="ignore"):
        norm = np.clip((flipped - vmin) / (vmax - vmin), 0.0, 1.0)
    norm_safe = np.nan_to_num(norm, nan=0.0, posinf=1.0, neginf=0.0)
    rgba = (matplotlib.colormaps[cmap_name](norm_safe) * 255).astype(np.uint8)
    if np.any(invalid_mask):
        rgba[invalid_mask] = np.array([0, 0, 0, 0], dtype=np.uint8)

    buffer = io.BytesIO()
    plt.imsave(buffer, rgba, format="png")
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return "data:image/png;base64," + encoded


def compute_detector_display_limits(frames, low_pct=5.0, high_pct=99.0):
    """Compute shared (vmin, vmax) display limits per detector from normalized frames.

    Using a detector-wide display stretch after background subtraction prevents
    partial-coverage frames (which may miss the bright central target) from
    stretching background noise to white during animation or blinking.
    """
    limits = {}
    detector_ids = sorted({int(frame["detector"]) for frame in frames})
    for det in detector_ids:
        det_frames = [f for f in frames if int(f["detector"]) == det]
        # Prefer full-coverage frames when estimating the 99th percentile peak
        full_frames = [
            f
            for f in det_frames
            if not f.get("validation", {})
            .get("alignment", {})
            .get("partial_coverage", False)
        ] or det_frames

        samples = []
        for f in full_frames:
            arr = np.asarray(f.get("normalized", f["aligned"]), dtype=float)
            finite = arr[np.isfinite(arr)]
            if finite.size > 0:
                samples.append(finite)

        if not samples:
            limits[det] = (-0.1, 1.0)
            continue

        pooled = np.concatenate(samples)
        vmin = float(np.percentile(pooled, low_pct))
        vmax = float(np.percentile(pooled, high_pct))
        if vmax <= vmin:
            center = float(np.median(pooled))
            scale = max(abs(center) * 1e-6, 1e-12)
            vmin, vmax = center - scale, center + scale
        limits[det] = (round(vmin, 6), round(vmax, 6))
    return limits


def _epoch_coverage_fraction(epoch):
    """Compute the fraction of finite pixels in the coadded epoch image."""
    if not epoch:
        return 0.0
    image = np.asarray(coadd(epoch), dtype=float)
    return round(float(np.isfinite(image).mean()), 6) if image.size > 0 else 0.0


def _epoch_summary(epoch, detector, epoch_index, display_limits=None):
    times = [float(frame["time"]) for frame in epoch]
    mean_mjd = float(np.mean(times))
    image = np.asarray(coadd(epoch), dtype=float)
    valid_pixel_fraction = (
        round(float(np.isfinite(image).mean()), 6) if image.size > 0 else 0.0
    )
    vmin, vmax = (
        display_limits.get(int(detector), (None, None))
        if display_limits
        else (None, None)
    )
    input_obs_ids = [
        f.get("metadata", {}).get("observation", {}).get("observation_id")
        for f in epoch
        if f.get("metadata", {}).get("observation", {}).get("observation_id")
    ]
    return {
        "epoch_id": epoch_index,
        "dataset_identity": "SPHEREx",
        "is_derived_image": True,
        "data_lineage": {
            "dataset_identity": "SPHEREx",
            "is_derived_product": True,
            "derivation_method": "median_coadd_of_tan_reprojected_zodi_sky_subtracted_frames",
            "input_observation_ids": input_obs_ids,
            "input_frame_count": len(epoch),
        },
        "detector": int(detector),
        "mjd_mean": round(mean_mjd, 4),
        "mjd_start": round(min(times), 4),
        "mjd_end": round(max(times), 4),
        "iso_date": Time(mean_mjd, format="mjd").isot[:10],
        "n_frames": len(epoch),
        "single_frame_epoch": len(epoch) == 1,
        "grouping_heuristic_gap_days": 30,
        "image_shape_yx": [int(image.shape[0]), int(image.shape[1])],
        "valid_pixel_fraction": valid_pixel_fraction,
        "display_limits_mjy_sr": [vmin, vmax] if vmin is not None else None,
        "invalid_pixel_rendering": "transparent_rgba_alpha_0",
        "png_b64": array_to_b64png(image, vmin=vmin, vmax=vmax, cmap_name="gray"),
    }


def _difference_for_detector(
    detector,
    detector_epochs,
    min_epoch_coverage=MIN_EPOCH_COVERAGE,
    min_common_coverage=MIN_COMMON_COVERAGE,
    cube_grid=None,
):
    """Build a coverage-aware first-versus-last epoch comparison for one detector.

    Enforces scientific comparison eligibility strictly: never silently lowers
    ``min_epoch_coverage`` by falling back to low-coverage partial epochs.
    """
    if len(detector_epochs) < 2:
        return {
            "detector": int(detector),
            "status": "insufficient_epochs",
            "comparison_eligible": False,
            "message": "At least two time-separated epochs are required.",
            "n_epochs": len(detector_epochs),
            "n_eligible_epochs": len(detector_epochs),
            "min_epoch_coverage_required": float(min_epoch_coverage),
            "png_b64": None,
        }

    epoch_coverages = [_epoch_coverage_fraction(ep) for ep in detector_epochs]
    well_covered_epochs = [
        epoch
        for epoch, cov in zip(detector_epochs, epoch_coverages)
        if cov >= float(min_epoch_coverage)
    ]

    if len(well_covered_epochs) < 2:
        return {
            "detector": int(detector),
            "status": "ineligible_low_coverage",
            "comparison_eligible": False,
            "message": (
                f"Fewer than two epochs meet the minimum spatial coverage threshold "
                f"({min_epoch_coverage:.0%}); found {len(well_covered_epochs)} eligible "
                f"of {len(detector_epochs)} epochs."
            ),
            "n_epochs": len(detector_epochs),
            "n_eligible_epochs": len(well_covered_epochs),
            "min_epoch_coverage_required": float(min_epoch_coverage),
            "epoch_coverage_fractions": epoch_coverages,
            "png_b64": None,
        }

    early_epoch = well_covered_epochs[0]
    late_epoch = well_covered_epochs[-1]
    early_mjd = float(np.mean([frame["time"] for frame in early_epoch]))
    late_mjd = float(np.mean([frame["time"] for frame in late_epoch]))
    comparison_mid_mjd = 0.5 * (early_mjd + late_mjd)

    try:
        diff, img_a_sub, img_b_sub, diagnostics = compute_difference_with_diagnostics(
            early_epoch,
            late_epoch,
            detector=int(detector),
            min_common_coverage=float(min_common_coverage),
        )
        finite = np.abs(diff[np.isfinite(diff)])
        if finite.size == 0:
            raise ValueError("Difference image contains no finite pixels")
        vlim = float(np.percentile(finite, 99))
        if not np.isfinite(vlim) or vlim <= 0:
            vlim = 1e-12
        png = array_to_b64png(
            diff,
            vmin=-vlim,
            vmax=vlim,
            cmap_name="RdBu_r",
        )
    except (ValueError, KeyError, TypeError) as exc:
        return {
            "detector": int(detector),
            "status": "processing_error",
            "comparison_eligible": False,
            "message": str(exc),
            "n_epochs": len(detector_epochs),
            "n_eligible_epochs": len(well_covered_epochs),
            "min_epoch_coverage_required": float(min_epoch_coverage),
            "min_common_coverage_required": float(min_common_coverage),
            "png_b64": None,
        }

    candidates = []
    if cube_grid is not None:
        try:
            candidates = detect_anomalies(
                diff,
                img_a_sub,
                img_b_sub,
                cube_grid,
                reference_mjd=comparison_mid_mjd,
            )
        except Exception:
            candidates = []

    early_count = len(early_epoch)
    late_count = len(late_epoch)
    warnings = []
    if early_count == 1:
        warnings.append("The early epoch contains only one observation.")
    if late_count == 1:
        warnings.append(
            "The late epoch contains only one observation; interpret the comparison cautiously."
        )

    common_valid_pixel_fraction = diagnostics["common_valid_pixel_fraction"]
    target_covered_in_difference = diagnostics["target_covered_in_difference"]
    difference_noise_mad_mjy_sr = diagnostics["difference_noise_mad_mjy_sr"]
    mean_wavelength_shift_um = diagnostics.get("mean_wavelength_shift_um") or 0.0

    if common_valid_pixel_fraction < 0.50:
        warnings.append(
            f"Common valid pixel coverage across compared epochs is low ({common_valid_pixel_fraction:.1%} < 50.0%)."
        )
    if not target_covered_in_difference:
        warnings.append(
            "Target cutout center is not covered in the joint difference mask."
        )
    if mean_wavelength_shift_um > 0.15:
        warnings.append(
            f"Mean wavelength shift between epochs is {mean_wavelength_shift_um:.3f} um (> 0.150 um); "
            "LVF bandpass shift may introduce chromatic residuals."
        )

    early_obs_ids = [
        f.get("metadata", {}).get("observation", {}).get("observation_id")
        for f in early_epoch
        if f.get("metadata", {}).get("observation", {}).get("observation_id")
    ]
    late_obs_ids = [
        f.get("metadata", {}).get("observation", {}).get("observation_id")
        for f in late_epoch
        if f.get("metadata", {}).get("observation", {}).get("observation_id")
    ]

    return {
        "detector": int(detector),
        "dataset_identity": "SPHEREx",
        "is_derived_image": True,
        "data_lineage": {
            "dataset_identity": "SPHEREx",
            "is_derived_product": True,
            "derivation_method": "late_epoch_median_coadd_minus_early_epoch_median_coadd",
            "early_input_observation_ids": early_obs_ids,
            "late_input_observation_ids": late_obs_ids,
        },
        "status": "available",
        "comparison_eligible": True,
        "min_epoch_coverage_required": float(min_epoch_coverage),
        "min_common_coverage_required": float(min_common_coverage),
        "mjd_early": round(early_mjd, 4),
        "mjd_late": round(late_mjd, 4),
        "baseline_days": round(late_mjd - early_mjd, 2),
        "n_epochs": len(detector_epochs),
        "n_eligible_epochs": len(well_covered_epochs),
        "n_frames_early": early_count,
        "n_frames_late": late_count,
        "image_shape_yx": [int(diff.shape[0]), int(diff.shape[1])],
        "common_valid_pixel_fraction": common_valid_pixel_fraction,
        "target_covered_in_difference": target_covered_in_difference,
        "difference_noise_mad_mjy_sr": difference_noise_mad_mjy_sr,
        "diagnostics": diagnostics,
        "candidates": candidates,
        "comparison_quality": "single_frame_epoch"
        if (early_count == 1 or late_count == 1)
        else "multi_frame_epochs",
        "warnings": warnings,
        "vlim_mjy_sr": round(vlim, 6),
        "invalid_pixel_rendering": "transparent_rgba_alpha_0",
        "png_b64": png,
    }


def package_cube_json(cube):
    """Serialize a cube while keeping timelines and differences detector-safe.

    The legacy top-level ``difference`` fields are retained for the frontend:
    they describe the first detector with an available comparison. The full,
    scientifically separated set is returned in ``differences``.
    """
    frames = cube.get("frames", [])
    if not frames:
        cube_status = cube.get("status", "empty")
        return {
            "status": cube_status if cube_status in {"empty", "error"} else "empty",
            "error_code": cube.get("error_code", "no_observations"),
            "error": cube.get("message", "No frames found"),
            "message": cube.get("message", "No frames found"),
            "frames": [],
            "epochs": [],
            "differences": [],
            "validation_summary": cube.get("validation_summary", {}),
        }

    ordered_frames = sorted(frames, key=lambda frame: float(frame["time"]))
    display_limits = compute_detector_display_limits(ordered_frames)

    frames_out = []
    for index, frame in enumerate(ordered_frames):
        mjd = float(frame["time"])
        det = int(frame["detector"])
        vmin, vmax = display_limits.get(det, (None, None))
        render_arr = np.asarray(frame.get("normalized", frame["aligned"]), dtype=float)
        cov_frac = (
            round(float(np.isfinite(render_arr).mean()), 6)
            if render_arr.size > 0
            else 0.0
        )
        obs_id = (
            frame.get("metadata", {}).get("observation", {}).get("observation_id")
        )
        frames_out.append(
            {
                "index": index,
                "dataset_identity": "SPHEREx",
                "is_derived_image": True,
                "data_lineage": {
                    "dataset_identity": "SPHEREx",
                    "original_observation_id": obs_id,
                    "original_fits_is_derived_product": False,
                    "rendered_image_is_derived_product": True,
                    "derivation_method": "icrs_tan_wcs_reprojection_plus_zodi_and_sky_pedestal_subtraction",
                },
                "mjd": round(mjd, 4),
                "iso_date": Time(mjd, format="mjd").isot[:10],
                "detector": det,
                "wave_um": round(float(frame.get("wave_um", 0.0)), 3),
                "collection": frame.get("collection"),
                "image_shape_yx": [int(render_arr.shape[0]), int(render_arr.shape[1])],
                "coverage_fraction": cov_frac,
                "display_limits_mjy_sr": [vmin, vmax] if vmin is not None else None,
                "invalid_pixel_rendering": "transparent_rgba_alpha_0",
                "metadata": frame.get("metadata", {}),
                "validation": frame.get("validation", {}),
                "png_b64": array_to_b64png(
                    render_arr, vmin=vmin, vmax=vmax, cmap_name="gray"
                ),
            }
        )

    grouped = group_epochs_by_detector(ordered_frames)
    cube_grid = cube.get("grid")
    epochs_out = []
    differences_out = []
    for detector, detector_epochs in grouped.items():
        for local_index, epoch in enumerate(detector_epochs):
            epochs_out.append(
                _epoch_summary(
                    epoch, detector, local_index, display_limits=display_limits
                )
            )
        differences_out.append(
            _difference_for_detector(detector, detector_epochs, cube_grid=cube_grid)
        )

    epochs_out.sort(key=lambda epoch: epoch["mjd_mean"])
    for global_index, epoch in enumerate(epochs_out):
        epoch["epoch_id"] = global_index

    available_differences = [
        item for item in differences_out if item.get("status") == "available"
    ]
    if available_differences:
        primary = dict(available_differences[0])
        difference_out = {**primary, "by_detector": differences_out}
    else:
        ineligible = [
            item
            for item in differences_out
            if item.get("status") == "ineligible_low_coverage"
        ]
        difference_out = {
            "status": "ineligible_low_coverage"
            if ineligible
            else "insufficient_observations",
            "comparison_eligible": False,
            "message": (
                ineligible[0]["message"]
                if ineligible
                else "No detector has two time-separated epochs."
            ),
            "png_b64": None,
            "by_detector": differences_out,
        }

    sky_health = compute_sky_health(cube)
    candidates_by_detector = {
        str(d["detector"]): d.get("candidates", [])
        for d in differences_out
        if d.get("status") == "available"
    }

    center = cube.get("center", (None, None))
    n_grid = cube.get("n")
    pix_scale = (
        round(float(cube["pix_arcsec"]), 3)
        if cube.get("pix_arcsec") is not None
        else None
    )
    return {
        "status": "ok",
        "dataset_identity": "SPHEREx",
        "archive_source": "NASA/IPAC IRSA",
        "official_archive_url": "https://irsa.ipac.caltech.edu/Missions/spherex.html",
        "center": {"ra": center[0], "dec": center[1]},
        "grid_size_px": n_grid,
        "pix_scale_arcsec": pix_scale,
        "wcs_summary": {
            "reference_frame": "ICRS",
            "projection": "TAN",
            "crval_ra_deg": center[0],
            "crval_dec_deg": center[1],
            "grid_shape_yx": [n_grid, n_grid] if n_grid is not None else None,
            "pix_scale_arcsec": pix_scale,
            "orientation": "North-up, East-left",
        },
        "processing_history": [
            "irsa_sia_query",
            "fits_structural_and_quality_validation",
            "icrs_tan_wcs_reprojection",
            "zodi_and_flag_aware_sky_pedestal_subtraction",
            "detector_separated_epoch_coadd",
            "joint_mask_difference_and_anomaly_detection",
        ],
        "display_limits_by_detector_mjy_sr": {
            str(k): list(v) for k, v in display_limits.items()
        },
        "epoch_grouping_policy": {
            "method": "time_gap_heuristic",
            "gap_days": 30,
            "official_mission_epoch": False,
            "note": (
                "Observations are grouped into visit windows separated by >30 days; "
                "these are heuristic visit clusters, not official NASA mission epoch boundaries."
            ),
        },
        "n_frames_total": len(frames_out),
        "n_epochs": len(epochs_out),
        "frames": frames_out,
        "epochs": epochs_out,
        "difference": difference_out,
        "differences": differences_out,
        "sky_health": sky_health,
        "experimental_analytics": {
            "experimental_heuristic": True,
            "sky_health": sky_health,
            "candidates_by_detector": candidates_by_detector,
        },
        "validation_summary": cube.get("validation_summary", {}),
    }
