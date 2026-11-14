import io
import base64

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time

from spherex_core import group_epochs_by_detector, coadd, difference


def array_to_b64png(arr, vmin=None, vmax=None, cmap_name="gray"):
    """Encode a 2-D array as a PNG data URI, safely handling NaN pixels."""
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
        # Avoid a zero-width display stretch for constant or nearly constant images.
        center = float(np.median(finite))
        scale = max(abs(center) * 1e-6, 1e-12)
        vmin, vmax = center - scale, center + scale

    flipped = np.flipud(image)
    norm = np.clip((flipped - vmin) / (vmax - vmin), 0.0, 1.0)
    norm = np.nan_to_num(norm, nan=0.0, posinf=1.0, neginf=0.0)
    rgba = (matplotlib.colormaps[cmap_name](norm) * 255).astype(np.uint8)

    buffer = io.BytesIO()
    plt.imsave(buffer, rgba, format="png")
    encoded = base64.b64encode(buffer.getvalue()).decode("ascii")
    return "data:image/png;base64," + encoded


def _epoch_summary(epoch, detector, epoch_index):
    times = [float(frame["time"]) for frame in epoch]
    mean_mjd = float(np.mean(times))
    image = coadd(epoch)
    return {
        "epoch_id": epoch_index,
        "detector": int(detector),
        "mjd_mean": round(mean_mjd, 4),
        "iso_date": Time(mean_mjd, format="mjd").isot[:10],
        "n_frames": len(epoch),
        "single_frame_epoch": len(epoch) == 1,
        "png_b64": array_to_b64png(image, cmap_name="gray"),
    }


def _difference_for_detector(detector, detector_epochs):
    """Build a first-versus-last epoch comparison for one detector."""
    if len(detector_epochs) < 2:
        return {
            "detector": int(detector),
            "status": "insufficient_epochs",
            "message": "At least two time-separated epochs are required.",
            "n_epochs": len(detector_epochs),
            "png_b64": None,
        }

    early_epoch = detector_epochs[0]
    late_epoch = detector_epochs[-1]
    early_mjd = float(np.mean([frame["time"] for frame in early_epoch]))
    late_mjd = float(np.mean([frame["time"] for frame in late_epoch]))

    try:
        diff, _, _ = difference(early_epoch, late_epoch, detector=int(detector))
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
            "message": str(exc),
            "n_epochs": len(detector_epochs),
            "png_b64": None,
        }

    early_count = len(early_epoch)
    late_count = len(late_epoch)
    warnings = []
    if early_count == 1:
        warnings.append("The early epoch contains only one observation.")
    if late_count == 1:
        warnings.append("The late epoch contains only one observation; interpret the comparison cautiously.")

    return {
        "detector": int(detector),
        "status": "available",
        "mjd_early": round(early_mjd, 4),
        "mjd_late": round(late_mjd, 4),
        "baseline_days": round(late_mjd - early_mjd, 2),
        "n_epochs": len(detector_epochs),
        "n_frames_early": early_count,
        "n_frames_late": late_count,
        "comparison_quality": "single_frame_epoch" if warnings else "multi_frame_epochs",
        "warnings": warnings,
        "vlim_mjy_sr": round(vlim, 6),
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
        return {"error": cube.get("message", "No frames found")}

    ordered_frames = sorted(frames, key=lambda frame: float(frame["time"]))

    frames_out = []
    for index, frame in enumerate(ordered_frames):
        mjd = float(frame["time"])
        frames_out.append({
            "index": index,
            "mjd": round(mjd, 4),
            "iso_date": Time(mjd, format="mjd").isot[:10],
            "detector": int(frame["detector"]),
            "wave_um": round(float(frame.get("wave_um", 0.0)), 3),
            "collection": frame.get("collection"),
            "png_b64": array_to_b64png(frame["aligned"], cmap_name="gray"),
        })

    grouped = group_epochs_by_detector(ordered_frames)
    epochs_out = []
    differences_out = []
    for detector, detector_epochs in grouped.items():
        for local_index, epoch in enumerate(detector_epochs):
            epochs_out.append(_epoch_summary(epoch, detector, local_index))
        differences_out.append(_difference_for_detector(detector, detector_epochs))

    # Put all detector-specific epochs on one chronological timeline.
    epochs_out.sort(key=lambda epoch: epoch["mjd_mean"])
    for global_index, epoch in enumerate(epochs_out):
        epoch["epoch_id"] = global_index

    available_differences = [
        item for item in differences_out if item.get("status") == "available"
    ]
    if available_differences:
        # Keep the existing singular response fields for basic frontend
        # compatibility, and identify which detector they represent.
        primary = dict(available_differences[0])
        difference_out = {**primary, "by_detector": differences_out}
    else:
        difference_out = {
            "status": "insufficient_observations",
            "message": "No detector has two time-separated epochs.",
            "png_b64": None,
            "by_detector": differences_out,
        }

    center = cube.get("center", (None, None))
    return {
        "center": {"ra": center[0], "dec": center[1]},
        "grid_size_px": cube.get("n"),
        "pix_scale_arcsec": round(float(cube["pix_arcsec"]), 3)
        if cube.get("pix_arcsec") is not None else None,
        "n_frames_total": len(frames_out),
        "n_epochs": len(epochs_out),
        "frames": frames_out,
        "epochs": epochs_out,
        "difference": difference_out,
        "differences": differences_out,
    }
