import numpy as np
from scipy.ndimage import binary_erosion, gaussian_filter, maximum_filter
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.ipac.irsa import Irsa
from config import TARGETS
from spherex_core import propagate_coords

_ASTROPHYSICAL_FLAG_NAMES = {"SOURCE", "TRANSIENT"}


def compute_sky_health(cube):
    """Transparent, formula-defined Sky Health Index (Blueprint Group II #3)."""
    frames = cube.get("frames", [])
    if not frames:
        return None

    aligned_stack = np.array([np.asarray(f["aligned"], dtype=float) for f in frames])
    valid_coverage = float(np.isfinite(aligned_stack).mean() * 100.0)

    per_frame_clean_pct = []
    for f in frames:
        defs = (f.get("validation") or {}).get("flag_definitions") or []
        if defs:
            defect_fraction = sum(
                float(d.get("flagged_fraction") or 0.0)
                for d in defs
                if str(d.get("name", "")).strip().upper() not in _ASTROPHYSICAL_FLAG_NAMES
            )
            defect_fraction = min(1.0, max(0.0, defect_fraction))
            per_frame_clean_pct.append((1.0 - defect_fraction) * 100.0)

    if per_frame_clean_pct:
        clean_detector_pct = float(np.mean(per_frame_clean_pct))
    else:
        flag_arrays = [
            np.asarray(f["flags"]).ravel()
            for f in frames
            if f.get("flags") is not None
        ]
        if flag_arrays:
            raw_flags = np.concatenate(flag_arrays)
            clean_detector_pct = float(((raw_flags & 31) == 0).mean() * 100.0)
        else:
            clean_detector_pct = 100.0

    mjd_span = float(frames[-1]["time"] - frames[0]["time"]) if len(frames) > 1 else 0.0
    time_baseline_score = min(100.0, (mjd_span / 180.0) * 100.0)

    health_score = (
        0.50 * valid_coverage + 0.30 * clean_detector_pct + 0.20 * time_baseline_score
    )
    status = (
        "Excellent"
        if health_score >= 85
        else ("Good" if health_score >= 70 else "Moderate")
    )

    return {
        "score_pct": round(health_score, 1),
        "status": status,
        "valid_coverage_pct": round(valid_coverage, 1),
        "clean_detector_pct": round(clean_detector_pct, 1),
        "clean_severe_flags_pct": round(clean_detector_pct, 1),
        "time_baseline_days": round(mjd_span, 1),
        "n_frames": len(frames),
        "formula": "0.50*valid_coverage_pct + 0.30*clean_detector_pct + 0.20*min(100, baseline_days/180*100)",
        "experimental_heuristic": True,
        "provenance_note": (
            "Derived heuristic index combining spatial coverage, non-source flag cleanliness, "
            "and temporal baseline; not an official NASA SPHEREx quality metric."
        ),
    }


def get_catalog_overlay(ra, dec, fov_deg, grid_wcs, n_grid, max_sources=50):
    """Query AllWISE catalog on IRSA and project onto our (n_grid x n_grid) pixel grid."""
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    try:
        tbl = Irsa.query_region(
            coordinates=coord, catalog="allwise_p3as_psd", radius=(fov_deg / 2.0) * u.deg
        )
    except Exception as e:
        print("Catalog query failed:", str(e)[:60])
        return []

    if len(tbl) == 0:
        return []

    tbl.sort("w2mpro")
    sources = []
    for row in tbl[:max_sources]:
        s_ra, s_dec = float(row["ra"]), float(row["dec"])
        sc = SkyCoord(s_ra * u.deg, s_dec * u.deg)
        x_px, y_px = grid_wcs.world_to_pixel(sc)
        if 0 <= x_px < n_grid and 0 <= y_px < n_grid:
            sources.append(
                {
                    "id": str(row["designation"]),
                    "catalog": "AllWISE",
                    "ra": round(s_ra, 6),
                    "dec": round(s_dec, 6),
                    "x_px": round(float(x_px), 2),
                    "y_px": round(float(y_px), 2),
                    "w1_mag": round(float(row["w1mpro"]), 2)
                    if np.isfinite(row["w1mpro"])
                    else None,
                    "w2_mag": round(float(row["w2mpro"]), 2)
                    if np.isfinite(row["w2mpro"])
                    else None,
                }
            )
    return sources


def _classify_candidate(cx, cy, shift_px, snr_combined, z_thresh, grid_wcs, cat_coords, catalog_sources):
    sky = grid_wcs.pixel_to_world(cx, cy)

    classification = "Anomaly Candidate"
    matched_id = None
    for t_info in TARGETS.values():
        tr_ra, tr_dec = propagate_coords(
            t_info["ra"],
            t_info["dec"],
            t_info["pm_ra_arcsec_yr"],
            t_info["pm_dec_arcsec_yr"],
            t_info["epoch"],
        )
        sep_arcsec = float(
            sky.separation(SkyCoord(tr_ra * u.deg, tr_dec * u.deg)).arcsec
        )
        if sep_arcsec < 12.0:
            classification = f"Confirmed Benchmark Object ({t_info['name']})"
            matched_id = t_info["name"]
            break

    if matched_id is None and cat_coords is not None and len(cat_coords) > 0:
        seps = sky.separation(cat_coords).arcsec
        min_idx = int(np.argmin(seps))
        if float(seps[min_idx]) < 6.5:
            matched_id = catalog_sources[min_idx]["id"]
            if shift_px == 0.0:
                classification = f"Known Object ({matched_id})"

    conf = round(
        float(min(0.95, max(0.15, (snr_combined - z_thresh) / 20.0 + 0.50))), 2
    )
    status = "High Priority Review" if snr_combined >= 8.0 else "Candidate - Review"
    return sky, classification, matched_id, conf, status


def detect_anomalies(
    diff_img,
    img_early,
    img_late,
    grid_wcs,
    catalog_sources=None,
    z_thresh=4.5,
    border=5,
):
    """
    Unsupervised Difference Blob Detector + Catalog Cross-Matcher + Uncertainty-Aware Ranking.
    """
    diff_arr = np.asarray(diff_img, dtype=float)
    early_arr = np.asarray(img_early, dtype=float)
    late_arr = np.asarray(img_late, dtype=float)

    valid_mask = np.isfinite(diff_arr) & np.isfinite(early_arr) & np.isfinite(late_arr)
    if int(np.count_nonzero(valid_mask)) < 25:
        return []

    eroded_mask = binary_erosion(valid_mask, iterations=2)
    if int(np.count_nonzero(eroded_mask)) < 9:
        eroded_mask = valid_mask

    raw_med = float(np.median(diff_arr[valid_mask]))
    filled_diff = np.where(valid_mask, diff_arr, raw_med)
    sm_diff = gaussian_filter(filled_diff, sigma=0.8)

    # Local sky noise strictly from finite pixels inside valid_mask
    finite_vals = sm_diff[eroded_mask]
    med = float(np.median(finite_vals))
    sigma_sky = float(1.4826 * np.median(np.abs(finite_vals - med))) + 1e-12

    early_clean = np.where(valid_mask, np.abs(early_arr), 0.0)
    late_clean = np.where(valid_mask, np.abs(late_arr), 0.0)
    mean_flux = 0.5 * (early_clean + late_clean)
    sigma_total = np.hypot(sigma_sky, 0.05 * gaussian_filter(mean_flux, sigma=1.0))

    z_map = np.zeros_like(sm_diff, dtype=float)
    z_map[valid_mask] = (sm_diff[valid_mask] - med) / sigma_total[valid_mask]

    ny_grid, nx_grid = z_map.shape
    border_mask = np.zeros_like(z_map, dtype=bool)
    if ny_grid > 2 * border and nx_grid > 2 * border:
        border_mask[border : ny_grid - border, border : nx_grid - border] = True
    else:
        border_mask[:, :] = True

    search_mask = eroded_mask & border_mask
    pos_max = (z_map == maximum_filter(z_map, size=5)) & (z_map >= z_thresh)
    neg_min = (-z_map == maximum_filter(-z_map, size=5)) & (-z_map >= z_thresh)

    pos_yx = np.argwhere(pos_max & search_mask)
    neg_yx = np.argwhere(neg_min & search_mask)

    cat_coords = (
        SkyCoord(
            [s["ra"] for s in catalog_sources] * u.deg,
            [s["dec"] for s in catalog_sources] * u.deg,
        )
        if catalog_sources
        else None
    )

    used_neg = set()
    candidates = []
    score_disclaimer = (
        "Heuristic SNR ranking score; does not represent the probability of a new astrophysical object."
    )

    for py, px in pos_yx:
        z_pos = float(z_map[py, px])
        best_neg_idx, best_dist = None, 999.0
        for idx, (ny, nx) in enumerate(neg_yx):
            if idx in used_neg:
                continue
            dist = float(np.hypot(py - ny, px - nx))
            if 0.8 <= dist <= 5.5 and dist < best_dist:
                best_dist = dist
                best_neg_idx = idx

        if best_neg_idx is not None:
            used_neg.add(best_neg_idx)
            ny, nx = neg_yx[best_neg_idx]
            z_neg = float(abs(z_map[ny, nx]))
            snr_combined = round(float(np.hypot(z_pos, z_neg)), 1)
            cy, cx = (py + ny) / 2.0, (px + nx) / 2.0
            cand_type = "Moving Object Candidate"
            flux_dir = "dipole"
            shift_px = round(best_dist, 2)
        else:
            snr_combined = round(z_pos, 1)
            cy, cx = float(py), float(px)
            cand_type = "Variable Source Candidate (Brightening)"
            flux_dir = "brightening"
            shift_px = 0.0

        sky, classification, matched_id, conf, status = _classify_candidate(
            cx, cy, shift_px, snr_combined, z_thresh, grid_wcs, cat_coords, catalog_sources
        )

        candidates.append(
            {
                "type": cand_type,
                "flux_change_direction": flux_dir,
                "classification": classification,
                "matched_catalog_id": matched_id,
                "x_px": round(float(cx), 2),
                "y_px": round(float(cy), 2),
                "ra": round(float(sky.ra.deg), 6),
                "dec": round(float(sky.dec.deg), 6),
                "snr_zscore": snr_combined,
                "shift_px": shift_px,
                "anomaly_score": conf,
                "confidence": conf,
                "experimental_heuristic": True,
                "score_disclaimer": score_disclaimer,
                "status": status,
            }
        )

    for idx, (ny, nx) in enumerate(neg_yx):
        if idx in used_neg:
            continue
        z_neg = round(float(abs(z_map[ny, nx])), 1)
        cy, cx = float(ny), float(nx)
        cand_type = "Variable Source Candidate (Fading)"
        shift_px = 0.0
        sky, classification, matched_id, conf, status = _classify_candidate(
            cx, cy, shift_px, z_neg, z_thresh, grid_wcs, cat_coords, catalog_sources
        )
        candidates.append(
            {
                "type": cand_type,
                "flux_change_direction": "fading",
                "classification": classification,
                "matched_catalog_id": matched_id,
                "x_px": round(float(cx), 2),
                "y_px": round(float(cy), 2),
                "ra": round(float(sky.ra.deg), 6),
                "dec": round(float(sky.dec.deg), 6),
                "snr_zscore": z_neg,
                "shift_px": shift_px,
                "anomaly_score": conf,
                "confidence": conf,
                "experimental_heuristic": True,
                "score_disclaimer": score_disclaimer,
                "status": status,
            }
        )

    candidates.sort(key=lambda d: d["snr_zscore"], reverse=True)
    return candidates[:8]
