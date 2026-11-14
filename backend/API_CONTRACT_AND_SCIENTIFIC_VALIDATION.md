# ASTRAEA Backend — API Contract (`INT01`) & Scientific Validation Report (`BE17`)

**System**: Project ASTRAEA SPHEREx Multi-Epoch Backend (Team 02 — Data Processing & API)
**Sample Handoff Fixture (`INT02`)**: `sample_barnard_response.json` (566 KB, generated from 30 validated SPHEREx QR2/DR1 FITS cutouts of Barnard's Star across Detectors 5 & 6).

---

## Part 1: Frontend–Backend API Contract (`INT01`)

All endpoints return JSON containing a top-level `"status"` field with one of three values:
- `"ok"` — Successful processing with valid scientific data.
- `"empty"` — Valid request parameters, but no SPHEREx observations covered the query or passed structural/alignment validation.
- `"error"` — Invalid request parameters (HTTP 400) or unresolvable object name (HTTP 404).

### Structured Error Schema (HTTP 400 / 404)
```json
{
  "status": "error",
  "error_code": "invalid_parameters | invalid_object_name | object_not_found",
  "detail": "Descriptive human-readable explanation",
  "message": "Descriptive human-readable explanation"
}
```

### 1. `GET /api/targets`
Returns the 5 verified NASA Space Apps 2026 benchmark targets, their exact spherical proper-motion-propagated coordinates (`ra_2026`, `dec_2026`), aliases, recommended detectors, and empirical IRSA observation counts (`irsa_audit`), along with the 6 SPHEREx detector bandpass ranges (`bands_um`).

### 2. `GET /api/search/coordinates`
Validates and normalizes sky coordinates before triggering a heavy cutout build.
- **Query Parameters**:
  - `ra` *(required)*: Decimal degrees (`0.0 <= ra < 360.0`) or sexagesimal string (`"08h55m10.831s"`).
  - `dec` *(required)*: Decimal degrees (`-90.0 <= dec <= 90.0`) or sexagesimal string (`"-7d14m42.53s"`).
  - `fov` *(optional, default `0.1`)*: Field of view in degrees (`0.03 <= fov <= 0.30`).
  - `detectors` *(optional, default `"5,6"`)*: Comma-separated detector IDs in `1..6`.
  - `max_frames` *(optional, default `30`)*: Integer in `1..100`.
- **Response (`200 OK`)**:
  `{"status": "ok", "valid": true, "ra": 133.795129, "dec": -7.245147, "hms_dms": "08h55m10.831s -07d14m42.530s", "fov_deg": 0.1, "detectors": [5, 6], "max_frames": 30, "nearby_benchmark_target": {...}}`

### 3. `GET /api/search/object` (Alias: `GET /api/resolve`)
Resolves an astronomical object name or benchmark alias (`"Barnard's Star"`, `"GJ 699"`, `"WISE 0855"`, `"eps Ind B"`, `"RW Aur"`, `"R Dra"`) offline against the verified NASA benchmark catalog with exact spherical proper-motion propagation, or falls back to Astropy Sesame when `allow_external=true`.
- **Query Parameters**:
  - `name` *(required)*: Object name or alias.
  - `to_mjd` *(optional, default `61000.0`)*: Target epoch MJD for proper-motion propagation.
  - `allow_external` *(optional, default `true`)*: Allow external Sesame resolution for non-benchmark targets.

### 4. `GET /api/frames`
Primary multi-epoch image, coadd, difference, and diagnostic cube endpoint (see `sample_barnard_response.json` for a complete real example).
- **Query Parameters**: `ra`, `dec`, `fov` (default `0.1`), `detectors` (default `"5,6"`), `max_frames` (default `30`).
- **Top-Level Response Keys (`status="ok"`)**:
  - `"status"`: `"ok"`
  - `"center"`: `{"ra": ..., "dec": ...}`
  - `"grid_size_px"`: Common ICRS North-up TAN grid dimension (`59` for `fov=0.1` at `6.155"/px`).
  - `"pix_scale_arcsec"`: Reference pixel scale (`6.155`).
  - `"display_limits_by_detector_mjy_sr"`: Shared `[vmin, vmax]` per detector so partial-footprint frames never flash white during animation.
  - `"epoch_grouping_policy"`: Explicit provenance record (`"method": "time_gap_heuristic"`, `"gap_days": 30`, `"official_mission_epoch": false`).
  - `"frames"`: Chronologically ordered validated frames, each containing `"mjd"`, `"iso_date"`, `"detector"`, `"wave_um"`, `"display_limits_mjy_sr"`, `"metadata"` (OBSID, timing, exposure, spatial/spectral provenance), `"validation"` (7-HDU check, decoded `flag_definitions`, `alignment`, `normalization`), and `"png_b64"`.
  - `"epochs"`: Detector-pure visit-group coadds with `"mjd_mean"`, `"mjd_start"`, `"mjd_end"`, `"n_frames"`, `"valid_pixel_fraction"`, `"grouping_heuristic_gap_days"`, and `"png_b64"`.
  - `"differences"`: Per-detector coverage-aware difference objects (`"detector"`, `"status"`, `"mjd_early"`, `"mjd_late"`, `"baseline_days"`, `"common_valid_pixel_fraction"`, `"target_covered_in_difference"`, `"difference_noise_mad_mjy_sr"`, `"diagnostics"`, `"candidates"`, `"warnings"`, `"vlim_mjy_sr"`, `"png_b64"`).
  - `"difference"`: Legacy primary difference object with `"by_detector": differences`.
  - `"sky_health"` & `"experimental_analytics"`: Derived Sky Health Index and candidate detections per detector, explicitly labeled with `"experimental_heuristic": true`.
  - `"validation_summary"`: End-to-end FITS validation, WCS registration, and ZODI+sky normalization audit counters.

### 5. `GET /api/spectrum`
Returns per-detector aperture surface-brightness excess measurements across Detectors 1–6 along with explicit LVF and `WCS-WAVE` scientific provenance.
- **Query Parameters**: `ra`, `dec`, `fov` (default `0.1`).
- **Response Keys**: `"status"`, `"center"`, `"fov_deg"`, `"spectrum"` (items include `"detector"`, `"n_frames"`, `"wave_um"`, `"flux_mjy_sr"`, `"err_mjy_sr"`, `"snr"`, `"band_range_um"`, `"sampled_wavelength_range_um"`, `"calibration_grade": false`, `"wavelength_source": "WCS-WAVE_approximate"`, `"observation_ids"`), and `"spectral_metadata"`.

### 6. `GET /api/candidates`
Lightweight endpoint returning `"candidates_by_detector"`, `"sky_health"`, and `"disclaimer"` for a sky position without base64 image payloads.

### 7. `GET /api/cache/status`
Returns local FITS cutout and SIA query cache counts and disk usage (`"fits_files_count"`, `"sia_queries_count"`, `"total_size_bytes"`, `"total_size_mb"`).

---

## Part 2: Scientific Validation Report (`BE17`)

### 1. Live IRSA Data Availability Across All 5 NASA Benchmark Objects (`BE12`)
All 5 benchmark targets from the NASA Space Apps 2026 specification were audited against live IRSA SPHEREx collections (`spherex_dr1`, `spherex_qr2`):
| Target Key | Object Name | Blueprint Test Case | Total FITS | Baseline (days) | D1–D6 Coverage | 30-Day Visit Groups |
| :--- | :--- | :--- | :---: | :---: | :--- | :--- |
| `barnard` | Barnard's Star | Proper Motion | 386 | 402.8 d | All 6 bands (62–68/det) | 3–4 epochs/det |
| `eps_ind_b` | Epsilon Indi B | Brown Dwarf Binary System | 486 | 525.6 d | All 6 bands (49–103/det) | 3–5 epochs/det |
| `wise0855` | WISE 0855-0714 | Fast-moving Brown Dwarf | 263 | 373.6 d | All 6 bands (37–49/det) | 3–4 epochs/det |
| `rw_aur` | RW Aurigae | Variable Star | 342 | 380.5 d | All 6 bands (57/det) | 3 epochs/det |
| `r_dra` | R Draconis | Variable Star | 578 | 425.9 d | All 6 bands (92–102/det) | 3 epochs/det |

### 2. FITS Structural & Flag Validation (`BE07`–`BE09`)
- Verified across 30 real SPHEREx Level-2 cutout FITS files (15 Detector 5, 15 Detector 6): all 30 contain the complete 7-HDU structure (`PRIMARY`, `IMAGE`, `FLAGS`, `VARIANCE`, `ZODI`, `PSF`, `WCS-WAVE`), unit `MJy / sr`, and 100% consistent non-finite masks between `IMAGE` and `VARIANCE`.
- Per-bit quality flags (`MSKN0000`–`MSKN0021`) are decoded directly from each FITS `FLAGS` header (`SUR_ERROR`, `NONFUNC`, `MISSING_DATA`, `HOT`, `COLD`, `OUTLIER`, `PERSIST`, `SOURCE`, `TRANSIENT`, etc.) without applying a destructive global mask or penalizing astrophysical `"SOURCE"`/`"TRANSIENT"` bits in `compute_sky_health`.

### 3. WCS Registration & Background Normalization (`BE10`–`BE11`)
- All 30 cached Barnard's Star frames reproject onto a common `59 x 59` North-up ICRS TAN grid (`6.154663"/px`) with sub-micro-pixel WCS round-trip error (`<= 5.18e-07 px`) and empirical centroid scatter of `0.067 px` (`0.41 arcsec`) within an epoch.
- Uncovered footprints are preserved strictly as `np.nan` (never `0.0`).
- Additive ZODI model subtraction + sigma-clipped (`3.0 sigma`, MAD) residual sky pedestal removal reduces inter-frame background median scatter from `2.5–3.1 MJy/sr` to `< 0.004 MJy/sr` while preserving physical `MJy / sr` surface-brightness units.

### 4. Coverage-Aware Difference Imaging & Proper-Motion Verification (`BE13`, `BE17`)
- On Barnard's Star (`ra=269.446294, dec=4.767881`), coverage-aware epoch selection (`min_epoch_coverage=0.50`) pairs well-covered multi-frame epochs across a ~1-year baseline (`365.88 d` in Detector 5 with `100.0%` joint coverage; `363.26 d` in Detector 6 with `96.1%` joint coverage), achieving background difference noise MAD of `0.021 MJy/sr` (D5) and `0.033 MJy/sr` (D6).
- Both Detector 5 and Detector 6 independently recover Barnard's Star at the cutout center (`x=29.0, y=29.0`) as a `"Moving Object Candidate"` (`Confirmed Benchmark Object (Barnard's Star)`) with `SNR = 48.1` (D5) and `SNR = 48.6` (D6) and a measured dipole separation of `4.0 px` (~24.6 arcsec, caused by proper motion plus annual parallax shift across the 1-year baseline).
