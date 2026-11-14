# ASTRAEA Backend API Contract

**Base URL (Development)**: `http://localhost:8000`  
**Entry Point**: [`backend/main.py`](../backend/main.py) (`uvicorn main:app --reload --port 8000`)  
**Version**: `2.0.0`

---

## 1. Error Response Contract

All error responses across all endpoints return a consistent JSON structure with an appropriate HTTP status code:

```json
{
  "status": "error",
  "error_code": "invalid_parameters",
  "detail": "RA must satisfy 0.0 <= ra < 360.0 deg; received 400.0",
  "message": "RA must satisfy 0.0 <= ra < 360.0 deg; received 400.0"
}
```

### Standard `error_code` Values & HTTP Status Codes

| HTTP Status | `error_code` | Meaning |
| :---: | :--- | :--- |
| `400` | `invalid_parameters` | Out-of-range `ra`, `dec`, `fov` ($[0.03, 0.30]$ deg), `detectors` ($\subseteq \{1..6\}$), or `max_frames` ($[1, 100]$). |
| `400` | `invalid_object_name` | Empty or whitespace-only `name` parameter on `/api/search/object`. |
| `404` | `object_not_found` | Object name could not be resolved in the benchmark catalog or CDS Sesame. |
| `503` | `upstream_service_unavailable` | All configured IRSA SIA collections failed (network/timeout/5xx error) and no valid complete cache exists. Never disguised as `no_observations_found`. |
| `503` | `partial_upstream_failure` | One or more IRSA SIA collections failed and no observations were returned for the requested detectors. |
| `200` (`status: "empty"`) | `no_observations_found` | Complete IRSA SIA query succeeded across all collections and returned 0 matching observations. |
| `200` (`status: "empty"`) | `validation_rejected` | Observations were found, but all downloaded/cached FITS files failed structural/WCS validation. |
| `200` (`status: "empty"`) | `alignment_failed` | Loaded FITS observations failed WCS reprojection onto the common TAN grid. |
| `500` | `internal_server_error` | Unhandled server exception; internal stack traces and host filesystem paths are never leaked. |

---

## 2. Endpoints

### 2.1 `GET /api/targets`
Returns the five NASA Space Apps 2026 benchmark targets, their catalog coordinates, exact `J2026.0` (`MJD 61041.0`) space-motion propagated coordinates, and provenance-backed local cache status.

- **Response (`200 OK`)**:
  - `status`: `"ok"`
  - `reference_epoch`: `{"jyear": 2026.0, "mjd": 61041.0, "iso_utc": "2026-01-01T00:00:00.000", "pm_convention": "pm_ra_cosdec"}`
  - `bands_um`: Wavelength passbands ($\mu\text{m}$) for Detectors `1..6`.
  - `targets`: Map of `barnard`, `eps_ind_b`, `wise0855`, `rw_aur`, `r_dra`. Each target includes:
    - `ra`, `dec`, `epoch`, `pm_ra_arcsec_yr`, `pm_dec_arcsec_yr`, `pm_convention` (`"pm_ra_cosdec"`)
    - `ra_2026`, `dec_2026`, `propagated_epoch_jyear` (`2026.0`), `propagated_to_mjd` (`61041.0`)
    - `coordinate_source`, `proper_motion_source`
    - `data_status`: `"verified_local_cache_subset"` (for `barnard`) or `"unverified_offline"` (for the other 4 targets).
    - `provenance_audit`: Honest local cache verification counts (`cached_fits_count`, `collections_verified_in_cache`, `by_detector_in_cache`).

---

### 2.2 `GET /api/targets/audit`
Runs a reproducible provenance and scientific readiness audit across all five benchmark targets (`audit_all_benchmark_targets()`).

- **Query Parameters**:
  - `query_live_irsa` (`bool`, default `false`): When `false`, performs a deterministic offline audit of local cache files and existing SIA caches without making network calls.
- **Response (`200 OK`)**:
  - `status`: `"ok"`
  - `checked_at_utc`: ISO-8601 UTC timestamp
  - `reference_mjd`: `61041.0`
  - `targets`: Per-target audit report (`verification_status`, `catalog_coordinates`, `propagated_coordinates`, `sia_query`, `local_cache`, `alignment_and_comparison`, `comparison_eligible`).

---

### 2.3 `GET /api/search/coordinates`
Validates and normalizes decimal-degree or sexagesimal sky coordinates and checks proximity ($\le 0.1^\circ$) to NASA benchmark targets.

> **Important Contract Note**: This endpoint **does not query IRSA** for observation availability (`irsa_queried: false`, `data_availability_checked: false`). A `200 OK` with `valid: true` means the coordinate parameters are syntactically and astronomically valid; call `/api/frames` to retrieve observations.

- **Query Parameters**:
  - `ra` (required): Decimal degrees (`[0, 360)`) or sexagesimal hourangle (e.g. `"08h55m10.831s"`).
  - `dec` (required): Decimal degrees (`[-90, 90]`) or sexagesimal degrees (e.g. `"-7d14m42.53s"`).
  - `fov` (`float`, default `0.1`): Field of view in degrees (`[0.03, 0.30]`).
  - `detectors` (`str`, default `"5,6"`): Comma-separated detector IDs in `1..6`.
  - `max_frames` (`int`, default `30`): Maximum frames to load (`[1, 100]`).

---

### 2.4 `GET /api/search/object` (Alias: `GET /api/resolve`)
Resolves an astronomical object name or alias against the NASA benchmark catalog first (applying exact space-motion propagation to `to_mjd`), falling back to Astropy CDS Sesame (`SkyCoord.from_name`) when `allow_external=true`.

- **Query Parameters**:
  - `name` (required): Object name or alias (e.g., `"Barnard's Star"`, `"GJ 699"`, `"WISE 0855"`).
  - `to_mjd` (`float`, default `61041.0`): Target MJD for proper-motion propagation (`61041.0` = `J2026.0`).
  - `allow_external` (`bool`, default `true`): Whether to query CDS Sesame if `name` is not a benchmark target.
- **Response (`200 OK`)**:
  - Includes `ra`, `dec`, `source` (`"nasa_benchmark_catalog"` or `"astropy_sesame"`), `benchmark_target` (`bool`), `proper_motion_applied` (`bool`), `propagated_to_mjd`, `propagated_epoch_jyear`, and `pm_convention`.

---

### 2.5 `GET /api/frames`
Queries IRSA SIA (or local SIA cache), loads and validates SPHEREx Level-2 FITS cutouts, reprojects them onto a common North-up ICRS TAN grid, subtracts `ZODI` and flag-aware sky pedestals, groups visits by detector (`gap_days=30`), computes coverage-gated epoch differences with variance propagation, and runs heuristic candidate detection and Sky Health scoring.

- **Query Parameters**: `ra`, `dec`, `fov` (default `0.1`), `detectors` (default `"5,6"`), `max_frames` (default `30`).
- **Response (`200 OK` when observations or valid empty state; `503` when upstream IRSA fails)**:
  - `status`: `"ok"` | `"empty"` | `"error"`
  - `center`: `{"ra": ..., "dec": ...}`
  - `grid_size_px`, `pix_scale_arcsec`, `wcs_summary`, `processing_history`
  - `epoch_grouping_policy`: `{"method": "time_gap_heuristic", "gap_days": 30, "official_mission_epoch": false, ...}`
  - `frames`: Ordered list of validated, aligned, normalized exposures. Each frame includes `mjd`, `iso_date`, `detector`, `wave_um`, `collection`, `image_shape_yx`, `coverage_fraction`, `invalid_pixel_rendering: "transparent_rgba_alpha_0"`, `metadata` (`observation_metadata.py` schema v1.1 with `sha256`), `validation` (`image_validation.py` report), and `png_b64` (transparent RGBA for `NaN` pixels).
  - `epochs`: Per-detector coadded visit summaries (`mjd_mean`, `mjd_start`, `mjd_end`, `n_frames`, `single_frame_epoch`, `valid_pixel_fraction`, `png_b64`).
  - `differences`: Per-detector comparison results (`status`: `"available"` | `"insufficient_epochs"` | `"ineligible_low_coverage"` | `"processing_error"`, `comparison_eligible`, `diagnostics` including `variance_propagated` and `median_propagated_sigma_mjy_sr`, `candidates`, `png_b64`).
  - `sky_health` & `experimental_analytics`: Explicitly flagged with `experimental_heuristic: true` and `official_nasa_metric: false`.
  - `validation_summary`: Counts of queried, valid, aligned, and rejected frames, plus `sia_query_diagnostics`, `alignment_summary.cross_frame_registration`, and `normalization_summary`.

---

### 2.6 `GET /api/candidates`
Returns the difference-image anomaly candidates and Sky Health summary for the requested sky position and detectors.

- **Candidate Object Contract**:
  - `type`: `"Moving Object Candidate"` | `"Variable Source Candidate (Brightening)"` | `"Variable Source Candidate (Fading)"`
  - `classification`: `"Benchmark Coordinate Match — Unconfirmed (<name>)"` | `"Catalog Coordinate Match — Unconfirmed (<id>)"` | `"Anomaly Candidate"` (never `"Confirmed"`)
  - `association_type`: `"benchmark_coordinate_match_unconfirmed"` | `"catalog_cross_match"` | `"unassociated_candidate"`
  - `independent_confirmation`: `false`
  - `human_reviewed`: `false`
  - `is_calibrated_probability`: `false`
  - `experimental_heuristic`: `true`
  - `snr_zscore`, `shift_px`, `confidence`, `anomaly_score`, `ra`, `dec`, `x_px`, `y_px`

---

### 2.7 `GET /api/spectrum`
Computes approximate 6-detector circular-aperture ($r = 2.5\text{ px}$) surface-brightness excess measurements (`MJy / sr`) across Detectors 1–6. Explicitly includes `spectral_metadata.calibration_grade = false` and documents SPHEREx LVF characteristics.

---

### 2.8 `GET /api/cache/status`
Returns local cache file counts and disk usage (`fits_files_count`, `sia_queries_count`, `corrupt_quarantined_count`, `total_size_bytes`, `total_size_mb`, and a sanitized basename `cache_dir: "cache"` that never leaks host paths).

---

### 9. `GET /api/sources`
Returns the machine-readable manifest of ASTRAEA's approved scientific data source ecosystem (`SPHEREx`, `Gaia`, `NEOSSat`, `WISE`, `Synthetic`), official SPHEREx collections (`spherex_dr1`, `spherex_dr1_deep`, `spherex_qr3`, `spherex_qr3_deep`), Gaia DR3 astrometric reference records for all 5 benchmark targets (including `source_id`, `parallax_mas`, `radial_velocity_km_s`, `ruwe`, and independent IR astrometric references for targets absent from Gaia DR3), dataset separation rules, and optional dynamic SPHEREx collection discovery (`?discover_live=false` by default).

