# ASTRAEA Backend — API Contract & Scientific Validation Summary

> [!NOTE]
> This file is maintained for backward compatibility with references inside `backend/`. The authoritative, comprehensive versions of these documents live in:
> - **[`docs/API_CONTRACT.md`](../docs/API_CONTRACT.md)** — Full HTTP API Contract
> - **[`docs/SCIENTIFIC_VALIDATION.md`](../docs/SCIENTIFIC_VALIDATION.md)** — Mathematical & Empirical Scientific Validation
> - **[`docs/DATA_PROVENANCE_AUDIT.md`](../docs/DATA_PROVENANCE_AUDIT.md)** — Complete Data Provenance & SHA-256 Ledger
> - **[`docs/AI_ORIGIN_AUDIT.md`](../docs/AI_ORIGIN_AUDIT.md)** — Forensic Audit of Removed Fabricated Claims & Reclassified Heuristics
> - **[`docs/KNOWN_LIMITATIONS.md`](../docs/KNOWN_LIMITATIONS.md)** — Explicit Scientific & Operational Limitations

---

## Part 1: Frontend–Backend API Contract Summary

All endpoints return JSON containing a top-level `"status"` field:
- `"ok"` — Successful processing with valid scientific data.
- `"empty"` — Valid request parameters and successful upstream/cache query, but no SPHEREx observations covered the query or passed validation.
- `"error"` — Invalid request parameters (`HTTP 400`), unresolvable object name (`HTTP 404`), upstream IRSA service failure (`HTTP 503`), or internal server error (`HTTP 500`).

### Endpoints

1. **`GET /api/targets`**
   - Returns the 5 NASA benchmark targets with coordinates propagated to exact `J2026.0` (`MJD 61041.0`, `2026-01-01T00:00:00.000` UTC), explicit `coordinate_epoch` metadata, and provenance-backed `data_status`:
     - `barnard`: `"verified_local_cache_subset"` (`30` verified cached FITS cutouts: `19 spherex_dr1`, `11 spherex_qr3`; `16` Detector 5, `14` Detector 6).
     - `eps_ind_b`, `wise0855`, `rw_aur`, `r_dra`: `"unverified_offline"` (`0` cached FITS files locally; live IRSA audit required via `/api/targets/audit?query_irsa_if_uncached=true`).
2. **`GET /api/targets/audit`**
   - Dynamically audits all 5 benchmark targets against local FITS cache (or live IRSA when `query_irsa_if_uncached=true`).
3. **`GET /api/search/coordinates`**
   - Validates and normalizes sky coordinates (`ra`, `dec`, `fov`, `detectors`, `max_frames`) and resolves nearby benchmark targets at `J2026.0`. Explicitly reports `"irsa_queried": false` and `"operation": "coordinate_validation_and_benchmark_matching_only"`.
4. **`GET /api/search/object`** (Alias: **`GET /api/resolve`**)
   - Resolves benchmark targets offline with exact spherical proper-motion propagation to `to_mjd` (default `61041.0`), or falls back to Astropy Sesame (`allow_external=true`) without double-applying proper motion.
5. **`GET /api/frames`**
   - Builds the validated 7-HDU SPHEREx cutout cube, reprojects `IMAGE` and `VARIANCE` onto a common North-up ICRS `TAN` grid (`6.154663"/px`), subtracts `ZODI` + flag-masked 3$\sigma$-clipped residual sky background (`MJy / sr`), computes detector-pure 30-day visit coadds (`min_epoch_coverage = 0.50`), and computes variance-propagated difference images (`median_propagated_sigma_mjy_sr`).
   - Encodes `NaN` / uncovered pixels as transparent RGBA (`alpha = 0`) in base64 PNGs.
6. **`GET /api/spectrum`**
   - Extracts per-detector aperture surface-brightness excess (`MJy / sr`) with `calibration_grade: false` and `wavelength_source: "WCS-WAVE_approximate"`.
7. **`GET /api/candidates`**
   - Returns heuristic candidate detections (`is_calibrated_probability: false`, `"Benchmark Coordinate Match — Unconfirmed (...)"`) and bitwise-union `sky_health` (`official_nasa_metric: false`).
8. **`GET /api/cache/status`**
   - Returns local cache counts and sizes with sanitized `"cache_dir": "cache"`.

---

## Part 2: Scientific Validation Summary (Verified Offline on 30 Barnard's Star Cutouts)

| Validation Stage | Verified Metric on 30 Cached Barnard's Star Cutouts |
| :--- | :--- |
| **FITS Provenance & Checksums** | 30/30 FITS files verified (`19 spherex_dr1`, `11 spherex_qr3`; `16` Detector 5, `14` Detector 6) with SHA-256 checksums in [`audit/data_provenance.json`](../audit/data_provenance.json). |
| **7-HDU Structural & Unit Check** | 30/30 contain `PRIMARY`, `IMAGE`, `FLAGS`, `VARIANCE`, `ZODI`, `PSF`, `WCS-WAVE` with `BUNIT = "MJy / sr"` and 100% matching finite masks between `IMAGE` and `VARIANCE`. |
| **WCS Round-Trip Accuracy** | Median numerical `pixel -> world -> pixel` RMS error of `4.90e-07 px` (`max <= 5.18e-07 px`) on the `59 x 59` North-up ICRS `TAN` grid (`6.154663"/px`). |
| **Empirical Cross-Frame Registration** | Across 24 full-coverage exposures in 6 multi-frame epochs, empirical flux-weighted source centroid dispersion is **`0.1268 px` (`0.7799 arcsec`) RMS**, max offset **`0.3793 px`**. |
| **Background Normalization** | Additive `ZODI` subtraction + source/defect-masked 3$\sigma$-clipped sky pedestal subtraction centers valid background medians at `0.000000 MJy/sr` while preserving `MJy / sr` units. |
| **Difference Imaging & Propagated Variance** | Detector 5 (`365.88 d` baseline, `100.0%` coverage): empirical MAD `0.020690 MJy/sr`, propagated $\sigma_{\text{diff}}$ `0.030369 MJy/sr`. Detector 6 (`363.26 d` baseline, `96.09%` coverage): empirical MAD `0.033373 MJy/sr`, propagated $\sigma_{\text{diff}}$ `0.048150 MJy/sr`. |
| **Proper-Motion + Parallax Dipole Recovery** | Both Detector 5 (`SNR = 48.07`) and Detector 6 (`SNR = 48.61`) detect Barnard's Star at `x = 29.0, y = 29.0` with a `4.0 px` (`24.62 arcsec`) dipole separation over ~1 year, labeled `"Benchmark Coordinate Match — Unconfirmed (Barnard's Star)"`. |
| **Bitwise-Union Sky Health** | Exact bitwise-OR union of non-astrophysical defect flags yields `96.8%` flag cleanliness (`score_pct = 92.1`, Grade `A`). |
