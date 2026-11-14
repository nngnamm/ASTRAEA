# Scientific Validation & Astronomical Methodology

**Repository**: `https://github.com/nngnamm/ASTRAEA`  
**Scope**: Team 02 — Backend & Astronomical Data Processing Pipeline

---

## 1. SPHEREx Instrument & Level-2 Data Product Architecture

NASA's **SPHEREx** (Spectro-Photometer for the History of the Universe, Epoch of Reionization and Ices Explorer) observes the sky across $0.75\text{–}5.00\,\mu\text{m}$ using six HgCdTe $2048 \times 2048$ detector arrays (Detectors 1–3 for $0.75\text{–}2.41\,\mu\text{m}$ and Detectors 4–6 for $2.42\text{–}5.00\,\mu\text{m}$) equipped with **Linear Variable Filters (LVFs)**.

Key physical properties enforced by the ASTRAEA backend:
1. **LVF Wavelength Gradient**: Unlike static broad-band filters, the effective central wavelength $\lambda(x, y)$ varies continuously across each detector's field of view. A single exposure does **not** produce 102 simultaneous full-field spectral images; rather, multiple dithered exposures step targets across the LVF gradient.
2. **Strict Detector Separation**: Each detector covers a distinct wavelength channel (`BANDS` in `backend/config.py`). Coaddition, display scaling, epoch grouping, and difference imaging **never mix observations from different detectors**.
3. **Physical Surface-Brightness Units**: SPHEREx Level-2 `IMAGE`, `ZODI`, and $\sqrt{\text{VARIANCE}}$ extensions are calibrated in **$\text{MJy/sr}$**. All background subtraction and coaddition operations are strictly additive/median-based so physical $\text{MJy/sr}$ amplitudes are preserved.

---

## 2. Coordinate Propagation & Reference Epoch (`propagate_coords`)

High-proper-motion benchmark objects (Barnard's Star at $10.39''/\text{yr}$, WISE 0855-0714 at $8.15''/\text{yr}$, Epsilon Indi B at $4.70''/\text{yr}$) move across dozens of SPHEREx pixels ($\sim 6.15''/\text{px}$) between their catalog epochs (`J2000.0` or `J2010.5`) and SPHEREx observation epochs (`2025–2026`).

- **Proper Motion Convention**: All catalog proper motions in `backend/config.py` use the IAU standard spherical tangent-plane convention:
  $$\mu_{\alpha^*} = \mu_\alpha \cos\delta \quad (\texttt{pm\_ra\_cosdec}), \qquad \mu_\delta \quad (\texttt{pm\_dec})$$
- **Exact Spherical Space Motion**: `propagate_coords()` in [`backend/spherex_core.py`](../backend/spherex_core.py) constructs an `astropy.coordinates.SkyCoord` in the `ICRS` frame at `obstime = Time(epoch0, format="jyear")` and applies `apply_space_motion(new_obstime=Time(to_mjd, format="mjd"))`.
- **Explicit Reference Epoch**:
  - Default propagation target is `DEFAULT_REFERENCE_MJD = 61041.0`, which is exact Julian year `J2026.0` (`2026-01-01T00:00:00.000` UTC).
  - During difference-image candidate cross-matching (`detect_anomalies`), benchmark coordinates are propagated to the actual midpoint MJD $\frac{1}{2}(\text{MJD}_{\text{early}} + \text{MJD}_{\text{late}})$ of the compared epochs.

---

## 3. FITS Validation & Quality-Flag Policy (`image_validation.py`)

`validate_image_payload()` separates validation into four explicit readiness dimensions:
1. **`structural_validity`**: Verifies 2-D numeric `IMAGE`, 2-D integer `FLAGS`, matching array dimensions across `IMAGE`/`FLAGS`/`VARIANCE`/`ZODI`, valid `DETECTOR` $\in \{1..6\}$, and finite positive `MJD-AVG`.
2. **`numerical_validity`**: Evaluates finite pixel fractions in `IMAGE`, `VARIANCE`, and `ZODI`, counts zero and negative variance pixels, and compares `IMAGE` vs. `VARIANCE` non-finite pixel masks (`image_variance_mask_consistency`).
3. **`alignment_validity`**: Verifies that the `IMAGE` header contains an invertible 2-D celestial WCS (`RA---TAN-SIP`, `DEC--TAN-SIP`) mapping the cutout center to finite sky coordinates.
4. **`scientific_operation_suitability`**: Reports operation-specific readiness booleans (`visual_display`, `background_estimation`, `epoch_coadd`, `difference_imaging`, `uncertainty_weighted_photometry`, `zodi_background_subtraction`, `psf_fitting`, `wavelength_lookup`).

### Header-Driven `FLAGS` Decoding & Bitwise Union Aggregation
- Each SPHEREx `FLAGS` HDU defines bit meanings via `MSKN0000..MSKN0021` and `MSKD0000..MSKD0021` header cards (e.g., Bit 0 `TRANSIENT`, Bit 1 `OVERFLOW`, Bit 6 `NONFUNC`, Bit 7 `OUTLIER`, Bit 9 `MISSING_DATA`, Bit 10 `HOT`, Bit 11 `COLD`, Bit 21 `SOURCE`).
- Because a single pixel can carry multiple bits simultaneously, `_decode_flag_header()` computes both per-bit counts and the exact **bitwise OR union** of non-astrophysical defect bits (`non_astrophysical_flagged_fraction`, excluding `SOURCE` and `TRANSIENT`).

---

## 4. WCS Registration & Empirical Cross-Frame Astrometric Verification

1. **Common Reference Grid (`make_grid`)**: Constructs a North-up, East-left `ICRS` tangential (`RA---TAN`, `DEC--TAN`) WCS grid centered on $(\alpha_0, \delta_0)$ with pixel scale equal to the median native pixel scale across loaded frames ($6.153''/\text{px}$ for the Barnard's Star dataset, $59 \times 59$ pixels for $\text{FOV}=0.1^\circ$).
2. **Reprojection & Mask Preservation (`align_frame_to_grid`)**: Reprojects `IMAGE` and `VARIANCE` onto the common grid via `reproject.reproject_interp`. Pixels outside the native detector footprint (`footprint <= 0` or non-finite) are strictly set to `np.nan` (never `0.0`).
3. **Empirical Cross-Frame Registration (`evaluate_cross_frame_registration`)**:
   - Single-frame WCS round-trip checks (`pixel_to_world` $\rightarrow$ `world_to_pixel`) only verify mathematical invertibility ($< 10^{-4}\text{ px}$), not physical alignment across exposures.
   - `evaluate_cross_frame_registration()` locates the primary compact source peak on each multi-frame visit coadd and measures flux-weighted first-moment centroids across individual aligned exposures within that visit epoch.
   - **Verified Result on the 30 Cached Barnard's Star Cutouts**: Across 24 full-coverage exposures in 6 multi-frame visit epochs, the empirical intra-epoch centroid RMS is **$0.1268\text{ px}$ ($0.780''$)**, with a maximum offset of **$0.3793\text{ px}$ ($2.33''$)**.
   - When no compact source with $\text{SNR} \ge 5$ is present in multi-frame epochs, the pipeline explicitly reports `status: "incomplete_no_reference_sources"` (`independent_astrometric_validation: False`).

---

## 5. Flag-Aware Background Normalization (`normalize_frame_background`)

SPHEREx exposures exhibit seasonal and roll-dependent zodiacal background variations ($0.2\text{–}1.5\,\text{MJy/sr}$ in Detectors 5 and 6).
1. **Zodiacal Model Subtraction**: When the `ZODI` HDU is present, it is reprojected onto the common grid and subtracted pixel-by-pixel (`aligned - aligned_zodi`).
2. **Flag-Aware Residual Sky Pedestal Removal**:
   - Constructs an operation-specific exclusion mask via `build_operation_flag_mask(..., operation="background_estimation")` that excludes pixels flagged as `SOURCE`, `TRANSIENT`, or detector defects (`OVERFLOW`, `NONFUNC`, `OUTLIER`, `HOT`, `COLD`, `MISSING_DATA`, `PERSIST`).
   - If excluding `SOURCE` leaves at least $\max(4, \min(25, 10\% \text{ of valid pixels}))$ unflagged sky pixels, the residual pedestal is estimated from those unflagged sky pixels using 3$\sigma$ MAD-clipped median (`sigma_clip(..., sigma=3.0, cenfunc="median", stdfunc=mad_std)`).
   - If a bright star's `SOURCE` halo covers almost the entire small cutout (e.g., `ebec6e19177daf652cc1bc959547e8c8.fits`, where `SOURCE` covers $97.8\%$ of pixels), it falls back to excluding `detector_defects_only` and then 3$\sigma$-clipped valid pixels, recording the exact `sky_estimation_mask_policy` and `flagged_pixels_excluded_from_sky` in `frame["validation"]["normalization"]`.

---

## 6. Coverage-Gated Difference Imaging & Variance Propagation

1. **Heuristic Visit Grouping (`group_epochs_by_detector`)**: Groups observations per detector into visit windows separated by $> 30\text{ days}$. Explicitly labeled as `method: "time_gap_heuristic"`, `official_mission_epoch: False`.
2. **Coverage Eligibility Gate (`_difference_for_detector`)**:
   - Requires at least 2 epochs with coadd spatial coverage $\ge 50\%$ (`MIN_EPOCH_COVERAGE = 0.50`).
   - **No Silent Fallback**: If fewer than 2 epochs meet the $50\%$ coverage threshold, the detector comparison returns `status: "ineligible_low_coverage"` (`comparison_eligible: False`) rather than silently falling back to low-coverage edge strips.
3. **Joint-Mask Pedestal & Difference (`compute_difference_with_diagnostics`)**:
   - Computes `common_mask = np.isfinite(img_early) & np.isfinite(img_late)` and requires `common_valid_pixel_fraction >= 0.10`.
   - Estimates and removes residual sky pedestals strictly over `common_mask` so non-overlapping edge regions cannot bias the difference baseline.
   - Propagates per-frame `VARIANCE` into epoch median-coadd variance maps ($\text{Var}(\text{median}) = \frac{\pi}{2 N^2}\sum_{i=1}^N \sigma_i^2$ for $N > 1$) and quadrature-sums $\sigma_{\text{diff}} = \sqrt{\sigma_{\text{early}}^2 + \sigma_{\text{late}}^2}$ over `common_mask`.
   - On the 30 cached Barnard's Star files:
     - **Detector 5**: `difference_noise_mad_mjy_sr = 0.02069`, `median_propagated_sigma_mjy_sr = 0.03037`, `common_valid_pixel_fraction = 1.0`.
     - **Detector 6**: `difference_noise_mad_mjy_sr = 0.03337`, `median_propagated_sigma_mjy_sr = 0.04815`, `common_valid_pixel_fraction = 0.9612`.
