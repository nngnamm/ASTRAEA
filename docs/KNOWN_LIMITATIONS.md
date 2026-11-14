# Known Scientific & Technical Limitations

**Repository**: `https://github.com/nngnamm/ASTRAEA`  
**Scope**: Team 02 — Backend, Data Pipeline & Astronomical Processing

In accordance with scientific software integrity standards, this document records the known physical, algorithmic, and operational limitations of the ASTRAEA backend.

---

## 1. Archive & Offline Data Availability Limitations

1. **Local Cache Scope (`barnard` Only)**:
   - Only **Barnard's Star (`barnard`)** has locally cached SPHEREx Level-2 FITS cutouts (`30` files across Detectors 5 and 6) in `backend/cache/`.
   - Even for Barnard's Star, the local cache is a **30-frame D5/D6 subset**, not the complete 6-detector SIA catalog.
   - The remaining four benchmark targets (`eps_ind_b`, `wise0855`, `rw_aur`, `r_dra`) have `"data_status": "unverified_offline"` in this repository and require live network connectivity to `https://irsa.ipac.caltech.edu` to query or download FITS products.
2. **Faint Brown Dwarf Sensitivity (`wise0855`)**:
   - Exploratory analysis in `notebooks/NASA_backend_exploration.ipynb` showed that WISE 0855-0714 ($T_{\text{eff}} \approx 250\text{ K}$, $W2 \approx 14.0\text{ mag}$) lies near the single-exposure surface-brightness noise floor of SPHEREx Level-2 cutouts. Detecting its proper-motion dipole requires multi-exposure coaddition and careful PSF-matched extraction rather than simple single-visit peak thresholding.

---

## 2. Instrument & Spectral Processing Limitations

1. **Linear Variable Filter (LVF) Chromatic Gradients**:
   - SPHEREx uses Linear Variable Filters where wavelength varies continuously across each detector. Two exposures of the same target taken at different spacecraft roll angles or field offsets sample slightly different wavelengths within the same detector band.
   - Consequently, a flux difference $\Delta I = I_{\text{late}} - I_{\text{early}}$ within a single detector can contain a chromatic component if the source has a steep spectral slope and $\Delta\lambda > 0.15\,\mu\text{m}$. The pipeline monitors `mean_wavelength_shift_um` and emits a warning when $\Delta\lambda > 0.15\,\mu\text{m}$, but does not perform spectral interpolation across the LVF gradient before differencing.
2. **Approximate Wavelength Lookup (`WCS-WAVE`)**:
   - The `WCS-WAVE` binary table extension provides a coarse $11 \times 11$ lookup grid of wavelength and bandpass values. As documented in `observation_metadata.py` and `/api/spectrum`, wavelengths extracted from `WCS-WAVE` (or fallback detector passband midpoints) carry `calibration_grade: False` and are for context only. Science-grade spectrophotometry requires the official SPHEREx Spectral WCS (`CWAVE`/`CBAND`) calibration products and PSF-fitted optimal extraction.
3. **Aperture Excess Spectrum (`compute_spectrum`)**:
   - `/api/spectrum` computes a fixed circular aperture ($r = 2.5\text{ px} \approx 15.4''$) surface-brightness excess in $\text{MJy/sr}$ above surrounding background apertures. It does **not** perform aperture-correction or PSF-fitted point-source flux density conversion (to Jy or AB magnitudes).

---

## 3. Image Alignment, Background & Difference Limitations

1. **Bilinear Interpolation During Reprojection**:
   - `align_frame_to_grid()` uses `reproject.reproject_interp` (bilinear interpolation in tangent-plane coordinates). While computationally fast and sufficient for visualization and coarse dipole detection, bilinear interpolation introduces slight spatial smoothing and correlated pixel noise compared to flux-conserving drizzle (`reproject_adaptive` or `reproject_exact`) or native-pixel forward modeling.
2. **No PSF Homogenization (ZOGY / Alard-Lupton) Before Differencing**:
   - `compute_difference_with_diagnostics()` subtracts background-normalized epoch coadds directly on the common TAN grid without convolving to a matched Point Spread Function (PSF) kernel. Asymmetric PSF wings rotated between spacecraft roll angles can produce false positive/negative residuals around very bright stars.
3. **Background Estimation in Crowded or Nebulous Cutouts**:
   - `normalize_frame_background()` subtracts the `ZODI` model and a flag-aware 3$\sigma$-clipped median pedestal. In small cutouts ($0.1^\circ \times 0.1^\circ$) dominated by extended reflection nebulosity (e.g., around T Tauri stars like RW Aurigae) or a bright star's extended halo (where the `SOURCE` flag bit 21 can cover $>90\%$ of pixels), local median sky estimation can partially absorb diffuse extended emission.

---

## 4. Heuristic Analytics Limitations (`analytics.py`)

1. **Epoch Grouping Heuristic (`gap_days=30`)**:
   - `group_epochs_by_detector()` clusters observations separated by $>30\text{ days}$. This is a simple temporal clustering heuristic (`official_mission_epoch: False`), not an official NASA SPHEREx survey-plan pass boundary.
2. **Sky Health Score (`compute_sky_health`)**:
   - The Sky Health score ($0.50 \times \text{coverage} + 0.30 \times \text{clean\_flags} + 0.20 \times \text{baseline}$) is an **experimental engineering heuristic** (`experimental_heuristic: True`, `official_nasa_metric: False`) designed for quick-look dataset triage. Its component weights and 180-day saturation scale are not derived from a physical instrument model.
3. **Candidate Detection & Ranking (`detect_anomalies`)**:
   - `detect_anomalies()` uses a Gaussian-smoothed local MAD $z$-score threshold (`z_thresh=4.5`) and pairs opposite-sign peaks within $0.8\text{–}5.5\text{ px}$ as moving-object dipole candidates.
   - The returned `confidence` / `anomaly_score` is a deterministic monotonic function of SNR (`is_calibrated_probability: False`), **not** a Bayesian or ML-calibrated probability of a genuine astrophysical transient.
   - Spatial matches within $12''$ of a benchmark target or $6.5''$ of an AllWISE source are reported strictly as `"Benchmark Coordinate Match — Unconfirmed"` or `"Catalog Coordinate Match — Unconfirmed"` (`independent_confirmation: False`, `human_reviewed: False`).
