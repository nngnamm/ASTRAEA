# ASTRAEA — COSMIC PULSE

**SPHEREx Multi-Epoch Infrared Time-Domain Inspection, Reprojection, Difference Imaging & Provenance-Verified Backend**

ASTRAEA — COSMIC PULSE is a scientific backend and interactive visual analytics system for inspecting multi-epoch infrared observations from NASA's **SPHEREx** (Spectro-Photometer for the History of the Universe, Epoch of Reionization, and Ices Explorer) mission hosted at the NASA/IPAC Infrared Science Archive (IRSA).

This repository has undergone a comprehensive **Data Provenance, Scientific Integrity, and Backend Engineering Audit (Team 02 — Backend + Data/Astronomy)** to ensure every scientific calculation, coordinate epoch, FITS validation check, uncertainty propagation step, and availability claim is reproducible and traceable to primary sources.

---

## Key Capabilities & Scientific Safeguards

1. **Multi-Collection IRSA SIA v2 Query & Provenance Tracking (`spherex_core.py`, `observation_metadata.py`)**
   - Queries IRSA SIA v2 collections (`spherex_qr2`, `spherex_dr1`) and differentiates `complete_with_observations`, `complete_empty`, `partial_failure`, and `upstream_failure`.
   - Incomplete multi-collection queries are never cached as complete results, and upstream outages return HTTP 503 (`upstream_service_unavailable`) rather than masquerading as empty sky.
   - Extracts full FITS provenance (`OBSID`, `MJD-BEG`, `MJD-AVG`, `MJD-END`, product `DATE`, `DETECTOR`, `BUNIT`, `SHA-256` checksum, and upstream `access_url`) and sorts time series strictly by `MJD-AVG`.

2. **Atomic FITS Download, Cache Verification & 7-HDU Validation (`image_validation.py`, `spherex_core.py`)**
   - Validates all 7 required SPHEREx Level-2 HDUs (`PRIMARY`, `IMAGE`, `FLAGS`, `VARIANCE`, `ZODI`, `PSF`, `WCS-WAVE`), array shapes, physical units (`MJy / sr`), and WCS metadata on a temporary `.tmp` file before atomic promotion into `backend/cache/`.
   - Automatically quarantines/evicts corrupted cached FITS files and retries download safely.
   - Separates validation diagnostics into `structural_validity`, `numerical_validity`, `alignment_validity`, and `scientific_operation_suitability`.

3. **Exact Spherical Proper-Motion Propagation (`spherex_core.py`, `config.py`)**
   - Propagates catalog positions using Astropy's exact spherical space-motion solver (`SkyCoord.apply_space_motion`) with explicit reference epochs (`DEFAULT_REFERENCE_MJD = 61041.0`, corresponding to `J2026.0` / `2026-01-01T00:00:00.000` UTC).
   - Propagates benchmark coordinates to the **actual midpoint MJD** of compared SPHEREx epochs when evaluating positional coincidence in difference images.

4. **WCS Reprojection, Variance Propagation & Empirical Registration Audit (`spherex_core.py`)**
   - Reprojects both `IMAGE` and `VARIANCE` onto a common North-up ICRS gnomonic (`TAN`) pixel grid (`6.154663"/px`), preserving uncovered pixels strictly as `np.nan` (and rendering them as transparent `alpha = 0` RGBA pixels in PNG previews).
   - Measures both numerical WCS round-trip error (`~4.9e-7 px` RMS) and **empirical cross-frame source centroid dispersion** (`0.1268 px` / `0.7799 arcsec` RMS across 24 full-coverage Barnard's Star exposures).

5. **Operation-Specific Flag Masking & ZODI + Residual Sky Normalization (`image_validation.py`, `spherex_core.py`)**
   - Subtracts the per-exposure `ZODI` zodiacal light model and estimates the residual sky pedestal using 3$\sigma$-clipped median background estimation after masking `SOURCE`, `TRANSIENT`, and detector defect flags (`build_operation_flag_mask(..., operation="background_estimation")`).
   - Preserves physical surface-brightness units (`MJy / sr`) without multiplicative rescaling.

6. **Coverage-Aware Epoch Coadds, Propagated Variance & Honest Heuristic Labeling (`serializers.py`, `analytics.py`)**
   - Enforces strict minimum epoch coverage (`min_epoch_coverage = 0.50`) and minimum overlap (`min_common_coverage = 0.10`); refuses to silently fall back to low-coverage edge cutouts (`status = "ineligible_low_coverage"`).
   - Propagates per-pixel `VARIANCE` into median-coadd variance maps ($\text{Var}(\text{median}) \approx \frac{\pi}{2 N^2}\sum \sigma_i^2$) and quadrature difference uncertainty (`median_propagated_sigma_mjy_sr`).
   - Clearly labels candidate detections (`"Benchmark Coordinate Match — Unconfirmed (...)"`, `is_calibrated_probability: false`) and the project-defined `Sky Health` index (`official_nasa_metric: false`, computed via exact bitwise-OR union of non-astrophysical defect flags).

---

## Verified Data Provenance Summary

Every dataset and benchmark target in this repository has been audited and recorded in [`docs/DATA_PROVENANCE_AUDIT.md`](docs/DATA_PROVENANCE_AUDIT.md) and [`audit/data_provenance.json`](audit/data_provenance.json):

| Benchmark Key | Object Name | Offline Audit Status | Verified Local FITS Cutouts | Detectors Verified Offline | Multi-Epoch Comparison Status |
| :--- | :--- | :--- | :---: | :---: | :--- |
| `barnard` | Barnard's Star | `verified_local_cache_subset` | **30** (`19 spherex_dr1`, `11 spherex_qr3`) | Detector 5 (**16**), Detector 6 (**14**) | **Eligible & Verified** (`365.88 d` D5, `363.26 d` D6) |
| `eps_ind_b` | Epsilon Indi B | `unverified_offline` | 0 cached locally | None offline | Requires live IRSA network access (`--live`) |
| `wise0855` | WISE 0855-0714 | `unverified_offline` | 0 cached locally | None offline | Requires live IRSA network access (`--live`) |
| `rw_aur` | RW Aurigae | `unverified_offline` | 0 cached locally | None offline | Requires live IRSA network access (`--live`) |
| `r_dra` | R Draconis | `unverified_offline` | 0 cached locally | None offline | Requires live IRSA network access (`--live`) |

---

## Repository Structure

```text
ASTRAEA/
├── audit/
│   └── data_provenance.json              # Machine-readable provenance ledger (SHA-256, URLs, MJDs)
├── backend/
│   ├── main.py                           # FastAPI application entry point & route handlers
│   ├── config.py                         # Configuration, epochs, & provenance-backed target catalog
│   ├── spherex_core.py                   # IRSA SIA client, FITS loader, WCS alignment, difference engine
│   ├── image_validation.py               # 7-HDU validator, bitwise flag masks, validation categories
│   ├── observation_metadata.py           # FITS header & SHA-256 provenance extractor (schema v1.1)
│   ├── serializers.py                    # Epoch grouping, RGBA PNG encoder, JSON response builder
│   ├── analytics.py                      # Heuristic candidate detection & Sky Health index
│   ├── diagnostics/
│   │   ├── audit_benchmark_targets.py    # Reproducible CLI audit for all 5 benchmark targets
│   │   └── regenerate_barnard_sample.py  # Regenerates sample_barnard_response.json from cached FITS
│   ├── cache/                            # 30 verified Barnard's Star FITS cutouts & SIA query cache
│   ├── tests/                            # 79 deterministic unit/provenance tests + 1 live IRSA test
│   └── requirements.txt                  # Pinned Python dependencies
├── docs/
│   ├── DATA_PROVENANCE_AUDIT.md          # Complete data provenance & FITS SHA-256 ledger
│   ├── AI_ORIGIN_AUDIT.md                # Forensic audit of hardcoded, heuristic, & prior fabricated claims
│   ├── SCIENTIFIC_VALIDATION.md          # Mathematical formulas, units, WCS & uncertainty validation
│   ├── KNOWN_LIMITATIONS.md              # Explicit scientific, coverage, & software limitations
│   └── API_CONTRACT.md                   # Complete HTTP API reference & response schemas
├── frontend/                             # Next.js frontend application (Team 01)
├── pytest.ini                            # Pytest configuration & markers
└── .github/workflows/ci.yml              # Automated GitHub Actions CI pipeline
```

---

## Quickstart (Backend & Scientific Audit Suite)

### 1. Create Environment & Install Dependencies

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r backend/requirements.txt
```

### 2. Run the Offline Test & Provenance Verification Suite

All 79 deterministic unit and real-FITS provenance tests run offline without network access:

```bash
pytest -v
```

To optionally run live IRSA network integration tests (requires internet connectivity to `irsa.ipac.caltech.edu`):

```bash
pytest -v -m live_irsa
```

### 3. Run the Reproducible Benchmark Target Audit

Inspect local cached FITS data across all 5 benchmark targets without claiming unverified availability:

```bash
python -m backend.diagnostics.audit_benchmark_targets --output /tmp/benchmark_audit.json
```

If live IRSA network access is available, run with `--live`:

```bash
python -m backend.diagnostics.audit_benchmark_targets --live --output /tmp/benchmark_audit_live.json
```

### 4. Start the FastAPI Backend Server

```bash
uvicorn backend.main:app --host 127.0.0.1 --port 8000 --reload
```

Interactive OpenAPI documentation is available at `http://127.0.0.1:8000/docs`.

---

## Documentation Index

- **[Data Provenance Audit (`docs/DATA_PROVENANCE_AUDIT.md`)](docs/DATA_PROVENANCE_AUDIT.md)** — Full SHA-256 checksums, observation IDs, `MJD-AVG` vs. product `DATE` timestamps, and upstream IRSA URLs for all 30 cached FITS cutouts.
- **[AI & Fabrication Forensic Audit (`docs/AI_ORIGIN_AUDIT.md`)](docs/AI_ORIGIN_AUDIT.md)** — Detailed accounting of removed fabricated `irsa_audit` numbers, corrected epoch labels, and reclassified heuristic scores.
- **[Scientific Validation Report (`docs/SCIENTIFIC_VALIDATION.md`)](docs/SCIENTIFIC_VALIDATION.md)** — Equations and empirical verification for spherical proper motion, WCS reprojection, background normalization, variance propagation, and difference imaging.
- **[Known Limitations (`docs/KNOWN_LIMITATIONS.md`)](docs/KNOWN_LIMITATIONS.md)** — Transparent disclosure of PSF differences, parallax + proper motion superposition, LVF approximate wavelengths, and offline cache scope.
- **[API Contract (`docs/API_CONTRACT.md`)](docs/API_CONTRACT.md)** — Full specification of all 8 backend HTTP endpoints, status codes, and response fields.
