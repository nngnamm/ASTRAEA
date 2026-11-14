# ASTRAEA Backend (`Team 02 — Backend + Data/Astronomy`)

This directory contains the FastAPI backend, SPHEREx Level-2 FITS validation pipeline, WCS gnomonic (`TAN`) reprojection engine, variance-propagated difference imaging pipeline, and automated scientific test suite for **ASTRAEA — COSMIC PULSE**.

---

## Architecture Overview

| Module | Responsibility |
| :--- | :--- |
| [`main.py`](main.py) | FastAPI application, CORS configuration, parameter validation, HTTP error mapping (`400`, `404`, `500`, `503`), and 8 REST endpoints. |
| [`config.py`](config.py) | Explicit epoch constants (`DEFAULT_REFERENCE_MJD = 61041.0` / `J2026.0`), detector wavelength bands (`BANDS_UM`), quality/coverage thresholds, and the provenance-backed `BENCHMARK_TARGETS` catalog. |
| [`spherex_core.py`](spherex_core.py) | Multi-collection IRSA SIA v2 queries (`find_images`), atomic FITS download & cache validation (`_load_frame_with_validation`), exact spherical proper-motion propagation (`propagate_coords_with_metadata`), WCS reprojection of `IMAGE` and `VARIANCE` (`align_frame_to_grid`), empirical source registration (`evaluate_cross_frame_registration`), flag-masked ZODI+sky normalization (`normalize_frame_background`), variance-propagated difference imaging (`compute_difference_with_diagnostics`), and reproducible benchmark auditing (`audit_benchmark_target`). |
| [`image_validation.py`](image_validation.py) | 7-HDU structural validator (`validate_spherex_fits`), bitwise flag decoding (`_decode_flag_header`), operation-specific pixel masks (`build_operation_flag_mask`), and four-category validation reporting (`structural_validity`, `numerical_validity`, `alignment_validity`, `scientific_operation_suitability`). |
| [`observation_metadata.py`](observation_metadata.py) | Extracts standardized schema `1.1` FITS metadata (`OBSID`, `MJD-BEG`, `MJD-AVG`, `MJD-END`, product `DATE`, `DETECTOR`, `BUNIT`, `SHA-256`, `access_url`) and enforces chronological ordering by `MJD-AVG`. |
| [`serializers.py`](serializers.py) | Detector-pure 30-day visit grouping (`group_frames_into_epochs`), strict coverage-gated epoch selection (`_difference_for_detector`), transparent-NaN RGBA PNG encoding (`array_to_b64png`), and JSON sanitization (`sanitize_for_json`). |
| [`analytics.py`](analytics.py) | Experimental heuristic candidate detection (`detect_anomalies`) with epoch-matched benchmark coordinate association and bitwise-union `Sky Health` scoring (`compute_sky_health`). |

---

## Setup & Installation

From the repository root (`/home/nngnamm/Projects/ASTRAEA`):

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r backend/requirements.txt
```

---

## Running the Backend Server

From the repository root:

```bash
uvicorn backend.main:app --host 127.0.0.1 --port 8000 --reload
```

Or from inside `backend/`:

```bash
uvicorn main:app --host 127.0.0.1 --port 8000 --reload
```

Open `http://127.0.0.1:8000/docs` to inspect the interactive OpenAPI schema.

---

## HTTP Endpoints

See [`docs/API_CONTRACT.md`](../docs/API_CONTRACT.md) for full request/response schemas:

- `GET /api/targets` — Returns the 5 NASA benchmark targets with exact `J2026.0` (`MJD 61041.0`) coordinates and honest provenance-backed availability status (`verified_local_cache_subset` for `barnard`, `unverified_offline` for the remaining 4 targets).
- `GET /api/targets/audit` — Runs a reproducible audit of all 5 benchmark targets against the local FITS cache (or live IRSA when `query_irsa_if_uncached=true`).
- `GET /api/search/coordinates` — Validates and normalizes sky coordinates (`irsa_queried: false` explicitly documented in response).
- `GET /api/search/object` (alias `GET /api/resolve`) — Resolves benchmark target aliases offline with exact spherical proper-motion propagation, or queries Astropy Sesame (`allow_external=true`).
- `GET /api/frames` — Builds the validated, aligned, normalized multi-epoch SPHEREx cutout cube, per-detector coadds, variance-propagated difference images, and candidate diagnostics. Returns HTTP `503` if upstream IRSA SIA fails when uncached.
- `GET /api/spectrum` — Computes per-detector aperture surface-brightness excess measurements (`calibration_grade: false`, `wavelength_source: "WCS-WAVE_approximate"`).
- `GET /api/candidates` — Lightweight endpoint returning per-detector candidate detections and `sky_health` without base64 PNG payloads.
- `GET /api/cache/status` — Reports local FITS cutout and SIA query cache counts and byte sizes (`cache_dir` sanitized to basename `"cache"`).

---

## Running Tests & Diagnostic Scripts

Run all 79 deterministic unit and real-FITS provenance tests from the repository root:

```bash
pytest -v
```

Run the reproducible benchmark audit CLI:

```bash
python -m backend.diagnostics.audit_benchmark_targets
```

Regenerate [`sample_barnard_response.json`](sample_barnard_response.json) directly from the 30 verified FITS cutouts in `backend/cache/`:

```bash
python -m backend.diagnostics.regenerate_barnard_sample
```
