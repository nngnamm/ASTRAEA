# ASTRAEA Backend Test Suite (`backend/tests/`)

The test suite contains **79 deterministic offline unit and real-FITS provenance tests** plus **1 opt-in live IRSA integration test**.

## Running Tests

From the repository root (`/home/nngnamm/Projects/ASTRAEA`):

```bash
# Run all deterministic offline tests (default):
.venv/bin/pytest -v

# Run opt-in live IRSA network tests:
.venv/bin/pytest -v -m live_irsa
```

## Test Modules

| File | Coverage Area |
| :--- | :--- |
| [`test_coordinates_and_resolution.py`](test_coordinates_and_resolution.py) | Exact spherical proper-motion propagation (`J2026.0` / `MJD 61041.0`), `pm_ra_cosdec` convention, offline benchmark alias resolution, Sesame fallback without double proper motion, and coordinate validation. |
| [`test_irsa_and_cache_integrity.py`](test_irsa_and_cache_integrity.py) | Multi-collection SIA v2 status differentiation (`complete_with_observations`, `complete_empty`, `partial_failure`, `upstream_failure`), refusal to cache incomplete queries, atomic `.tmp` download validation, and corrupted cache eviction + retry. Plus 1 `@pytest.mark.live_irsa` test. |
| [`test_fits_validation_and_metadata.py`](test_fits_validation_and_metadata.py) | 7-HDU structural checks, `BUNIT` enforcement, NaN/Inf/negative variance detection, bitwise-OR union flag counts (`non_astrophysical_flagged_pixels`), operation-specific flag masks (`background_estimation` vs. `detector_defects_only`), and `MJD-AVG` vs. product `DATE` separation. |
| [`test_alignment_normalization_and_epochs.py`](test_alignment_normalization_and_epochs.py) | WCS reprojection of `IMAGE` and `VARIANCE`, `np.nan` preservation outside footprint, empirical cross-frame registration (`evaluate_cross_frame_registration`), ZODI + source-masked background normalization, and detector-pure 30-day epoch grouping. |
| [`test_difference_and_analytics.py`](test_difference_and_analytics.py) | Strict `min_epoch_coverage = 0.50` enforcement (`ineligible_low_coverage`), propagated coadd/difference variance, dipole/brightening/fading candidate classification, `"Benchmark Coordinate Match — Unconfirmed (...)"` labeling, and bitwise-union `compute_sky_health` scoring. |
| [`test_api_endpoints.py`](test_api_endpoints.py) | End-to-end FastAPI `TestClient` tests for all 8 routes (`/api/targets`, `/api/targets/audit`, `/api/search/coordinates`, `/api/search/object`, `/api/frames`, `/api/spectrum`, `/api/candidates`, `/api/cache/status`), HTTP 400/404/503 error contracts, and transparent-NaN RGBA PNG rendering. |
| [`test_real_cached_fits_provenance.py`](test_real_cached_fits_provenance.py) | Verifies all 30 real cached Barnard's Star FITS cutouts against [`audit/data_provenance.json`](../../audit/data_provenance.json) (SHA-256 checksums, `19 spherex_dr1` + `11 spherex_qr3`, `16` Detector 5 + `14` Detector 6), empirical registration RMS (`0.1268 px`), and regenerated [`sample_barnard_response.json`](../sample_barnard_response.json). |
| [`test_backend_scientific_suite.py`](test_backend_scientific_suite.py) | Comprehensive integration checks across the full Barnard's Star cutout cube, spectrum extraction, and benchmark audit contract. |
| [`test_propagate_coords.py`](test_propagate_coords.py) | Focused regression tests for spherical coordinate propagation. |
