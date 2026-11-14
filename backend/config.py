import os
from pathlib import Path

from approved_data_sources import (
    GAIA_DR3_BENCHMARK_ASTROMETRY,
    OFFICIAL_SPHEREX_COLLECTIONS,
)

COLLECTIONS = ["spherex_dr1", "spherex_qr3"]
ALL_OFFICIAL_SPHEREX_COLLECTIONS = list(OFFICIAL_SPHEREX_COLLECTIONS)

_DEFAULT_CACHE_DIR = str(Path(__file__).resolve().parent / "cache")
CACHE_DIR = os.environ.get("ASTRAEA_CACHE_DIR", _DEFAULT_CACHE_DIR)
os.makedirs(CACHE_DIR, exist_ok=True)

# Coordinate propagation defaults:
# J2026.0 in Julian years corresponds exactly to MJD 61041.0 (2026-01-01T00:00:00 UTC).
# Previously an undocumented MJD 61000.0 (2025-11-21, J2025.8877) was used while
# labeling outputs as ra_2026 / dec_2026.
DEFAULT_REFERENCE_EPOCH_JYEAR = 2026.0
DEFAULT_REFERENCE_MJD = 61041.0
LEGACY_CUTOUT_QUERY_MJD = 61000.0

# Scientific comparison eligibility thresholds
MIN_EPOCH_COVERAGE = 0.50
MIN_COMMON_COVERAGE = 0.10

# Network and API configuration
HTTP_TIMEOUT_SECONDS = float(os.environ.get("ASTRAEA_HTTP_TIMEOUT", "60.0"))
_raw_cors_origins = os.environ.get(
    "ASTRAEA_CORS_ORIGINS",
    "http://localhost:5173,http://127.0.0.1:5173,http://localhost:3000,http://127.0.0.1:3000",
)
CORS_ALLOWED_ORIGINS = [
    origin.strip() for origin in _raw_cors_origins.split(",") if origin.strip()
]

# Approximate wavelength passbands (um) per SPHEREx detector (Linear Variable Filter channels)
BANDS = {
    1: (0.75, 1.09),
    2: (1.10, 1.62),
    3: (1.63, 2.41),
    4: (2.42, 3.82),
    5: (3.83, 4.41),
    6: (4.42, 5.00),
}

# NASA Space Apps 2026 Benchmark Targets.
# Provenance note:
# - Only "barnard" has cryptographically verified SPHEREx Level-2 FITS cutouts cached
#   locally in backend/cache/ (30 files across Detectors 5 and 6).
# - Static full-mission observation counts that previously appeared under "irsa_audit"
#   were removed after provenance auditing proved they contained unverified/fabricated
#   per-detector breakdowns. Use run_benchmark_audit() or
#   backend/diagnostics/audit_benchmark_targets.py for live IRSA verification.
TARGETS = {
    "barnard": {
        "name": "Barnard's Star",
        "category": "High Proper Motion Star",
        "test_case": "Proper Motion",
        "ra": 269.452077,
        "dec": 4.693365,
        "epoch": 2000.0,
        "pm_ra_arcsec_yr": -0.8015,
        "pm_dec_arcsec_yr": 10.3624,
        "pm_convention": "pm_ra_cosdec",
        "coordinate_source": "SIMBAD / Hipparcos / Gaia DR3 ICRS propagated to J2000.0 (GJ 699 / HIP 87937)",
        "proper_motion_source": "Gaia DR3 (pmra=-801.551 mas/yr, pmdec=+10362.394 mas/yr)",
        "aliases": [
            "barnard",
            "barnard's star",
            "barnards star",
            "gj 699",
            "hip 87937",
        ],
        "data_status": "verified_local_cache_subset",
        "recommended_detectors": [5, 6],
        "cached_cutout_query_center": {
            "ra_deg": 269.446294,
            "dec_deg": 4.767881,
            "propagated_to_mjd": 61000.0,
            "propagated_epoch_jyear": 2025.887748,
            "fov_deg": 0.1,
        },
        "provenance_audit": {
            "verification_status": "verified_local_cache_subset",
            "cached_fits_count": 30,
            "complete_sia_catalog_cached": False,
            "collections_verified_in_cache": {"spherex_dr1": 19, "spherex_qr3": 11},
            "by_detector_in_cache": {5: 16, 6: 14},
            "mjd_min": 60902.865381,
            "mjd_max": 61305.628658,
            "baseline_days": 402.76,
            "visit_groups_30d_in_cache": {"D5": 3, "D6": 4},
            "note": (
                "30 NASA/IPAC IRSA SPHEREx Level-2 FITS cutouts (Detectors 5 & 6) are cached "
                "and SHA-256 verified in backend/cache/. Detectors 1-4 and full SIA row counts "
                "require a live IRSA query."
            ),
        },
    },
    "eps_ind_b": {
        "name": "Epsilon Indi B",
        "category": "Binary Brown Dwarf System",
        "test_case": "Brown Dwarf Binary System",
        "ra": 331.044147,
        "dec": -56.782844,
        "epoch": 2000.0,
        "pm_ra_arcsec_yr": 3.9614,
        "pm_dec_arcsec_yr": -2.5383,
        "pm_convention": "pm_ra_cosdec",
        "coordinate_source": "SIMBAD / Scholz et al. (2003) ICRS J2000.0 (HD 209100 B)",
        "proper_motion_source": "Scholz et al. (2003, A&A, 398, L29) / Gaia DR3",
        "aliases": [
            "eps_ind_b",
            "epsilon indi b",
            "eps ind b",
            "hd 209100 b",
        ],
        "data_status": "unverified_offline",
        "recommended_detectors": [5, 6],
        "provenance_audit": {
            "verification_status": "unverified_offline",
            "cached_fits_count": 0,
            "complete_sia_catalog_cached": False,
            "note": (
                "No FITS cutouts or SIA query results are cached in this repository. "
                "Run backend/diagnostics/audit_benchmark_targets.py with live IRSA access to verify."
            ),
        },
    },
    "wise0855": {
        "name": "WISE 0855-0714",
        "category": "Fast-Moving Y-type Brown Dwarf (~250 K)",
        "test_case": "Fast-moving Brown Dwarf",
        "ra": 133.795129,
        "dec": -7.245147,
        "epoch": 2010.5,
        "pm_ra_arcsec_yr": -8.118,
        "pm_dec_arcsec_yr": 0.680,
        "pm_convention": "pm_ra_cosdec",
        "coordinate_source": "Luhman (2014, ApJ, 786, L18) WISE discovery epoch J2010.5 (08h55m10.831s -07d14m42.53s)",
        "proper_motion_source": "Luhman & Esplin (2016, AJ, 152, 78) / Kirkpatrick et al. (2021)",
        "aliases": [
            "wise0855",
            "wise 0855",
            "wise 0855-0714",
            "wise j085510.83-071442.5",
            "wise j0855-0714",
        ],
        "data_status": "unverified_offline",
        "recommended_detectors": [5, 6],
        "provenance_audit": {
            "verification_status": "unverified_offline",
            "cached_fits_count": 0,
            "complete_sia_catalog_cached": False,
            "historical_notebook_reference": {
                "source": "notebooks/NASA_backend_exploration.ipynb (Cells 8, 31, 36)",
                "collection_queried": "spherex_dr1",
                "total_sia_rows_reported": 263,
                "d5_d6_frames_reported": 94,
                "baseline_days_reported": 373.6,
                "verified_in_local_cache": False,
            },
            "note": (
                "Exploratory notebook recorded 263 spherex_dr1 SIA rows during an earlier session, "
                "but no FITS cutouts or SIA JSON files for WISE 0855-0714 are cached in backend/cache/."
            ),
        },
    },
    "rw_aur": {
        "name": "RW Aurigae",
        "category": "T Tauri Variable Star",
        "test_case": "Variable Star",
        "ra": 76.956526,
        "dec": 30.401438,
        "epoch": 2000.0,
        "pm_ra_arcsec_yr": 0.003,
        "pm_dec_arcsec_yr": -0.024,
        "pm_convention": "pm_ra_cosdec",
        "coordinate_source": "SIMBAD / Gaia DR3 ICRS J2000.0 (HD 240764)",
        "proper_motion_source": "Gaia DR3 (pmra=+3.34 mas/yr, pmdec=-24.16 mas/yr)",
        "aliases": [
            "rw_aur",
            "rw aurigae",
            "rw aur",
            "hd 240764",
        ],
        "data_status": "unverified_offline",
        "recommended_detectors": [1, 2, 3, 4, 5, 6],
        "provenance_audit": {
            "verification_status": "unverified_offline",
            "cached_fits_count": 0,
            "complete_sia_catalog_cached": False,
            "note": (
                "No FITS cutouts or SIA query results are cached in this repository. "
                "Previous static uniform counts ({1: 57, ..., 6: 57}) were removed as unverified."
            ),
        },
    },
    "r_dra": {
        "name": "R Draconis",
        "category": "Mira Cet Variable Star",
        "test_case": "Variable Star",
        "ra": 248.167611,
        "dec": 66.754979,
        "epoch": 2000.0,
        "pm_ra_arcsec_yr": -0.006,
        "pm_dec_arcsec_yr": -0.033,
        "pm_convention": "pm_ra_cosdec",
        "coordinate_source": "SIMBAD / Gaia DR3 ICRS J2000.0 (HD 150072)",
        "proper_motion_source": "Gaia DR3 (pmra=-5.92 mas/yr, pmdec=-33.14 mas/yr)",
        "aliases": [
            "r_dra",
            "r draconis",
            "r dra",
            "hd 150072",
        ],
        "data_status": "unverified_offline",
        "recommended_detectors": [1, 2, 3, 4, 5, 6],
        "provenance_audit": {
            "verification_status": "unverified_offline",
            "cached_fits_count": 0,
            "complete_sia_catalog_cached": False,
            "note": (
                "No FITS cutouts or SIA query results are cached in this repository. "
                "Run backend/diagnostics/audit_benchmark_targets.py with live IRSA access to verify."
            ),
        },
    },
}

for _target_key, _target_cfg in TARGETS.items():
    if _target_key in GAIA_DR3_BENCHMARK_ASTROMETRY:
        _target_cfg["astrometric_reference"] = dict(
            GAIA_DR3_BENCHMARK_ASTROMETRY[_target_key]
        )
        _target_cfg["parallax_mas"] = GAIA_DR3_BENCHMARK_ASTROMETRY[_target_key].get(
            "parallax_mas"
        )
        _target_cfg["radial_velocity_km_s"] = GAIA_DR3_BENCHMARK_ASTROMETRY[
            _target_key
        ].get("radial_velocity_km_s")
        _target_cfg["gaia_dr3_source_id"] = GAIA_DR3_BENCHMARK_ASTROMETRY[
            _target_key
        ].get("source_id")
