"""Approved Scientific Data Sources & Dataset Separation Registry for ASTRAEA.

Enforces strict dataset separation, collection provenance verification,
Gaia DR3 astrometric reference metadata (including parallax, radial velocity,
and independent reference fallback for targets absent from Gaia DR3), and an
instrument-specific adapter for supplementary CSA/CADC NEOSSat observations.

Approved Source Ecosystem:
1. Primary Data Source: NASA/IPAC SPHEREx
   - Archive: https://irsa.ipac.caltech.edu/Missions/spherex.html
   - Tutorials: https://caltech-ipac.github.io/irsa-tutorials/spherex/
   - Official Collections: spherex_dr1, spherex_dr1_deep, spherex_qr3, spherex_qr3_deep
2. Astrometric Reference: Gaia DR3
   - Reference: https://archive.stsci.edu/missions-and-data/gaia
3. Supplementary Mission: NEOSSat (CSA)
   - Dataset: https://data.asc-csa.gc.ca/dataset/neossat
   - Tutorial: https://github.com/asc-csa/NEOSSAT_Tutorial
4. Supplementary Archive: Canadian Astronomical Data Centre (CADC)
   - Reference: https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/
5. Development & Visualization Workflow Inspiration: WiseView
   - Reference: http://byw.tools/wiseview (never used as a substitute for SPHEREx data)
"""
from __future__ import annotations

import logging
import os
import re
import warnings
from typing import Any, Mapping, Sequence

import astropy.units as u
import numpy as np
from astropy.coordinates import Distance, SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

try:
    from erfa import ErfaWarning
except ImportError:
    ErfaWarning = Warning

logger = logging.getLogger(__name__)

OFFICIAL_SPHEREX_COLLECTIONS: tuple[str, ...] = (
    "spherex_dr1",
    "spherex_dr1_deep",
    "spherex_qr3",
    "spherex_qr3_deep",
)

DATASET_IDENTITIES: dict[str, dict[str, Any]] = {
    "SPHEREx": {
        "dataset_identity": "SPHEREx",
        "role": "primary_observational_dataset",
        "provider": "NASA/IPAC Infrared Science Archive (IRSA)",
        "official_archive_url": "https://irsa.ipac.caltech.edu/Missions/spherex.html",
        "official_tutorials_url": "https://caltech-ipac.github.io/irsa-tutorials/spherex/",
        "official_collections": list(OFFICIAL_SPHEREX_COLLECTIONS),
        "wavelength_regime": "Near- to Mid-Infrared (0.75 - 5.00 um, 6 LVF detectors)",
        "native_surface_brightness_unit": "MJy / sr",
        "nominal_pixel_scale_arcsec": 6.155,
        "is_real_observation": True,
        "allowed_in_production_responses": True,
    },
    "Gaia": {
        "dataset_identity": "Gaia",
        "role": "astrometric_reference_catalog",
        "provider": "ESA Gaia Data Release 3 / STScI / CDS",
        "official_archive_url": "https://archive.stsci.edu/missions-and-data/gaia",
        "wavelength_regime": "Optical astrometry (G band, 330 - 1050 nm)",
        "is_spherex_measurement": False,
        "is_real_observation": True,
        "allowed_in_production_responses": True,
        "note": (
            "Used strictly as an astrometric reference catalog (source_id, RA/Dec, "
            "reference epoch, proper motion, parallax, radial velocity, RUWE). "
            "Gaia measurements are never claimed to be SPHEREx measurements."
        ),
    },
    "NEOSSat": {
        "dataset_identity": "NEOSSat",
        "role": "supplementary_observational_mission",
        "provider": "Canadian Space Agency (CSA) / Canadian Astronomical Data Centre (CADC)",
        "official_dataset_url": "https://data.asc-csa.gc.ca/dataset/neossat",
        "official_tutorials_url": "https://github.com/asc-csa/NEOSSAT_Tutorial",
        "cadc_archive_url": "https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/",
        "cadc_collection": "NEOSSAT",
        "wavelength_regime": "Optical broadband CCD (350 - 1050 nm, unfiltered)",
        "native_unit_examples": ["ADU", "count", "e-/s"],
        "nominal_pixel_scale_arcsec": 3.0,
        "is_spherex_measurement": False,
        "direct_pixel_combination_with_spherex_allowed": False,
        "is_real_observation": True,
        "allowed_in_production_responses": True,
        "citation_requirement": (
            "Contains information licensed under the Open Government Licence - Canada "
            "(Canadian Space Agency NEOSSat) and distributed via CADC."
        ),
    },
    "WISE": {
        "dataset_identity": "WISE",
        "role": "development_and_visualization_workflow_inspiration",
        "provider": "WiseView (byw.tools) / NASA WISE & NEOWISE",
        "reference_url": "http://byw.tools/wiseview",
        "is_spherex_measurement": False,
        "is_substitute_for_spherex": False,
        "is_real_observation": False,
        "allowed_in_production_responses": False,
        "note": (
            "WiseView is a reference for multi-epoch blinking and motion visualization "
            "workflows, never a source of SPHEREx observations or scientific validation."
        ),
    },
    "Synthetic": {
        "dataset_identity": "Synthetic",
        "role": "deterministic_unit_test_fixtures_only",
        "provider": "ASTRAEA automated unit test suite",
        "is_spherex_measurement": False,
        "is_real_observation": False,
        "allowed_in_production_responses": False,
        "note": (
            "Artificially generated arrays and FITS structures used solely in offline "
            "unit tests. Never substituted for unavailable real observations."
        ),
    },
}


# Authoritative Astrometric Reference Table for the 5 NASA Benchmark Targets.
# Preserves Gaia DR3 source_id, reference epoch, proper motion (pm_ra_cosdec, pm_dec),
# parallax, radial velocity, and astrometric quality indicators where available.
# For targets absent from Gaia DR3 (e.g., WISE 0855-0714), documents the independent
# peer-reviewed infrared astrometric reference and confidence explicitly.
GAIA_DR3_BENCHMARK_ASTROMETRY: dict[str, dict[str, Any]] = {
    "barnard": {
        "dataset_identity": "Gaia",
        "is_spherex_measurement": False,
        "gaia_dr3_present": True,
        "catalog_name": "Gaia DR3",
        "reference_url": "https://archive.stsci.edu/missions-and-data/gaia",
        "source_id": "4472832130942575872",
        "designation": "Gaia DR3 4472832130942575872",
        "ra_deg": 269.4485025254,
        "dec_deg": 4.7394200511,
        "reference_epoch_jyear": 2016.0,
        "reference_epoch_mjd": 57388.5,
        "pm_ra_cosdec_mas_yr": -801.551,
        "pm_dec_mas_yr": 10362.394,
        "pm_convention": "pm_ra_cosdec",
        "parallax_mas": 546.9759,
        "parallax_error_mas": 0.0401,
        "radial_velocity_km_s": -110.51,
        "radial_velocity_error_km_s": 0.12,
        "astrometric_quality": {
            "ra_error_mas": 0.0318,
            "dec_error_mas": 0.0268,
            "pmra_error_mas_yr": 0.0318,
            "pmdec_error_mas_yr": 0.0247,
            "ruwe": 1.08,
            "astrometric_excess_noise_mas": 0.21,
        },
        "reference_confidence": "authoritative_gaia_dr3",
        "citation": "Gaia Collaboration et al. (2023, A&A, 674, A1) — Gaia Data Release 3",
    },
    "eps_ind_b": {
        "dataset_identity": "Gaia",
        "is_spherex_measurement": False,
        "gaia_dr3_present": False,
        "catalog_name": "Scholz et al. (2003) / Chen et al. (2022) + Gaia DR3 Comoving Host (eps Ind A)",
        "reference_url": "https://archive.stsci.edu/missions-and-data/gaia",
        "source_id": None,
        "comoving_primary_gaia_dr3_source_id": "6412595290592307840",
        "designation": "HD 209100 B (eps Ind Ba + Bb)",
        "absence_reason": (
            "Epsilon Indi B is a close T1+T6 brown dwarf binary (Ba/Bb) resolved in "
            "near-IR astrometry; the comoving primary host Epsilon Indi A is Gaia DR3 "
            "6412595290592307840."
        ),
        "ra_deg": 331.044147,
        "dec_deg": -56.782844,
        "reference_epoch_jyear": 2000.0,
        "reference_epoch_mjd": 51544.5,
        "pm_ra_cosdec_mas_yr": 3961.4,
        "pm_dec_mas_yr": -2538.3,
        "pm_convention": "pm_ra_cosdec",
        "parallax_mas": 274.84,
        "parallax_error_mas": 0.10,
        "radial_velocity_km_s": -40.4,
        "radial_velocity_error_km_s": 0.5,
        "astrometric_quality": {
            "host_star_gaia_dr3_source_id": "6412595290592307840",
            "host_star_parallax_mas": 274.8431,
            "host_star_parallax_error_mas": 0.0956,
            "host_star_ruwe": 0.99,
        },
        "independent_reference_source": (
            "Scholz et al. (2003, A&A, 398, L29); Dieterich et al. (2018, ApJ, 865, 28); "
            "Chen et al. (2022, ApJ, 934, 178); Gaia DR3 6412595290592307840 (host eps Ind A)"
        ),
        "reference_confidence": "high_peer_reviewed_ir_astrometry_plus_gaia_dr3_host",
        "citation": "Scholz et al. (2003, A&A, 398, L29) & Gaia DR3 (host star 6412595290592307840)",
    },
    "wise0855": {
        "dataset_identity": "Gaia",
        "is_spherex_measurement": False,
        "gaia_dr3_present": False,
        "catalog_name": "Independent Mid-IR Astrometry (WISE + Spitzer + HST)",
        "reference_url": "https://ui.adsabs.harvard.edu/abs/2021ApJS..253....7K",
        "source_id": None,
        "designation": "WISE J085510.83-071442.5",
        "absence_reason": (
            "Substellar Y4 brown dwarf (~250 K) with optical magnitude G > 25 mag, "
            "invisible in the Gaia optical passband (330 - 1050 nm, limit G ~ 21)."
        ),
        "ra_deg": 133.795129,
        "dec_deg": -7.245147,
        "reference_epoch_jyear": 2010.5,
        "reference_epoch_mjd": 55378.625,
        "pm_ra_cosdec_mas_yr": -8118.0,
        "pm_dec_mas_yr": 680.0,
        "pm_convention": "pm_ra_cosdec",
        "parallax_mas": 439.0,
        "parallax_error_mas": 2.4,
        "radial_velocity_km_s": None,
        "radial_velocity_error_km_s": None,
        "astrometric_quality": {
            "pmra_error_mas_yr": 8.0,
            "pmdec_error_mas_yr": 8.0,
            "parallax_error_mas": 2.4,
            "astrometric_facilities": ["WISE/NEOWISE", "Spitzer/IRAC", "HST"],
        },
        "independent_reference_source": (
            "Luhman (2014, ApJ, 786, L18); Luhman & Esplin (2016, AJ, 152, 78); "
            "Kirkpatrick et al. (2021, ApJS, 253, 7)"
        ),
        "reference_confidence": "high_peer_reviewed_ir_astrometry",
        "citation": "Luhman & Esplin (2016, AJ, 152, 78) & Kirkpatrick et al. (2021, ApJS, 253, 7)",
    },
    "rw_aur": {
        "dataset_identity": "Gaia",
        "is_spherex_measurement": False,
        "gaia_dr3_present": True,
        "catalog_name": "Gaia DR3",
        "reference_url": "https://archive.stsci.edu/missions-and-data/gaia",
        "source_id": "156430822114424576",
        "designation": "Gaia DR3 156430822114424576",
        "ra_deg": 76.9565399,
        "dec_deg": 30.4013895,
        "reference_epoch_jyear": 2016.0,
        "reference_epoch_mjd": 57388.5,
        "pm_ra_cosdec_mas_yr": 3.34,
        "pm_dec_mas_yr": -24.16,
        "pm_convention": "pm_ra_cosdec",
        "parallax_mas": 6.115,
        "parallax_error_mas": 0.068,
        "radial_velocity_km_s": 15.87,
        "radial_velocity_error_km_s": 1.45,
        "astrometric_quality": {
            "ra_error_mas": 0.062,
            "dec_error_mas": 0.044,
            "pmra_error_mas_yr": 0.078,
            "pmdec_error_mas_yr": 0.056,
            "ruwe": 1.42,
            "astrometric_excess_noise_mas": 0.34,
        },
        "reference_confidence": "authoritative_gaia_dr3",
        "citation": "Gaia Collaboration et al. (2023, A&A, 674, A1) — Gaia Data Release 3",
    },
    "r_dra": {
        "dataset_identity": "Gaia",
        "is_spherex_measurement": False,
        "gaia_dr3_present": True,
        "catalog_name": "Gaia DR3",
        "reference_url": "https://archive.stsci.edu/missions-and-data/gaia",
        "source_id": "1646028373578916608",
        "designation": "Gaia DR3 1646028373578916608",
        "ra_deg": 248.1675509,
        "dec_deg": 66.7548317,
        "reference_epoch_jyear": 2016.0,
        "reference_epoch_mjd": 57388.5,
        "pm_ra_cosdec_mas_yr": -5.92,
        "pm_dec_mas_yr": -33.14,
        "pm_convention": "pm_ra_cosdec",
        "parallax_mas": 2.302,
        "parallax_error_mas": 0.056,
        "radial_velocity_km_s": -133.0,
        "radial_velocity_error_km_s": 0.8,
        "astrometric_quality": {
            "ra_error_mas": 0.049,
            "dec_error_mas": 0.053,
            "pmra_error_mas_yr": 0.061,
            "pmdec_error_mas_yr": 0.066,
            "ruwe": 1.15,
            "astrometric_excess_noise_mas": 0.25,
        },
        "reference_confidence": "authoritative_gaia_dr3",
        "citation": "Gaia Collaboration et al. (2023, A&A, 674, A1) — Gaia Data Release 3",
    },
}


def get_astrometric_reference(target_key: str) -> dict[str, Any] | None:
    """Return a copy of the astrometric reference record for a benchmark target."""
    ref = GAIA_DR3_BENCHMARK_ASTROMETRY.get(str(target_key).strip().lower())
    return dict(ref) if ref is not None else None


def infer_spherex_collection_from_url(
    access_url: str | None, fallback_collection: str | None = None
) -> str | None:
    """Verify SPHEREx collection provenance directly from an IRSA access URL or URI.

    Recognizes official collections (`spherex_dr1_deep`, `spherex_dr1`,
    `spherex_qr3_deep`, `spherex_qr3`) so that a cross-collection SIA response
    never mislabels an observation's actual collection.
    """
    if access_url:
        text = str(access_url).lower()
        for col_token, canonical in (
            ("spherex_dr1_deep", "spherex_dr1_deep"),
            ("/spherex/dr1_deep/", "spherex_dr1_deep"),
            ("spherex_dr1", "spherex_dr1"),
            ("/spherex/dr1/", "spherex_dr1"),
            ("spherex_qr3_deep", "spherex_qr3_deep"),
            ("/spherex/qr3_deep/", "spherex_qr3_deep"),
            ("spherex_qr3", "spherex_qr3"),
            ("/spherex/qr3/", "spherex_qr3"),
            ("spherex_qr2", "spherex_qr2"),
            ("/spherex/qr2/", "spherex_qr2"),
        ):
            if col_token in text:
                return canonical
    return fallback_collection


def discover_spherex_collections(
    coord: SkyCoord | None = None,
    fov_deg: float = 0.1,
    candidate_collections: Sequence[str] = OFFICIAL_SPHEREX_COLLECTIONS,
    *,
    allow_network: bool = False,
    query_fn: Any = None,
) -> dict[str, Any]:
    """Dynamically discover and probe official SPHEREx collections at IRSA.

    Never assumes that all collections contain observations for every target.
    Records which collections were queried, per-collection status, and observation
    counts. When ``allow_network=False``, reports offline status cleanly without
    inventing query results.
    """
    candidates = [str(c) for c in candidate_collections]
    if not allow_network:
        return {
            "dataset_identity": "SPHEREx",
            "discovery_mode": "offline_configuration_only",
            "network_queried": False,
            "official_collections_investigated": candidates,
            "collections_queried": [],
            "collections_with_observations": [],
            "collection_status": {
                col: {
                    "queried": False,
                    "status": "unverified_offline",
                    "observation_count": None,
                    "error": None,
                }
                for col in candidates
            },
            "note": (
                "Live IRSA collection probing disabled (allow_network=False). "
                "Enable allow_network=True to dynamically query spherex_dr1, "
                "spherex_dr1_deep, spherex_qr3, and spherex_qr3_deep."
            ),
        }

    if coord is None:
        raise ValueError("coord (SkyCoord) is required when allow_network=True")

    if query_fn is None:
        from astroquery.ipac.irsa import Irsa

        query_fn = Irsa.query_sia

    collection_status: dict[str, dict[str, Any]] = {}
    collections_with_obs: list[str] = []
    n_succeeded = 0
    n_failed = 0

    for col in candidates:
        try:
            res = query_fn(pos=(coord, (float(fov_deg) / 2.0) * u.deg), collection=col)
            table = res.to_table() if hasattr(res, "to_table") else res
            count = len(table)
            n_succeeded += 1
            if count > 0:
                collections_with_obs.append(col)
            collection_status[col] = {
                "queried": True,
                "status": "ok_with_observations" if count > 0 else "ok_empty",
                "observation_count": int(count),
                "error": None,
            }
        except Exception as exc:
            n_failed += 1
            err_msg = f"{type(exc).__name__}: {exc}"
            logger.warning("SPHEREx collection discovery failed for %s: %s", col, err_msg)
            collection_status[col] = {
                "queried": True,
                "status": "error",
                "observation_count": 0,
                "error": err_msg,
            }

    return {
        "dataset_identity": "SPHEREx",
        "discovery_mode": "live_irsa_sia_probe",
        "network_queried": True,
        "official_collections_investigated": candidates,
        "collections_queried": candidates,
        "collections_with_observations": collections_with_obs,
        "complete": bool(len(candidates) > 0 and n_failed == 0),
        "collections_succeeded_count": n_succeeded,
        "collections_failed_count": n_failed,
        "collection_status": collection_status,
    }


def propagate_astrometry_with_parallax(
    ra_deg: float,
    dec_deg: float,
    *,
    pm_ra_cosdec_mas_yr: float = 0.0,
    pm_dec_mas_yr: float = 0.0,
    reference_epoch_jyear: float = 2016.0,
    to_mjd: float = 61041.0,
    parallax_mas: float | None = None,
    radial_velocity_km_s: float | None = None,
) -> dict[str, Any]:
    """Propagate astrometric reference coordinates accounting for proper motion,
    parallax distance, and radial velocity where available.

    Returns both the space-motion propagated barycentric ICRS position and the
    annual trigonometric parallax ellipse amplitude (`parallax_amplitude_arcsec`)
    so multi-epoch difference analysis explicitly accounts for annual parallax.
    """
    warnings.filterwarnings("ignore", category=ErfaWarning)
    t0 = Time(float(reference_epoch_jyear), format="jyear")
    t1 = Time(float(to_mjd), format="mjd")

    kwargs: dict[str, Any] = {
        "ra": float(ra_deg) * u.deg,
        "dec": float(dec_deg) * u.deg,
        "pm_ra_cosdec": float(pm_ra_cosdec_mas_yr) * u.mas / u.yr,
        "pm_dec": float(pm_dec_mas_yr) * u.mas / u.yr,
        "obstime": t0,
        "frame": "icrs",
    }
    has_parallax = (
        parallax_mas is not None
        and np.isfinite(float(parallax_mas))
        and float(parallax_mas) > 0.0
    )
    if has_parallax:
        kwargs["distance"] = Distance(parallax=float(parallax_mas) * u.mas)
    if radial_velocity_km_s is not None and np.isfinite(float(radial_velocity_km_s)):
        kwargs["radial_velocity"] = float(radial_velocity_km_s) * u.km / u.s

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        c0 = SkyCoord(**kwargs)
        c1 = c0.apply_space_motion(new_obstime=t1)

    parallax_arcsec = float(parallax_mas) / 1000.0 if has_parallax else None
    return {
        "dataset_identity": "Gaia",
        "is_spherex_measurement": False,
        "ra_deg": round(float(c1.ra.deg), 6),
        "dec_deg": round(float(c1.dec.deg), 6),
        "catalog_ra_deg": float(ra_deg),
        "catalog_dec_deg": float(dec_deg),
        "reference_epoch_jyear": float(reference_epoch_jyear),
        "propagated_to_mjd": round(float(t1.mjd), 6),
        "propagated_epoch_jyear": round(float(t1.jyear), 6),
        "pm_ra_cosdec_mas_yr": float(pm_ra_cosdec_mas_yr),
        "pm_dec_mas_yr": float(pm_dec_mas_yr),
        "pm_convention": "pm_ra_cosdec",
        "parallax_mas": float(parallax_mas) if has_parallax else None,
        "annual_parallax_amplitude_arcsec": (
            round(parallax_arcsec, 4) if parallax_arcsec is not None else None
        ),
        "max_peak_to_peak_parallax_shift_arcsec": (
            round(2.0 * parallax_arcsec, 4) if parallax_arcsec is not None else None
        ),
        "radial_velocity_km_s": (
            float(radial_velocity_km_s)
            if radial_velocity_km_s is not None and np.isfinite(float(radial_velocity_km_s))
            else None
        ),
        "parallax_accounted_in_space_motion": bool(has_parallax),
    }


def validate_dataset_separation(items: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Ensure observations from different missions or synthetic fixtures are never
    silently merged into a single observational image cube or coadd.

    Raises ``ValueError`` if multiple incompatible ``dataset_identity`` values are
    present or if non-observational datasets (e.g., ``Gaia``, ``WISE``, ``Synthetic``)
    are mixed with primary ``SPHEREx`` observations.
    """
    identities: set[str] = set()
    for item in items:
        ident = (
            item.get("dataset_identity")
            or item.get("metadata", {}).get("dataset_identity")
            or "SPHEREx"
        )
        identities.add(str(ident))

    if len(identities) > 1:
        raise ValueError(
            f"Dataset separation violation: cannot merge multiple dataset identities "
            f"{sorted(identities)} into a single observational collection."
        )
    if identities and next(iter(identities)) not in {"SPHEREx", "NEOSSat"}:
        only_ident = next(iter(identities))
        raise ValueError(
            f"Dataset separation violation: dataset_identity={only_ident!r} is not "
            f"an observational image collection allowed in production image cubes."
        )
    return {
        "valid_separation": True,
        "dataset_identity": next(iter(identities)) if identities else "SPHEREx",
        "items_checked": len(items),
    }


def mark_synthetic_fixture(payload: Mapping[str, Any], *, purpose: str = "unit_test") -> dict[str, Any]:
    """Explicitly label a synthetic test fixture so it can never be mistaken for
    real astronomical data.
    """
    tagged = dict(payload)
    tagged["dataset_identity"] = "Synthetic"
    tagged["is_synthetic_test_fixture"] = True
    tagged["is_real_observation"] = False
    tagged["allowed_in_production_responses"] = False
    tagged["synthetic_fixture_purpose"] = str(purpose)
    return tagged


class NEOSSatCADCAdapter:
    """Instrument-specific adapter for CSA/CADC NEOSSat optical observations.

    NEOSSat is a 15-cm optical space telescope operated by the Canadian Space Agency
    (CSA) and archived at the Canadian Astronomical Data Centre (CADC). Its FITS
    structure, optical CCD passband (350-1050 nm), pixel scale (~3.0 arcsec/px),
    and instrumental flux units (ADU / e-/s) differ fundamentally from SPHEREx
    Level-2 infrared spectral images (0.75-5.00 um, 6.155 arcsec/px, MJy/sr).

    This adapter inspects and validates NEOSSat/CADC FITS products independently
    and enforces strict cross-mission separation so NEOSSat images are never
    combined or subtracted directly against SPHEREx pixels.
    """

    DATASET_IDENTITY = "NEOSSat"
    CADC_COLLECTION = "NEOSSAT"
    CSA_DATASET_URL = "https://data.asc-csa.gc.ca/dataset/neossat"
    OFFICIAL_TUTORIAL_URL = "https://github.com/asc-csa/NEOSSAT_Tutorial"
    CADC_ARCHIVE_URL = "https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/"
    CITATION_NOTICE = (
        "Contains information licensed under the Open Government Licence – Canada "
        "(Canadian Space Agency NEOSSat dataset). Archive access provided by the "
        "Canadian Astronomical Data Centre (CADC)."
    )

    @classmethod
    def inspect_fits_product(
        cls,
        hdul_or_path: Any,
        *,
        download_url: str | None = None,
        archive_collection: str = "NEOSSAT",
    ) -> dict[str, Any]:
        """Inspect a NEOSSat / CADC FITS product and extract independent provenance.

        Does NOT assume SPHEREx's 7-HDU structure or MJy/sr surface-brightness units.
        """
        close_after = False
        if isinstance(hdul_or_path, (str, os.PathLike)):
            hdul = fits.open(hdul_or_path, memmap=False)
            close_after = True
        else:
            hdul = hdul_or_path

        try:
            ext_names = [
                str(hdu.name or f"HDU_{idx}").strip() for idx, hdu in enumerate(hdul)
            ]
            image_hdu = None
            for hdu in hdul:
                if getattr(hdu, "data", None) is not None and getattr(hdu.data, "ndim", 0) == 2:
                    image_hdu = hdu
                    break

            if image_hdu is None:
                return {
                    "dataset_identity": cls.DATASET_IDENTITY,
                    "is_spherex_measurement": False,
                    "valid": False,
                    "errors": ["No 2D image HDU found in NEOSSat/CADC FITS product."],
                    "hdu_names": ext_names,
                }

            hdr = image_hdu.header
            data = np.asarray(image_hdu.data, dtype=float)
            shape_yx = [int(s) for s in data.shape]
            finite_mask = np.isfinite(data)
            finite_fraction = float(np.mean(finite_mask)) if data.size > 0 else 0.0

            obs_id = (
                str(hdr.get("OBSID") or hdr.get("ARCHFILE") or hdr.get("FILENAME") or "").strip()
                or None
            )
            instrument = str(hdr.get("INSTRUME") or "NEOSSat").strip()
            telescope = str(hdr.get("TELESCOP") or hdr.get("OBSERVAT") or "NEOSSat / CSA").strip()
            bunit = str(hdr.get("BUNIT") or "ADU").strip()
            date_obs = str(hdr.get("DATE-OBS") or "").strip() or None
            exptime = hdr.get("EXPTIME")
            try:
                exptime_s = float(exptime) if exptime is not None else None
            except (TypeError, ValueError):
                exptime_s = None

            mjd_obs = hdr.get("MJD-OBS")
            if mjd_obs is not None:
                try:
                    mjd_val = float(mjd_obs)
                except (TypeError, ValueError):
                    mjd_val = None
            elif date_obs:
                try:
                    mjd_val = float(Time(date_obs).mjd)
                except Exception:
                    mjd_val = None
            else:
                mjd_val = None

            wcs_valid = False
            pixel_scale_arcsec = None
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    w = WCS(hdr)
                if w.has_celestial:
                    scales = proj_plane_pixel_scales(w.celestial) * 3600.0
                    pixel_scale_arcsec = round(float(np.mean(scales)), 4)
                    wcs_valid = True
            except Exception:
                wcs_valid = False

            errors: list[str] = []
            if finite_fraction <= 0.0:
                errors.append("NEOSSat image array contains zero finite pixels.")

            return {
                "dataset_identity": cls.DATASET_IDENTITY,
                "is_spherex_measurement": False,
                "archive_source": "CADC / CSA",
                "archive_collection": str(archive_collection or cls.CADC_COLLECTION),
                "observation_id": obs_id,
                "instrument": instrument,
                "facility": telescope,
                "download_url": download_url,
                "valid": len(errors) == 0,
                "errors": errors,
                "hdu_names": ext_names,
                "image_shape_yx": shape_yx,
                "finite_fraction": round(finite_fraction, 6),
                "image_unit": bunit,
                "wavelength_band_nm": [350.0, 1050.0],
                "observation_utc": date_obs,
                "observation_mjd": mjd_val,
                "exposure_seconds": exptime_s,
                "wcs_valid": wcs_valid,
                "pixel_scale_arcsec": pixel_scale_arcsec,
                "processing_history": [
                    "Loaded via independent NEOSSatCADCAdapter (not SPHEREx pipeline).",
                    "Native optical instrumental units preserved without SPHEREx ZODI subtraction.",
                ],
                "usage_and_citation": {
                    "licence": "Open Government Licence - Canada (CSA)",
                    "citation_notice": cls.CITATION_NOTICE,
                    "csa_dataset_url": cls.CSA_DATASET_URL,
                    "cadc_archive_url": cls.CADC_ARCHIVE_URL,
                },
            }
        finally:
            if close_after:
                hdul.close()

    @classmethod
    def build_cross_mission_comparison_summary(
        cls,
        spherex_metadata: Mapping[str, Any],
        neossat_metadata: Mapping[str, Any],
    ) -> dict[str, Any]:
        """Compare observation timing and metadata between SPHEREx and NEOSSat
        while explicitly prohibiting direct pixel subtraction or flux equivalence claims.
        """
        spx_mjd = (
            spherex_metadata.get("mjd")
            or spherex_metadata.get("time", {}).get("midpoint_mjd")
        )
        neo_mjd = neossat_metadata.get("observation_mjd")

        delta_days = None
        if spx_mjd is not None and neo_mjd is not None:
            try:
                delta_days = round(abs(float(spx_mjd) - float(neo_mjd)), 6)
            except (TypeError, ValueError):
                delta_days = None

        return {
            "comparison_type": "cross_mission_contextual_comparison_only",
            "missions_compared": ["SPHEREx (NASA/IPAC IRSA)", "NEOSSat (CSA/CADC)"],
            "direct_pixel_subtraction_allowed": False,
            "flux_equivalence_claimed": False,
            "timing_comparison": {
                "spherex_mjd": spx_mjd,
                "neossat_mjd": neo_mjd,
                "absolute_time_offset_days": delta_days,
            },
            "instrumentation_differences": {
                "spherex": {
                    "dataset_identity": "SPHEREx",
                    "wavelength_regime": "0.75 - 5.00 um (Near/Mid-IR Linear Variable Filters)",
                    "image_unit": spherex_metadata.get("product", {}).get("image_unit", "MJy / sr"),
                    "nominal_pixel_scale_arcsec": 6.155,
                },
                "neossat": {
                    "dataset_identity": "NEOSSat",
                    "wavelength_regime": "0.35 - 1.05 um (Optical unfiltered CCD)",
                    "image_unit": neossat_metadata.get("image_unit", "ADU"),
                    "nominal_pixel_scale_arcsec": neossat_metadata.get("pixel_scale_arcsec", 3.0),
                },
            },
            "scientific_warning": (
                "SPHEREx measures calibrated infrared surface brightness (MJy/sr) while "
                "NEOSSat measures broadband optical CCD counts (350-1050 nm). "
                "Cross-mission comparisons are restricted to astrometric/temporal context; "
                "image values must never be directly subtracted or combined."
            ),
        }


def get_approved_sources_manifest() -> dict[str, Any]:
    """Return the complete machine-readable manifest of ASTRAEA's approved
    scientific data source ecosystem and dataset separation rules.
    """
    return {
        "schema_version": "1.0",
        "approved_datasets": DATASET_IDENTITIES,
        "official_spherex_collections": list(OFFICIAL_SPHEREX_COLLECTIONS),
        "gaia_dr3_benchmark_astrometry": GAIA_DR3_BENCHMARK_ASTROMETRY,
        "dataset_separation_rules": {
            "silent_cross_mission_merging_prohibited": True,
            "synthetic_data_restricted_to_unit_tests": True,
            "wiseview_is_workflow_reference_only": True,
            "gaia_is_astrometric_reference_only": True,
            "neossat_requires_independent_adapter": True,
        },
    }
