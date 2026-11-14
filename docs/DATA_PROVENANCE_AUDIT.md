# Data Provenance & Scientific Integrity Audit

**Repository**: `https://github.com/nngnamm/ASTRAEA`  
**Scope**: Team 02 — Backend, Data Pipeline & Astronomical Processing  
**Machine-Readable Ledger**: [`audit/data_provenance.json`](../audit/data_provenance.json)

---

## 1. Provenance Classification Taxonomy

Every dataset, coordinate, header value, derived image, score, and static claim in the repository has been audited and assigned to one of eight strict provenance categories:

1. **Verified real data (`verified_real_data`)** — Traced to an identifiable authoritative upstream source (NASA/IPAC IRSA SPHEREx Level-2 FITS archive) with cryptographic MD5/SHA-256 and FITS header verification.
2. **Derived real data (`derived_real_data`)** — Calculated deterministically from verified observations using documented astronomical processing (ICRS TAN WCS reprojection, `ZODI` subtraction, flag-aware sigma-clipped sky pedestal removal, 2.5$\sigma$-clipped median epoch coaddition, and joint-mask difference imaging).
3. **Heuristic estimate (`heuristic_estimate`)** — Generated using a documented approximation or engineering formula (`compute_sky_health`, `detect_anomalies` SNR confidence mapping, 30-day visit-gap epoch clustering).
4. **Synthetic data (`synthetic_data`)** — Artificially generated arrays or FITS structures used strictly inside `backend/tests/` for deterministic unit testing.
5. **Mock or placeholder data (`mock_or_placeholder_data`)** — Intentionally substituted values when upstream metadata was not persisted (e.g., the prior offline generator labeling all 30 Barnard's Star frames as `spherex_qr2`; remediated in this audit).
6. **Hardcoded configuration (`hardcoded_configuration`)** — Manually supplied catalog coordinates, proper motions, and wavelength passband bounds in `backend/config.py` with documented literature/catalog citations.
7. **Unverified data (`unverified_data`)** — Claims or live archive states whose origin or current availability cannot be established offline (`eps_ind_b`, `wise0855`, `rw_aur`, and `r_dra` live FITS availability without network access).
8. **Unsupported or potentially fabricated claim (`unsupported_or_fabricated_claim`)** — Presented as factual without sufficient evidence (e.g., the static `irsa_audit` blocks and `"Confirmed Benchmark Object"` labels; remediated in this audit).

---

## 2. Audit of the 30 Cached SPHEREx FITS Products (`verified_real_data`)

All 30 `.fits` files in `backend/cache/` are **genuine NASA/IPAC Infrared Science Archive (IRSA) SPHEREx Level-2 spectral image cutouts** (`VERSION = '7.0.5'`, SSDC pipeline tags `l2b-v27-2026-178` through `l2b-v27-2026-278`).

### 2.1 Upstream Cryptographic & Header Verification

- **Cutout Query Parameters**: All 30 files were retrieved from the IRSA Image Backend (`https://irsa.ipac.caltech.edu/ibe/data/spherex/...`) for Barnard's Star at:
  - `center = (269.44629, 4.76788)` deg (Barnard's Star propagated to MJD `61000.0`)
  - `size = 0.1` deg ($6' \times 6'$ cutout)
- **MD5 Preimage Verification**: Every filename `<md5>.fits` in `backend/cache/` was cryptographically verified against:
  $$\text{MD5}\left(\texttt{access\_url} \;\|\; \texttt{"269.44629"} \;\|\; \texttt{"4.76788"} \;\|\; \texttt{"0.1"}\right)$$
  - **19 files** match `spherex_dr1` (`https://irsa.ipac.caltech.edu/ibe/data/spherex/dr1/level2/...`)
  - **11 files** match `spherex_qr3` (`https://irsa.ipac.caltech.edu/ibe/data/spherex/qr3/level2/...`)
- **Detector Distribution**:
  - **Detector 5** ($3.83\text{–}4.41\,\mu\text{m}$): **16 files**
  - **Detector 6** ($4.42\text{–}5.00\,\mu\text{m}$): **14 files**
- **HDU Structure (all 30 files)**:
  - `0: PRIMARY`
  - `1: IMAGE` (`float32`, `BUNIT = 'MJy / sr'`, `RADESYS = 'ICRS'`, `CTYPE1 = 'RA---TAN-SIP'`, `CTYPE2 = 'DEC--TAN-SIP'`)
  - `2: FLAGS` (`int32`, 32-bit quality bitmap with `MSKN0000..MSKN0021` header definitions)
  - `3: VARIANCE` (`float32`, surface-brightness variance in $(\text{MJy/sr})^2$)
  - `4: ZODI` (`float32`, zodiacal light model in $\text{MJy/sr}$)
  - `5: EPSF` (`BinTableHDU`, 441 rows, $33 \times 33$ empirical PSF grid)
  - `6: WCS-WAVE` (`BinTableHDU`, 1 row, $11 \times 11 \times 2$ approximate wavelength/bandpass lookup grid)
- **Observation Time vs. Product Creation Time**:
  - Astronomical observation midpoints (`MJD-AVG` / `DATE-AVG`) span **MJD `60902.865381`** (`2025-08-15T20:46:08.896` UTC) to **MJD `61305.628658`** (`2026-09-22T15:05:16.058` UTC) — a **402.76-day** temporal baseline.
  - Archive product creation timestamps (`DATE`) are pipeline processing dates in August–October 2026 (`2026-08-15` to `2026-10-05`) and are strictly separated from `MJD-AVG` in `observation_metadata.py` and `spherex_core.py`.

### 2.2 Complete Ledger of the 30 Verified Cached FITS Files

| Cache File (`backend/cache/`) | Collection | Det | `OBSID` | `EXPIDN` | `MJD-AVG` | `DATE-AVG` (UTC) | `DATE` Created (UTC) | Shape $(Y,X)$ | SHA-256 Checksum (first 16 hex) |
| :--- | :--- | :---: | :--- | :---: | :---: | :--- | :--- | :---: | :--- |
| `0140ea2816e7d2d22dbbca47d8ed45e0.fits` | `spherex_dr1` | 5 | `2025W34_2B_0557_1` | `202534202055715` | `60911.491271` | `2025-08-24T11:47:25` | `2026-08-17T22:11:36` | `58x59` | `8403ffbd35c978d6` |
| `0e7458b5b765209bc3dbadbc10517f8c.fits` | `spherex_dr1` | 6 | `2026W12_1A_0598_3` | `202612101059836` | `61120.728945` | `2026-03-21T17:29:40` | `2026-06-27T11:25:58` | `58x59` | `19cc6deb00b899c8` |
| `2cf2d3ce43550860e4343e2f935b41bb.fits` | `spherex_qr3` | 5 | `2026W36_2A_0276_3` | `202636201027635` | `61291.771252` | `2026-09-08T18:30:36` | `2026-09-21T20:20:32` | `58x58` | `0da82dc983ca5a0f` |
| `370d9894e8a5513d595417a44a3e08f0.fits` | `spherex_dr1` | 5 | `2025W36_2C_0135_2` | `202536203013525` | `60925.568165` | `2025-09-07T13:38:09` | `2026-08-24T22:16:37` | `58x58` | `ec9934c69a67a541` |
| `3cec8479d10e8b116c737521e5309b06.fits` | `spherex_dr1` | 6 | `2025W39_1A_0568_2` | `202539101056826` | `60946.060454` | `2025-09-28T01:27:03` | `2026-08-31T21:16:40` | `58x58` | `fd732ff7481c8549` |
| `426bb54bae875f163c9b637e479b378d.fits` | `spherex_dr1` | 6 | `2025W37_1A_0632_2` | `202537101063226` | `60932.273837` | `2025-09-14T06:34:19` | `2026-08-29T00:54:51` | `10x58` | `4c6e03dc131e02cd` |
| `4322681cb4f8653db9a29466f2ede97d.fits` | `spherex_dr1` | 6 | `2025W33_2C_0386_3` | `202533203038636` | `60902.863992` | `2025-08-15T20:44:08` | `2026-08-15T13:49:46` | `58x59` | `b7346c04eb8a30dc` |
| `46046aeb3c058ecb486ee1fa9fbe2e7c.fits` | `spherex_qr3` | 6 | `2026W35_2A_0275_2` | `202635201027526` | `61285.910212` | `2026-09-02T21:50:42` | `2026-09-20T05:19:45` | `58x58` | `69c35f143b5ab6d4` |
| `49e43b144df52f4649a7d8d4c66a9d31.fits` | `spherex_qr3` | 5 | `2026W34_2B_0328_3` | `202634202032835` | `61276.647218` | `2026-08-24T15:31:59` | `2026-09-16T04:56:41` | `58x59` | `07f6a77e122d86df` |
| `5be173f068aef88913ae954374edda98.fits` | `spherex_dr1` | 6 | `2026W16_1A_0648_2` | `202616101064826` | `61144.307507` | `2026-04-14T07:22:48` | `2026-06-29T18:10:04` | `58x58` | `55fd619426ff1e3a` |
| `636a3c537b5e7887aaf5d862388b8786.fits` | `spherex_qr3` | 6 | `2026W38_1A_0260_1` | `202638101026016` | `61299.166530` | `2026-09-16T03:59:48` | `2026-10-03T18:56:59` | `58x59` | `8ff15ff3b14d28b5` |
| `6e7567d875e6c55ae4ddd93d94b03b4f.fits` | `spherex_dr1` | 5 | `2026W16_1A_0391_3` | `202616101039135` | `61141.737742` | `2026-04-11T17:42:20` | `2026-06-29T14:43:08` | `59x29` | `157e08819e35d5c6` |
| `6ed482a6acdc2bf5d92d6dbabd51653d.fits` | `spherex_dr1` | 5 | `2025W37_2A_0630_3` | `202537201063035` | `60934.398216` | `2025-09-16T09:33:25` | `2026-08-29T05:19:01` | `58x59` | `f664cfbd8e7969e8` |
| `71aee69864d5981d539a05bf61a5470a.fits` | `spherex_dr1` | 5 | `2025W38_1B_0261_3` | `202538102026135` | `60939.265417` | `2025-09-21T06:22:12` | `2026-08-29T05:49:10` | `21x58` | `0dc786a5e4f5be7e` |
| `78e0dd2be6e5bf934e73f926856c5d5b.fits` | `spherex_dr1` | 5 | `2026W12_1A_0664_1` | `202612101066415` | `61121.391167` | `2026-03-22T09:23:16` | `2026-06-27T11:46:47` | `59x58` | `f524c75645cb36a7` |
| `7c3ffc77f3cfaad1f0a364ca0870048a.fits` | `spherex_qr3` | 5 | `2026W33_2A_0277_3` | `202633201027735` | `61268.629577` | `2026-08-16T15:06:35` | `2026-09-09T04:17:04` | `58x58` | `2aa8e97eb50dc7c9` |
| `84d6482ff47b3660b9cc6e748843aa4d.fits` | `spherex_qr3` | 5 | `2026W35_2A_0341_1` | `202635201034115` | `61286.634873` | `2026-09-03T15:14:13` | `2026-09-20T05:47:13` | `58x58` | `a2ab7e140fd9de90` |
| `94373761ba6ee0cc6cc94093b34f3755.fits` | `spherex_qr3` | 6 | `2026W39_1A_0448_2` | `202639101044826` | `61305.628658` | `2026-09-22T15:05:16` | `2026-10-05T07:53:42` | `58x58` | `e40886bb24dc4bd9` |
| `94d62774fd06a9cb0c527bf1a8e89376.fits` | `spherex_dr1` | 6 | `2025W33_2C_0386_4` | `202533203038646` | `60902.865381` | `2025-08-15T20:46:08` | `2026-08-15T13:49:47` | `58x58` | `ef3f2a8b173c8502` |
| `95b8371b8d74707410276b70fa06c435.fits` | `spherex_dr1` | 5 | `2026W14_2A_0330_3` | `202614201033035` | `61131.676177` | `2026-04-01T16:13:41` | `2026-06-28T13:02:13` | `59x58` | `bb60179bc0f729f2` |
| `b8951e640e7056445ad8656ebc484b4e.fits` | `spherex_dr1` | 5 | `2026W19_1A_0448_3` | `202619101044835` | `61164.133758` | `2026-05-04T03:12:36` | `2026-07-02T01:40:00` | `59x26` | `6f4aa8798a93f944` |
| `ba89f804b7948f98d9e55ec98e969e64.fits` | `spherex_qr3` | 5 | `2026W37_1C_0276_3` | `202637103027635` | `61296.453037` | `2026-09-13T10:52:22` | `2026-09-25T20:32:54` | `58x59` | `159fa4308dc4d9cf` |
| `bbcb28c5e4300752b75e2f97e9c6df24.fits` | `spherex_dr1` | 6 | `2025W35_2B_0619_1` | `202535202061916` | `60918.563727` | `2025-08-31T13:31:46` | `2026-08-21T20:51:04` | `58x59` | `58cd00bb54b9fd4c` |
| `be12b58d91e7db2abb2e27b014b9df70.fits` | `spherex_qr3` | 6 | `2026W36_2A_0276_2` | `202636201027626` | `61291.769713` | `2026-09-08T18:28:23` | `2026-09-21T20:20:32` | `58x58` | `78535d8259dc9cf0` |
| `bed9acfa8fb34a3a4a8d70559abaeb7b.fits` | `spherex_dr1` | 6 | `2026W14_2A_0330_2` | `202614201033026` | `61131.674789` | `2026-04-01T16:11:41` | `2026-06-28T13:02:13` | `59x58` | `1e8355bc14a162f6` |
| `c22793f0b31233b7ad9b86094d1545f2.fits` | `spherex_qr3` | 6 | `2026W33_2A_0277_2` | `202633201027726` | `61268.628188` | `2026-08-16T15:04:35` | `2026-09-09T04:17:04` | `58x59` | `6814217eb04fc323` |
| `e394055ec44abbe35608738c87340c13.fits` | `spherex_dr1` | 5 | `2026W15_2A_0391_3` | `202615201039135` | `61137.456764` | `2026-04-07T10:57:44` | `2026-06-29T00:05:19` | `58x58` | `b8dc43949c0ea585` |
| `eb3f217b688e66154f4d8679394bde6e.fits` | `spherex_dr1` | 6 | `2025W36_2C_0135_1` | `202536203013516` | `60925.566776` | `2025-09-07T13:36:09` | `2026-08-24T22:16:36` | `58x58` | `e5268f7f7d3a0808` |
| `ebec6e19177daf652cc1bc959547e8c8.fits` | `spherex_dr1` | 6 | `2026W18_2B_0503_2` | `202618202050326` | `61159.455774` | `2026-04-29T10:56:18` | `2026-07-01T14:50:00` | `58x14` | `e97ce33fa3e59f15` |
| `ffe72d5befb7f1da0f159dafa42cb49e.fits` | `spherex_qr3` | 5 | `2026W39_1A_0514_2` | `202639101051425` | `61305.893506` | `2026-09-22T21:26:38` | `2026-10-05T08:11:36` | `58x59` | `1da60a6c3089cb88` |

---

## 3. Independent Audit of All Five NASA Benchmark Targets

Each of the five benchmark targets configured in `backend/config.py` was audited independently for coordinate provenance, reference epoch, proper motion convention, local cache availability, alignment behavior, and comparison eligibility:

### 3.1 `barnard` — Barnard's Star (GJ 699 / HIP 87937)
- **Catalog Coordinates**: $\alpha = 269.452077^\circ$, $\delta = +4.693365^\circ$, Epoch `J2000.0` (ICRS; SIMBAD / Hipparcos / Gaia DR3).
- **Proper Motion**: $\mu_{\alpha^*} = -0.8015''/\text{yr}$, $\mu_\delta = +10.3624''/\text{yr}$ (`pm_ra_cosdec`; Gaia DR3 reports $-801.551\text{ mas/yr}$, $+10362.394\text{ mas/yr}$).
- **Propagated Position (`J2026.0`, MJD `61041.0`)**: $\alpha = 269.446204^\circ$, $\delta = +4.768205^\circ$ (total motion $\sim 270.6''$ since J2000.0; proper motion is critically relevant).
- **Local Cache Status**: **`verified_local_cache_subset`** — 30 verified FITS cutouts in `backend/cache/` (19 `spherex_dr1`, 11 `spherex_qr3`; 16 Detector 5, 14 Detector 6).
- **Alignment & Coverage**: All 30 files pass structural, numerical, and WCS validation and reproject onto a $59 \times 59$ pixel ICRS TAN grid ($6.153''/\text{px}$). 25 of 30 cutouts cover the central target (median coverage $97.3\%$); 5 cutouts (`426bb54b`, `6e7567d8`, `71aee698`, `b8951e64`, `ebec6e19`) are partial-footprint edge strips ($16.5\%\text{–}45.0\%$ coverage).
- **Usable Epochs & Comparison Eligibility**:
  - **Detector 5**: 16 frames across 3 visit groups (`MJD 60911–60939`, `61121–61164`, `61268–61305`), all 3 epochs have $100\%$ coadd coverage ($\ge 50\%$). **Eligible for comparison (`comparison_eligible: True`)**.
  - **Detector 6**: 14 frames across 4 visit groups (`MJD 60902–60946`, `61120–61144`, `61159` single partial frame at $23.1\%$ coverage, `61268–61305`). 3 of 4 epochs have $\ge 96\%$ coverage; the low-coverage single-frame epoch is excluded by the $50\%$ coverage gate. **Eligible for comparison (`comparison_eligible: True`)**.

### 3.2 `eps_ind_b` — Epsilon Indi B (HD 209100 B)
- **Catalog Coordinates**: $\alpha = 331.044147^\circ$, $\delta = -56.782844^\circ$, Epoch `J2000.0` (ICRS; Scholz et al. 2003 / SIMBAD).
- **Proper Motion**: $\mu_{\alpha^*} = +3.9614''/\text{yr}$, $\mu_\delta = -2.5383''/\text{yr}$ (`pm_ra_cosdec`).
- **Propagated Position (`J2026.0`, MJD `61041.0`)**: $\alpha = 331.096402^\circ$, $\delta = -56.801162^\circ$ (proper motion highly relevant: $\sim 122.4''$ over 26 yr).
- **Local Cache Status**: **`unverified_offline`** — 0 cached FITS files and 0 cached SIA JSON files in `backend/cache/`.
- **Limitations**: Requires live network access to `irsa.ipac.caltech.edu` (`python backend/diagnostics/audit_benchmark_targets.py --live`) to retrieve SIA metadata and FITS cutouts.

### 3.3 `wise0855` — WISE J085510.83-071442.5 (WISE 0855-0714)
- **Catalog Coordinates**: $\alpha = 133.795129^\circ$ (`08h55m10.831s`), $\delta = -7.245147^\circ$ (`-07d14m42.53s`), Epoch `J2010.5` (WISE discovery epoch; Luhman 2014, ApJ, 786, L18). Note: SIMBAD's back-propagated J2000 position is `(133.818908, -7.247048)`.
- **Proper Motion**: $\mu_{\alpha^*} = -8.118''/\text{yr}$, $\mu_\delta = +0.680''/\text{yr}$ (`pm_ra_cosdec`; Luhman & Esplin 2016).
- **Propagated Position (`J2026.0`, MJD `61041.0`)**: $\alpha = 133.759898^\circ$, $\delta = -7.242219^\circ$ (proper motion critically relevant: $8.15''/\text{yr}$).
- **Local Cache Status**: **`unverified_offline`** — 0 cached FITS files in `backend/cache/`.
- **Historical Notebook Evidence**: In `notebooks/NASA_backend_exploration.ipynb` (Cells 8, 31, 36), a live query against `spherex_dr1` returned 263 total SIA metadata rows (94 in D5+D6, spanning 373.6 days across 2 D5 visit groups and 3 D6 visit groups), and noted that WISE 0855-0714 ($T_{\text{eff}} \approx 250\text{ K}$, $W2 \approx 14\text{ mag}$) is near the single-visit surface-brightness noise floor of SPHEREx Level-2 cutouts. Because those FITS files were not saved in `backend/cache/`, its current repository status is `unverified_offline`.

### 3.4 `rw_aur` — RW Aurigae (HD 240764)
- **Catalog Coordinates**: $\alpha = 76.956526^\circ$, $\delta = +30.401438^\circ$, Epoch `J2000.0` (ICRS; Gaia DR3 / SIMBAD).
- **Proper Motion**: $\mu_{\alpha^*} = +0.003''/\text{yr}$, $\mu_\delta = -0.024''/\text{yr}$ (negligible over 26 yr: $< 0.7''$).
- **Propagated Position (`J2026.0`, MJD `61041.0`)**: $\alpha = 76.956551^\circ$, $\delta = +30.401265^\circ$.
- **Local Cache Status**: **`unverified_offline`** — 0 cached FITS files and 0 cached SIA JSON files in `backend/cache/`. Previous static counts (`total_observations: 342`, `by_detector: {1: 57, ..., 6: 57}`) were identified as fabricated and removed.

### 3.5 `r_dra` — R Draconis (HD 150072)
- **Catalog Coordinates**: $\alpha = 248.167611^\circ$, $\delta = +66.754979^\circ$, Epoch `J2000.0` (ICRS; Gaia DR3 / SIMBAD).
- **Proper Motion**: $\mu_{\alpha^*} = -0.006''/\text{yr}$, $\mu_\delta = -0.033''/\text{yr}$ (small: $\sim 0.9''$ over 26 yr).
- **Propagated Position (`J2026.0`, MJD `61041.0`)**: $\alpha = 248.167501^\circ$, $\delta = +66.754741^\circ$.
- **Local Cache Status**: **`unverified_offline`** — 0 cached FITS files and 0 cached SIA JSON files in `backend/cache/`.

---

## 4. Summary of Detected & Remediated Data Fabrications

Full forensic details of all AI-originated / unverified claims found during the audit are documented in [`docs/AI_ORIGIN_AUDIT.md`](AI_ORIGIN_AUDIT.md):

1. **Fabricated `irsa_audit` Observation Counts & `verified_available` Labels (`config.py`, P0)** — Removed and replaced with honest `verified_local_cache_subset` (Barnard's Star) and `unverified_offline` (other 4 targets) plus `audit_all_benchmark_targets()`.
2. **Unsupported `"Confirmed Benchmark Object"` Candidate Classification (`analytics.py`, P0)** — Replaced with `"Benchmark Coordinate Match — Unconfirmed (<name>)"` (`independent_confirmation=False`, `human_reviewed=False`, `is_calibrated_probability=False`).
3. **Synthetic `collection="spherex_qr2"` in `sample_barnard_response.json` (P1)** — Regenerated using the cryptographically verified `spherex_dr1` (19) and `spherex_qr3` (11) provenance records.
4. **Epoch Mismatch (`to_mjd=61000.0` Labeled as `ra_2026`/`dec_2026`, P1)** — Corrected to exact `J2026.0` (`MJD 61041.0`) with explicit epoch metadata.
5. **Silent Conversion of Upstream IRSA Failure into `no_observations_found` (`spherex_core.py`, P0)** — Fixed to track per-collection status, prevent caching incomplete queries, and return HTTP 503 `upstream_service_unavailable`.

---

## 7. Approved Scientific Data Source Ecosystem & Dataset Separation (`backend/approved_data_sources.py`)

ASTRAEA enforces strict provenance separation across five distinct dataset identities (`validate_dataset_separation()`):

| Dataset Identity | Role | Authoritative Source / Archive | Key Rules & Safeguards |
| :--- | :--- | :--- | :--- |
| **`SPHEREx`** | Primary Observational Dataset | NASA/IPAC IRSA (`https://irsa.ipac.caltech.edu/Missions/spherex.html`) | Official collections investigated dynamically: `spherex_dr1`, `spherex_dr1_deep`, `spherex_qr3`, `spherex_qr3_deep`. Collection provenance verified from `access_url` (`infer_spherex_collection_from_url()`). |
| **`Gaia`** | Astrometric Reference Catalog | ESA Gaia DR3 / STScI (`https://archive.stsci.edu/missions-and-data/gaia`) | Preserves `source_id`, `RA/Dec`, `reference_epoch_jyear`, `pm_ra_cosdec_mas_yr`, `pm_dec_mas_yr`, `parallax_mas`, `radial_velocity_km_s`, and `ruwe`. Explicitly marked `is_spherex_measurement: false`. Targets absent from Gaia DR3 (`wise0855`, `eps_ind_b`) document independent peer-reviewed IR astrometry and confidence. |
| **`NEOSSat`** | Supplementary Mission | CSA (`https://data.asc-csa.gc.ca/dataset/neossat`) & CADC (`https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/`) | Handled via independent `NEOSSatCADCAdapter`. Preserves optical broadband (`350-1050 nm`) units (`ADU`), pixel scale (`~3.0"/px`), CADC collection (`NEOSSAT`), and Open Government Licence — Canada citation. Direct pixel subtraction/merging with SPHEREx is prohibited (`direct_pixel_subtraction_allowed: false`). |
| **`WISE`** | Workflow & Visualization Inspiration | WiseView (`http://byw.tools/wiseview`) | Used strictly as a product design and blinking workflow reference (`is_substitute_for_spherex: false`). |
| **`Synthetic`** | Unit Test Fixtures Only | `backend/tests/` (`mark_synthetic_fixture()`) | Explicitly tagged `dataset_identity: "Synthetic"`, `allowed_in_production_responses: false`. Never substituted for unavailable archive data. |

