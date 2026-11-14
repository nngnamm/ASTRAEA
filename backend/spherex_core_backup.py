import os, re, time, hashlib, logging, warnings
import numpy as np
import requests
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.stats import sigma_clip, mad_std
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astroquery.ipac.irsa import Irsa
from reproject import reproject_interp
from concurrent.futures import ThreadPoolExecutor

from config import COLLECTIONS, CACHE_DIR, BANDS

logging.getLogger("astropy").setLevel(logging.ERROR)

def propagate_coords(ra0, dec0, pm_ra_arcsec_yr=0.0, pm_dec_arcsec_yr=0.0, epoch0=2000.0, to_mjd=61000.0):
    if pm_ra_arcsec_yr == 0.0 and pm_dec_arcsec_yr == 0.0:
        return ra0, dec0
    year_now = Time(to_mjd, format="mjd").decimalyear
    dt = year_now - epoch0
    ra1 = ra0 + (pm_ra_arcsec_yr * dt / 3600.0) / np.cos(np.radians(dec0))
    dec1 = dec0 + (pm_dec_arcsec_yr * dt / 3600.0)
    return float(ra1), float(dec1)

def find_images(coord, fov_deg):
    found = []
    for col in COLLECTIONS:
        try:
            r = Irsa.query_sia(pos=(coord, (fov_deg / 2) * u.deg), collection=col)
            t = r.to_table() if hasattr(r, "to_table") else r
        except Exception as e:
            print(col, "query failed:", str(e)[:60])
            continue
        for row in t:
            url = str(row["access_url"])
            m = re.search(r"_\d+D(\d)_spx", url)
            found.append({
                "url": url,
                "mjd": float(row["t_min"]),
                "detector": int(m.group(1)) if m else None,
                "collection": col
            })
    found.sort(key=lambda d: d["mjd"])
    return found

def load_frame(meta, coord, fov_deg):
    key = hashlib.md5(f"{meta['url']}{coord.ra.deg:.5f}{coord.dec.deg:.5f}{fov_deg}".encode()).hexdigest()
    path = f"{CACHE_DIR}/{key}.fits"
    if not os.path.exists(path):
        url = f"{meta['url']}?center={coord.ra.deg},{coord.dec.deg}&size={fov_deg}"
        ok = False
        for attempt in range(4):
            try:
                r = requests.get(url, timeout=120)
                if r.status_code == 503:
                    time.sleep(2 ** attempt)
                    continue
                r.raise_for_status()
                with open(path, "wb") as f_out:
                    f_out.write(r.content)
                ok = True
                break
            except Exception:
                time.sleep(2 ** attempt)
        if not ok:
            return None

    try:
        with fits.open(path) as hdul:
            hdr = hdul["IMAGE"].header.copy()
            det = int(hdr["DETECTOR"])
            x, y = WCS(hdr).world_to_pixel(coord)
            try:
                swcs = WCS(hdr, fobj=hdul, key="W")
                swcs.sip = None
                wave_val = float(swcs.pixel_to_world(x, y)[0].to(u.um).value)
            except Exception:
                wave_val = np.nan
            if not np.isfinite(wave_val) and det in BANDS:
                wave_val = float(np.mean(BANDS[det]))
            return {
                "time": float(hdr["MJD-AVG"]),
                "detector": det,
                "wave_um": wave_val,
                "image": hdul["IMAGE"].data.astype("float32"),
                "flags": hdul["FLAGS"].data.copy(),
                "header": hdr,
                "collection": meta["collection"]
            }
    except Exception as e:
        print("read failed:", str(e)[:60])
        return None

def make_grid(coord, fov_deg, pix_arcsec):
    n = int(round(fov_deg * 3600 / pix_arcsec))
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.crval = [coord.ra.deg, coord.dec.deg]
    w.wcs.crpix = [(n + 1) / 2, (n + 1) / 2]
    w.wcs.cdelt = [-pix_arcsec / 3600, pix_arcsec / 3600]
    return w, n

def build_cube(ra, dec, fov_deg=0.1, detectors=(5, 6), max_frames=40):
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    imgs = [m for m in find_images(coord, fov_deg) if m["detector"] in detectors]
    if not imgs:
        return {"frames": [], "message": "no SPHEREx images cover this position"}
    if len(imgs) > max_frames:
        imgs = [imgs[i] for i in np.linspace(0, len(imgs) - 1, max_frames).astype(int)]
    with ThreadPoolExecutor(max_workers=3) as ex:
        loaded = [f for f in ex.map(lambda m: load_frame(m, coord, fov_deg), imgs) if f]
    if not loaded:
        return {"frames": [], "message": "downloads failed"}

    scales = proj_plane_pixel_scales(WCS(loaded[0]["header"])) * 3600
    pix_arcsec = float(np.mean(scales))
    grid, n = make_grid(coord, fov_deg, pix_arcsec)
    for f in loaded:
        f["aligned"], _ = reproject_interp((f["image"], WCS(f["header"])), grid, shape_out=(n, n))
    loaded.sort(key=lambda f: f["time"])
    return {"frames": loaded, "grid": grid, "n": n, "pix_arcsec": pix_arcsec, "center": (ra, dec)}

def group_epochs(frames, gap_days=30):
    """Group one detector's observations into epochs separated by large gaps.

    Frames are sorted by observation time inside this function. Empty input
    returns an empty list instead of raising IndexError.
    """
    if not frames:
        return []

    ordered = sorted(frames, key=lambda frame: float(frame["time"]))
    epochs = []
    current = [ordered[0]]

    for previous, following in zip(ordered[:-1], ordered[1:]):
        gap = float(following["time"]) - float(previous["time"])
        if gap > gap_days:
            epochs.append(current)
            current = []
        current.append(following)

    epochs.append(current)
    return epochs


def group_epochs_by_detector(frames, gap_days=30):
    """Group observations independently for each detector.

    SPHEREx detectors represent different wavelength channels. Keeping their
    timelines separate prevents a temporal epoch or coadd from mixing bands.
    Returns ``{detector_number: [epoch, ...]}`` in detector-number order.
    """
    by_detector = {}
    detector_ids = sorted({int(frame["detector"]) for frame in frames})

    for detector in detector_ids:
        detector_frames = [
            frame for frame in frames
            if int(frame["detector"]) == detector
        ]
        by_detector[detector] = group_epochs(detector_frames, gap_days=gap_days)

    return by_detector

def coadd(frames):
    """Robustly combine aligned frames that share a detector and grid."""
    if not frames:
        raise ValueError("Cannot coadd an empty observation list")

    arrays = [np.asarray(frame["aligned"], dtype=float) for frame in frames]
    reference_shape = arrays[0].shape
    if any(array.shape != reference_shape for array in arrays):
        raise ValueError("All aligned images in a coadd must have the same shape")

    # A one-frame epoch is still a usable image product, but the caller should
    # report that its comparison has less redundancy than a multi-frame epoch.
    if len(arrays) == 1:
        return arrays[0].copy()

    stack = np.stack(arrays, axis=0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        clipped = sigma_clip(
            stack,
            sigma=2.5,
            cenfunc="median",
            stdfunc=mad_std,
            axis=0,
            masked=True,
        )
        with np.errstate(all="ignore"):
            return np.nanmedian(clipped.filled(np.nan), axis=0)


def difference(frames_a, frames_b, detector=None):
    """Create a background-subtracted difference for compatible observations.

    When ``detector`` is supplied, both epochs are filtered to that detector.
    Empty epochs and incompatible image shapes fail explicitly rather than
    producing a misleading or all-NaN difference product.
    """
    if detector is not None:
        frames_a = [f for f in frames_a if int(f["detector"]) == int(detector)]
        frames_b = [f for f in frames_b if int(f["detector"]) == int(detector)]

    if not frames_a or not frames_b:
        raise ValueError("Both comparison epochs must contain observations")

    img_a = coadd(frames_a)
    img_b = coadd(frames_b)
    if img_a.shape != img_b.shape:
        raise ValueError("Comparison epochs must use the same aligned image grid")

    valid_a = img_a[np.isfinite(img_a)]
    valid_b = img_b[np.isfinite(img_b)]
    if valid_a.size == 0 or valid_b.size == 0:
        raise ValueError("Comparison epoch contains no finite image pixels")

    img_a_sub = img_a - np.median(valid_a)
    img_b_sub = img_b - np.median(valid_b)
    diff = img_b_sub - img_a_sub

    if not np.isfinite(diff).any():
        raise ValueError("Difference image contains no finite pixels")

    return diff, img_a_sub, img_b_sub

def compute_spectrum(ra, dec, fov_deg=0.1, max_per_det=8, r_ap=2.5):
    coord = SkyCoord(ra * u.deg, dec * u.deg)
    all_imgs = find_images(coord, fov_deg)
    if not all_imgs:
        return []

    picked = []
    for d in range(1, 7):
        d_imgs = [m for m in all_imgs if m["detector"] == d]
        if len(d_imgs) > max_per_det:
            d_imgs = [d_imgs[i] for i in np.linspace(0, len(d_imgs) - 1, max_per_det).astype(int)]
        picked.extend(d_imgs)

    with ThreadPoolExecutor(max_workers=3) as ex:
        loaded = [f for f in ex.map(lambda m: load_frame(m, coord, fov_deg), picked) if f]
    if not loaded:
        return []

    scales = proj_plane_pixel_scales(WCS(loaded[0]["header"])) * 3600
    grid, n = make_grid(coord, fov_deg, float(np.mean(scales)))
    c = (n - 1) / 2.0
    yy, xx = np.mgrid[0:n, 0:n]

    def aper_mean(img, cx, cy):
        m = np.hypot(xx - cx, yy - cy) <= r_ap
        v = img[m]
        return float(np.nanmean(v)) if np.isfinite(v).sum() > 0.5 * m.sum() else np.nan

    spectrum = []
    for d in range(1, 7):
        g = [f for f in loaded if f["detector"] == d]
        if len(g) < 2:
            continue
        for f in g:
            f["aligned"], _ = reproject_interp((f["image"], WCS(f["header"])), grid, shape_out=(n, n))
        img = coadd(g)

        vals = [
            aper_mean(img, x, y)
            for x in range(8, n - 8, 3)
            for y in range(8, n - 8, 3)
            if np.hypot(x - c, y - c) > 8
        ]
        vals = np.array([v for v in vals if np.isfinite(v)])
        center_flux = aper_mean(img, c, c)
        if len(vals) == 0 or not np.isfinite(center_flux):
            continue

        med = float(np.median(vals))
        sig = float(1.4826 * np.median(np.abs(vals - med)))
        if not np.isfinite(sig) or sig <= 0:
            continue

        excess = float(center_flux - med)
        wave = float(np.nanmean([f["wave_um"] for f in g]))
        if not np.isfinite(wave):
            wave = float(np.mean(BANDS[d]))

        spectrum.append({
            "detector": d,
            "n_frames": len(g),
            "wave_um": round(wave, 3),
            "flux_mjy_sr": round(excess, 4),
            "err_mjy_sr": round(sig, 4),
            "snr": round(excess / sig, 1)
        })
    return spectrum
