import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import astropy.units as u
from astropy.coordinates import SkyCoord

from config import TARGETS
from spherex_core import find_images, group_epochs, propagate_coords

print("========== NASA BENCHMARK TARGETS IRSA AUDIT ==========")
for key, t in TARGETS.items():
    name = t["name"]
    ra_now, dec_now = propagate_coords(
        t["ra"], t["dec"], t["pm_ra_arcsec_yr"], t["pm_dec_arcsec_yr"], t["epoch"]
    )
    coord = SkyCoord(ra_now * u.deg, dec_now * u.deg, frame="icrs")
    imgs = find_images(coord, fov_deg=0.1)
    det_counts = dict(sorted(Counter(m["detector"] for m in imgs).items()))

    ep_counts = {}
    for d in sorted(det_counts):
        d_frames = [{"time": m["mjd"]} for m in imgs if m["detector"] == d]
        ep_counts[f"D{d}"] = len(group_epochs(d_frames, gap_days=30))

    span_days = round(imgs[-1]["mjd"] - imgs[0]["mjd"], 1) if len(imgs) > 1 else 0.0
    print(
        f"{key:10s} ({name:16s}) | Total={len(imgs):3d} | "
        f"Span={span_days:5.1f} d | By det={det_counts} | Epochs={ep_counts}"
    )
