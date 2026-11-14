from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

CACHE_DIR = Path("cache")

TARGET = SkyCoord(ra=269.446294 * u.deg, dec=4.767881 * u.deg, frame="icrs")

files = sorted(CACHE_DIR.glob("*.fits"))

if not files:
    raise RuntimeError("No cached FITS files found.")

results = []

print(f"Inspecting {len(files)} FITS files")


for path in files:
    try:
        with fits.open(path, memmap=False) as hdul:
            image_hdu = hdul["IMAGE"]
            image = image_hdu.data
            header = image_hdu.header

            ny, nx = image.shape

            detector = header.get("DETECTOR")
            obsid = header.get("OBSID")
            mjd = header.get("MJD-AVG")

            wcs = WCS(header).celestial

            if not wcs.has_celestial:
                raise ValueError("No celestial WCS found")

            # Check pixel scale
            pixel_scales = proj_plane_pixel_scales(wcs) * 3600

            # Image center
            cx = (nx - 1) / 2
            cy = (ny - 1) / 2

            center_sky = wcs.pixel_to_world(cx, cy)

            # Convert requested target to pixel coordinates
            tx, ty = wcs.world_to_pixel(TARGET)

            target_inside = -0.5 <= tx < nx - 0.5 and -0.5 <= ty < ny - 0.5

            # WCS round-trip test
            rx, ry = wcs.world_to_pixel(center_sky)

            roundtrip_error = np.hypot(rx - cx, ry - cy)

            # Check image corners in sky coordinates
            corner_pixels = [
                (-0.5, -0.5),
                (nx - 0.5, -0.5),
                (-0.5, ny - 0.5),
                (nx - 0.5, ny - 0.5),
            ]

            corner_coordinates = [wcs.pixel_to_world(x, y) for x, y in corner_pixels]

            results.append(
                {
                    "file": path.name,
                    "detector": detector,
                    "obsid": obsid,
                    "mjd": mjd,
                    "shape": image.shape,
                    "pixel_scale_arcsec": pixel_scales.tolist(),
                    "target_pixel": (float(tx), float(ty)),
                    "target_inside": target_inside,
                    "center_separation_arcsec": float(
                        center_sky.separation(TARGET).arcsec
                    ),
                    "roundtrip_error_pixels": float(roundtrip_error),
                    "corners": corner_coordinates,
                }
            )

    except Exception as exc:
        results.append({"file": path.name, "error": str(exc)})


# -------------------------------------------------
# SUMMARY
# -------------------------------------------------

valid = [r for r in results if "error" not in r]

print("\n========== WCS SUMMARY ==========")
print("Total files:", len(files))
print("Valid celestial WCS:", len(valid))
print("WCS failures:", len(results) - len(valid))

for detector in sorted({r["detector"] for r in valid}):
    subset = [r for r in valid if r["detector"] == detector]

    inside = sum(r["target_inside"] for r in subset)

    print(f"\nDetector {detector}")
    print("Observations:", len(subset))
    print("Target inside footprint:", inside)

    print(
        "Median pixel scale:",
        np.median([np.mean(r["pixel_scale_arcsec"]) for r in subset]),
        "arcsec/pixel",
    )

    print(
        "Maximum WCS round-trip error:",
        max(r["roundtrip_error_pixels"] for r in subset),
        "pixels",
    )


# -------------------------------------------------
# INDIVIDUAL OBSERVATIONS
# -------------------------------------------------

print("\n========== OBSERVATION DETAILS ==========")

for r in sorted(valid, key=lambda x: (x["detector"], x["mjd"])):
    print(
        f"D{r['detector']} | "
        f"{r['obsid']} | "
        f"MJD={r['mjd']} | "
        f"Shape={r['shape']} | "
        f"Target inside={r['target_inside']} | "
        f"Target pixel={r['target_pixel']} | "
        f"Scale={r['pixel_scale_arcsec']}"
    )


# -------------------------------------------------
# ERRORS
# -------------------------------------------------

print("\n========== WCS ERRORS ==========")

for r in results:
    if "error" in r:
        print(r["file"], r["error"])
