from pathlib import Path

import numpy as np
from astropy.io import fits

CACHE_DIR = Path("cache")

files = sorted(CACHE_DIR.glob("*.fits"))

if not files:
    raise RuntimeError("No cached FITS files found.")

samples = {}

print(f"Inspecting {len(files)} FITS files")


# --------------------------------------------------
# 1. Numerical integrity across all observations
# --------------------------------------------------

print("\n========== NUMERICAL INTEGRITY ==========")

for path in files:
    with fits.open(path, memmap=False) as hdul:
        image = np.asarray(hdul["IMAGE"].data)
        flags = np.asarray(hdul["FLAGS"].data)
        variance = np.asarray(hdul["VARIANCE"].data)
        zodi = np.asarray(hdul["ZODI"].data)

        detector = int(hdul["IMAGE"].header["DETECTOR"])

        samples.setdefault(detector, path)

        variance_finite = np.isfinite(variance)

        negative_variance = np.count_nonzero(variance[variance_finite] < 0)

        print(
            f"{path.name[:8]} | "
            f"D{detector} | "
            f"Shape: {image.shape} | "
            f"Image finite: {np.isfinite(image).mean():.4f} | "
            f"Variance finite: {variance_finite.mean():.4f} | "
            f"Negative variance: {negative_variance} | "
            f"ZODI finite: {np.isfinite(zodi).mean():.4f}"
        )


# --------------------------------------------------
# 2. Inspect one example per detector
# --------------------------------------------------

for detector, path in sorted(samples.items()):
    print(f"\n========== DETECTOR {detector} ==========")
    print("File:", path.name)

    with fits.open(path, memmap=False) as hdul:
        # FLAGS
        flags_hdu = hdul["FLAGS"]
        flags = np.asarray(flags_hdu.data).astype(np.uint32, copy=False)

        print("\nFLAG DEFINITIONS:")

        for key in flags_hdu.header.keys():
            if key.startswith("MSKN"):
                suffix = key[4:]
                name = flags_hdu.header[key]

                description = flags_hdu.header.get("MSKD" + suffix, "No description")

                bit = int(suffix)
                mask = np.uint32(1 << bit)

                count = np.count_nonzero(flags & mask)

                print(
                    f"Bit {bit}: {name} | Pixels: {count} | Description: {description}"
                )

        # EPSF
        epsf = hdul["EPSF"].data

        print("\nEPSF:")
        print("Rows:", len(epsf))
        print("Columns:", hdul["EPSF"].columns.names)
        print("EPSF array shape:", epsf["EPSF"].shape)

        # WCS-WAVE
        wave = hdul["WCS-WAVE"].data

        print("\nWCS-WAVE:")
        print("Rows:", len(wave))
        print("Columns:", hdul["WCS-WAVE"].columns.names)

        row = wave[0]

        x = np.asarray(row["X"])
        y = np.asarray(row["Y"])
        values = np.asarray(row["VALUES"])

        print("X shape:", x.shape)
        print("Y shape:", y.shape)
        print("VALUES shape:", values.shape)
        print("All values finite:", np.isfinite(values).all())

        if values.ndim == 3 and values.shape[-1] == 2:
            wavelength = values[..., 0]
            bandwidth = values[..., 1]

            print("Wavelength range:", np.nanmin(wavelength), np.nanmax(wavelength))

            print("Bandwidth range:", np.nanmin(bandwidth), np.nanmax(bandwidth))

        # Important image metadata
        header = hdul["IMAGE"].header

        print("\nIMAGE METADATA:")

        for key in ["BUNIT", "DETECTOR", "MJD-AVG", "OBSID", "DATE"]:
            print(f"{key}: {header.get(key, 'Unavailable')}")
