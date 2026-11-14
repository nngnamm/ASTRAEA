from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

CACHE_DIR = Path("cache")

files = sorted(CACHE_DIR.glob("*.fits"))

if not files:
    raise RuntimeError("No FITS files found.")


flag_stats = defaultdict(
    lambda: {
        "flagged_pixels": 0,
        "files_with_flag": 0,
        "definitions_found": 0,
    }
)

detector_stats = defaultdict(
    lambda: {
        "files": 0,
        "total_pixels": 0,
        "image_nonfinite": 0,
        "variance_nonfinite": 0,
        "variance_zero": 0,
        "variance_negative": 0,
    }
)

date_records = []

# -------------------------------------------------
# Inspect every FITS file
# -------------------------------------------------

for path in files:
    with fits.open(path, memmap=False) as hdul:
        image_hdu = hdul["IMAGE"]
        flags_hdu = hdul["FLAGS"]
        variance_hdu = hdul["VARIANCE"]

        detector = int(image_hdu.header["DETECTOR"])

        image = np.asarray(image_hdu.data)
        flags = np.asarray(flags_hdu.data).astype(np.uint64)
        variance = np.asarray(variance_hdu.data)

        stats = detector_stats[detector]

        stats["files"] += 1
        stats["total_pixels"] += image.size
        stats["image_nonfinite"] += np.count_nonzero(~np.isfinite(image))

        finite_variance = np.isfinite(variance)

        stats["variance_nonfinite"] += np.count_nonzero(~finite_variance)

        stats["variance_zero"] += np.count_nonzero(finite_variance & (variance == 0))

        stats["variance_negative"] += np.count_nonzero(finite_variance & (variance < 0))

        # Observation date: use MJD-AVG, not DATE.
        header = image_hdu.header

        mjd = header.get("MJD-AVG")
        obsid = header.get("OBSID")
        product_date = header.get("DATE")

        observation_date = (
            Time(mjd, format="mjd", scale="utc").isot[:10] if mjd is not None else None
        )

        date_records.append(
            {
                "file": path.name,
                "detector": detector,
                "obsid": obsid,
                "mjd": mjd,
                "observation_date": observation_date,
                "product_date": product_date,
            }
        )

        # Dynamically discover flag definitions.
        definitions = {}

        for key in flags_hdu.header:
            if not key.startswith("MSKN"):
                continue

            suffix = key[4:]

            if not suffix.isdigit():
                continue

            bit = int(suffix)

            name = str(flags_hdu.header[key])

            description = str(flags_hdu.header.get("MSKD" + suffix, ""))

            definitions[bit] = (name, description)

        # Count each flag independently.
        for bit, (name, description) in definitions.items():
            mask = np.uint64(1) << np.uint64(bit)

            flagged = (flags & mask) != 0

            count = int(np.count_nonzero(flagged))

            record = flag_stats[(detector, bit, name)]

            record["flagged_pixels"] += count
            record["files_with_flag"] += int(count > 0)
            record["definitions_found"] += 1

            record["description"] = description


# -------------------------------------------------
# 1. Numerical integrity summary
# -------------------------------------------------

print("\n========== NUMERICAL INTEGRITY ==========")

for detector, stats in sorted(detector_stats.items()):
    total = stats["total_pixels"]

    print(f"\nDetector {detector}")
    print("Files:", stats["files"])

    print("Image finite:", f"{100 * (1 - stats['image_nonfinite'] / total):.4f}%")

    print("Variance finite:", f"{100 * (1 - stats['variance_nonfinite'] / total):.4f}%")

    print("Zero variance pixels:", stats["variance_zero"])
    print("Negative variance pixels:", stats["variance_negative"])


# -------------------------------------------------
# 2. Aggregate quality flag distribution
# -------------------------------------------------

print("\n========== AGGREGATE QUALITY FLAGS ==========")

for (detector, bit, name), record in sorted(flag_stats.items()):
    total = detector_stats[detector]["total_pixels"]

    percentage = 100 * record["flagged_pixels"] / total if total else 0

    print(
        f"D{detector} | "
        f"Bit {bit:2d} | "
        f"{name:12s} | "
        f"Pixels: {record['flagged_pixels']:7d} | "
        f"Rate: {percentage:8.4f}% | "
        f"Observations affected: "
        f"{record['files_with_flag']}/{detector_stats[detector]['files']}"
    )


# -------------------------------------------------
# 3. Observation and product timestamps
# -------------------------------------------------

print("\n========== OBSERVATION PROVENANCE ==========")

for record in sorted(
    date_records,
    key=lambda x: (x["mjd"] is None, x["mjd"] if x["mjd"] is not None else 0),
):
    print(
        f"{record['obsid']} | "
        f"MJD: {record['mjd']} | "
        f"Observation: {record['observation_date']} | "
        f"Product DATE: {record['product_date']}"
    )
