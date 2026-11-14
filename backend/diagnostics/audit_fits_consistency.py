from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

CACHE_DIR = Path("cache")
files = sorted(CACHE_DIR.glob("*.fits"))

if not files:
    raise RuntimeError("No FITS files found.")

definitions = defaultdict(lambda: defaultdict(set))
definition_presence = defaultdict(int)
detector_stats = defaultdict(lambda: defaultdict(int))
shape_variants = defaultdict(set)
unknown_active_bits = []

print(f"Auditing {len(files)} FITS products")


for path in files:
    with fits.open(path, memmap=False) as hdul:
        image = np.asarray(hdul["IMAGE"].data)
        flags = np.asarray(hdul["FLAGS"].data)
        variance = np.asarray(hdul["VARIANCE"].data)
        zodi = np.asarray(hdul["ZODI"].data)

        detector = int(hdul["IMAGE"].header["DETECTOR"])
        stats = detector_stats[detector]

        stats["files"] += 1
        stats["pixels"] += image.size

        # ------------------------------------------
        # 1. IMAGE / VARIANCE mask consistency
        # ------------------------------------------

        if image.shape != variance.shape:
            stats["shape_mismatches"] += 1
        else:
            image_bad = ~np.isfinite(image)
            variance_bad = ~np.isfinite(variance)

            stats["image_only_bad"] += int(np.count_nonzero(image_bad & ~variance_bad))

            stats["variance_only_bad"] += int(
                np.count_nonzero(variance_bad & ~image_bad)
            )

            stats["both_bad"] += int(np.count_nonzero(image_bad & variance_bad))

            stats["mask_disagreements"] += int(
                np.count_nonzero(image_bad != variance_bad)
            )

        # ------------------------------------------
        # 2. Extension integrity
        # ------------------------------------------

        for name, array in [
            ("FLAGS", flags),
            ("VARIANCE", variance),
            ("ZODI", zodi),
        ]:
            if array.shape != image.shape:
                stats[f"{name}_shape_mismatches"] += 1

        stats["zodi_nonfinite"] += int(np.count_nonzero(~np.isfinite(zodi)))

        # ------------------------------------------
        # 3. EPSF integrity
        # ------------------------------------------

        epsf_data = hdul["EPSF"].data
        epsf_values = np.asarray(epsf_data["EPSF"])

        shape_variants[detector].add(epsf_values.shape)

        if epsf_values.shape != (441, 33, 33):
            stats["epsf_shape_issues"] += 1

        stats["epsf_nonfinite"] += int(np.count_nonzero(~np.isfinite(epsf_values)))

        # ------------------------------------------
        # 4. WCS-WAVE integrity
        # ------------------------------------------

        wave_data = hdul["WCS-WAVE"].data

        if len(wave_data) != 1:
            stats["wcs_wave_row_issues"] += 1
        else:
            row = wave_data[0]

            x = np.asarray(row["X"])
            y = np.asarray(row["Y"])
            values = np.asarray(row["VALUES"])

            expected_shape = (len(y), len(x), 2)

            if values.shape != expected_shape:
                stats["wcs_wave_shape_issues"] += 1

            if not np.isfinite(values).all():
                stats["wcs_wave_nonfinite"] += 1

            if values.shape[-1] == 2:
                wavelength = values[..., 0]
                bandwidth = values[..., 1]

                stats["invalid_wavelength_values"] += int(
                    np.count_nonzero(np.isfinite(wavelength) & (wavelength <= 0))
                )

                stats["invalid_bandwidth_values"] += int(
                    np.count_nonzero(np.isfinite(bandwidth) & (bandwidth <= 0))
                )

        # ------------------------------------------
        # 5. Dynamic flag definitions
        # ------------------------------------------

        flag_header = hdul["FLAGS"].header
        file_definitions = {}

        for key in flag_header:
            if not key.startswith("MSKN"):
                continue

            suffix = key[4:]

            if not suffix.isdigit():
                continue

            bit = int(suffix)
            name = str(flag_header[key]).strip()

            description = str(flag_header.get("MSKD" + suffix, "")).strip()

            file_definitions[bit] = (name, description)

            definitions[detector][bit].add((name, description))

            definition_presence[(detector, bit)] += 1

        # ------------------------------------------
        # 6. Find active bits without definitions
        # ------------------------------------------

        flags_u = flags.astype(np.uint32, copy=False)

        active_bits = []

        for bit in range(32):
            mask = np.uint32(1 << bit)

            if np.any((flags_u & mask) != 0):
                active_bits.append(bit)

        missing = [bit for bit in active_bits if bit not in file_definitions]

        if missing:
            unknown_active_bits.append(
                {
                    "file": path.name,
                    "detector": detector,
                    "bits": missing,
                }
            )


# ==================================================
# REPORT
# ==================================================

print("\n========== IMAGE / VARIANCE MASK CONSISTENCY ==========")

for detector, stats in sorted(detector_stats.items()):
    print(f"\nDetector {detector}")

    for key in [
        "image_only_bad",
        "variance_only_bad",
        "both_bad",
        "mask_disagreements",
    ]:
        print(f"{key}: {stats[key]}")


print("\n========== EXTENSION INTEGRITY ==========")

for detector, stats in sorted(detector_stats.items()):
    print(f"\nDetector {detector}")

    for key in [
        "FLAGS_shape_mismatches",
        "VARIANCE_shape_mismatches",
        "ZODI_shape_mismatches",
        "zodi_nonfinite",
        "epsf_shape_issues",
        "epsf_nonfinite",
        "wcs_wave_row_issues",
        "wcs_wave_shape_issues",
        "wcs_wave_nonfinite",
        "invalid_wavelength_values",
        "invalid_bandwidth_values",
    ]:
        print(f"{key}: {stats[key]}")


print("\n========== FLAG DEFINITION CONSISTENCY ==========")

for detector, bit_map in sorted(definitions.items()):
    print(f"\nDetector {detector}")

    for bit, variants in sorted(bit_map.items()):
        print(
            f"Bit {bit}: "
            f"Definitions: {len(variants)} | "
            f"Found in {definition_presence[(detector, bit)]} files"
        )

        for name, description in sorted(variants):
            print(f"  {name}: {description}")


print("\n========== ACTIVE BITS WITHOUT DEFINITIONS ==========")

if unknown_active_bits:
    for item in unknown_active_bits:
        print(item)
else:
    print("No undefined active bits found.")
