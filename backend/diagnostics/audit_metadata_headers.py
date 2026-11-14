from collections import defaultdict
from pathlib import Path

from astropy.io import fits

files = sorted(Path("cache").glob("*.fits"))

header_keys = defaultdict(lambda: defaultdict(int))
samples = defaultdict(dict)

for path in files:
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            extension = hdu.name

            for key in hdu.header:
                if key in ("", "COMMENT", "HISTORY"):
                    continue

                header_keys[extension][key] += 1

                samples[extension][key] = hdu.header[key]


print("TOTAL FILES:", len(files))

for extension in sorted(header_keys):
    print(f"\n========== {extension} ==========")

    for key, count in sorted(header_keys[extension].items()):
        value = samples[extension][key]

        print(f"{key:20s} ({count}/{len(files)} files): {value}")
