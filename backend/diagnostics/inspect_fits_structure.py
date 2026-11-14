from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.time import Time

cache_dir = Path("cache")
files = sorted(cache_dir.glob("*.fits"))

print("Cached FITS files:", len(files))

for path in files:
    try:
        with fits.open(path, memmap=False) as hdul:
            image_header = hdul["IMAGE"].header

            detector = image_header.get("DETECTOR")
            mjd = image_header.get("MJD-AVG")

            date = Time(mjd, format="mjd").isot[:10] if mjd is not None else "Unknown"

            print("\n" + "=" * 65)
            print("FILE:", path.name)
            print("DATE:", date)
            print("DETECTOR:", detector)

            for index, hdu in enumerate(hdul):
                data = hdu.data
                name = hdu.name

                print(f"\nHDU {index}: {name}")

                if data is None:
                    print("  No data")
                    continue

                print("  Type:", type(data).__name__)

                if isinstance(data, np.ndarray):
                    print("  Shape:", data.shape)
                    print("  Dtype:", data.dtype)

                elif hasattr(data, "names"):
                    print("  Columns:", data.names)

    except Exception as exc:
        print("\nFAILED:", path.name, exc)

    print("-" * 65)
