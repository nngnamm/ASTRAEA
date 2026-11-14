import json
from collections import defaultdict
from statistics import median

with open("/tmp/astraea_be09.json") as f:
    data = json.load(f)

groups = defaultdict(list)

for frame in data["frames"]:
    validation = frame.get("validation", {})

    groups[frame["detector"]].append(
        {
            "date": frame["iso_date"],
            "shape": validation.get("image_shape"),
            "finite": validation.get("finite_pixel_fraction"),
            "coverage": validation.get("reprojection_covered_fraction"),
        }
    )

for detector, observations in sorted(groups.items()):
    coverages = [x["coverage"] for x in observations if x["coverage"] is not None]

    print(f"\n========== DETECTOR {detector} ==========")
    print("Observations:", len(observations))
    print("Minimum coverage:", min(coverages))
    print("Median coverage:", median(coverages))
    print("Maximum coverage:", max(coverages))

    print("\nLowest coverage observations:")

    for item in sorted(observations, key=lambda x: x["coverage"] or 0)[:5]:
        print(item)
