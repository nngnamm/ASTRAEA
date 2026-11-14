"""Reproducible audit script for the five NASA Space Apps 2026 benchmark targets.

Usage:
    # Offline local-cache audit (deterministic, no network required):
    python backend/diagnostics/audit_benchmark_targets.py

    # Live IRSA SIA + local-cache audit (requires internet access to irsa.ipac.caltech.edu):
    python backend/diagnostics/audit_benchmark_targets.py --live --json audit_output.json
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from spherex_core import audit_all_benchmark_targets


def main():
    parser = argparse.ArgumentParser(
        description="Audit NASA SPHEREx benchmark target coordinates, SIA status, and local cache."
    )
    parser.add_argument(
        "--live",
        action="store_true",
        help="Query live IRSA SIA collections (spherex_dr1, spherex_qr3) in addition to local cache.",
    )
    parser.add_argument(
        "--fov",
        type=float,
        default=0.1,
        help="Field of view in degrees (default: 0.1).",
    )
    parser.add_argument(
        "--json",
        type=str,
        default=None,
        help="Optional path to write full JSON audit report.",
    )
    args = parser.parse_args()

    report = audit_all_benchmark_targets(
        fov_deg=args.fov,
        query_live_irsa=args.live,
    )

    print("========== NASA BENCHMARK TARGETS PROVENANCE & READINESS AUDIT ==========")
    print(f"Checked at (UTC) : {report['checked_at_utc']}")
    print(f"Live IRSA query  : {report['query_live_irsa']}")
    print(f"Reference MJD    : {report['reference_mjd']} (J2026.0)")
    print("-" * 73)

    for key, t in report["targets"].items():
        local = t["local_cache"]
        sia = t["sia_query"]
        prop = t["propagated_coordinates"]
        print(
            f"{key:10s} | {t['name']:16s} | status={t['verification_status']:27s} | "
            f"J2026=({prop['ra']:.5f}, {prop['dec']:.5f}) | "
            f"SIA={sia['query_status']} (obs={sia['total_observations']}) | "
            f"LocalFITS={local['valid_fits_count']:2d} (span={local['baseline_days']:6.2f}d, "
            f"det={local['by_detector']}, col={local['by_collection']}) | "
            f"CompareEligible={t['comparison_eligible']}"
        )

    if args.json:
        out_path = Path(args.json)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"\nSaved full JSON audit report to {out_path}")


if __name__ == "__main__":
    main()
