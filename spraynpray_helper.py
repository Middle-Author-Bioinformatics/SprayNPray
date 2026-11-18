#!/usr/bin/env python3
"""Refresh a completed comparison: concise nearest references and expanded labels.

Requires saved ani_pairs.csv, network_nodes.csv and layout_diagnostics.csv.
No Prodigal, DIAMOND, downloads or ANI calculations are run.
Keep this helper beside the updated spraynpray_compare.py.
"""

import argparse
import sys
from pathlib import Path

from spraynpray_compare import refresh_existing_analysis


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, required=True, help="Existing spraynpray_compare output folder")
    parser.add_argument("--summary-only", action="store_true", help="Only write nearest_references.csv; no plotting dependencies needed")
    args = parser.parse_args(argv)
    try:
        result = refresh_existing_analysis(args.analysis_dir, args.summary_only)
    except (ValueError, OSError, ImportError, KeyError) as exc:
        parser.exit(1, str(exc) + "\n")
    print("Wrote nearest_references.csv for %s samples" % len(result))
    if not args.summary_only:
        print("Updated HTML/SVG reference labels using the saved genome coordinates")
    print("Existing ANI measurements were reused; no analyses were rerun")
    return 0


if __name__ == "__main__":
    sys.exit(main())
