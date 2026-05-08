#!/usr/bin/env python3

"""
Description: Script to create empty JSON for single reads info to be used downstream.

Usage: ./create_empty_fastp_json.py -n sample_name [-V]

Modules required: None

Created by Jill Hagey (qpk9@cdc.gov)
Converted to Python: (05/06/2025)
"""

import argparse
import json
import sys

__version__ = "2.0.0"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create an empty fastp JSON for single reads."
    )
    p.add_argument("-n", "--sample-name", required=True,
                   help="Sample name of isolate")
    p.add_argument("-V", "--version", action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    empty = {
        "summary": {
            "after_filtering": {
                "total_reads":       0,
                "total_bases":       0,
                "q20_bases":         0,
                "q30_bases":         0,
                "q20_rate":          0,
                "q30_rate":          0,
                "read1_mean_length": 0,
                "gc_content":        0,
            }
        }
    }

    out = f"{args.sample_name}_singles.fastp.json"
    with open(out, "w") as fh:
        json.dump(empty, fh, indent="\t")
        fh.write("\n")


if __name__ == "__main__":
    main()