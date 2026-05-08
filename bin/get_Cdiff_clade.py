#!/usr/bin/env python3

"""
Description: Pull out the clade of a C. diff sample from the pubmlst profile lists.

Usage: ./get_Cdiff_clade.py -m mlst_combined_file -r path_to_mlst_database_folder [-V]

Output: sample_name_cdifficile_clade.tsv

Modules required: None

v1.1 (09/27/2023)
Converted to Python: (05/06/2025)

Created by Nick Vlachos (nvx4@cdc.gov)
"""

import argparse
import sys
from pathlib import Path

__version__ = "1.1"

# Types that indicate an ambiguous / unresolved MLST result
UNDEFINED_TYPES = {"A-SUB", "P-SUB", "novel_allele", "novel_profile"}


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Determine C. difficile clade from MLST combined file."
    )
    p.add_argument("-m", "--input",    required=True,
                   help="MLST combined file")
    p.add_argument("-r", "--mlst-db",  required=True,
                   help="Path to MLST database folder")
    p.add_argument("-V", "--version",  action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Type resolution
# ---------------------------------------------------------------------------

def is_undefined(t: str) -> bool:
    return t in UNDEFINED_TYPES


def resolve_type_single(raw_type: str) -> str:
    """Handle the simple single-data-line case."""
    return "MLST_type_not_defined" if is_undefined(raw_type) else raw_type


# Source priority: assembly > reads  (higher = better)
SOURCE_PRIORITY = {"assembly": 2, "reads": 1, "none": 0}


def resolve_type_multi(lines: list[str]) -> str:
    """
    Select the best MLST type from multiple data lines.

    Priority rules (mirrors original bash logic):
      1. A definitive assembly type wins immediately — break on first find.
      2. A definitive reads type beats an undefined assembly type.
      3. A definitive reads type does NOT overwrite a definitive reads type
         already found (first-found wins among equals).
      4. An undefined type from any source sets "MLST_type_not_defined" only
         if nothing better has been found.
    """
    best_type   = "MLST_type_not_found"
    best_source = "none"

    for line in lines:
        fields      = line.split("\t")
        src_raw     = fields[1].strip() if len(fields) > 1 else ""
        type_temp   = fields[4].strip() if len(fields) > 4 else ""

        # Normalise source to "assembly" or "reads"
        if src_raw.startswith("assembly"):
            src = "assembly"
        elif src_raw == "reads":
            src = "reads"
        else:
            continue

        definitive = not is_undefined(type_temp)

        if definitive:
            if src == "assembly":
                # A definitive assembly answer is the best possible — take it and stop
                return type_temp
            elif src == "reads":
                # Definitive reads answer: only update if we don't already have
                # a definitive reads answer or any assembly answer
                if best_source == "none" or (
                    best_source == "reads" and is_undefined(best_type)
                ) or (
                    best_source == "assembly" and best_type == "MLST_type_not_defined"
                ):
                    best_type   = type_temp
                    best_source = "reads"
        else:
            # Undefined type — only record if we have nothing yet
            if best_source == "none":
                best_type   = "MLST_type_not_defined"
                best_source = src
            elif best_source == src == "assembly" and best_type == "MLST_type_not_defined":
                pass   # already undefined from assembly — no change
            elif best_source == src == "reads" and best_type not in UNDEFINED_TYPES:
                pass   # already have a real reads answer — keep it

    return best_type


def select_mlst_type(input_path: Path) -> str:
    """
    Read the MLST file and return the best resolved type string.
    Fixed: original used `${line_Count}` (capital C) — always unset — so the
    <2 branch never triggered. Corrected here.
    """
    lines = [l.rstrip("\n") for l in input_path.read_text().splitlines() if l.strip()]

    # lines[0] is the header; data starts at lines[1]
    data_lines = lines[1:]
    line_count  = len(lines)   # total including header

    if line_count == 2:
        # Exactly one data line
        fields = data_lines[0].split("\t")
        return resolve_type_single(fields[4].strip() if len(fields) > 4 else "")

    elif line_count < 2:
        # Header only or empty
        return "MLST_type_not_found"

    else:
        return resolve_type_multi(data_lines)


# ---------------------------------------------------------------------------
# Clade lookup
# ---------------------------------------------------------------------------

def lookup_clade(mlst_type: str, mlst_db: str) -> str:
    """
    Look up the clade for a given MLST type in the pubmlst cdifficile profile.
    Fixed: original used '%{profile_line}' (% instead of $) so the existence
    check always passed regardless of whether a match was found.
    """
    profile_path = Path(mlst_db) / "pubmlst" / "cdifficile" / "cdifficile.txt"
    if not profile_path.is_file():
        return "Undefined"

    with profile_path.open() as fh:
        for line in fh:
            if line.startswith(f"{mlst_type}\t"):
                fields = line.split("\t")
                clade  = fields[8].strip() if len(fields) > 8 else ""
                return clade if clade else "Undefined"

    return "Undefined"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    if not input_path.is_file():
        print("Empty sample path supplied (-m), exiting", file=sys.stderr)
        sys.exit(1)

    # Derive sample name: reverse path, take first component, drop leading field
    # Mirrors: rev | cut -d'/' -f1 | cut -d'_' -f2- | rev
    stem        = input_path.name
    parts       = stem.split("_")
    sample_name = "_".join(parts[1:]) if len(parts) > 1 else stem

    mlst_type = select_mlst_type(input_path)
    print(f"Type: {mlst_type}")

    clade = lookup_clade(mlst_type, args.mlst_db)
    print(f"Clade: {clade}")

    out = Path(f"{sample_name}_cdifficile_clade.tsv")
    with out.open("w") as fh:
        fh.write("Sample name\tClade\tMLST\n")
        fh.write(f"{sample_name}\t{clade}\t{mlst_type}\n")


if __name__ == "__main__":
    main()