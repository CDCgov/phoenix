#!/usr/bin/env python3

"""
Description: Quick and dirty way to run blat with dnaX mode AND to sort and
             trim the output psl files.

Usage: ./blat_toxinotypes.py -i input_assembly -o sample_name -d AA_database -t toxinotype_def_file [-p terra] [-V]

Output location: Varies on contents

Modules required: blat (external binary)

v1.0 (09/08/2020) - original shell script
Converted to Python: (05/06/2025)

Created by Nick Vlachos (nvx4@cdc.gov)
"""
import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Dict, Tuple

__version__ = "1.0"

TERRA_BLAT_PATH = "/opt/conda/envs/phoenix/bin/blat"
DEFAULT_BLAT_PATH = "blat"

HEADER = "ID\tToxinotype\tToxin\tsub-type\tContig\tStart\tStop"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run blat in dnaX mode and parse toxinotype results."
    )
    parser.add_argument("-i", "--input",       required=True, help="Input assembly file")
    parser.add_argument("-o", "--sample-name", required=True, help="Output/sample name prefix")
    parser.add_argument("-d", "--db",          required=True, help="AA database (.fa) to compare against")
    parser.add_argument("-t", "--tox-def-file",required=True, help="Toxinotype definition file")
    parser.add_argument("-p", "--terra",       default="",    help="Set to 'terra' if running on Terra")
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s: {__version__}")
    return parser.parse_args()


def run_blat(blat_bin: str, input_file: Path, db: Path, psl_out: Path) -> None:
    """Run blat in dnaX mode against the AA database."""
    cmd = [
        blat_bin,
        "-q=prot",
        "-t=dnax",
        str(input_file),
        str(db),
        "-minIdentity=100",
        "-noHead",
        str(psl_out),
    ]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"blat exited with code {result.returncode}", file=sys.stderr)
        sys.exit(result.returncode)


def classify_subtype(sub_type: str) -> Optional[str]:
    """
    Return 'Toxin-A', 'Toxin-B', or None based on the sub_type string.
    Mirrors the original bash logic including sordellii group handling.
    """
    first = sub_type[0].upper() if sub_type else ""

    if first == "A":
        return "Toxin-A"
    elif first == "B":
        return "Toxin-B"
    elif first == "S":
        if "sordellii_group" in sub_type or "sordellii_TcsL" in sub_type:
            return "Toxin-B"
        elif "sordelii_TcsH" in sub_type:
            return "Toxin-A"
    return None


def parse_psl(psl_path: Path, sample_name: str) -> Tuple[List[Dict], Dict[str, str]]:
    """
    Parse the blat PSL output. Only records where matches == total_length are kept.

    Returns:
        hits       - list of dicts with keys: type, sub_type, contig, start, stop
        subtypes   - dict with keys 'A' and 'B', values are the last matched sub_type
                     strings (or empty string if not found)
    """
    hits = []
    subtypes = {"A": "", "B": ""}

    if not psl_path.is_file():
        print(f"PSL file not found: {psl_path}", file=sys.stderr)
        return hits, subtypes

    with psl_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 17:
                continue

            try:
                matches      = int(fields[0])
                total_length = int(fields[10])
            except ValueError:
                continue

            print(f"M={matches}:TL-{total_length}")

            if matches != total_length:
                continue

            sub_type = fields[9]
            contig   = fields[13]
            start    = fields[15]
            stop     = fields[16]

            tox_type = classify_subtype(sub_type)
            if tox_type is None:
                continue

            hits.append({
                "type":     tox_type,
                "sub_type": sub_type,
                "contig":   contig,
                "start":    start,
                "stop":     stop,
            })

            if tox_type == "Toxin-A":
                subtypes["A"] = sub_type
            elif tox_type == "Toxin-B":
                subtypes["B"] = sub_type

    return hits, subtypes


def lookup_toxinotype(a_sub: str, b_sub: str, tox_def_path: Path) -> str:
    """
    Look up the toxinotype from the definition file.

    Definition file columns (tab-delimited):
        0: toxinotype
        4: sub-type A
        5: sub-type B

    Returns the toxinotype string, 'Undefined' if not found.
    """
    tmp_a = a_sub if a_sub else "-"
    tmp_b = b_sub if b_sub else "-"

    with tox_def_path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            file_ttype = fields[0]
            file_sub_a = fields[4]
            file_sub_b = fields[5]
            print(f"|{tmp_a}| = |{file_sub_a}| : |{tmp_b}| = |{file_sub_b}|")
            if tmp_a == file_sub_a and tmp_b == file_sub_b:
                return file_ttype

    return "Undefined"


def write_output(
    sample_name: str,
    db_name: str,
    hits: List[Dict],
    a_found: bool,
    b_found: bool,
    toxinotype: str,
) -> None:
    """Write the sorted .tox output file and clean up temporaries."""
    tmx_path = Path(f"{sample_name}_{db_name}.tmx")
    tsx_path = Path(f"{sample_name}_{db_name}.tsx")
    tox_path = Path(f"{sample_name}_{db_name}.tox")

    # Write unsorted hits
    with tmx_path.open("w") as fh:
        for h in hits:
            fh.write(
                f"{sample_name}\t{h['type']}\t{h['sub_type']}"
                f"\t{h['contig']}\t{h['start']}\t{h['stop']}\n"
            )
        if not a_found:
            fh.write(f"{sample_name}\tToxin-A\tN/A\tN/A\tN/A\tN/A\n")
        if not b_found:
            fh.write(f"{sample_name}\tToxin-B\tN/A\tN/A\tN/A\tN/A\n")

    # Sort by column 2 (toxin type) — matches `sort -k2`
    with tmx_path.open() as fh:
        rows = [line for line in fh if line.strip()]
    rows.sort(key=lambda r: r.split("\t")[1] if len(r.split("\t")) > 1 else "")
    with tsx_path.open("w") as fh:
        fh.writelines(rows)

    # Write final .tox file: header + sorted hits + toxinotype
    with tox_path.open("w") as fh:
        fh.write(HEADER + "\n")
        with tsx_path.open() as src:
            fh.write(src.read())
        fh.write(f"Toxinotype:\t{toxinotype}\n")

    # Cleanup temporaries
    if tmx_path.exists():
        tmx_path.unlink()

    if tsx_path.exists():
        tsx_path.unlink()


def main() -> None:
    args = parse_args()

    input_file   = Path(args.input)
    db_file      = Path(args.db)
    tox_def_file = Path(args.tox_def_file)
    sample_name  = args.sample_name

    if not input_file.is_file():
        print("Assembly empty or non-existent, exiting", file=sys.stderr)
        sys.exit(1)
    if not db_file.is_file():
        print("Database empty or non-existent, exiting", file=sys.stderr)
        sys.exit(1)

    blat_bin = TERRA_BLAT_PATH if args.terra == "terra" else DEFAULT_BLAT_PATH
    db_name  = db_file.stem
    psl_out  = Path(f"{sample_name}_{db_name}.psl")

    run_blat(blat_bin, input_file, db_file, psl_out)

    hits, subtypes = parse_psl(psl_out, sample_name)

    a_found = any(h["type"] == "Toxin-A" for h in hits)
    b_found = any(h["type"] == "Toxin-B" for h in hits)
    a_count = sum(1 for h in hits if h["type"] == "Toxin-A")
    b_count = sum(1 for h in hits if h["type"] == "Toxin-B")
    total   = a_count + b_count

    print(f"Count check {a_count},{b_count},{total}")

    if total == 0:
        toxinotype = "No subtypes found"
    elif total <= 2:
        toxinotype = lookup_toxinotype(subtypes["A"], subtypes["B"], tox_def_file)
    else:
        toxinotype = "Too_many_subtypes"

    write_output(sample_name, db_name, hits, a_found, b_found, toxinotype)

    global_end_time = datetime.now().strftime("%m-%d-%Y @ %Hh_%Mm_%Ss")
    print(f"All isolates completed - {global_end_time}")


if __name__ == "__main__":
    main()