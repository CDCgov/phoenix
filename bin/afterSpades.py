#!/usr/bin/env python3

"""
Description: Script to clean up the output of SPAdes.

Usage: ./afterspades.py [-V]

Modules required: None

V1.0 (07/03/2022)
Converted to Python: (05/06/2025)

Created by Jill Hagey (qpk9@cdc.gov)
"""

import argparse
import gzip
import shutil
import sys
from pathlib import Path

__version__ = "2.0"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Clean up SPAdes output files.")
    p.add_argument("-V", "--version", action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


def gz_no_timestamp(src: Path) -> None:
    """
    Gzip a file in-place without embedding a timestamp,
    matching `gzip -n` behaviour.
    """
    dest = src.with_suffix(src.suffix + ".gz")
    with src.open("rb") as f_in, gzip.GzipFile(dest, "wb", mtime=0) as f_out:
        shutil.copyfileobj(f_in, f_out)
    src.unlink()


def append_outcome(csv_path: Path, value: str) -> None:
    """Append ,value with no trailing newline — mirrors `echo ,x | tr -d '\n' >>`."""
    with csv_path.open("a") as fh:
        fh.write(f",{value}")


def main() -> None:
    args = parse_args()

    # Locate the .spades.log file and derive the sample prefix
    logs = list(Path(".").glob("*.spades.log"))
    if not logs:
        print("No .spades.log file found.", file=sys.stderr)
        sys.exit(1)
    prefix  = logs[0].stem.removesuffix(".spades")
    csv_out = Path(f"{prefix}_spades_outcome.csv")

    # scaffolds.fasta
    if Path("scaffolds.fasta").is_file():
        renamed = Path(f"{prefix}.scaffolds.fa")
        Path("scaffolds.fasta").rename(renamed)
        gz_no_timestamp(renamed)
        append_outcome(csv_out, "scaffolds_created")
    else:
        append_outcome(csv_out, "no_scaffolds")
        old = Path(f"{prefix}_summary_old_3.txt")
        if old.is_file():
            old.rename(f"{prefix}_trimstats_summary.txt")

    # contigs.fasta
    if Path("contigs.fasta").is_file():
        renamed = Path(f"{prefix}.contigs.fa")
        Path("contigs.fasta").rename(renamed)
        gz_no_timestamp(renamed)
        append_outcome(csv_out, "contigs_created")
    else:
        append_outcome(csv_out, "no_contigs")
        old = Path(f"{prefix}_summary_old_3.txt")
        if old.is_file():
            old.rename(f"{prefix}_trimstats_summary.txt")

    # assembly graph (optional — no outcome logged if missing)
    gfa = Path("assembly_graph_with_scaffolds.gfa")
    if gfa.is_file():
        renamed = Path(f"{prefix}.assembly.gfa")
        gfa.rename(renamed)
        gz_no_timestamp(renamed)


if __name__ == "__main__":
    main()