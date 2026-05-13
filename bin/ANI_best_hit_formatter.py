#!/usr/bin/env python3

"""
Description: Script to format ANI output to include more information on the top line.

Usage: ./ANI_best_hit_formatter.py -a ani_file -n sample_name [-d db_name] [-t terra] [-V]

Output location: /sample_name/fastANI/

Modules required: None

Created by Nick Vlachos (nvx4@cdc.gov)
Converted to Python: (05/06/2025)
"""

import argparse
import re
import sys
from datetime import datetime
from pathlib import Path

__version__ = "2.0"
ANI_THRESHOLD = 80.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Format ANI output to include more information on the top line."
    )
    parser.add_argument("-a", "--ani-file",   required=True,  help="ANI input file")
    parser.add_argument("-n", "--sample-name", required=True, help="Sample name")
    parser.add_argument("-d", "--db-name",    default="",     help="Database name")
    parser.add_argument("-t", "--terra",      default="",     help="Set to 'terra' if running on Terra")
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s: {__version__}")
    return parser.parse_args()


def parse_genus_species(filename: str) -> tuple[str, str]:
    """
    Extract genus and species from an ANI reference filename.

    Handles three filename formats:
      1. Uncultured_Genus_species[_more]_GCF...
      2. Genus_species[_more]_GCF...
      3. Fallback: Genus_species (no GCF token present)

    Returns:
        (genus, species) — species may be multi-word with spaces.
    """
    # Format 1: Uncultured_Genus_species_GCF...
    m = re.match(r"^Uncultured_(.+?)_GCF", filename)
    if m:
        after = m.group(1)
        parts = after.split("_")
        genus = parts[0]
        species = " ".join(parts[1:]) if len(parts) > 1 else ""
        return genus, species

    # Format 2: Genus_species_GCF...
    m = re.match(r"^([^_]+)_(.+?)_GCF", filename)
    if m:
        genus = m.group(1)
        species = m.group(2).replace("_", " ")
        return genus, species

    # Format 3: Fallback — no GCF token
    parts = filename.split("_")
    genus = parts[0] if len(parts) > 0 else filename
    species = parts[1] if len(parts) > 1 else ""
    return genus, species


def load_ani_file(ani_path: Path) -> list[dict]:
    """
    Parse tab-delimited ANI file into a list of hit dicts, sorted by ANI value descending.

    Expected columns (0-indexed):
        0: query assembly path
        1: reference assembly path
        2: ANI value
        3: matching fragment count
        4: total fragment count
    """
    hits = []
    with ani_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            hits.append({
                "query":            fields[0],
                "reference":        fields[1],
                "ani":              float(fields[2]),
                "fragment_matches": int(fields[3]),
                "total_fragments":  int(fields[4]),
            })

    hits.sort(key=lambda h: h["ani"], reverse=True)
    return hits


def write_sorted_ani(ani_path: Path, hits: list[dict]) -> Path:
    """Write sorted ANI hits back to a companion .sorted.txt file."""
    sorted_path = ani_path.with_suffix("").with_suffix(".sorted.txt")
    with sorted_path.open("w") as fh:
        for h in hits:
            fh.write(
                f"{h['query']}\t{h['reference']}\t{h['ani']}"
                f"\t{h['fragment_matches']}\t{h['total_fragments']}\n"
            )
    return sorted_path


def write_output(out_path: Path, best: dict, ref_filename: str, genus: str, species: str) -> None:
    """Write the formatted FastANI summary file."""
    best_percent  = round(best["ani"], 2)
    best_coverage = round(100 * best["fragment_matches"] / best["total_fragments"], 2)
    organism      = f"{genus} {species}".strip()

    with out_path.open("w") as fh:
        fh.write("% ID\t% Coverage\tOrganism\tSource File\n")
        fh.write(f"{best_percent}\t{best_coverage}\t{organism}\t{ref_filename}\n")


def main() -> None:
    args = parse_args()

    ani_path = Path(args.ani_file)
    if not ani_path.is_file():
        print("ani file does not exist, exiting")
        sys.exit(1)

    out_filename = f"{args.sample_name}_{args.db_name}.fastANI_initial.txt"
    out_path     = Path(out_filename)

    # Load and check the top hit against the ANI threshold
    hits = load_ani_file(ani_path)
    if not hits or hits[0]["ani"] < ANI_THRESHOLD:
        out_path.write_text(
            "Mash/FastANI Error: No hits above an ANI value >=80%\n"
        )
        print("Mash/FastANI Error: No hits above an ANI value >=80%")
        sys.exit(0)

    # Write sorted ANI file
    write_sorted_ani(ani_path, hits)

    # Parse the best hit
    best         = hits[0]
    ref_filename = Path(best["reference"]).name
    print(ref_filename)

    genus, species = parse_genus_species(ref_filename)

    write_output(out_path, best, ref_filename, genus, species)

    end_time = datetime.now().strftime("%m-%d-%Y_at_%Hh_%Mm_%Ss")
    print(f"ENDed ANI at {end_time}")


if __name__ == "__main__":
    main()