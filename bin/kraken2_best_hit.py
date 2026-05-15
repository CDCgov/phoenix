#!/usr/bin/env python3

"""
Description: Grabs the best species match based on %/read hits from the kraken
             tool run. Simplified for nextflow inclusion.

Usage: ./kraken_best_hit.py -i path_to_list_file [-q count_file] -n sample_name [-V]

Output location: same directory as input list file

Modules required: None

Created by Nick Vlachos (nvx4@cdc.gov)
Converted to Python: (05/06/2025)
"""

import argparse
import sys
from pathlib import Path

__version__ = "2.0"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Grab best kraken2 species hit from a summary list file."
    )
    p.add_argument("-i", "--list-file",    required=True, help="Path to kraken2 .summary.txt list file")
    p.add_argument("-q", "--read-file",    default="", help="Count file (reads or contigs)")
    p.add_argument("-n", "--sample-name",  required=True)
    p.add_argument("-V", "--version", action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


def parse_kraken_list(list_path: Path) -> dict:
    """
    Parse a kraken2 summary list file, tracking the best (highest-read)
    representative at each taxonomic level.

    Species selection is intentionally restricted to entries that fall under
    the top genus — matching the original bash logic.

    Returns a dict of all extracted values.
    """
    result = dict(
        unclass_percent=0.0,   unclass_reads=0,
        root_percent=0.0,      classified_reads=0,
        domain="N/A",          domain_percent=0.0,   domain_reads=0,
        phylum="N/A",          phylum_percent=0.0,   phylum_reads=0,
        klass="N/A",           class_percent=0.0,    class_reads=0,
        order="N/A",           order_percent=0.0,    order_reads=0,
        family="N/A",          family_percent=0.0,   family_reads=0,
        top_genus="N/A",       genus_percent=0.0,    genus_reads=0,
        species="N/A",         species_percent=0.0,  species_reads=0,
    )

    current_genus = "N/A"   # tracks genus of the current block (may not be top)

    with list_path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue

            parts = line.split()
            if len(parts) < 4:
                continue

            percent        = float(parts[0])
            reads          = int(parts[1])
            classification = parts[3]
            # Description is everything from index 5 onward (same as bash ${arrLine[*]:5})
            description    = " ".join(parts[5:]).strip() if len(parts) > 5 else ""
            first_desc     = description.split()[0] if description else ""
            desc_cap = description[0].upper() + description[1:] if description else ""

            if classification == "U":
                result["unclass_percent"] = percent
                result["unclass_reads"]   = reads

            elif classification == "R" or (classification == "-" and first_desc == "root"):
                result["classified_reads"] = reads
                result["root_percent"]     = percent

            elif classification == "R2" and reads > result["domain_reads"]:
                result["domain"]         = desc_cap
                result["domain_percent"] = percent
                result["domain_reads"]   = reads

            elif classification == "P" and reads > result["phylum_reads"]:
                result["phylum"]         = desc_cap
                result["phylum_percent"] = percent
                result["phylum_reads"]   = reads

            elif classification == "C" and reads > result["class_reads"]:
                result["klass"]          = desc_cap
                result["class_percent"]  = percent
                result["class_reads"]    = reads

            elif classification == "O" and reads > result["order_reads"]:
                result["order"]          = desc_cap
                result["order_percent"]  = percent
                result["order_reads"]    = reads

            elif classification == "F" and reads > result["family_reads"]:
                result["family"]         = desc_cap
                result["family_percent"] = percent
                result["family_reads"]   = reads

            elif classification == "G" and reads > result["genus_reads"]:
                result["top_genus"]      = desc_cap
                result["genus_percent"]  = percent
                result["genus_reads"]    = reads
                current_genus            = desc_cap

            elif classification == "G":
                # Track current genus even if not the top — needed for species gating
                current_genus = desc_cap

            elif (classification == "S"
                  and reads > result["species_reads"]
                  and current_genus == result["top_genus"]):
                # Species name in kraken output includes genus as first word — strip it
                words = desc_cap.split()
                species_name = " ".join(words[1:]) if len(words) > 1 else desc_cap
                print(f"Old: {result['species']}-{result['species_reads']}")
                result["species"]         = species_name
                result["species_percent"] = percent
                result["species_reads"]   = reads
                print(f"New: {result['species']}-{result['species_reads']}")

    return result


def recalc_weighted_percents(r: dict) -> dict:
    """
    For weighted assembly reports, recalculate all percentages relative to
    (unclassified + root) total — replacing bc/awk with native Python floats.
    """
    total = r["unclass_percent"] + r["root_percent"]
    if total == 0:
        return r

    def pct(val: float) -> float:
        return round((val * 100) / total, 2)

    r["unclass_percent"] = pct(r["unclass_percent"])
    r["root_percent"]    = pct(r["root_percent"])
    r["domain_percent"]  = pct(r["domain_percent"])
    r["phylum_percent"]  = pct(r["phylum_percent"])
    r["class_percent"]   = pct(r["class_percent"])
    r["order_percent"]   = pct(r["order_percent"])
    r["family_percent"]  = pct(r["family_percent"])
    r["genus_percent"]   = pct(r["genus_percent"])
    r["species_percent"] = pct(r["species_percent"])
    return r


def write_summary(sample_name: str, r: dict) -> None:
    out = Path(f"{sample_name}.summary.txt")
    with out.open("w") as fh:
        fh.write("Taxon level\tMatch percentage\tTaxa\n")
        fh.write(f"U: {r['unclass_percent']} unclassified\n")
        fh.write(f"D: {r['domain_percent']} {r['domain']}\n")
        fh.write(f"P: {r['phylum_percent']} {r['phylum']}\n")
        fh.write(f"C: {r['class_percent']} {r['klass']}\n")
        fh.write(f"O: {r['order_percent']} {r['order']}\n")
        fh.write(f"F: {r['family_percent']} {r['family']}\n")
        fh.write(f"G: {r['genus_percent']} {r['top_genus']}\n")
        fh.write(f"s: {r['species_percent']} {r['species']}\n")


def main() -> None:
    args = parse_args()

    list_path = Path(args.list_file)
    if not list_path.is_file():
        print(f"List file not found: {list_path}", file=sys.stderr)
        sys.exit(1)

    r = parse_kraken_list(list_path)

    name = list_path.name
    if "kraken2_trimd.summary.txt" in name:
        print("doing trimd")
    elif "kraken2_asmbld.summary.txt" in name:
        print("doing asmbld")
    elif "kraken2_wtasmbld.summary.txt" in name:
        print("doing weighted")
        r = recalc_weighted_percents(r)

    write_summary(args.sample_name, r)


if __name__ == "__main__":
    main()