#!/usr/bin/env python3

"""
Description: Compare local assembly to expected assembly size based upon
             taxonomy file or directly given genus/species.

Usage: ./calculate_assembly_ratio.py -d db_file -q report.tsv -x tax_file
           -s sample_name [-f "genus species"] [-t terra] [-V]

Output: sample_name_Assembly_ratio_DATE.txt
        sample_name_GC_content_DATE.txt

Modules required: None

Created by Nick Vlachos (nvx4@cdc.gov)
Converted to Python: (05/06/2025)
"""

import argparse
import sys
from pathlib import Path

__version__ = "2.0"

MIN_REFERENCES = 10


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Calculate assembly ratio vs. NCBI expected length."
    )
    p.add_argument("-d", "--db-path",      required=True,
                   help="Path to sorted NCBI ratio database file")
    p.add_argument("-q", "--quast-report", required=True,
                   help="QUAST report.tsv file")
    p.add_argument("-x", "--tax-file",     default="",
                   help="Tax file from determine_taxID.sh")
    p.add_argument("-f", "--force",        default="",
                   metavar="\"genus species\"",
                   help="Manually specify taxonomy to compare against")
    p.add_argument("-s", "--sample-name",  required=True)
    p.add_argument("-t", "--terra",        default="")
    p.add_argument("-V", "--version",      action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


# ---------------------------------------------------------------------------
# File parsers
# ---------------------------------------------------------------------------

def parse_quast(quast_report: str) -> tuple[str, str]:
    """
    Extract assembly length (line 16) and GC% (line 17) from a QUAST tsv.
    Returns (assembly_length_str, gc_percent_str) — "NA" on failure.
    """
    try:
        lines = Path(quast_report).read_text().splitlines()
        def field(n):   # 1-based, tab-delimited, third field
            try:
                return lines[n - 1].split("\t")[1].strip()
            except IndexError:
                return "NA"
        return field(16), field(17)
    except (OSError, IndexError):
        return "NA", "NA"


def parse_tax_file(tax_file: str) -> tuple[str, str]:
    """
    Extract genus (line 7) and species (line 8) from a .tax file.
    Returns (genus, species).
    """
    lines = Path(tax_file).read_text().splitlines()
    def taxfield(n):
        try:
            val = lines[n - 1].split("\t")[1].strip()
            return val if val else None
        except IndexError:
            return None
    genus   = taxfield(7) or "No genus found"
    species = taxfield(8) or "No species found"
    return genus, species


def parse_db_line(line: str) -> list[str]:
    """
    Parse one tab-delimited database line, uppercasing the first character
    of the organism name and stripping brackets — replacing the sed calls
    that wrote db_path_update.txt in the original.
    """
    fields = line.rstrip("\n").split("\t")
    if fields:
        name = fields[0].lstrip()
        name = name.replace("[", "").replace("]", "")
        if name:
            name = name[0].upper() + name[1:]
        fields[0] = name
    return fields


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_ratio(path: str, total_tax: str, taxid: str, stdev,
                stdevs, assembly_length, expected_length, ratio) -> None:
    with open(path, "w") as fh:
        fh.write(f"Tax: {total_tax}\n")
        fh.write(f"NCBI_TAXID: {taxid}\n")
        fh.write(f"Species_St.Dev: {stdev}\n")
        fh.write(f"Isolate_St.Devs: {stdevs}\n")
        fh.write(f"Actual_length: {assembly_length}\n")
        fh.write(f"Expected_length: {expected_length}\n")
        fh.write(f"Ratio: {ratio}\n")


def write_gc(path: str, total_tax: str, taxid: str,
             gc_stdev, gc_min, gc_max, gc_mean,
             gc_count, sample_gc: str) -> None:
    with open(path, "w") as fh:
        fh.write(f"Tax: {total_tax}\n")
        fh.write(f"NCBI_TAXID: {taxid}\n")
        fh.write(f"Species_GC_StDev: {gc_stdev}\n")
        fh.write(f"Species_GC_Min: {gc_min}\n")
        fh.write(f"Species_GC_Max: {gc_max}\n")
        fh.write(f"Species_GC_Mean: {gc_mean}\n")
        fh.write(f"Species_GC_Count: {gc_count}\n")
        fh.write(f"Sample_GC_Percent: {sample_gc}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    # Default values
    taxid           = "NA"
    stdev           = "NA"
    stdevs          = "NA"
    assembly_length = "NA"
    expected_length = "NA"
    total_tax       = "NA"
    sample_gc       = "NA"
    gc_stdev = gc_min = gc_max = gc_mean = gc_count = "NA"

    # Database date from filename
    db_path = Path(args.db_path)
    if db_path.is_file():
        ncbi_ratio_date = db_path.stem.rsplit("_", 1)[-1] if "_" in db_path.stem else "19991231"
    else:
        ncbi_ratio_date = "19991231"

    ratio_out = f"{args.sample_name}_Assembly_ratio_{ncbi_ratio_date}.txt"
    gc_out    = f"{args.sample_name}_GC_content_{ncbi_ratio_date}.txt"

    # QUAST
    print(f"Checking if quast Assembly_stats exists: {args.quast_report}")
    if Path(args.quast_report).is_file():
        assembly_length, sample_gc = parse_quast(args.quast_report)
    else:
        print("No quast exists, cannot continue")
        write_ratio(ratio_out, total_tax, taxid, stdev, stdevs,
                    assembly_length, expected_length, -2)
        write_gc(gc_out, "No genus Found\tNo species found", "No Match Found",
                 "No Match Found", "No Match Found", "No Match Found",
                 "No Match Found", "No Match Found", "No Match Found")
        sys.exit(0)

    print(f"Assembly length: {assembly_length}\nGC%: {sample_gc}")

    # Taxonomy — forced or from file
    if args.force:
        parts   = args.force.split(maxsplit=1)
        genus   = parts[0].capitalize()
        species = parts[1].lower() if len(parts) > 1 else ""
        total_tax = f"{genus} {species}\t(selected manually)"
    else:
        print(f"Checking if Tax summary exists: {args.tax_file}")
        if not args.tax_file or not Path(args.tax_file).is_file():
            print("No Tax file to find accession for lookup, exiting")
            sys.exit(1)
        genus, species = parse_tax_file(args.tax_file)
        total_tax = f"{genus} {species}"

    # Database lookup
    if not db_path.is_file():
        print("No ratio DB file found, exiting")
        sys.exit(1)

    with db_path.open() as fh:
        first_line = True
        for raw_line in fh:
            if first_line:         # skip header
                first_line = False
                continue

            arr = parse_db_line(raw_line)
            if not arr or not arr[0]:
                continue

            db_name = arr[0]
            print(f"|{genus} {species}| vs |{db_name}|")

            # Case-insensitive match
            if f"{genus} {species}".lower() == db_name.lower():
                # taxID
                taxid = arr[19] if len(arr) > 19 else "NA"
                if taxid == "-2":
                    taxid = "No mode available when determining tax id"
                elif taxid == "-1":
                    taxid = "No tax id given or empty when making lookup"

                # Assembly length and stdev (values in DB are in Mbp — scale to bp)
                expected_length = int(float(arr[4]) * 1_000_000) if len(arr) > 4 else 0
                reference_count = int(arr[6]) if len(arr) > 6 else 0
                stdev_raw       = int(float(arr[5]) * 1_000_000) if len(arr) > 5 else 0

                print(f"{arr} - {expected_length} - {stdev_raw} - {reference_count}")

                if reference_count < MIN_REFERENCES:
                    stdev  = "Not calculated on species with n<10 references"
                    stdevs = "NA"
                else:
                    stdev  = stdev_raw
                    bigger  = max(int(assembly_length), expected_length)
                    smaller = min(int(assembly_length), expected_length)
                    stdevs  = round((bigger - smaller) / stdev_raw, 4) if stdev_raw else "NA"

                # GC content
                gc_min   = arr[7]  if len(arr) > 7  else "NA"
                gc_max   = arr[8]  if len(arr) > 8  else "NA"
                gc_mean  = arr[10] if len(arr) > 10 else "NA"
                gc_count = arr[12] if len(arr) > 12 else "NA"
                gc_stdev = (arr[11] if len(arr) > 11 else "NA") \
                    if (gc_count.lstrip('-').isdigit() and int(gc_count) >= MIN_REFERENCES) \
                    else "Not calculated on species with n<10 references"

                write_gc(gc_out, total_tax, taxid, gc_stdev,
                         gc_min, gc_max, gc_mean, gc_count, sample_gc)
                break

            # Early exit: DB is sorted alphabetically — stop if we've passed the genus
            elif genus[0].lower() < db_name[0].lower():
                break

    print(f"{expected_length}-{assembly_length}")

    # Handle no-match cases
    if expected_length == "NA" or not str(expected_length).strip():
        print("No expected length was found to compare to")
        write_ratio(ratio_out, total_tax, taxid, "NA", "NA",
                    assembly_length, "NA", -1)
        write_gc(gc_out, total_tax, taxid,
                 "No Match Found", "No Match Found", "No Match Found",
                 "No Match Found", "No Match Found", sample_gc)
        sys.exit(0)

    if assembly_length == "NA" or not str(assembly_length).strip():
        print("No assembly length was found to compare with")
        write_ratio(ratio_out, total_tax, taxid, stdev, "NA",
                    "NA", expected_length, -2)
        write_gc(gc_out, total_tax, taxid, gc_stdev,
                 gc_min, gc_max, gc_mean, gc_count, "NA")
        sys.exit(0)

    ratio = round(int(assembly_length) / int(expected_length), 4)

    print(f"Actual - {assembly_length}\nExpected - {expected_length}\n"
          f"Ratio - {ratio}\nSpecies_St.Devs - {stdev}\nIsolate_St.Dev - {stdevs}")

    write_ratio(ratio_out, total_tax, taxid, stdev, stdevs,
                assembly_length, expected_length, ratio)


if __name__ == "__main__":
    main()