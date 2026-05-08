#!/usr/bin/env python3

"""
Description: Script to check for FASTQ file integrity and log errors.

Usage: ./fairy_proc.py -f fastq_file -p prefix [-r input_read] [-b] [-v phx_version] [-V]

v1.0.0 (07/13/2023)
v1.1.0 (10/13/2023)
v2.0   (11/15/2023)
v2.0.1 - JH fixed R1/R2 parsing
Converted to Python: (05/06/2025)

Created by Maria Diaz (lex0@cdc.gov), additions by Jill Hagey (qpk9@cdc.gov)
"""

import argparse
import gzip
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

__version__ = "2.0.1"

# ---------------------------------------------------------------------------
# Synopsis failure entries — (label, detail) driven by a list so there are
# no 30 individual printf calls. Entries marked True get the sample_name
# interpolated into the detail string.
# ---------------------------------------------------------------------------

SYNOPSIS_FAILURES = [
    ("FASTQs",                    "{s}_raw_read_counts reads QC file does not exist -- FASTQs couldn't be unzipped!"),
    ("RAW_READ_COUNTS",           "{s}_raw_read_counts reads QC file does not exist -- FASTQs couldn't be unzipped!"),
    ("RAW_Q30_R1%",               "{s}_raw_read_counts.txt not found -- FASTQs couldn't be unzipped!"),
    ("RAW_Q30_R2%",               "{s}_raw_read_counts.txt not found -- FASTQs couldn't be unzipped!"),
    ("TRIMMED_R1",                "No trimming performed -- FASTQs couldn't be unzipped!"),
    ("TRIMMED_R2",                "No trimming performed -- FASTQs couldn't be unzipped!"),
    ("TRIMMED_FASTQs",            "Trimmed reads QC file does not exist -- FASTQs couldn't be unzipped!"),
    ("TRIMMED_READ_COUNTS",       "Trimmed reads QC file does not exist -- FASTQs couldn't be unzipped!"),
    ("TRIMMED_Q30_R1%",           "No trimming performed -- FASTQs couldn't be unzipped!"),
    ("TRIMMED_Q30_R2%",           "No trimming performed -- FASTQs couldn't be unzipped!"),
    ("KRAKEN2_READS",             "{s}.kraken2_trimd.report.txt not found -- FASTQs couldn't be unzipped!"),
    ("KRONA_READS",               "kraken2 reads do not exist -- FASTQs couldn't be unzipped!"),
    ("KRAKEN2_CLASSIFY_READS",    "There are no classified reads -- FASTQs couldn't be unzipped!"),
    ("QC_COUNTS",                 "FASTQs couldn't be unzipped!"),
    ("Q30_STATS",                 "FASTQs couldn't be unzipped!"),
    ("BBDUK",                     "FASTQs couldn't be unzipped!"),
    ("TRIMMING",                  "FASTQs couldn't be unzipped!"),
    ("KROFAILED_READS",           "FASTQs couldn't be unzipped!"),
    ("ASSEMBLY",                  "{s}.scaffolds.fa.gz not found -- FASTQs couldn't be unzipped!"),
    ("SRST2",                     "SPAdes not performed -- FASTQs couldn't be unzipped!"),
    ("SCAFFOLD_TRIM",             "{s}.filtered.scaffolds.fa.gz not found -- FASTQs couldn't be unzipped!"),
    ("KRAKEN2_ASMBLD",            "{s}.kraken2_asmbld.report.txt not found -- FASTQs couldn't be unzipped!"),
    ("KRONA_ASMBLD",              "kraken2 unweighted not performed -- FASTQs couldn't be unzipped!"),
    ("KRAKEN2_CLASSIFY_ASMBLD",   "kraken2 assembly not performed -- FASTQs couldn't be unzipped!"),
    ("KRAKEN2_ASMBLD_CONTAM",     "kraken2 assembly not performed -- FASTQs couldn't be unzipped!"),
    ("KRAKEN2_WEIGHTED",          "kraken2 weighted not performed -- FASTQs couldn't be unzipped!"),
    ("KRONA_WEIGHTED",            "kraken2 weighted not performed -- FASTQs couldn't be unzipped!"),
    ("KRAKEN2_CLASSIFY_WEIGHTED", "kraken2 weighted not performed -- FASTQs couldn't be unzipped!"),
    ("KRAKEN2_WEIGHTED_CONTAM",   "kraken2 weighted not performed -- FASTQs couldn't be unzipped!"),
    ("QUAST",                     "QUAST not performed -- FASTQs couldn't be unzipped!"),
    ("TAXA-UNK",                  "No Taxa File found -- FASTQs couldn't be unzipped!"),
    ("ASSEMBLY_RATIO(SD)",        "No Ratio File exists -- FASTQs couldn't be unzipped!"),
    ("COVERAGE",                  "No trimmed reads to review for coverage -- FASTQs couldn't be unzipped!"),
    ("BUSCO",                     "BUSCO was not performed -- FASTQs couldn't be unzipped!"),
    ("FASTANI_REFSEQ",            "FASTANI was not performed -- FASTQs couldn't be unzipped!"),
    ("MLST",                      "MLST was not performed -- FASTQs couldn't be unzipped!"),
    ("GAMMA_AR",                  "GAMMA_AR was not performed -- FASTQs couldn't be unzipped!"),
    ("PLASMID_REPLICONS",         "GAMMA was not performed -- FASTQs couldn't be unzipped!"),
    ("HYPERVIRULENCE",            "GAMMA was not performed -- FASTQs couldn't be unzipped!"),
    ("Auto PASS/FAIL",            "FASTQs couldn't be unzipped!"),
]

# TSV column headers differ only by the BUSCO columns
SUMMARY_HEADERS_BUSCO = [
    "WGS_ID", "PHX_Version", "Auto_QC_Outcome", "Warning_Count",
    "Estimated_Coverage", "Genome_Length", "Assembly_Ratio_(STDev)",
    "#_of_Scaffolds_>500bp", "GC_%", "BUSCO", "BUSCO_DB",
    "Final_Taxa_ID", "Taxa_Source", "FastANI_Organism", "FastANI_%ID",
    "FastANI_%Coverage", "ShigaPass_Organism", "Kraken2_Trimd",
    "Kraken2_Weighted", "MLST_Scheme_1", "MLST_1", "MLST_Scheme_2",
    "MLST_2", "GAMMA_Beta_Lactam_Resistance_Genes", "GAMMA_Other_AR_Genes",
    "AMRFinder_Point_Mutations", "Hypervirulence_Genes",
    "Plasmid_Incompatibility_Replicons", "Auto_QC_Failure_Reason",
]
SUMMARY_HEADERS_NO_BUSCO = [h for h in SUMMARY_HEADERS_BUSCO
                             if h not in ("BUSCO", "BUSCO_DB")]


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Check FASTQ file integrity and log errors for PHoeNIx."
    )
    p.add_argument("-f", "--fname",        required=True,  help="FASTQ file to check")
    p.add_argument("-p", "--prefix",       required=True,  help="Sample prefix (meta.id)")
    p.add_argument("-r", "--input-read",   default="",     help="Fallback read orientation (R1/R2)")
    p.add_argument("-b", "--busco",        action="store_true", help="Include BUSCO columns")
    p.add_argument("-v", "--phx-version",  default="",     help="PHoeNIx version string")
    p.add_argument("-V", "--version",      action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


# ---------------------------------------------------------------------------
# gzip integrity check  (replaces `gzip -t`, no subprocess needed)
# ---------------------------------------------------------------------------

def check_gzip_integrity(fname: str) -> tuple[bool, str]:
    """
    Attempt to fully read the gzip file to verify integrity.
    Returns (is_corrupt, error_message).
    Mirrors `gzip -t` behaviour — reads and discards all decompressed bytes.
    """
    try:
        with gzip.open(fname, "rb") as fh:
            while fh.read(1 << 20):   # read in 1 MB chunks
                pass
        return False, ""
    except (gzip.BadGzipFile, EOFError, OSError) as e:
        return True, str(e)


# ---------------------------------------------------------------------------
# R1/R2 detection  (three-tier fallback, mirrors original bash logic)
# ---------------------------------------------------------------------------

def detect_read_orientation(fname: str, full_name: str,
                            input_read: str) -> str:
    """
    Try to determine R1 or R2 in priority order:
      1. Parse the first FASTQ read header for the Illumina '1:N:' / '2:N:' pattern.
      2. Extract R1/R2 from the filename (removing SRR prefix first).
      3. Extract a bare 1/2 from the filename.
      4. Fall back to input_read.
    """
    # Tier 1: read header
    try:
        with gzip.open(fname, "rt") as fh:
            header = fh.readline()
        m = re.search(r"([12]):[NY]:", header)
        if m:
            read = f"R{m.group(1)}"
            if read in ("R1", "R2"):
                return read
    except Exception:
        pass

    print("read orientation not captured trying another method.")

    # Tier 2: R1/R2 token in filename (strip SRR prefix first)
    clean_name = re.sub(r"SRR", "", full_name)
    m = re.search(r"(?<![a-zA-Z0-9])(R[12])(?![a-zA-Z0-9])", clean_name)
    if m and m.group(1) in ("R1", "R2"):
        return m.group(1)

    # Tier 3: bare 1 or 2 in full path
    m = re.search(r"(?<![a-zA-Z0-9])([12])(?![a-zA-Z0-9])", fname)
    if m:
        read = f"R{m.group(1)}"
        if read in ("R1", "R2"):
            return read

    print("read orientation could not be captured correctly.")
    return input_read


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_synopsis_failure(prefix: str, sample_name: str,
                           read: str, status: str = "FAILED") -> None:
    out = Path(f"{sample_name}.synopsis")
    today = datetime.now().strftime("%c")

    with out.open("w") as fh:
        fh.write(f"----------Checking {sample_name} for successful completion on ----------\n")
        fh.write(f"{'Summarized':<30}: {'SUCCESS':<8} : {today}\n")

        for label, detail in SYNOPSIS_FAILURES:
            detail = detail.format(s=sample_name)
            fh.write(f"{label:<30}: {'FAILED':<8} : {detail}\n")

        fh.write(f"---------- {sample_name} completed as {status} ----------\n")
        fh.write("File corruption detected in one or more isolate FASTQs. "
                 "Please check the FAIry result folder for details.\n")
        fh.write("PHoeNIx will only analyze FASTQ files that can be unzipped without error.\n")
        fh.write("\n")
        fh.write("WARNINGS: out of line with what is expected and MAY cause problems downstream.\n")
        fh.write("ALERT: something to note, does not mean it is a poor-quality assembly.\n")


def write_summaryline_failure(prefix: str, phx_version: str,
                              read: str, busco: bool) -> None:
    """Write the TSV summary line for a failed/corrupt sample."""
    headers = SUMMARY_HEADERS_BUSCO if busco else SUMMARY_HEADERS_NO_BUSCO
    n_unknown = len(headers) - 2   # all columns except WGS_ID and failure reason
    unknowns  = "\t".join(["Unknown"] * (n_unknown - 1))
    failure_reason = f"CANNOT UNZIP FASTQ {prefix}_{read} FILE. CHECK FASTQ FILE(S) FOR CORRUPTION!"

    out = Path(f"{prefix}_summaryline.tsv")
    with out.open("w") as fh:
        fh.write("\t".join(headers) + "\n")
        # WGS_ID + PHX_Version + FAIL + 0 + unknowns + failure_reason
        row = f"{prefix}\t{phx_version}\tFAIL\t0\t{unknowns}\t{failure_reason}"
        fh.write(row)   # tr -d '\n' in original — no trailing newline


def write_corruption_summary(prefix: str, read: str, corrupt: bool,
                             err_msg: str = "") -> None:
    out = Path(f"{prefix}_corruption_summary.txt")
    with out.open("a") as fh:
        if corrupt:
            fh.write(
                f"FAILED CORRUPTION CHECK! CANNOT UNZIP FASTQ FILE. "
                f"CHECK FASTQ FILE {prefix}_{read} FOR CORRUPTION!\n"
            )
            if err_msg:
                fh.write(f"  Error: {err_msg}\n")
        else:
            fh.write(f"PASSED: File {prefix}_{read} is not corrupt.\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    fname  = args.fname
    prefix = args.prefix

    # Derive full_name (strip .fastq.gz or .fq.gz)
    full_name = Path(fname).name
    for sfx in (".fastq.gz", ".fq.gz"):
        if full_name.endswith(sfx):
            full_name = full_name[: -len(sfx)]
            break

    # Detect read orientation
    read = detect_read_orientation(fname, full_name, args.input_read)

    # Check gzip integrity (replaces `gzip -t` + grep for "error"/"unexpected")
    corrupt, err_msg = check_gzip_integrity(fname)

    if corrupt:
        write_corruption_summary(prefix, read, corrupt=True, err_msg=err_msg)
        write_summaryline_failure(prefix, args.phx_version, read, args.busco)
        write_synopsis_failure(prefix, prefix, read)
    else:
        write_corruption_summary(prefix, read, corrupt=False)
        print(f"PASSED: File {prefix}_{read} is not corrupt.")


if __name__ == "__main__":
    main()