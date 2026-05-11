#!/usr/bin/env python3

"""
Description: Checks sample output folders for correct files and tests thresholds
             for passability. Edited for use in PHoeNIx.

Usage: ./pipeline_stats_writer.py -d sample_name [options]

Output location: results/ID/ID.synopsis

Modules required: None

Created by Nick Vlachos (nvx4@cdc.gov), edits by Jill Hagey (qpk9@cdc.gov)
Converted to Python: (05/06/2025)
"""

import argparse
import gzip
import sys
from datetime import datetime
from functools import lru_cache
from pathlib import Path

__version__ = "2.0"

# Thresholds
KRAKEN2_UNCLASS_FLAG          = 30.0
KRAKEN2_TOP_HIT_THRESHOLD     = 70.0
KRAKEN2_CONTAMINATION_THRESH  = 25
ANI_COVERAGE_THRESHOLD        = 90

# Coverage thresholds (may be overridden by --coverage)
READS_MIN_DEFAULT  = 30
READS_LOW          = 40
READS_HIGH         = 150

GAMMA_HEADER = "Gene\tContig\tStart\tStop\tMatch_Type"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Pipeline QC stats writer for PHoeNIx.")
    p.add_argument("-a", "--raw-read-counts")
    p.add_argument("-b", "--total-read-counts")
    p.add_argument("-c", "--gc-content-file")
    p.add_argument("-d", "--sample-name",         required=True)
    p.add_argument("-e", "--kraken2-trimd-report")
    p.add_argument("-f", "--kraken2-trimd-summary")
    p.add_argument("-g", "--krona-trimd")
    p.add_argument("-i", "--filtered-assembly")
    p.add_argument("-H", "--assembly",      dest="spades_assembly")
    p.add_argument("-j", "--kraken2-asmbld-report")
    p.add_argument("-k", "--kraken2-asmbled-summary")
    p.add_argument("-l", "--krona-asmbld")
    p.add_argument("-m", "--kraken2-weighted-report")
    p.add_argument("-n", "--kraken2-weighted-summary")
    p.add_argument("-o", "--krona-weighted")
    p.add_argument("-p", "--quast-report")
    p.add_argument("-q", "--taxid-file")
    p.add_argument("-r", "--assembly-ratio-file")
    p.add_argument("-s", "--busco-summary")
    p.add_argument("-t", "--formatted-fastani")
    p.add_argument("-u", "--gamma-ar")
    p.add_argument("-v", "--gamma-replicon")
    p.add_argument("-w", "--gamma-hv")
    p.add_argument("-x", "--srst2-file")
    p.add_argument("-y", "--mlst-file")
    p.add_argument("-z", "--assembly-only",        action="store_true")
    p.add_argument("-1", "--amr-file")
    p.add_argument("-5", "--coverage",             type=int, default=0)
    p.add_argument("-3", "--cdc-phoenix-mode",     action="store_true")
    p.add_argument("-V", "--version",              action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Core I/O helpers
# ---------------------------------------------------------------------------

def nonempty(path_str: str | None) -> bool:
    """Return True if path_str points to a non-empty file."""
    if not path_str:
        return False
    p = Path(path_str)
    return p.is_file() and p.stat().st_size > 0


@lru_cache(maxsize=32)
def read_tsv_last_line(path_str: str) -> list[str]:
    """
    Return the last non-empty line of a TSV file as a list of fields.
    Result is cached — the same file is never read more than once per run,
    even when multiple checkers need fields from it.
    """
    last: list[str] = []
    with open(path_str) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.strip():
                last = line.split("\t")
    return last


@lru_cache(maxsize=32)
def read_all_lines(path_str: str) -> list[str]:
    """
    Return all non-empty lines of a plain-text file.
    Cached so the same file is never read twice.
    """
    with open(path_str) as fh:
        return [l.rstrip("\n") for l in fh if l.strip()]


def gz_count_headers(path_str: str) -> int:
    """
    Count '>' lines in a (possibly gzipped) FASTA file efficiently.
    Reads raw bytes and checks only the first byte of each line —
    avoids full UTF-8 decoding of sequence data.
    """
    path  = Path(path_str)
    count = 0
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rb") as fh:
        for line in fh:
            if line[0:1] == b">":
                count += 1
    return count


def safe_float(val: str) -> float:
    try:
        return float(val)
    except (ValueError, TypeError):
        return 0.0


def safe_int(val: str | int) -> int:
    try:
        return int(str(val).split(".")[0])
    except (ValueError, TypeError):
        return 0


def status_priority(current: str, new: str) -> str:
    order = {"FAILED": 4, "WARNING": 3, "ALERT": 2, "SUCCESS": 1, "NA": 0}
    return new if order.get(new, 0) > order.get(current, 0) else current


# ---------------------------------------------------------------------------
# Synopsis writer
# ---------------------------------------------------------------------------

class Synopsis:
    def __init__(self, path: Path):
        self._fh   = path.open("w")
        self.status = "SUCCESS"

    def write(self, label: str, result: str, detail: str) -> None:
        self._fh.write(f"{label:<30}: {result:<8} : {detail}\n")

    def update_status(self, new: str) -> None:
        self.status = status_priority(self.status, new)

    def record(self, label: str, result: str, detail: str,
               escalate: bool = True) -> None:
        self.write(label, result, detail)
        if escalate and result in ("FAILED", "WARNING", "ALERT"):
            self.update_status(result)

    def close(self) -> None:
        self._fh.close()


# ---------------------------------------------------------------------------
# Kraken2 summary parser  (cached — called up to 3× on the same files)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=16)
def parse_kraken2_summary(path_str: str) -> dict:
    """
    Parse a kraken2 summary file and return a dict with:
        unclass, unclass_str, domain, genus, species, species_pct, genus_pct.
    Result is cached so repeated calls with the same path are free.
    """
    result = dict(unclass=0.0, unclass_str="0",
                  domain=0.0, genus="", species="",
                  species_pct=0.0, genus_pct=0.0)
    if not nonempty(path_str):
        return result
    try:
        lines = read_all_lines(path_str)   # uses shared cache
        def field(n, col):
            try:
                return lines[n - 1].split()[col - 1]
            except IndexError:
                return "0"
        unclass_raw = field(2, 2)
        if unclass_raw == "UNK":
            result["unclass"]     = 0.0
            result["unclass_str"] = "UNKNOWN"
        else:
            result["unclass"]     = safe_float(unclass_raw)
            result["unclass_str"] = unclass_raw
        result["domain"]      = safe_float(field(3, 2))
        result["genus"]       = field(8, 3)
        result["genus_pct"]   = safe_float(field(8, 2))
        result["species"]     = " ".join(lines[8].split()[2:]) if len(lines) >= 9 else ""
        result["species_pct"] = safe_float(field(9, 2))
    except Exception:
        pass
    return result


# ---------------------------------------------------------------------------
# Kraken2 contamination checker  (cached — report files read only once each)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=16)
def check_kraken2_contamination(report_path: str,
                                weighted: bool = False) -> tuple[int, float]:
    """
    Count genera above the contamination threshold.
    Cached — the report file is read at most once per path+mode combination.
    """
    if not nonempty(report_path):
        return 0, 0.0

    lines        = read_all_lines(report_path)   # uses shared cache
    n_genera     = 0
    total_pct    = 0.0
    root = unclass = 0.0

    for line in lines:
        cols = line.split()
        if len(cols) < 4:
            continue
        level   = cols[3]
        raw_pct = safe_float(cols[0])

        if weighted:
            if level == "U":
                unclass = raw_pct
            elif level == "R":
                root      = raw_pct
                total_pct = unclass + root
            elif level == "G" and total_pct > 0:
                if int((raw_pct * 100) / total_pct) > KRAKEN2_CONTAMINATION_THRESH:
                    n_genera += 1
        else:
            if level == "G":
                if int(str(raw_pct).split(".")[0]) > KRAKEN2_CONTAMINATION_THRESH:
                    n_genera += 1

    return n_genera, total_pct


# ---------------------------------------------------------------------------
# GAMMA file checker
# ---------------------------------------------------------------------------

def count_gamma_hits(gamma_path: str) -> int:
    """Count non-header data lines. Returns -1 if file missing/empty."""
    if not nonempty(gamma_path):
        return -1
    return sum(
        1 for l in read_all_lines(gamma_path)
        if not l.startswith(GAMMA_HEADER) and l.strip()
    )


def check_gamma(syn: Synopsis, label: str, gamma_path: str,
                gene_type: str, sample_name: str, db_tag: str) -> None:
    count = count_gamma_hits(gamma_path)
    if count < 0:
        syn.record(label, "FAILED",
                   f"{sample_name}_{db_tag}.gamma does not exist")
    elif count == 0:
        syn.record(label, "SUCCESS",
                   f"No {gene_type} genes were found", escalate=False)
    else:
        syn.record(label, "SUCCESS",
                   f"{count} {gene_type} gene(s) found", escalate=False)


# ---------------------------------------------------------------------------
# Section checkers
# ---------------------------------------------------------------------------

def check_raw_reads(syn: Synopsis, args) -> tuple[int, int, bool]:
    if not nonempty(args.raw_read_counts):
        syn.record("FASTQs",          "FAILED",
                   f"{args.raw_read_counts} reads QC file does not exist")
        syn.record("RAW_READ_COUNTS", "FAILED",
                   f"{args.raw_read_counts} reads QC file does not exist")
        return 0, 0, False

    last = read_tsv_last_line(args.raw_read_counts)
    r1   = safe_int(last[2])  if len(last) > 2 else 0
    r2   = safe_int(last[4])  if len(last) > 4 else 0

    if r1 > 0 and r2 > 0:
        syn.record("FASTQs", "SUCCESS", f"R1: {r1}bps R2: {r2}bps", escalate=False)
    else:
        if r1 <= 0:
            syn.record("FASTQs_R1", "FAILED",
                       f"R1 not represented correctly in {args.raw_read_counts}")
            r1 = 0
        else:
            syn.record("FASTQs_R1", "SUCCESS", f"{r1}bps", escalate=False)
        if r2 <= 0:
            syn.record("FASTQs_R2", "FAILED",
                       f"R2 not represented correctly in {args.raw_read_counts}")
            r2 = 0
        else:
            syn.record("FASTQs_R2", "SUCCESS", f"{r2}bps", escalate=False)

    return r1, r2, True


def check_raw_counts(syn: Synopsis, args, raw_exists: bool) -> None:
    if not raw_exists:
        for lbl in ("RAW_READ_COUNTS", "RAW_Q30_R1%", "RAW_Q30_R2%"):
            syn.record(lbl, "FAILED",
                       f"{args.sample_name}_raw_read_counts.txt not found")
        return

    last      = read_tsv_last_line(args.raw_read_counts)   # cached — free
    raw_reads = safe_int(last[16]) if len(last) > 16 else 0
    raw_pairs = raw_reads // 2
    q30_r1_s  = last[13] if len(last) > 13 else "0"
    q30_r2_s  = last[14] if len(last) > 14 else "0"
    q30_r1    = safe_int(q30_r1_s.split(".")[1][:2]) if "." in q30_r1_s else 0
    q30_r2    = safe_int(q30_r2_s.split(".")[1][:2]) if "." in q30_r2_s else 0

    if raw_reads <= 0:
        syn.record("RAW_READ_COUNTS", "FAILED",
                   f"No individual read count before trimming: {raw_reads} ({raw_pairs} paired)")
    elif raw_reads <= 1_000_000:
        syn.record("RAW_READ_COUNTS", "WARNING",
                   f"Low individual read count before trimming: {raw_reads} ({raw_pairs} paired)")
    else:
        syn.record("RAW_READ_COUNTS", "SUCCESS",
                   f"{raw_reads} individual reads ({raw_pairs} paired)", escalate=False)

    for pct, thresh, label in ((q30_r1, 90, "RAW_Q30_R1%"), (q30_r2, 70, "RAW_Q30_R2%")):
        if pct < thresh:
            syn.record(label, "WARNING", f"Q30 at {pct}% (Threshold is {thresh}%)")
        else:
            syn.record(label, "SUCCESS",  f"Q30 at {pct}% (Threshold is {thresh}%)",
                       escalate=False)


def check_trimmed_reads(syn: Synopsis, args) -> tuple[int, int, int, bool]:
    if not nonempty(args.total_read_counts):
        syn.record("TRIMMED_FASTQs",      "FAILED",
                   f"{args.total_read_counts} reads QC file does not exist")
        syn.record("TRIMMED_READ_COUNTS", "FAILED",
                   f"{args.total_read_counts} reads QC file does not exist")
        return 0, 0, 0, False

    last     = read_tsv_last_line(args.total_read_counts)
    r1       = safe_int(last[2])  if len(last) > 2 else 0
    r2       = safe_int(last[4])  if len(last) > 4 else 0
    unpaired = safe_int(last[6])  if len(last) > 6 else 0

    if r1 > 0 and r2 > 0 and unpaired > 0:
        syn.record("TRIMMED_BPS", "SUCCESS",
                   f"R1: {r1}bps R2: {r2}bps Unpaired: {unpaired}bps", escalate=False)
    else:
        for val, lbl in ((r1, "TRIMMED_R1"), (r2, "TRIMMED_R2")):
            if val <= 0:
                syn.record(lbl, "FAILED",
                           f"{lbl[-2:]} not represented correctly in {args.total_read_counts}")
            else:
                syn.record(lbl, "SUCCESS", f"{val}bps", escalate=False)
        if unpaired > 0:
            syn.record("TRIMMED_UNPAIRED", "SUCCESS", f"{unpaired}bps", escalate=False)
        else:
            syn.record("TRIMMED_UNPAIRED", "ALERT", "No orphaned reads were found")
            unpaired = 0

    r1 = max(r1, 0)
    r2 = max(r2, 0)
    return r1, r2, unpaired, True


def check_trimmed_counts(syn: Synopsis, args, total_exists: bool) -> tuple[int, int]:
    if not total_exists:
        for lbl in ("TRIMMED_READ_COUNTS", "TRIMMED_Q30_R1%", "TRIMMED_Q30_R2%"):
            syn.record(lbl, "FAILED",
                       f"{args.sample_name}_trimmed_read_counts.txt not found")
        return 0, 0

    last     = read_tsv_last_line(args.total_read_counts)   # cached — free
    total    = safe_int(last[23]) if len(last) > 23 else 0
    orphaned = safe_int(last[5])  if len(last) > 5  else 0
    trimmed  = total - orphaned
    paired   = trimmed // 2
    bps_all  = safe_int(last[21]) if len(last) > 21 else 0

    q30_r1_s = last[18] if len(last) > 18 else "0"
    q30_r2_s = last[19] if len(last) > 19 else "0"
    q30_r1   = safe_int(q30_r1_s.split(".")[1][:2]) if "." in q30_r1_s else 0
    q30_r2   = safe_int(q30_r2_s.split(".")[1][:2]) if "." in q30_r2_s else 0

    if total <= 0:
        syn.record("TRIMMED_READ_COUNTS", "FAILED",
                   f"No individual read count after trimming: {trimmed} ({paired} paired, {orphaned} singled)")
    elif total <= 1_000_000:
        syn.record("TRIMMED_READ_COUNTS", "WARNING",
                   f"Low individual read count after trimming: {total} ({paired} paired, {orphaned} singled)")
    else:
        syn.record("TRIMMED_READ_COUNTS", "SUCCESS",
                   f"{total} individual reads ({paired} paired, {orphaned} singled)", escalate=False)

    for pct, thresh, label in (
        (q30_r1, 90, "TRIMMED_Q30_R1%"),
        (q30_r2, 70, "TRIMMED_Q30_R2%"),
    ):
        if pct < thresh:
            syn.record(label, "WARNING", f"Q30 at {pct}% (Threshold is {thresh}%)")
        else:
            syn.record(label, "SUCCESS",  f"Q30 at {pct}% (Threshold is {thresh}%)",
                       escalate=False)

    return total, bps_all


def _kraken2_contamination_report(syn: Synopsis, report_path: str,
                                  label: str, weighted: bool = False) -> None:
    """Shared contamination reporting logic for all three kraken2 contexts."""
    n_genera, _ = check_kraken2_contamination(report_path, weighted)
    if n_genera > 1:
        syn.record(label, "WARNING",
                   f"{n_genera} genera found above {KRAKEN2_CONTAMINATION_THRESH}% threshold")
    elif n_genera == 1:
        syn.record(label, "SUCCESS",
                   f"Only one genus found above {KRAKEN2_CONTAMINATION_THRESH}% threshold",
                   escalate=False)
    else:
        result = "ALERT" if "ASMBLD" in label else "WARNING"
        syn.record(label, result,
                   f"No genera found above {KRAKEN2_CONTAMINATION_THRESH}% threshold")


def check_kraken2_reads(syn: Synopsis, args) -> None:
    pre_success = nonempty(args.kraken2_trimd_report)
    if not pre_success:
        syn.record("KRAKEN2_READS", "FAILED",
                   f"{args.sample_name}.kraken2_trimd.summary.txt not found")

    if pre_success and not nonempty(args.krona_trimd):
        syn.record("KRONA_READS", "FAILED",
                   f"{args.sample_name}_trimd.html not found")
    elif not pre_success:
        syn.record("KRONA_READS", "FAILED",
                   "kraken2 reads did not complete successfully")

    if nonempty(args.kraken2_trimd_summary):
        k = parse_kraken2_summary(args.kraken2_trimd_summary)
        if k["domain"] <= 0:
            syn.record("KRAKEN2_CLASSIFY_READS", "FAILED",
                       "There are no classified reads" if pre_success
                       else "KRAKEN2_READS did not complete successfully")
        elif k["unclass"] > KRAKEN2_UNCLASS_FLAG:
            syn.record("KRAKEN2_CLASSIFY_READS", "WARNING",
                       f"unclassified reads comprise {k['unclass_str']}% of total")
        else:
            syn.record("KRAKEN2_CLASSIFY_READS", "SUCCESS",
                       f"{k['species_pct']}% {k['genus']} {k['species']} "
                       f"with {k['unclass_str']}% unclassified reads", escalate=False)
    else:
        syn.record("KRAKEN2_CLASSIFY_READS", "FAILED",
                   f"{args.sample_name}.kraken2_trimd.classifiedreads.txt not found")

    if pre_success:
        _kraken2_contamination_report(syn, args.kraken2_trimd_report,
                                      "KRAKEN2_READS_CONTAM")


def check_kraken2_assembly(syn: Synopsis, args) -> None:
    ok = nonempty(args.kraken2_asmbld_report)
    if not ok:
        syn.record("KRAKEN2_ASMBLD", "FAILED",
                   f"{args.sample_name}.kraken2_asmbld.summary.txt not found")

    if ok and not nonempty(args.krona_asmbld):
        syn.record("KRONA_ASMBLD", "FAILED",
                   f"{args.sample_name}_asmbld.html not found")
    elif not ok:
        syn.record("KRONA_ASMBLD", "FAILED",
                   "kraken2 unweighted did not complete successfully")

    if nonempty(args.kraken2_asmbled_summary):
        k = parse_kraken2_summary(args.kraken2_asmbled_summary)
        if k["domain"] <= 0:
            syn.record("KRAKEN2_CLASSIFY_ASMBLD", "FAILED",
                       "There are no classified reads (Did post assembly kraken2 fail too?)"
                       if ok else "kraken2 assembly did not complete successfully")
        elif k["unclass"] > KRAKEN2_UNCLASS_FLAG:
            syn.record("KRAKEN2_CLASSIFY_ASMBLD", "WARNING",
                       f"unclassified scaffolds comprise {k['unclass_str']}% of total")
        elif k["genus_pct"] < 70:
            syn.record("KRAKEN2_CLASSIFY_ASMBLD", "WARNING",
                       f"Genus-{k['genus']}({k['genus_pct']}%) under 70% "
                       f"(species {k['species']} ({k['species_pct']}%)), possibly contaminated")
        else:
            syn.record("KRAKEN2_CLASSIFY_ASMBLD", "SUCCESS",
                       f"{k['genus']}({k['genus_pct']}%) {k['species']}({k['species_pct']}%) "
                       f"with {k['unclass_str']}% unclassified scaffolds", escalate=False)
    else:
        syn.record("KRAKEN2_CLASSIFY_ASMBLD", "FAILED",
                   f"{args.sample_name}.kraken2_asmbld.classifiedreads.txt not found")

    # Fixed: original had 'lassification' typo that silently broke this check
    if ok:
        _kraken2_contamination_report(syn, args.kraken2_asmbld_report,
                                      "KRAKEN2_ASMBLD_CONTAM")


def check_kraken2_weighted(syn: Synopsis, args) -> None:
    ok = nonempty(args.kraken2_weighted_report)
    if not ok:
        syn.record("KRAKEN2_WEIGHTED", "FAILED",
                   f"{args.sample_name}.kraken2_wtasmbld.summary.txt not found")

    if ok and not nonempty(args.krona_weighted):
        syn.record("KRONA_WEIGHTED", "FAILED",
                   f"{args.sample_name}_weighted.html not found")
    elif not ok:
        syn.record("KRONA_WEIGHTED", "FAILED",
                   "kraken2 weighted did not complete successfully")

    if nonempty(args.kraken2_weighted_summary):
        k = parse_kraken2_summary(args.kraken2_weighted_summary)
        if k["domain"] <= 0:
            syn.record("KRAKEN2_CLASSIFY_WEIGHTED", "FAILED",
                       "There are no classified reads" if ok
                       else "Kraken2 weighted did not complete successfully")
        elif k["unclass"] > KRAKEN2_UNCLASS_FLAG:
            syn.record("KRAKEN2_CLASSIFY_WEIGHTED", "WARNING",
                       f"unclassified reads comprise {k['unclass_str']}% of total")
        elif k["genus_pct"] < 70:
            syn.record("KRAKEN2_CLASSIFY_WEIGHTED", "FAILED",
                       f"Genus-{k['genus']} under 70% "
                       f"(species-{k['species']} {k['species_pct']}%), likely contaminated")
        else:
            syn.record("KRAKEN2_CLASSIFY_WEIGHTED", "SUCCESS",
                       f"{k['genus']}({k['genus_pct']}%) {k['species']}({k['species_pct']}%) "
                       f"with {k['unclass_str']}% unclassified scaffolds", escalate=False)
    else:
        syn.record("KRAKEN2_CLASSIFY_WEIGHTED", "FAILED",
                   f"{args.sample_name}.kraken2_wtasmbld.classifiedreads.txt not found")

    if ok:
        _kraken2_contamination_report(syn, args.kraken2_weighted_report,
                                      "KRAKEN2_WEIGHTED_CONTAM", weighted=True)


def check_assembly(syn: Synopsis, args) -> tuple[int, int]:
    full_scaffolds = 0
    full_longies   = 0

    if nonempty(args.spades_assembly):
        full_scaffolds = gz_count_headers(args.spades_assembly)
        syn.record("ASSEMBLY", "SUCCESS",
                   f"{full_scaffolds} scaffolds found", escalate=False)
    else:
        syn.record("ASSEMBLY", "FAILED",
                   f"{args.sample_name}.scaffolds.fa.gz not found")

    if nonempty(args.trimmed_assembly):
        full_longies  = gz_count_headers(args.trimmed_assembly)
        full_shorties = full_scaffolds - full_longies
        if full_longies <= 200:
            syn.record("SCAFFOLD_TRIM", "SUCCESS",
                       f"{full_longies} scaffolds remain. {full_shorties} removed (too short)",
                       escalate=False)
        elif full_longies <= 500:
            syn.record("SCAFFOLD_TRIM", "WARNING",
                       f"{full_longies} scaffolds remain (high). {full_shorties} removed (too short)")
        else:
            syn.record("SCAFFOLD_TRIM", "FAILED",
                       f"{full_longies} scaffolds remain (too high). {full_shorties} removed (too short)")
    else:
        syn.record("SCAFFOLD_TRIM", "FAILED",
                   f"{args.sample_name}.filtered.scaffolds.fa.gz not found")

    return full_scaffolds, full_longies


def check_quast(syn: Synopsis, args) -> tuple[int, str]:
    if not nonempty(args.quast_report):
        syn.record("QUAST",          "FAILED",
                   f"{args.sample_name}_report.tsv does not exist")
        syn.record("QUAST_GC_Content","FAILED",
                   f"{args.sample_name}_report.tsv does not exist")
        return 0, ""

    lines = read_all_lines(args.quast_report)   # cached

    def qfield(row, col, sep="\t"):
        try:
            return lines[row - 1].replace("\t", " ").split()[col - 1]
        except IndexError:
            return ""

    contig_num      = qfield(14, 3)
    assembly_length = safe_int(qfield(16, 3))
    n50             = qfield(18, 2)
    gc_con          = qfield(17, 3)

    syn.record("QUAST", "SUCCESS",
               f"#-{contig_num} length-{assembly_length} n50-{n50} %GC-{gc_con}",
               escalate=False)

    if nonempty(args.gc_content_file):
        gc_lines    = read_all_lines(args.gc_content_file)  # cached
        gc_stdev_s  = gc_lines[2].split()[1] if len(gc_lines) > 2 else "Not"
        if gc_stdev_s.startswith("Not"):
            syn.record("QUAST_GC_Content", "ALERT",
                       f"Low References for STDev - {gc_con}x({gc_stdev_s}-SD)")
        else:
            gc_stdev = safe_float(gc_stdev_s)
            gc_mean  = safe_float(gc_lines[5].split()[1]) if len(gc_lines) > 5 else 0.0
            devs     = 2.58 * gc_stdev
            gc_left  = gc_mean - devs
            gc_right = gc_mean + devs
            gc_val   = safe_float(gc_con)
            if gc_val < gc_left:
                syn.record("QUAST_GC_Content", "WARNING",
                           f"%GC-{gc_con} below {gc_left:.5f} "
                           f"(2.58*{gc_stdev:.5f}SD from mean {gc_mean:.5f})")
            elif gc_val > gc_right:
                syn.record("QUAST_GC_Content", "WARNING",
                           f"%GC-{gc_con} above {gc_right:.5f} "
                           f"(2.58*{gc_stdev:.5f}SD from mean {gc_mean:.5f})")
            else:
                syn.record("QUAST_GC_Content", "SUCCESS",
                           f"%GC-{gc_con} within {gc_left:.5f}-{gc_right:.5f} "
                           f"(2.58*{gc_stdev:.5f}SD from mean {gc_mean:.5f})",
                           escalate=False)
    else:
        syn.record("QUAST_GC_Content", "FAILED",
                   f"%GC-{gc_con}, but GC content file does not exist.")

    return assembly_length, gc_con


def check_taxa(syn: Synopsis, args) -> tuple[str, str]:
    dec_genus   = ""
    dec_species = ""
    tax_source  = "UNK"

    if not nonempty(args.taxid_file):
        syn.record(f"TAXA-{tax_source}", "FAILED", "No Taxa File found")
        return dec_genus, dec_species

    lines = read_all_lines(args.taxid_file)
    if lines:
        tax_source = lines[0].split("\t")[0].strip()
    for line in lines:
        if line.startswith("G:"):
            dec_genus   = line.split("\t")[1].strip()
        elif line.startswith("s:"):
            dec_species = line.split("\t")[1].strip()

    if dec_genus != "Not_assigned" and dec_species != "Not_assigned":
        syn.record(f"TAXA-{tax_source}", "SUCCESS",
                   f"{dec_genus} {dec_species}", escalate=False)
    elif dec_genus != "Not_assigned":
        syn.record("TAXA", "WARNING", "No Species was able to be determined")
    else:
        syn.record("TAXA", "FAILED",
                   "None of the classifiers completed successfully")

    return dec_genus, dec_species


def check_assembly_ratio(syn: Synopsis, args,
                         dec_genus: str, dec_species: str) -> str:
    genus_initial = dec_genus[0] if dec_genus else "?"
    assembly_id   = f"{genus_initial}.{dec_species}"
    qc_fail       = ""

    if not args.assembly_ratio_file or not Path(args.assembly_ratio_file).is_file():
        syn.record("ASSEMBLY_RATIO(SD)", "FAILED", "No Ratio File exists")
        return qc_fail

    lines = read_all_lines(args.assembly_ratio_file)
    assembly_ratio   = safe_float(lines[-1].split()[1]) if lines else 0.0
    stdev_line       = lines[3] if len(lines) > 3 else ""
    species_stdev_ln = lines[2] if len(lines) > 2 else ""
    st_dev_str       = (stdev_line.split()[1]
                        if stdev_line and "N/A" not in stdev_line else "N/A")
    st_dev           = safe_float(st_dev_str) if st_dev_str != "N/A" else 0.0

    if assembly_ratio < 0:
        syn.record("ASSEMBLY_RATIO(SD)", "WARNING",
                   f"No Reference - {assembly_ratio}x({st_dev_str}-SD) against {assembly_id}")
    elif ("Not calculated on species with n<10 references" in species_stdev_ln
          or st_dev_str == "N/A"):
        syn.record("ASSEMBLY_RATIO(SD)", "ALERT",
                   f"Low References for STDev - {assembly_ratio}x({st_dev_str}-SD) against {assembly_id}")
    elif st_dev > 2.58:
        syn.record("ASSEMBLY_RATIO(SD)", "FAILED",
                   f"St. dev. too large - {assembly_ratio}x({st_dev}-SD) against {assembly_id}")
        qc_fail = f"STDev_above_2.58({st_dev})-"
    else:
        syn.record("ASSEMBLY_RATIO(SD)", "SUCCESS",
                   f"{assembly_ratio}x({st_dev}-SD) against {assembly_id}", escalate=False)

    return qc_fail


def check_coverage(syn: Synopsis, args, bps_post_all: int,
                   assembly_length: int, reads_min: int) -> str:
    if not (assembly_length > 0 and bps_post_all > 0):
        return ""

    avg_cov = round(bps_post_all / assembly_length, 2)

    if READS_LOW <= avg_cov < READS_HIGH:
        syn.record("COVERAGE", "SUCCESS",
                   f"{avg_cov}x coverage (Target:40x, Cutoff:{reads_min}x)",
                   escalate=False)
    elif avg_cov >= READS_HIGH:
        syn.record("COVERAGE", "ALERT",
                   f"{avg_cov}x coverage (Target:<150x)")
    elif avg_cov > reads_min:
        syn.record("COVERAGE", "ALERT",
                   f"{avg_cov}x coverage (Target:40x, Cutoff:{reads_min}x)")
    else:
        syn.record("COVERAGE", "FAILED",
                   f"{avg_cov}x coverage (Min:30x)")
        return f"coverage_below_30({avg_cov})-"

    return ""


def check_busco(syn: Synopsis, args) -> None:
    if not nonempty(args.busco_summary):
        syn.record("BUSCO", "FAILED",
                   f"short_summary.*.{args.sample_name}.scaffolds.fa.txt not found")
        return

    found = total = 0
    db    = ""
    for line in read_all_lines(args.busco_summary):
        if "Complete BUSCOs (C)" in line:
            found = safe_int(line.split()[0])
        elif "Total BUSCO groups searched" in line:
            total = safe_int(line.split()[0])
        elif "The lineage dataset is:" in line:
            db = line.split()[5]

    organism = db.split("_odb")[0] if "_odb" in db else db
    pct      = (found * 100 // total) if total else 0
    label    = f"BUSCO_{db.upper()}"

    if pct > 97:
        syn.record(label, "SUCCESS",
                   f"{pct}% core genes for {organism} found ({found}/{total}) (Target:90%)",
                   escalate=False)
    else:
        syn.record(label, "WARNING",
                   f"only {pct}% core genes for {organism} found ({found}/{total}) (Target:90%)")


def check_fastani(syn: Synopsis, args) -> None:
    if not nonempty(args.formatted_fastani):
        syn.record("FASTANI_REFSEQ", "FAILED",
                   f"No {args.sample_name}_REFSEQ_*.fastANI.txt file")
        return

    lines  = read_all_lines(args.formatted_fastani)
    info   = lines[1] if len(lines) > 1 else ""
    fields = info.split("\t")
    pct_s  = fields[0] if fields else "0"
    cov_s  = fields[1] if len(fields) > 1 else "0"
    org    = fields[2] if len(fields) > 2 else ""
    ref    = fields[3] if len(fields) > 3 else ""
    pct    = safe_int(pct_s.split(".")[0])
    cov    = safe_int(cov_s.split(".")[0])

    if pct_s.startswith("0."):
        syn.record("FASTANI_REFSEQ", "FAILED", "No assembly file to work with")
        return

    if pct >= 95 and cov >= ANI_COVERAGE_THRESHOLD:
        syn.record("FASTANI_REFSEQ", "SUCCESS",
                   f"{pct}%ID {cov}%cov  tax={org}  ref={ref}", escalate=False)
    elif pct < 95 and cov < ANI_COVERAGE_THRESHOLD:
        syn.record("FASTANI_REFSEQ", "WARNING",
                   f"% Identity({pct}%) and % coverage({cov}%) too low. {info}")
    elif pct < 95:
        syn.record("FASTANI_REFSEQ", "WARNING",
                   f"% Identity({pct}%) is too low: {info}")
    else:
        syn.record("FASTANI_REFSEQ", "WARNING",
                   f"% coverage too low ({cov}%). {info}")


def check_mlst(syn: Synopsis, args) -> None:
    if not nonempty(args.mlst_file):
        with open(args.mlst_file) as fh:
            main_line = fh.readline()
        result = ("FAILED",
                  f"No MLST entries in {args.sample_name}.tsv"
                  if main_line.startswith("Sample  Source")
                  else f"{args.sample_name}.tsv does not exist")
        syn.record("MLST", *result)
        return

    lines = read_all_lines(args.mlst_file)
    if len(lines) < 2:
        syn.record("MLST", "FAILED",
                   f"No MLST entries in {args.sample_name}.tsv")
        return

    for line in lines[1:]:
        f          = line.split("\t")
        mlst_db    = f[3] if len(f) > 3 else "-"
        mlst_type  = f[4] if len(f) > 4 else "-"
        mlst_src   = f[1] if len(f) > 1 else ""
        if mlst_db == "-":
            syn.record("MLST", "FAILED", "No scheme identified")
        elif mlst_type == "-":
            syn.record(f"MLST-{mlst_db.upper()}", "FAILED",
                       f"No type identified, scheme={mlst_db} via {mlst_src}")
        else:
            syn.record(f"MLST-{mlst_db.upper()}", "SUCCESS",
                       f"ST{mlst_type} via {mlst_src}", escalate=False)


def check_srst2(syn: Synopsis, args) -> None:
    if not nonempty(args.srst2_file):
        syn.record("SRST2", "FAILED",
                   f"{args.sample_name}__fullgenes__*.txt file does not exist")
        return
    lines    = read_all_lines(args.srst2_file)
    cols     = lines[0].split("\t") if lines else []
    narc_num = len(cols) - 1
    db       = "_".join(Path(args.srst2_file).stem.split("_")[3:5])
    if narc_num == 0:
        syn.record("SRST2", "ALERT",
                   f"Completed, but NO KNOWN AMR genes present from {db}")
    else:
        syn.record("SRST2", "SUCCESS",
                   f"{narc_num} gene(s) found from {db}", escalate=False)


def check_amrfinder(syn: Synopsis, args) -> None:
    if not nonempty(args.amr_file):
        syn.record("AMRFINDER", "FAILED",
                   f"{args.sample_name}_all_genes.tsv does not exist")
        return
    count = sum(1 for l in read_all_lines(args.amr_file) if "POINT" in l)
    syn.record("AMRFINDER", "SUCCESS",
               f"No point mutations were found" if count == 0
               else f"{count} point mutation(s) found", escalate=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    reads_min = max(args.coverage, READS_MIN_DEFAULT) if args.coverage >= 30 else READS_MIN_DEFAULT

    syn  = Synopsis(Path(f"{args.sample_name}.synopsis"))
    today = datetime.now().strftime("%c")
    syn._fh.write(
        f"---------- Checking {args.sample_name} for successful completion"
        f" on {today} ----------\n"
    )
    syn.record("Summarized", "SUCCESS", today, escalate=False)

    run_type     = "assembly-only" if args.assembly_only else "all"
    qc_fail      = ""
    bps_post_all = 0
    assembly_length = 0

    if run_type == "all":
        _, _, raw_exists = check_raw_reads(syn, args)
        check_raw_counts(syn, args, raw_exists)
        _, _, _, total_exists = check_trimmed_reads(syn, args)
        _, bps_post_all = check_trimmed_counts(syn, args, total_exists)
        check_kraken2_reads(syn, args)
    else:
        for lbl in ("QC_COUNTS", "Q30_STATS", "BBDUK", "TRIMMING",
                    "KRAKEN2_READS", "KRONA_READS"):
            syn.write(lbl, "NA", "Assembly only isolate")

    check_assembly(syn, args)

    if args.cdc_phoenix_mode:
        check_kraken2_assembly(syn, args)

    check_kraken2_weighted(syn, args)
    assembly_length, _ = check_quast(syn, args)

    if assembly_length < 1_000_000 and assembly_length > 0:
        qc_fail += f"smaller_than_1000000_bps({assembly_length})-"

    dec_genus, dec_species = check_taxa(syn, args)

    if args.assembly_ratio_file:
        qc_fail += check_assembly_ratio(syn, args, dec_genus, dec_species)

    qc_fail += check_coverage(syn, args, bps_post_all, assembly_length, reads_min)

    if args.cdc_phoenix_mode:
        check_busco(syn, args)

    check_fastani(syn, args)

    if args.mlst_file:
        check_mlst(syn, args)

    if args.gamma_ar:
        db_tag = "_".join(Path(args.gamma_ar).stem.split("_")[-3:-1])
        check_gamma(syn, "GAMMA_AR", args.gamma_ar, "AR", args.sample_name, db_tag)

    if args.amr_file:
        check_amrfinder(syn, args)

    if args.cdc_phoenix_mode:
        if run_type == "all":
            check_srst2(syn, args)
        else:
            syn.write("SRST2", "NA", "Assembly only isolate")

    if args.gamma_replicon:
        db_tag = "_".join(Path(args.gamma_replicon).stem.split("_")[-3:-1])
        check_gamma(syn, "PLASMID_REPLICONS", args.gamma_replicon,
                    "replicon", args.sample_name, db_tag)

    if args.gamma_hv:
        db_tag = "_".join(Path(args.gamma_hv).stem.split("_")[-3:-1])
        check_gamma(syn, "HYPERVIRULENCE", args.gamma_hv,
                    "hypervirulence", args.sample_name, db_tag)

    if qc_fail:
        syn.record("Auto Pass/FAIL", "FAIL", qc_fail.rstrip("-"))
    else:
        syn.write("Auto Pass/FAIL", "PASS",
                  "Minimum Requirements met for coverage(30x)/ratio_stdev(<2.58)"
                  "/min_length(>1000000) to pass auto QC filtering")

    syn._fh.write(
        f"---------- {args.sample_name} completed as {syn.status} ----------\n"
        "WARNINGS: out of line with what is expected and MAY cause problems downstream.\n"
        "ALERT: something to note, does not mean it is a poor-quality assembly.\n"
    )

    if args.cdc_phoenix_mode:
        syn._fh.write(
            "\n*BUSCO defines core genes as single-copy orthologs that should be "
            "highly conserved among the closely related species.\n"
        )

    syn.close()
    print(f"---------- {args.sample_name} completed as {syn.status} ----------")


if __name__ == "__main__":
    main()