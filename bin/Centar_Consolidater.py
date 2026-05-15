#!/usr/bin/env python3

"""
Description: Extract GOIs for C. diff samples from a PHoeNIx source isolate.

Usage: ./Centar_consolidater.py -t tox_file -c clade_file -o output_file
           -y toxinotype_file -r ribotype_file -p plasmid_file
           -s sample_name -a gamma_ar_file -n gamma_nt_file -x crosswalk_file

v1.0.0 (07/15/2024)
Converted to Python: (05/06/2025)

Created by Nick Vlachos (nvx4@cdc.gov)
"""

import argparse
import sys
from pathlib import Path

__version__ = "2.0.0"

# ---------------------------------------------------------------------------
# Known mutation tables (direct port of bash associative arrays)
# ---------------------------------------------------------------------------

AA_POINTS: dict[str, list[str]] = {
    "dacS":  ["E238D", "V183A"],
    "fur":   ["E41K"],
    "gdpP":  ["E328Stop", "truncation at codon 328 (of 665 codons)"],
    "glyC":  ["A229T"],
    "murG":  ["P109L"],
    "nifJ":  ["Q803R", "G423E"],
    "rpoB":  ["V1143D", "V1143F", "V1143G", "V1143L", "Q1074R", "Q1074H", "Q1074K"],
    "rpoC":  ["D245Y", "D1127E", "D237Y", "Q781R", "R89G"],
    "thiH":  ["S328F"],
    "vanR":  ["T115A"],
    "vanS":  ["G319D", "R314L", "R314H", "S313F", "T349I"],
    "gyrA":  ["A117S", "A118S", "A118T", "A384D", "A92E", "D103N", "D71G", "D71V",
              "D81N", "E123K", "L345I", "P116A", "R90K", "T82A", "T82I", "T82V", "V43D"],
    "gyrB":  ["D426N", "D426V", "E466K", "E466V", "I139R", "L444F", "Q434K", "R377G",
              "R447K", "R447L", "S366V", "S464T", "V130I"],
    "feoB":  ["117DelA", "1 bp Deletion at 120"],
    "hemN":  ["Y214Stop", "1 bp Deletion at 642"],
    "hsmA":  ["372DelA", "1 bp Deletion at 372"],
    "lscR":  ["V76A", "153DelA", "1 bp Deletion at 153"],
    "marR":  ["349DelT", "1 bp Deletion at 356"],
    "sdaB":  ["883DelGCA", "879DelACG", "3 bp Deletion at 887"],
}

NT_POINTS: dict[str, list[str]] = {
    "PNimB": ["T115G"],
}

# ---------------------------------------------------------------------------
# Gene configs for GAMMA tox processing.
# Each entry describes how to parse and format a gene's GAMMA hits.
#
# Fields:
#   pattern        : string to grep for in the tox file
#   has_allele     : whether to extract allele field (col 2 of '__'-split name)
#   has_muts       : whether to extract mutation description (col 6)
#   trunc_100_is_1 : if True, Trunc match with IDA==100 counts as present (=1)
#   skip_trunc     : if True, skip the codon-trunc check entirely (cdtAB1/2)
#   extra_cols     : number of extra tab fields before stats  (affects _set format)
# ---------------------------------------------------------------------------

GENE_CONFIGS: dict[str, dict] = {
    # PaLoc toxins
    "tcdA":         {"pattern": "tcdA",        "has_allele": False, "has_muts": False, "trunc_100_is_1": False, "skip_trunc": False},
    "tcdB":         {"pattern": "tcdB",        "has_allele": False, "has_muts": False, "trunc_100_is_1": False, "skip_trunc": False},
    "tcdC":         {"pattern": "tcdC",        "has_allele": True,  "has_muts": True,  "trunc_100_is_1": False, "skip_trunc": False},
    "tcdD":         {"pattern": "tcdD",        "has_allele": False, "has_muts": False, "trunc_100_is_1": False, "skip_trunc": False},
    "tcdE":         {"pattern": "tcdE",        "has_allele": False, "has_muts": False, "trunc_100_is_1": False, "skip_trunc": False},
    "PaLoc_NonTox": {"pattern": "PaLoc_NonTox","has_allele": True,  "has_muts": True,  "trunc_100_is_1": True,  "skip_trunc": False},
    # CDT locus
    "cdtA":         {"pattern": "cdtA_",       "has_allele": False, "has_muts": False, "trunc_100_is_1": False, "skip_trunc": False},
    "cdtB":         {"pattern": "cdtB_",       "has_allele": False, "has_muts": False, "trunc_100_is_1": False, "skip_trunc": False},
    "cdtR":         {"pattern": "cdtR",        "has_allele": True,  "has_muts": True,  "trunc_100_is_1": True,  "skip_trunc": False},
    "cdtAB1":       {"pattern": "cdtAB1",      "has_allele": False, "has_muts": False, "trunc_100_is_1": False, "skip_trunc": True},
    "cdtAB2":       {"pattern": "cdtAB2",      "has_allele": False, "has_muts": False, "trunc_100_is_1": False, "skip_trunc": True},
    "nontox":       {"pattern": "cdtNonTox",   "has_allele": False, "has_muts": False, "trunc_100_is_1": True,  "skip_trunc": False},
}

OUTPUT_HEADER = "\t".join([
    "WGS_ID", "MLST Clade", "Diffbase_Toxinotype",
    "tcdA_presence", "Diffbase_Toxin-A_sub-type", "tcdA [%Nuc_Identity | %AA_Identity | %Coverage]",
    "tcdB_presence", "Diffbase_Toxin-B_sub-type", "tcdB [%Nuc_Identity | %AA_Identity | %Coverage]",
    "tcdC_presence", "tcdC_Variant", "tcdC other mutations", "tcdC [%Nuc_Identity | %AA_Identity | %Coverage]",
    "tcdR_presence", "tcdR [%Nuc_Identity | %AA_Identity | %Coverage]",
    "tcdE_presence", "tcdE [%Nuc_Identity | %AA_Identity | %Coverage]",
    "PaLoc_NonTox_presence", "PaLoc_NonTox_Variant", "PaLoc_NonTox other mutations", "PaLoc_NonTox [%Nuc_Identity | %AA_Identity | %Coverage]",
    "cdtA_presence", "cdtA [%Nuc_Identity | %AA_Identity | %Coverage]",
    "cdtB_presence", "cdtB [%Nuc_Identity | %AA_Identity | %Coverage]",
    "cdtR_presence", "cdtR_Variant", "cdtR other mutations", "cdtR [%Nuc_Identity | %AA_Identity | %Coverage]",
    "cdtAB1_presence", "cdtAB1 [%Nuc_Identity | %AA_Identity | %Coverage]",
    "cdtAB2_presence", "cdtAB2 [%Nuc_Identity | %AA_Identity | %Coverage]",
    "cdt_NonTox_presence", "cdt_NonTox[%Nuc_Identity | %AA_Identity | %Coverage]",
    "gyrA known mutations", "gyrA other mutations", "gyrA [%Nuc_Identity | %AA_Identity | %Coverage]",
    "gyrB known mutations", "gyrB other mutations", "gyrB [%Nuc_Identity | %AA_Identity | %Coverage]",
    "dacS known mutations", "dacS other mutations", "dacS [%Nuc_Identity | %AA_Identity | %Coverage]",
    "feoB known mutations", "feoB other mutations", "feoB [%Nuc_Identity | %AA_Identity | %Coverage]",
    "fur known mutations",  "fur other mutations",  "fur [%Nuc_Identity | %AA_Identity | %Coverage]",
    "gdpP known mutations", "gdpP other mutations", "gdpP [%Nuc_Identity | %AA_Identity | %Coverage]",
    "glyC known mutations", "glyC other mutations", "glyC [%Nuc_Identity | %AA_Identity | %Coverage]",
    "hemN known mutations", "hemN other mutations", "hemN [%Nuc_Identity | %AA_Identity | %Coverage]",
    "hsmA known mutations", "hsmA other mutations", "hsmA [%Nuc_Identity | %AA_Identity | %Coverage]",
    "lscR known mutations", "lscR other mutations", "lscR [%Nuc_Identity | %AA_Identity | %Coverage]",
    "marR known mutations", "marR other mutations", "marR [%Nuc_Identity | %AA_Identity | %Coverage]",
    "murG known mutations", "murG other mutations", "murG [%Nuc_Identity | %AA_Identity | %Coverage]",
    "nifJ known mutations", "nifJ other mutations", "nifJ [%Nuc_Identity | %AA_Identity | %Coverage]",
    "PNimB known mutations","PNimB other mutations","PNimB [%Nuc_Identity | %Coverage]",
    "rpoB known mutations", "rpoB other mutations", "rpoB [%Nuc_Identity | %AA_Identity | %Coverage]",
    "rpoC known mutations", "rpoC other mutations", "rpoC [%Nuc_Identity | %AA_Identity | %Coverage]",
    "sdaB known mutations", "sdaB other mutations", "sdaB [%Nuc_Identity | %AA_Identity | %Coverage]",
    "thiH known mutations", "thiH other mutations", "thiH [%Nuc_Identity | %AA_Identity | %Coverage]",
    "vanR known mutations", "vanR other mutations", "vanR [%Nuc_Identity | %AA_Identity | %Coverage]",
    "vanS known mutations", "vanS other mutations", "vanS [%Nuc_Identity | %AA_Identity | %Coverage]",
    "CEMB RT Crosswalk", "Inferred RT", "Probability", "ML Method", "ML Note",
])


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="CENTAR consolidator for C. diff PHoeNIx isolates.")
    p.add_argument("-t", "--tox-input",    default="", help="GAMMA tox file")
    p.add_argument("-c", "--clade-input",  default="", help="Clade/MLST file")
    p.add_argument("-o", "--output",       required=True)
    p.add_argument("-y", "--ttype-file",   default="", help="Toxinotype definition file")
    p.add_argument("-a", "--aa-mut-file",  default="", help="GAMMA AA mutation file")
    p.add_argument("-n", "--nt-mut-file",  default="", help="GAMMA NT mutation file")
    p.add_argument("-r", "--rt-file",      default="", help="Ribotype file")
    p.add_argument("-s", "--sample-name",  required=True)
    p.add_argument("-p", "--plasmid-file", default="", help="Plasmid file (future use)")
    p.add_argument("-x", "--xwalk-rt-file",default="", help="MLST→RT crosswalk file")
    p.add_argument("-V", "--version",      action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


# ---------------------------------------------------------------------------
# GAMMA file helpers
# ---------------------------------------------------------------------------

def pct(raw: str) -> int:
    """Convert a 0–1 float string to integer percent, truncating decimals."""
    try:
        return int(float(raw) * 100)
    except (ValueError, TypeError):
        return 0


def extract_allele(line: str) -> str:
    """Extract allele name: split gene name on '__', take 3rd field (index 2)."""
    gene_field = line.split("\t")[0] if "\t" in line else line.split()[0]
    parts = gene_field.split("__")
    return parts[2] if len(parts) > 2 else ""


def clean_muts(line: str) -> str:
    """Extract col 6 (0-based index 5) mutation description, strip trailing comma."""
    fields = line.split("\t")
    desc = fields[5].strip() if len(fields) > 5 else ""
    return desc.rstrip(",")


def parse_gamma_fields(line: str) -> tuple[int, int, int, str]:
    """Return (IDA%, IDN%, length%, matchtype) from a GAMMA output line."""
    f = line.split("\t")
    ida       = pct(f[9])  if len(f) > 9  else 0
    idn       = pct(f[10]) if len(f) > 10 else 0
    length    = pct(f[11]) if len(f) > 11 else 0
    matchtype = f[4].strip() if len(f) > 4 else ""
    return ida, idn, length, matchtype


def classify_hit(ida: int, length: int, matchtype: str,
                 trunc_100_is_1: bool, skip_trunc: bool) -> str:
    """
    Return the presence string for a single GAMMA hit:
      "1", "TRUNC:CODON[NAA]", or "TRUNC:LENGTH<90[NCOV]"
    """
    if not skip_trunc and "Trunc" in matchtype:
        if trunc_100_is_1 and ida == 100:
            return "1"
        return f"TRUNC:CODON[{ida}AA]"
    if length < 90:
        return f"TRUNC:LENGTH<90[{length}COV]"
    return "1"


def stats_str(idn: int, ida: int, length: int) -> str:
    return f"[{idn}NT|{ida}AA|{length}COV]"


def process_gamma_gene(gene: str, cfg: dict,
                       tox_lines: list[str],
                       subtype: str = "") -> str:
    """
    Process all GAMMA hits for one gene and return the formatted set string.

    For genes with subtype (tcdA, tcdB): set is  presence \\t subtype \\t stats
    For genes with allele+muts (tcdC, PaLoc_NonTox, cdtR): presence \\t allele \\t muts \\t stats
    For simple genes (tcdD, tcdE, cdtA, cdtB, cdtAB1, cdtAB2, nontox): presence \\t stats
    """
    pattern        = cfg["pattern"]
    has_allele     = cfg["has_allele"]
    has_muts       = cfg["has_muts"]
    trunc_100_is_1 = cfg["trunc_100_is_1"]
    skip_trunc     = cfg["skip_trunc"]

    hits = [l for l in tox_lines if pattern in l.split("\t")[0]]

    # ---- zero hits ----
    if not hits:
        if has_allele and has_muts:
            return "0\tNOT_FOUND\tNOT_FOUND\tNOT_FOUND"
        elif subtype:
            return "0\tNOT_FOUND\tNOT_FOUND"
        else:
            return "0\tNOT_FOUND"

    presence_parts: list[str] = []
    allele_parts:   list[str] = []
    muts_parts:     list[str] = []
    stats_parts:    list[str] = []
    presence = ""

    for line in hits:
        ida, idn, length, matchtype = parse_gamma_fields(line)
        hit_str = classify_hit(ida, length, matchtype, trunc_100_is_1, skip_trunc)

        if presence == "":
            presence = hit_str
        else:
            presence = f"{presence}-{hit_str}"

        if has_allele:
            allele_parts.append(extract_allele(line))
        if has_muts:
            m = clean_muts(line)
            muts_parts.append("-" if m in ("No coding mutations", "") else m)

        stats_parts.append(stats_str(idn, ida, length))

    stats = "-".join(stats_parts)

    if has_allele and has_muts:
        alleles = "-".join(allele_parts)
        muts    = "-".join(muts_parts)
        return f"{presence}\t{alleles}\t{muts}\t{stats}"
    elif subtype:
        return f"{presence}\t{subtype}\t{stats}"
    else:
        return f"{presence}\t{stats}"


# ---------------------------------------------------------------------------
# Ribotype file parser
# ---------------------------------------------------------------------------

def parse_rt_file(rt_path: str) -> str:
    """Parse the ribotype file and return the ML_RT string."""
    try:
        lines = [l.rstrip("\n") for l in Path(rt_path).read_text().splitlines() if l.strip()]
        line  = lines[-1] if lines else ""
        f     = line.split("\t")
        raw_prob = f[4].strip() if len(f) > 4 else "0"
        prob     = f"{float(raw_prob) * 100:.2f}"
        rt       = f[7].strip() if len(f) > 7 else ""
        source   = f[5].strip() if len(f) > 5 else ""
        note     = f[6].strip() if len(f) > 6 else ""
        if not note:
            return f"{rt}\t{prob}\t{source}"
        elif note == "-":
            return f"{rt}\t{prob}\t{source}\tNo Notes"
        else:
            return f"{rt}\t{prob}\t{source}\t{note}"
    except Exception:
        return "NO_RT_FILE"


# ---------------------------------------------------------------------------
# Toxinotype file parser
# ---------------------------------------------------------------------------

def parse_ttype_file(ttype_path: str) -> tuple[str, str, str]:
    """Return (toxinotype, subtype_A, subtype_B)."""
    toxinotype = ""
    subtype_a: list[str] = []
    subtype_b: list[str] = []

    try:
        lines = Path(ttype_path).read_text().splitlines()
        for i, line in enumerate(lines):
            if i == 0:
                continue
            line = line.rstrip("\n")
            if line.startswith("Toxinotype:"):
                toxinotype = line.split("\t")[1].strip() if "\t" in line else ""
            else:
                f = line.split("\t")
                tox_type    = f[1].strip() if len(f) > 1 else ""
                tox_subtype = f[2].strip() if len(f) > 2 else ""
                if tox_type == "Toxin-A":
                    subtype_a.append(tox_subtype)
                elif tox_type == "Toxin-B":
                    subtype_b.append(tox_subtype)
    except Exception:
        return "NO_ttype_file", "NO_ttype_file", "NO_ttype_file"

    return (
        toxinotype,
        ",".join(subtype_a) if subtype_a else "Unknown",
        ",".join(subtype_b) if subtype_b else "Unknown",
    )


# ---------------------------------------------------------------------------
# Clade/MLST file parser
# ---------------------------------------------------------------------------

def parse_clade_file(clade_path: str, xwalk_path: str) -> tuple[str, str, str]:
    """Return (clade, mlst_type, xrt)."""
    try:
        lines = [l.rstrip("\n") for l in Path(clade_path).read_text().splitlines() if l.strip()]
        if len(lines) != 2:
            return "Clade_file_incorrect", "", "Clade_file_incorrect"
        f     = lines[1].split("\t")
        clade = f[1].strip() if len(f) > 1 else ""
        mlst  = f[2].strip() if len(f) > 2 else ""
    except Exception:
        return "No_clade/MLST_file", "", ""

    xrt = lookup_xwalk(mlst, xwalk_path) if xwalk_path and Path(xwalk_path).is_file() else "No_lookup_file"
    if xrt == "unset":
        xrt = f"{mlst} has no crosswalk match"
    return clade, mlst, xrt


def lookup_xwalk(mlst: str, xwalk_path: str) -> str:
    """Look up MLST → RT crosswalk, zero-padding RT numbers."""
    try:
        with open(xwalk_path) as fh:
            for line in fh:
                f = line.split("\t")
                if f[0].strip() == mlst:
                    raw_rts = f[1].strip() if len(f) > 1 else ""
                    parts   = [r.strip() for r in raw_rts.split(",")]
                    padded  = []
                    for r in parts:
                        if len(r) == 1:
                            padded.append(f"00{r}")
                        elif len(r) == 2:
                            padded.append(f"0{r}")
                        else:
                            padded.append(r)
                    result = ",".join(padded).lstrip(",")
                    return result if result else "unset"
    except Exception:
        pass
    return "unset"


# ---------------------------------------------------------------------------
# Mutation file processing
# ---------------------------------------------------------------------------

def unescape_muts(s: str) -> str:
    """Replace __ with space and strip whitespace — mirrors bash ${var//__/ }."""
    return s.replace("__", " ").strip()


def classify_mutations(all_muts_str: str,
                       known_set: list[str],
                       ida: int, length: int,
                       matchtype: str,
                       is_nt: bool = False,
                       max_other: int = 10) -> tuple[str, str]:
    """
    Split a comma-separated mutation string into known and other mutations.
    Returns (known_muts_string, other_muts_string).
    """
    known:  list[str] = []
    other:  list[str] = []

    # Structural truncation flags
    if not is_nt:
        if "Truncation" in matchtype or "Indel Truncation" in matchtype:
            other.append(f"TRUNC:CODON[{ida}AA]")
        elif length <= 90:
            other.append(f"TRUNC:LENGTH<90[{length}COV]")
    else:
        if length <= 90:
            other.append(f"TRUNC:LENGTH<90[{length}COV]")

    if all_muts_str and all_muts_str not in ("No coding mutations", "Exact match"):
        for raw_mut in all_muts_str.split(","):
            mut     = unescape_muts(raw_mut)
            matched = False
            for known_mut in known_set:
                if mut == unescape_muts(known_mut):
                    known.append(mut)
                    matched = True
                    break
            if not matched and mut not in ("No mutations", ""):
                other.append(mut)

    known_str = ",".join(known) if known else "-"
    if len(other) > max_other:
        other_str = f"{len(other)} other mutations"
    else:
        other_str = ",".join(other) if other else "-"

    return known_str, other_str


def process_mut_file(mut_path: str, points: dict[str, list[str]],
                     is_nt: bool = False) -> dict[str, str]:
    """
    Process a GAMMA mutation file against a points dict.
    Returns {gene_name: set_string}.
    The file is read once; each gene's lines are matched in one pass.
    """
    results: dict[str, str] = {g: f"NOT_FOUND\tNOT_FOUND\tNOT_FOUND" for g in points}

    try:
        lines = [l.rstrip("\n") for l in Path(mut_path).read_text().splitlines() if l.strip()]
    except Exception:
        return {g: f"NO_mut_file\tNO_mut_file\tNO_mut_file" for g in points}

    for gene, known_set in points.items():
        for line in lines:
            gene_section = line.split()[0] if line.split() else ""
            if gene not in gene_section:
                continue

            all_muts = clean_muts(line)
            f        = line.split("\t")

            if is_nt:
                # Fixed: original overwrote line_IDNT with line_raw_IDA on second assignment
                idn    = pct(f[13]) if len(f) > 13 else 0
                length = pct(f[14]) if len(f) > 14 else 0
                ida    = 0   # not used for NT
                matchtype = ""
                stats  = f"[{idn}NT|{length}COV]"
            else:
                ida, idn, length, matchtype = parse_gamma_fields(line)
                stats = stats_str(idn, ida, length)

            known_str, other_str = classify_mutations(
                all_muts, known_set, ida, length, matchtype, is_nt=is_nt
            )
            results[gene] = f"{known_str}\t{other_str}\t{stats}"
            break   # first match per gene

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    # ---- Ribotype ----
    ml_rt = parse_rt_file(args.rt_file) if args.rt_file and Path(args.rt_file).is_file() \
            else "NO_RT_FILE"

    # ---- Toxinotype ----
    if args.ttype_file and Path(args.ttype_file).is_file():
        toxinotype, subtype_a, subtype_b = parse_ttype_file(args.ttype_file)
    else:
        print("No toxinotype input file")
        toxinotype = subtype_a = subtype_b = "NO_ttype_file"

    # ---- Tox gene hits — read file once, pass lines to each gene processor ----
    if args.tox_input and Path(args.tox_input).is_file():
        tox_lines = [
            l.rstrip("\n")
            for l in Path(args.tox_input).read_text().splitlines()
            if l.strip()
        ]
        # subtypes only apply to tcdA and tcdB
        gene_subtypes = {"tcdA": subtype_a, "tcdB": subtype_b}
        tox_sets: dict[str, str] = {}
        for gene, cfg in GENE_CONFIGS.items():
            tox_sets[gene] = process_gamma_gene(
                gene, cfg, tox_lines,
                subtype=gene_subtypes.get(gene, "")
            )
    else:
        print("No tox file found")
        tox_sets = {
            "tcdA":         "No_Tox_file\tNo_Tox_file\tNo_Tox_file",
            "tcdB":         "No_Tox_file\tNo_Tox_file\tNo_Tox_file",
            "tcdC":         "No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file",
            "tcdD":         "No_Tox_file\tNo_Tox_file",
            "tcdE":         "No_Tox_file\tNo_Tox_file",
            "PaLoc_NonTox": "No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file",
            "cdtA":         "No_Tox_file\tNo_Tox_file",
            "cdtB":         "No_Tox_file\tNo_Tox_file",
            "cdtR":         "No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file",
            "cdtAB1":       "No_Tox_file\tNo_Tox_file",
            "cdtAB2":       "No_Tox_file\tNo_Tox_file",
            "nontox":       "No_Tox_file\tNo_Tox_file",
        }

    # ---- Clade / MLST ----
    if args.clade_input and Path(args.clade_input).is_file():
        clade, mlst, xrt = parse_clade_file(args.clade_input, args.xwalk_rt_file)
    else:
        print("No clade/mlst file, carry on")
        clade = "No_clade/MLST_file"
        mlst  = ""
        xrt   = ""

    # ---- AA mutations — file read once, all genes processed in one pass ----
    if args.aa_mut_file and Path(args.aa_mut_file).is_file():
        aa_sets = process_mut_file(args.aa_mut_file, AA_POINTS, is_nt=False)
    else:
        print("NO_AA_mut_file")
        aa_sets = {g: "NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file" for g in AA_POINTS}

    # ---- NT mutations ----
    if args.nt_mut_file and Path(args.nt_mut_file).is_file():
        nt_sets = process_mut_file(args.nt_mut_file, NT_POINTS, is_nt=True)
    else:
        print("NO_NT_mut_file")
        nt_sets = {g: "NO_NT_mut_file\tNO_NT_mut_file\tNO_NT_mut_file" for g in NT_POINTS}

    # ---- Write output ----
    out_path = Path(args.output)
    if not out_path.is_file():
        out_path.write_text(OUTPUT_HEADER + "\n")

    row = "\t".join([
        args.sample_name,
        clade,
        toxinotype,
        tox_sets["tcdA"],
        tox_sets["tcdB"],
        tox_sets["tcdC"],
        tox_sets["tcdD"],
        tox_sets["tcdE"],
        tox_sets["PaLoc_NonTox"],
        tox_sets["cdtA"],
        tox_sets["cdtB"],
        tox_sets["cdtR"],
        tox_sets["cdtAB1"],
        tox_sets["cdtAB2"],
        tox_sets["nontox"],
        aa_sets.get("gyrA",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("gyrB",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("dacS",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("feoB",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("fur",   "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("gdpP",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("glyC",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("hemN",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("hsmA",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("lscR",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("marR",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("murG",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("nifJ",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        nt_sets.get("PNimB", "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("rpoB",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("rpoC",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("sdaB",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("thiH",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("vanR",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        aa_sets.get("vanS",  "NOT_FOUND\tNOT_FOUND\tNOT_FOUND"),
        xrt,
        ml_rt,
    ])

    with out_path.open("a") as fh:
        fh.write(row + "\n")


if __name__ == "__main__":
    main()