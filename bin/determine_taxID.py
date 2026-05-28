#!/usr/bin/env python3

"""
Description: Creates a single file that attempts to pull the best taxonomic
             information from the isolate. Operates in priority order:
             1. ANI  2. kraken2 weighted assembly  3. kraken2 trimmed reads

Usage: ./determine_taxID.py -k weighted_kraken_report -s sample_name
           -f formatted_fastani_file -d nodes_X.dmp -m names_X.dmp
           [-r reads_kraken_report] [-V]

Output: sample_name.tax

Modules required: None (taxonomy dumps may be gzipped)

Created by Nick Vlachos (nvx4@cdc.gov)
Converted to Python: (05/06/2025)
"""

import argparse
import gzip
import re
import sys
from pathlib import Path

__version__ = "2.0.1"

TAXA_LEVELS = ("kingdom", "phylum", "class", "order", "family", "genus", "species")
NOT_ASSIGNED = "Not_assigned"

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Determine best taxonomic ID from ANI / kraken2 sources."
    )
    p.add_argument("-k", "--weighted-kraken",  default="", help="Weighted kraken2 report")
    p.add_argument("-s", "--sample-name",      required=True)
    p.add_argument("-f", "--fastani-file",     default="", help="Formatted FastANI output file")
    p.add_argument("-d", "--nodes",            required=True, help="NCBI nodes.dmp (plain or gzipped)")
    p.add_argument("-m", "--names",            required=True, help="NCBI names.dmp (plain or gzipped)")
    p.add_argument("-r", "--trimmed-kraken",   default="", help="Trimmed-reads kraken2 report (optional)")
    p.add_argument("-V", "--version",          action="version", version=f"%(prog)s: {__version__}")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Taxonomy dump loaders
# ---------------------------------------------------------------------------

def open_possibly_gzipped(path: str):
    p = Path(path)
    return gzip.open(p, "rt") if p.suffix == ".gz" else open(p, "rt")


def load_names(names_path: str) -> tuple[dict, dict]:
    """
    Parse names.dmp into two lookup structures:
        by_taxid  : taxid -> scientific_name  (str)
        by_name   : lowercase_name -> taxid   (int)

    names.dmp columns (tab-pipe-tab delimited):
        0: taxid | 2: name_txt | 4: unique_name | 6: name_class
    """
    by_taxid: dict[int, str] = {}
    by_name:  dict[str, int] = {}

    with open_possibly_gzipped(names_path) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 4:
                continue
            taxid      = int(parts[0].strip())
            name_txt   = parts[1].strip()
            name_class = parts[3].strip().rstrip("|").strip()
            if name_class == "scientific name":
                by_taxid[taxid]           = name_txt
                by_name[name_txt.lower()] = taxid

    return by_taxid, by_name


def load_nodes(nodes_path: str) -> dict[int, tuple[int, str]]:
    """
    Parse nodes.dmp into:
        taxid -> (parent_taxid, rank)
    """
    nodes: dict[int, tuple[int, str]] = {}
    with open_possibly_gzipped(nodes_path) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 3:
                continue
            taxid  = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank   = parts[2].strip()
            nodes[taxid] = (parent, rank)
    return nodes


# ---------------------------------------------------------------------------
# Species string cleanup — mirrors bash logic exactly, same priority order:
#   1. complex sp.  2. strain  3. plain complex  4. everything else
# ---------------------------------------------------------------------------

def clean_species(raw: str) -> tuple[str, str]:
    """
    Normalize the species string extracted from ANI / kraken2.
    Returns (cleaned_species, strain_suffix).

    Priority order matches bash:
      1. "complex sp." substring  -> extract strain, truncate species at complex sp.
      2. "strain" substring       -> extract strain, truncate species at strain
      3. plain "complex"          -> blank both (bash sets species="" strain="")
      4. everything else          -> fix sp.<nospace> and trailing -chromosome only
    """
    species = raw.strip()
    strain  = ""

    if "complex sp." in species:
        # strain = "complex sp." and everything after
        m_strain = re.search(r"complex sp\.", species)
        strain   = species[m_strain.start():] if m_strain else ""
        # species = everything before "complex sp." with trailing "- " stripped,
        # then further truncated at the first "-"
        m_pre = re.search(r"^(.*?)[- ]*complex sp\.", species)
        if m_pre:
            species = re.sub(r"-.*$", "", m_pre.group(1)).strip()
        else:
            species = ""

    elif "strain" in species:
        # strain = "strain" and everything after
        m_strain = re.search(r"strain", species)
        strain   = species[m_strain.start():] if m_strain else ""
        # species = everything before "strain" with trailing "- " stripped,
        # then further truncated at the first "-"
        m_pre = re.search(r"^(.*?)[- ]*strain", species)
        if m_pre:
            species = re.sub(r"-.*$", "", m_pre.group(1)).strip()
        else:
            species = ""

    elif "complex" in species:
        # bash: plain complex (no "sp.") blanks both
        species = ""
        strain  = ""

    # Fix "sp.<non-space>" -> "sp. <rest>"  (applied in all branches)
    if species:
        species = re.sub(r"sp\.(\S)", r"sp. \1", species)
        # Drop trailing "-chromosome"
        species = re.sub(r"-chromosome$", "", species)

    return species, strain


# ---------------------------------------------------------------------------
# Source readers
# ---------------------------------------------------------------------------

def read_ani(fastani_file: str) -> tuple[str, str, str, str]:
    """
    Returns (genus, species, confidence_index, source_file).
    Raises ValueError if the file is unusable.
    """
    path = Path(fastani_file)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"FastANI file missing or empty: {fastani_file}")

    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]

    header = lines[0] if lines else ""
    bad_prefixes = (
        "No matching ANI database found for",
        "Mash/FastANI Error:",
        "0.00%",
    )
    if any(header.startswith(pfx) for pfx in bad_prefixes):
        raise ValueError(f"FastANI header indicates no usable result: {header}")

    info   = lines[-1] if len(lines) > 1 else ""
    fields = info.split("\t")
    organism    = fields[2] if len(fields) > 2 else ""
    genus       = organism.split()[0] if organism else ""
    # Strip surrounding brackets, matching bash: sed 's/[][]//g'
    species_raw = " ".join(organism.split()[1:]).replace("[", "").replace("]", "") if organism else ""
    confidence  = fields[0] if fields else "0"

    return genus, species_raw, confidence, fastani_file


def read_kraken(kraken_path: str, source_label: str) -> tuple[str, str, str, str]:
    """
    Returns (genus, species, confidence_index, source_file).
    Raises ValueError if file is missing/empty.
    """
    path = Path(kraken_path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Kraken file missing or empty: {kraken_path}")

    genus   = ""
    species = ""

    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]

    for line in lines:
        first = line[0] if line else ""
        parts = line.split()
        if first == "s" and len(parts) >= 3:
            species = parts[2]
        elif first == "G" and len(parts) >= 3:
            genus   = parts[2]

    confidence = lines[-1].split()[1] if lines and len(lines[-1].split()) > 1 else "0"
    return genus, species, confidence, kraken_path


# ---------------------------------------------------------------------------
# Source priority waterfall
# ---------------------------------------------------------------------------

def determine_source(args) -> tuple[str, str, str, str, str]:
    """
    Try sources in priority order: ANI -> weighted kraken -> trimmed kraken.
    Returns (source_label, genus, species, confidence, source_file).
    """
    # 1. ANI
    if args.fastani_file and "fastANI" in args.fastani_file and args.fastani_file.endswith(".txt"):
        try:
            genus, species, conf, src = read_ani(args.fastani_file)
            return "ANI_REFSEQ", genus, species, conf, src
        except ValueError:
            pass

    # 2. Weighted kraken (assembly)
    if args.weighted_kraken:
        try:
            genus, species, conf, src = read_kraken(args.weighted_kraken, "kraken2_wtasmbld")
            return "kraken2_wtasmbld", genus, species, conf, src
        except ValueError:
            pass

    # 3. Trimmed reads kraken
    if args.trimmed_kraken:
        try:
            genus, species, conf, src = read_kraken(args.trimmed_kraken, "kraken2_trimmed")
            return "kraken2_trimmed", genus, species, conf, src
        except ValueError:
            pass

    return NOT_ASSIGNED, "", "", "0", NOT_ASSIGNED


# ---------------------------------------------------------------------------
# TaxID lookup
# ---------------------------------------------------------------------------

def find_species_taxid(genus: str, species: str, by_name: dict) -> int:
    """
    Look up species taxID from names index.
    First tries "Genus species", then replaces hyphens with spaces and retries.
    Returns 0 if not found.
    """
    key   = f"{genus.capitalize()} {species}".lower()
    taxid = by_name.get(key, 0)
    if taxid:
        return taxid
    # Retry with hyphens -> spaces, matching bash: species="${species//-/ }"
    key2 = key.replace("-", " ")
    return by_name.get(key2, 0)

def print_tax_map(taxid_map: dict) -> None:
    print("TaxID map:")
    for lvl in TAXA_LEVELS:
        tid = taxid_map.get(lvl, 0)
        print(f"  {lvl}: {tid}")


def find_genus_taxid(genus: str, by_name: dict) -> int:
    return by_name.get(genus.lower(), 0)


def walk_taxonomy(start_taxid: int, nodes: dict, max_levels: int) -> dict:
    """
    Walk up the taxonomy tree from start_taxid, filling in levels.
    Matches bash while loop behavior: advances taxid every iteration
    regardless of whether the rank matched, counter only increments on match.
    """
    taxid_map = {lvl: 0 for lvl in TAXA_LEVELS}
    taxid     = start_taxid
    counter   = 0
    visited   = set()

    while counter < max_levels:
        if taxid not in nodes or taxid in visited:
            break
        visited.add(taxid)
        parent, rank = nodes[taxid]
        index        = max_levels - 1 - counter
        level_name   = TAXA_LEVELS[index] if 0 <= index < len(TAXA_LEVELS) else ""

        if rank == level_name or rank == "superkingdom":
            taxid_map[level_name] = taxid
            counter += 1

        taxid = parent

    return taxid_map


def taxa_indices_at(index: int) -> str:
    if 0 <= index < len(TAXA_LEVELS):
        return TAXA_LEVELS[index]
    return ""


def resolve_names(taxid_map: dict, by_taxid: dict) -> dict:
    """Convert taxID map to scientific name map."""
    name_map: dict[str, str] = {}
    for lvl in TAXA_LEVELS:
        tid  = taxid_map.get(lvl, 0)
        name = by_taxid.get(tid, "NA") if tid else "NA"
        name_map[lvl] = name.capitalize() if name != "NA" else "NA"
    return name_map


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------

def write_unknown_tax(sample_name: str, source: str,
                      confidence: str, source_file: str) -> None:
    out = Path(f"{sample_name}.tax")
    with out.open("w") as fh:
        fh.write(f"{source}\t{confidence}\t{source_file}\n")
        for prefix in ("K", "P", "C", "O", "F", "G", "s"):
            fh.write(f"{prefix}:\tUnknown\n")


def write_tax(sample_name: str, source: str, confidence: str, source_file: str,
              taxid_map: dict, name_map: dict, species: str) -> None:
    """Write the final .tax file."""
    lvl_prefixes = {
        "kingdom": "K", "phylum": "P", "class": "C", "order": "O",
        "family": "F", "genus": "G", "species": "s",
    }

    # Species display: lowercase full name then cut field 2+, matching bash:
    #   tax_name_list[species]=$(echo ${tax_name_list[species],,} | cut -d' ' -f2-)
    # Exception: if "sp." is in the name, skip the lowercase and just cut field 2+
#    if sp_name not in ("NA", "Unknown"):
#        if "sp." in sp_name:
#            # bash: cut -d' ' -f2- without lowercasing
#            sp_name = " ".join(sp_name.split()[1:])
#        else:
#            # bash: lowercase entire string first, then cut field 2+
    sp_name = name_map.get("species", "NA")
    genus_name = name_map.get("genus", "NA")
    parts = sp_name.split()

    if len(parts) > 1:
        # Strip genus prefix first if something meaningful remains after
        if parts[0].lower() == genus_name.lower():
            remaining = parts[1:]
            if not (len(remaining) == 1 and remaining[0].lower() == genus_name.lower()):
                parts = remaining

        # Now check what we're left with and format accordingly
        if parts[0].lower() == "sp.":
            sp_name = "sp. " + "_".join(parts[1:])
        else:
            sp_name = "_".join(parts)
    else:
        sp_name = sp_name.lower()

    out = Path(f"{sample_name}.tax")
    with out.open("w") as fh:
        fh.write(f"{source}\t{confidence}\t{source_file}\n")
        for lvl in TAXA_LEVELS:
            prefix = lvl_prefixes[lvl]
            tid    = taxid_map.get(lvl, 0) or "NA"
            name   = sp_name if lvl == "species" else name_map.get(lvl, "NA")
            fh.write(f"{prefix}:{tid}\t{name}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    source, genus, species_raw, confidence, source_file = determine_source(args)

    if not species_raw:
        genus    = re.sub(r"[\s\[\]]", "", genus)
        genus_ok = bool(genus)
    else:
        genus_ok = True

    # Clean species — priority order: complex sp. > strain > complex > other
    species = ""
    strain  = ""
    if species_raw:
        species, strain = clean_species(species_raw)

    genus = re.sub(r"[\s\[\]]", "", genus)

    # Load taxonomy dumps once
    by_taxid, by_name = load_names(args.names)
    nodes = load_nodes(args.nodes)

    # Look up taxIDs
    species_taxid = 0
    genus_taxid   = 0

    if species:
        species_taxid = find_species_taxid(genus, species, by_name)

    if not species_taxid:
        genus_taxid = find_genus_taxid(genus, by_name)

    # Bail out if nothing found — matches bash unknown output exactly
    if not species_taxid and not genus_taxid:
        write_unknown_tax(args.sample_name, source, confidence, source_file)
        print("No ACCEPTABLE source found to determine taxonomy")
        sys.exit(0)


    # Walk taxonomy tree
    taxid_map = {lvl: 0 for lvl in TAXA_LEVELS}

    if species_taxid:
        max_counter = 7
        start_taxid = species_taxid
    else:
        max_counter = 6
        start_taxid = genus_taxid
        # Species name falls back to cleaned species string or Unknown, matching bash
        taxid_map["species"] = 0

    filled = walk_taxonomy(start_taxid, nodes, max_counter)
    taxid_map.update(filled)

    name_map = resolve_names(taxid_map, by_taxid)

    # If no species taxid, use cleaned species string as name, matching bash behavior
    if not species_taxid:
        name_map["species"] = species if species else "Unknown"

    write_tax(args.sample_name, source, confidence, source_file,
              taxid_map, name_map, species)

if __name__ == "__main__":
    main()