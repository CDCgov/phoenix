#!/usr/bin/env python3

import copy
import collections
import argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Process MLST results from assembly and read-based typing'
    )
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='MLST output file from mlst tool (assembly-based)'
    )
    parser.add_argument(
        '-s', '--srst2',
        help='SRST2 output file (read-based, optional)'
    )
    parser.add_argument(
        '-t', '--taxonomy',
        required=True,
        help='Taxonomy file (.tax format)'
    )
    parser.add_argument(
        '-d', '--mlst_database',
        required=True,
        help='Path to MLST database directory'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file prefix (default: derived from taxonomy file)'
    )
    parser.add_argument(
        '-v','--version', 
        action='version',
        version=f'%(prog)s {get_version()}'
    )
    
    return parser.parse_args()


def read_taxonomy(taxonomy_file: str) -> Tuple[str, str]:
    """Extract genus and species from taxonomy file."""
    genus = "Unknown"
    species = "Unknown"
    
    try:
        with open(taxonomy_file, 'r') as f:
            for line in f:
                units = line.strip().replace("\t", "|").split("|")
                if len(units) < 2:
                    continue
                    
                tier, info = units[0], units[1]
                if tier == "G:":
                    genus = info
                elif tier == "s:":
                    species = info
                    print(f"Taxonomy: {genus} {species}")
                    break
    except FileNotFoundError:
        print(f"Warning: Taxonomy file not found: {taxonomy_file}")
    
    return genus, species

# Function to get the script version
def get_version():
	return "4.0.0"

def get_pull_date(mlst_db_path: str) -> str:
    """Read database version/pull date."""
    try:
        with open(f"{mlst_db_path}/db_version", 'r') as f:
            return f.readline().strip()
    except FileNotFoundError:
        print(f"Warning: db_version file not found at {mlst_db_path}")
        return "Unknown"


def parse_alleles(items: List[str], start_col: int) -> Tuple[List[str], List[List[str]], bool]:
    """
    Parse allele information from MLST output.
    
    Returns:
        allele_names: List of locus names
        allele_options: List of possible alleles for each locus
        has_multiple: Whether any locus has multiple alleles
    """
    allele_names = []
    allele_options = []
    has_multiple = False
    
    for col in range(start_col, len(items)):
        if "(" not in items[col]:
            continue
            
        locus_name = items[col].split("(")[0]
        alleles = items[col].split("(")[1].split(")")[0].split(",")
        
        if len(alleles) > 1:
            has_multiple = True
        
        allele_names.append(locus_name)
        allele_options.append(alleles)
    
    return allele_names, allele_options, has_multiple


def expand_allele_combinations(template: List, allele_options: List[List[str]]) -> List[List]:
    """
    Generate all possible combinations of alleles (Cartesian product).
    
    Args:
        template: Base profile structure [sample, db, type, count, names, [], source, date]
        allele_options: List of allele choices for each locus
        
    Returns:
        List of all possible allele combinations
    """
    if not allele_options:
        return [template]
    
    expanded = [copy.deepcopy(template)]
    
    for locus_alleles in allele_options:
        new_expanded = []
        for allele_value in locus_alleles:
            # Copy all current combinations
            current_combos = copy.deepcopy(expanded)
            # Add this allele to each combination
            for combo in current_combos:
                combo[5].append(allele_value)
            new_expanded.extend(current_combos)
        expanded = new_expanded
    
    return expanded


def parse_mlst_line(items: List[str], source_type: str, pull_date: str) -> List[List]:
    """
    Parse a single MLST line and expand into all possible profiles.
    
    Args:
        items: Split line from MLST output
        source_type: "assembly" or "reads"
        pull_date: Database pull date
        
    Returns:
        List of expanded profiles
    """
    sample = items[0]
    db_name = items[1]
    
    # Determine starting column based on source type
    if source_type == "assembly":
        st_type = items[2]
        start_col = 3
    else:  # reads/srst2
        if db_name == "No match found":
            return []
        st_type = items[2]
        start_col = 7
    
    # Check if we have alleles
    if len(items) <= start_col:
        return [[sample, db_name, "-", 1, ["-"], ["-"], source_type, pull_date]]
    

    # Parse alleles
    allele_names, allele_options, has_multiple = parse_alleles(items, start_col)

    print("allele names:", allele_names, "options:", allele_options, "multiple:", has_multiple)
    
    if not allele_names:
        return [[sample, db_name, "-", 1, ["-"], ["-"], source_type, pull_date]]
    
    # Create template
    template = [sample, db_name, st_type, len(allele_names), allele_names, [], source_type, pull_date]
    
    # Mark as novel profile if multiple alleles detected
    if has_multiple:
        template[2] = "-"
    
    # Expand combinations
    return expand_allele_combinations(template, allele_options)


def merge_duplicate_profiles(profiles: List[List]) -> List[List]:
    """
    Find and merge duplicate profiles from different sources.
    
    Profiles with identical alleles are merged, with source updated to
    "assembly/reads" if they come from different sources.
    """
    seen_indices = set()
    merged = []
    
    for i in range(len(profiles)):
        if i in seen_indices:
            continue
            
        current = profiles[i]
        current_alleles = current[5]
        current_source = current[6]
        
        # Look for duplicates
        for j in range(i + 1, len(profiles)):
            if j in seen_indices:
                continue
                
            other = profiles[j]
            other_alleles = other[5]
            other_source = other[6]
            
            # Check if alleles match (order-independent)
            if collections.Counter(current_alleles) != collections.Counter(other_alleles):
                continue
            
            # Same database?
            if current[1] != other[1]:
                continue
            
            # Merge sources
            if current_source != other_source:
                sources = {current_source, other_source}
                if sources == {"assembly", "reads"}:
                    current[6] = "assembly/reads"
                elif "assembly/reads" in sources:
                    current[6] = "assembly/reads"
            
            seen_indices.add(j)
        
        merged.append(current)
    
    return merged


def classify_profile_type(profile: List, mlst_db_path: str) -> str:
    """
    Determine the ST type by checking alleles against database.
    
    Returns:
        ST number or status (Novel_allele, Missing_allele, Novel_profile, etc.)
    """
    database = profile[1]
    st_type = profile[2]
    alleles = profile[5]
    source = profile[6]

    print("classifying profile:", profile, database, st_type, alleles)
    
    # Check for special characters indicating issues, this should only be occuring in srst2 profiles, so imma edit somee
#    bad_markers = ["*", "?", "~", "NF"]
#    if any(marker in st_type for marker in bad_markers):
#        return "Novel_profile"

    # Handle dash
    if (database == "-" and int(profile[3]) == 1) or (database.startswith("No match found for") and int(profile[3]) == 1):
#        return "Novel_profile" if profile[3] > 1 else "-"
        return "-" 
    
    if st_type == "failed":
        return "Failed" if profile[3] > 1 else "failed"
    
    # Check if ST is already a valid number
    try:
        int(st_type)
        return st_type  # Already validated
    except ValueError:
        pass
    
    # Check alleles for issues
    has_novel = False
    has_missing = False

    print(alleles)
    
    for allele in alleles:
        if any(marker in allele for marker in ['*', '?', '~']):
            has_novel = True
        elif '-' in allele:
            has_missing = True
    
    if has_missing:
        return "Missing_allele"
    if has_novel:
        return "Novel_allele"
    
    # Check if all alleles are dashes
    if all(a == "-" for a in alleles):
        return "-"
    
    # Look up in database
    db_name = profile[1].split("(")[0]  # Remove any suffix
    profile_file = Path(mlst_db_path) / "pubmlst" / db_name / f"{db_name}.txt"
    
    if not profile_file.exists():
        print(f"Warning: Profile file not found: {profile_file}")
        return "File_not_Found"
    
    try:
        with open(profile_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                    
                parts = line.split("\t")
                if len(parts) < len(alleles) + 1:
                    continue
                
                db_st = parts[0]
                db_alleles = parts[1:len(alleles) + 1]
                
                if db_alleles == alleles:
                    return db_st
    except Exception as e:
        print(f"Error reading profile file: {e}")
    
    return "Novel_profile"


def standardize_database_name(db_name: str) -> str:
    """Standardize database names to consistent format."""
    name_map = {
        "abaumannii": "abaumannii(Oxford)",
        "abaumannii(Oxford)": "abaumannii(Oxford)",
        "Acinetobacter_baumannii#1": "abaumannii(Oxford)",
        "abaumannii_2": "abaumannii(Pasteur)",
        "Acinetobacter_baumannii#2": "abaumannii(Pasteur)",
        "ecoli": "ecoli(Achtman)",
        "ecoli_achtman_4": "ecoli(Achtman)",
        "ecoli_2": "ecoli_2(Pasteur)",
        "Escherichia_coli#2": "ecoli_2(Pasteur)",
        "aparagallinarum_Ghanem": "aparagallinarum(Ghanem)",
        "aparagallinarum_Guo": "aparagallinarum(Guo)",
        # None for 'regular' efaecium
        "efaecium_Bezdicek": "efaecium(Bezdicek)",
        "salmonella_Oxford": "salmonella(Oxford)",
        "salmonella_Achtman": "salmonella(Achtman)",
        "mgallisepticum_Ghanem": "mgallisepticum(Ghanem)",
        "mgallisepticum_Beko": "mgallisepticum(Beko)",
        "pmultocida_multihost": "pmultocida(multihost)",
        "pmultocida_rirdc": "pmultocida(rirdc)",
        # None for 'regular' mbovis
        "mbovis_legacy": "mbovis(legacy)",
        "smutans_Do": "smutans(Do)",
        "smutans_Kakano": "smutans(Kakano)",
        "tpallidum_Grillova": "tpallidum(Grillova)",
        "tpallidum_Pla-Diaz": "tpallidum(Pla-Diaz)"
    }
    
    return name_map.get(db_name, db_name)


def detect_paralogs(profile: List) -> bool:
    """Check for known paralog markers in A. baumannii."""
    db_name = profile[1]
    alleles = profile[5]
    
    # Only check A. baumannii Oxford scheme
    if not any(x in db_name for x in ["abaumannii", "abaumannii(Oxford)", "Acinetobacter_baumannii#1"]):
        return False
    
    # Check third allele for known paralog values
    if len(alleles) > 2:
        if str(alleles[2]) in ["182", "189"]:
            return True
    
    return False


def consolidate_novel_alleles(profiles: List[List]) -> List[List]:
    """
    Consolidate novel allele profiles that can be combined.
    
    Novel allele profiles from assembly and reads that have compatible
    alleles (same or complementary mismatches) are merged.
    """
    # Separate profiles by type
    complete = []
    novel = []
    
    for profile in profiles:
        if profile[2] not in ["Novel_allele", "Missing_allele"]:
            complete.append(profile)
        else:
            novel.append(profile)
    
    # If 0 or 1 novel profiles, nothing to consolidate
    if len(novel) <= 1:
        return complete + novel
    
    # Try to consolidate novel profiles
    consolidated = complete
    processed = set()
    
    for i in range(len(novel)):
        if i in processed:
            continue
        
        primary = novel[i]
        primary_db = primary[1]
        primary_source = primary[6]
        primary_alleles = primary[5]
        
        # Classify primary alleles
        primary_status = []
        for allele in primary_alleles:
            if any(m in allele for m in ['*', '?', '~']):
                clean = allele.replace('*', '').replace('?', '').replace('~', '')
                primary_status.append(('mismatch', clean))
            elif '-' in allele:
                primary_status.append(('missing', None))
            else:
                primary_status.append(('complete', allele))
        
        # Look for compatible profiles
        merged = False
        for j in range(i + 1, len(novel)):
            if j in processed:
                continue
            
            secondary = novel[j]
            
            # Must be same database
            if secondary[1] != primary_db:
                continue
            
            secondary_alleles = secondary[5]
            
            # Classify secondary alleles
            secondary_status = []
            for allele in secondary_alleles:
                if any(m in allele for m in ['*', '?', '~']):
                    clean = allele.replace('*', '').replace('?', '').replace('~', '')
                    secondary_status.append(('mismatch', clean))
                elif '-' in allele:
                    secondary_status.append(('missing', None))
                else:
                    secondary_status.append(('complete', allele))
            
            # Check if compatible
            if not are_alleles_compatible(primary_status, secondary_status):
                continue
            
            # Merge
            new_alleles = []
            for k in range(len(primary_alleles)):
                p_status, p_val = primary_status[k]
                s_status, s_val = secondary_status[k]
                
                if p_status == 'complete' or s_status == 'complete':
                    # Use the complete one
                    new_alleles.append(p_val if p_status == 'complete' else s_val)
                elif p_status == 'mismatch' and s_status == 'mismatch':
                    if p_val == s_val:
                        new_alleles.append(primary_alleles[k])
                    else:
                        # Different mismatches - annotate
                        new_alleles.append(f"{primary_alleles[k]},{secondary_alleles[k]}")
                else:
                    # One is missing - use the other
                    new_alleles.append(primary_alleles[k] if p_status != 'missing' else secondary_alleles[k])
            
            # Determine new source
            sources = {primary_source, secondary[6]}
            if sources == {"assembly", "reads"}:
                new_source = "assembly/reads"
            else:
                new_source = list(sources)[0]
            
            # Create merged profile
            merged_profile = [
                primary[0], primary[1], primary[2], primary[3],
                primary[4], new_alleles, new_source, primary[7]
            ]
            consolidated.append(merged_profile)
            
            processed.add(i)
            processed.add(j)
            merged = True
            break
        
        if not merged:
            consolidated.append(primary)
            processed.add(i)
    
    return consolidated


def are_alleles_compatible(status1: List[Tuple], status2: List[Tuple]) -> bool:
    """Check if two allele status lists can be merged."""
    if len(status1) != len(status2):
        return False
    
    for (s1, v1), (s2, v2) in zip(status1, status2):
        # Complete alleles must match
        if s1 == 'complete' and s2 == 'complete':
            if v1 != v2:
                return False
        # Mismatches must have same base value or be complementary
        elif s1 == 'mismatch' and s2 == 'mismatch':
            # Allow different markers on same base
            continue
        # Complete + mismatch is compatible
        elif (s1 == 'complete' and s2 == 'mismatch') or (s1 == 'mismatch' and s2 == 'complete'):
            continue
        # Missing can be filled by anything
        elif 'missing' in (s1, s2):
            continue
        else:
            return False
    
    return True


def write_output_file(profiles: List[List], isolate_name: str, genus: str, species: str):
    """Write final MLST results to TSV file."""
    outfile = f"{isolate_name}_combined.tsv"
    all_complete = "unknown"
    
    with open(outfile, 'w') as f:
        # Write header
        f.write("WGS_ID\tSource\tPulled_on\tDatabase\tST\t")
        f.write("\t".join([f"locus_{i+1}" for i in range(10)]))
        f.write("\n")
        
        if not profiles:
            f.write(f"{isolate_name}\tNone-{genus} {species}\t-\t-\t-\t")
            f.write("\t".join(["-"] * 10) + "\n")
            all_complete = "False"
        else:
            # Track if all STs are valid
            all_complete = "True"
            
            for profile in profiles:
                sample, db_name, st_type, num_loci, loci_names, alleles, source, date = profile
                
                # Build allele section
                if num_loci == 1 and loci_names[0] == "-":
                    allele_section = "-"
                else:
                    allele_parts = []
                    for locus, allele in zip(loci_names, alleles):
                        # Strip prefix if present
                        locus_name = locus.split("_")[1] if "_" in locus else locus
                        allele_parts.append(f"{locus_name}({allele})")
                    allele_section = "\t".join(allele_parts)
                
                # Check if ST is valid
                invalid_markers = ["Novel_allele", "Missing_allele", "Novel_profile", "-", "failed", "Failed"]
                if st_type in invalid_markers or (not st_type.isnumeric() and "-PARALOG" not in st_type):
                    all_complete = "False"
                
                # Write line
                f.write(f"{sample}\t{source}\t{date}\t{db_name}\t{st_type}\t{allele_section}\n")
    
    # Write status file
    status_file = f"{isolate_name}_status.txt"
    with open(status_file, 'w') as f:
        f.write(f"{all_complete}\n")
    
    print(f"Output written to: {outfile}")
    return outfile, all_complete


def read_mlst_file(filepath: str) -> List[Tuple[str, str]]:
    """Read MLST output file and return list of (line, type) tuples."""
    results = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    results.append((line, "mlst"))
    except FileNotFoundError:
        print(f"Warning: MLST file not found: {filepath}")
    
    return results


def read_srst2_file(filepath: str) -> List[Tuple[str, str]]:
    """Read SRST2 output file and return list of (line, type) tuples."""
    results = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and not line.startswith("Sample\tdatabase\tST"):
                    results.append((line, "srst2"))
    except FileNotFoundError:
        print(f"Warning: SRST2 file not found: {filepath}")
    
    return results


def do_MLST_check(input_MLST_line_tuples: List[Tuple], taxonomy_file: str, mlst_db_path: str, output_prefix: Optional[str] = None):
    """
    Main function to process MLST results and generate combined output.
    
    Args:
        input_MLST_line_tuples: List of (line, type) tuples from MLST tools
        taxonomy_file: Path to taxonomy file
        mlst_db_path: Path to MLST database directory
        output_prefix: Optional output file prefix
    """
    # Extract sample name
    if output_prefix:
        isolate_name = output_prefix
    else:
        isolate_name = ".".join(Path(taxonomy_file).stem.split(".")[:-1]) if "." in Path(taxonomy_file).stem else Path(taxonomy_file).stem
    
    # Read taxonomy
    genus, species = read_taxonomy(taxonomy_file)
    
    # Get database pull date
    pull_date = get_pull_date(mlst_db_path)
    
    # Parse all MLST lines and expand combinations
    all_profiles = []
    
    for mlst_line, file_type in input_MLST_line_tuples:
        line = mlst_line.strip()
        print(line)
        if not line or line == "source_file  Database  ST  locus_1 locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 locus_8 locus_9  locus_10" or line == "source_file  Database  ST  locus_1 locus_2 locus_3 locus_4 locus_5 locus_6 locus_7 locus_8 lous_9  locus_10":
            continue
        
        items = line.split("\t")
        
        if file_type == "mlst":
            profiles = parse_mlst_line(items, "assembly", pull_date)
        elif file_type == "srst2":
            profiles = parse_mlst_line(items, "reads", pull_date)
        else:
            print(f"Unknown MLST filetype: {file_type}")
            continue
        
        all_profiles.extend(profiles)
    
    if not all_profiles:
        print("No MLST profiles found")
        write_output_file([], isolate_name, genus, species)
        return
    
    print(f"Found {len(all_profiles)} profile combinations")
    
    # Merge duplicates from different sources
    merged_profiles = merge_duplicate_profiles(all_profiles)
    print(f"After merging duplicates: {len(merged_profiles)} profiles")
    
    # Classify each profile
    for profile in merged_profiles:
        print("about to classify profile", profile)
        profile[2] = classify_profile_type(profile, mlst_db_path)
        
        # Check for paralogs
        if detect_paralogs(profile):
            profile[2] = f"{profile[2]}-PARALOG"
        
        print("classified profile as", profile[2])

        # Standardize database name
        profile[1] = standardize_database_name(profile[1])
    
    # Consolidate novel alleles
    final_profiles = consolidate_novel_alleles(merged_profiles)
    print(f"Final profiles: {len(final_profiles)}")
    
    # Sort by ST type (reverse to put numbers first)
    final_profiles.sort(key=lambda x: str(x[2]), reverse=True)
    
    # Write output
    write_output_file(final_profiles, isolate_name, genus, species)


def main():
    """Main entry point for command line usage."""
    args = parse_args()

    print(f"MLST Profile Processor")
    print(f"=" * 60)
    print(f"MLST file: {args.input}")
    if args.srst2:
        print(f"SRST2 file: {args.srst2}")
    print(f"Taxonomy: {args.taxonomy}")
    print(f"Database: {args.mlst_database}")
    print(f"=" * 60)
    print()
    
    # Read input files
    input_lines = []
    
    # Read MLST file
    mlst_lines = read_mlst_file(args.input)
    input_lines.extend(mlst_lines)
    print(f"Read {len(mlst_lines)} lines from MLST file")
    
    # Read SRST2 file if provided
    if args.srst2:
        srst2_lines = read_srst2_file(args.srst2)
        input_lines.extend(srst2_lines)
        print(f"Read {len(srst2_lines)} lines from SRST2 file")
    
    if not input_lines:
        print("Error: No MLST data found in input files")
        return 1
    
    print()
    
    # Process MLST data
    do_MLST_check(input_lines, args.taxonomy, args.mlst_database, args.output)
    
    print()
    print("Processing complete!")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())