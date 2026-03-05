import pandas as pd
import sys
import os
import glob
import re

def is_non_wildtype_ftsi(gene_symbol):
    """
    Check if a ftsI gene symbol represents a non-wild type mutation.
    Non-wild type means the amino acid before the number differs from the amino acid(s) after.
    
    Examples:
        ftsI_N337NYRIN -> "N" vs "NYRIN" -> different -> non-wild type (True)
        ftsI_I336IKYRI -> "I" vs "IKYRI" -> different -> non-wild type (True)
        ftsI_I336I -> "I" vs "I" -> same -> wild type (False)
        ftsI_N337N -> "N" vs "N" -> same -> wild type (False)
    """
    if not gene_symbol or not gene_symbol.startswith('ftsI_'):
        return False
    
    # Extract the mutation part after "ftsI_"
    mutation_part = gene_symbol[5:]  # Remove "ftsI_" prefix
    
    # Pattern: one or more letters, followed by numbers, followed by one or more letters
    # e.g., "N337NYRIN" -> before="N", number="337", after="NYRIN"
    match = re.match(r'^([A-Za-z]+)(\d+)([A-Za-z]+)$', mutation_part)
    
    if match:
        before = match.group(1)  # Amino acid(s) before the number
        after = match.group(3)   # Amino acid(s) after the number
        
        # Non-wild type if before != after
        return before != after
    
    return False

def get_ftsi_mutations(project_dir, griphin_summary_file, identity_threshold, coverage_threshold):
    """
    Extract specific ftsI mutation data (I336I, N337N) from AMRFinder all_mutations files.
    Only E. coli samples will have ftsI mutations processed; other species will have empty values.
    
    Args:
        project_dir: Path to the project directory containing WGS_ID folders
        griphin_summary_file: Path to GRiPHin_Summary.tsv file
        identity_threshold: Minimum identity percentage (e.g., 90)
        coverage_threshold: Minimum coverage percentage (e.g., 90)
    
    Returns:
        DataFrame with ID and ftsI_mutations columns.
        Non-E.coli samples will have empty ftsI_mutations.
    """
    # Read the GRiPHin summary to get WGS_ID and FastANI_Organism
    griphin_df = pd.read_csv(griphin_summary_file, sep='\t')
    
    # Remove extra notes at the end of the file
    # Find first index where the entire row is NaN
    nan_indices = griphin_df.index[griphin_df.isna().all(axis=1)]
    if len(nan_indices) > 0:
        idx = nan_indices[0]
        # Keep everything *before* that row
        griphin_df = griphin_df.loc[:idx-1]
    
    # Also remove rows where WGS_ID is NaN or empty
    griphin_df = griphin_df[griphin_df['WGS_ID'].notna()]
    griphin_df = griphin_df[griphin_df['WGS_ID'] != '']
    
    # Select only WGS_ID and FastANI_Organism columns, rename WGS_ID to ID
    wgs_data = griphin_df[['WGS_ID', 'FastANI_Organism']].copy()
    wgs_data = wgs_data.rename(columns={'WGS_ID': 'ID'})
    
    # Initialize ftsI_mutations column as empty for all samples
    wgs_data['ftsI_mutations'] = ''
    
    # Identify E. coli samples (case-insensitive check)
    ecoli_mask = wgs_data['FastANI_Organism'].str.contains('Escherichia coli', case=False, na=False)
    
    # Process only E. coli samples
    for idx, row in wgs_data[ecoli_mask].iterrows():
        wgs_id = row['ID']
        
        # Find the mutations file - use glob to handle wildcards
        mutations_pattern = os.path.join(project_dir, wgs_id, 'AMRFinder', f'{wgs_id}_all_mutations*.tsv')
        mutations_files = glob.glob(mutations_pattern)
        
        if mutations_files:
            mutations_file = mutations_files[0]  # Take the first match
            
            try:
                # Read the mutations file
                mutations_df = pd.read_csv(mutations_file, sep='\t')
                
                # Filter for ftsI gene mutations that are non-wild type
                if 'Gene symbol' in mutations_df.columns and 'Sequence name' in mutations_df.columns:
                    # Filter for ftsI genes first, then check for non-wild type
                    ftsi_mutations = mutations_df[
                        mutations_df['Gene symbol'].str.startswith('ftsI_', na=False) &
                        mutations_df['Gene symbol'].apply(is_non_wildtype_ftsi)
                    ].copy()
                    
                    if len(ftsi_mutations) > 0:
                        # Apply identity and coverage thresholds
                        filtered_mutations = ftsi_mutations[
                            (ftsi_mutations['% Identity to reference sequence'].astype(float) >= identity_threshold) &
                            (ftsi_mutations['% Coverage of reference sequence'].astype(float) >= coverage_threshold)
                        ]
                        
                        if len(filtered_mutations) > 0:
                            # Extract mutation information - only gene_symbol
                            mutation_list = []
                            seen_mutations = set()  # Track which mutations we've already added
                            
                            for _, mut_row in filtered_mutations.iterrows():
                                gene_symbol = mut_row.get('Gene symbol', '')
                                
                                if pd.notna(gene_symbol):
                                    # Only add if we haven't seen this gene_symbol before
                                    if gene_symbol not in seen_mutations:
                                        mutation_list.append(gene_symbol)
                                        seen_mutations.add(gene_symbol)
                            
                            # Store results
                            if mutation_list:
                                wgs_data.at[idx, 'ftsI_mutations'] = ', '.join(mutation_list)
                        else:
                            print(f"{wgs_id}: ftsI non-wild type mutations found but below threshold (ID>={identity_threshold}%, Cov>={coverage_threshold}%)")
                    else:
                        print(f"{wgs_id}: No ftsI non-wild type mutations found")
                else:
                    print(f"{wgs_id}: 'Gene symbol' or 'Sequence name' column not found in mutations file")
                    
            except Exception as e:
                print(f"Error processing {mutations_file}: {e}")
        else:
            print(f"Mutations file not found for {wgs_id}: {mutations_pattern}")
    
    # Return only ID and ftsI_mutations columns (all samples, non-E.coli will have empty ftsI_mutations)
    return wgs_data[['ID', 'ftsI_mutations']].reset_index(drop=True)

if __name__ == "__main__":
    if len(sys.argv) == 5:
        project_dir = sys.argv[1]
        griphin_summary_file = sys.argv[2]
        identity_threshold = float(sys.argv[3])
        coverage_threshold = float(sys.argv[4])
        
        # Run the get_ftsi_mutations function
        wgs_data = get_ftsi_mutations(project_dir, griphin_summary_file, identity_threshold, coverage_threshold)
        
        # Save to CSV
        output_csv_file = f"ftsI_mutations_ID{identity_threshold}_Cov{coverage_threshold}.csv"
        wgs_data.to_csv(output_csv_file, index=False)
        print(f"\nResults saved to: {output_csv_file}")
        print(f"Total samples: {len(wgs_data)}")
        print(f"Samples with ftsI non-wild type mutations: {(wgs_data['ftsI_mutations'] != '').sum()}")
    else:
        print("Usage: python AMRFinder_ftsI-mutations_ID_Cov.py <project_dir> <griphin_summary_file> <identity_threshold> <coverage_threshold>")
        print("Example: python AMRFinder_ftsI-mutations_ID_Cov.py /path/to/project GRiPHin_Summary.tsv 90 90")
