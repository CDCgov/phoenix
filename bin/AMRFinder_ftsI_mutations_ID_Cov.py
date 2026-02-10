import pandas as pd
import sys
import os
import glob

def get_ftsi_mutations(project_dir, griphin_summary_file, identity_threshold, coverage_threshold):
    """
    Extract specific ftsI mutation data (I336I, N337N) from AMRFinder all_mutations files for E. coli samples
    
    Args:
        project_dir: Path to the project directory containing WGS_ID folders
        griphin_summary_file: Path to GRiPHin_Summary.tsv file
        identity_threshold: Minimum identity percentage (e.g., 90)
        coverage_threshold: Minimum coverage percentage (e.g., 90)
    
    Returns:
        DataFrame with WGS_ID, FastANI_Organism, and ftsI mutation information formatted as:
        ftsI_I336I (Sequence name), ftsI_N337N (Sequence name)
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
    
    # Filter for E. coli samples (case-insensitive check for "coli")
    ecoli_df = griphin_df[griphin_df['FastANI_Organism'].str.contains('Escherichia coli', case=False, na=False)].copy()
    
    # Select and rename columns
    wgs_data = ecoli_df[['WGS_ID', 'FastANI_Organism']].copy()
    wgs_data = wgs_data.rename(columns={'WGS_ID': 'ID', 'FastANI_Organism': 'FastANI Organism'})
    
    # Initialize result columns
    wgs_data['ftsI_mutations'] = ''
    
    # Target mutations to look for
    target_mutations = ['ftsI_I336IKYRI', 'ftsI_N337NYRIN']
    
    # Process each E. coli WGS_ID
    for idx, row in wgs_data.iterrows():
        wgs_id = row['ID']
        
        # Find the mutations file - use glob to handle wildcards
        mutations_pattern = os.path.join(project_dir, wgs_id, 'AMRFinder', f'{wgs_id}_all_mutations*.tsv')
        mutations_files = glob.glob(mutations_pattern)
        
        if mutations_files:
            mutations_file = mutations_files[0]  # Take the first match
            
            try:
                # Read the mutations file
                mutations_df = pd.read_csv(mutations_file, sep='\t')
                
                # Filter for specific ftsI gene mutations (I336I and N337N)
                if 'Gene symbol' in mutations_df.columns and 'Sequence name' in mutations_df.columns:
                    # Filter for the target mutations
                    ftsi_mutations = mutations_df[mutations_df['Gene symbol'].isin(target_mutations)].copy()
                    
                    if len(ftsi_mutations) > 0:
                        # Apply identity and coverage thresholds
                        filtered_mutations = ftsi_mutations[
                            (ftsi_mutations['% Identity to reference sequence'].astype(float) >= identity_threshold) &
                            (ftsi_mutations['% Coverage of reference sequence'].astype(float) >= coverage_threshold)
                        ]
                        
                        if len(filtered_mutations) > 0:
                            # Extract mutation information formatted as: ftsI_I336I (Sequence name)
                            mutation_list = []
                            seen_mutations = set()  # Track which mutations we've already added
                            
                            for _, mut_row in filtered_mutations.iterrows():
                                gene_symbol = mut_row.get('Gene symbol', '')
                                sequence_name = mut_row.get('Sequence name', '')
                                
                                if pd.notna(gene_symbol) and pd.notna(sequence_name):
                                    # Only add if we haven't seen this gene_symbol before
                                    if gene_symbol not in seen_mutations:
                                        # Format as: ftsI_I336I (penicillin-binding protein PBP3 FtsI [WILDTYPE])
                                        mutation_list.append(f"{gene_symbol} ({sequence_name})")
                                        seen_mutations.add(gene_symbol)
                            
                            # Store results
                            if mutation_list:
                                wgs_data.at[idx, 'ftsI_mutations'] = ', '.join(mutation_list)
                        else:
                            print(f"{wgs_id}: ftsI mutations (I336I/N337N) found but below threshold (ID>={identity_threshold}%, Cov>={coverage_threshold}%)")
                    else:
                        print(f"{wgs_id}: No ftsI I336I or N337N mutations found")
                else:
                    print(f"{wgs_id}: 'Gene symbol' or 'Sequence name' column not found in mutations file")
                    
            except Exception as e:
                print(f"Error processing {mutations_file}: {e}")
        else:
            print(f"Mutations file not found for {wgs_id}: {mutations_pattern}")
    
    # Filter out samples with no ftsI mutations
    wgs_data = wgs_data[wgs_data['ftsI_mutations'] != ''].reset_index(drop=True)
    
    return wgs_data

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
        print(f"Total E. coli samples processed: {len(wgs_data)}")
        print(f"Samples with ftsI mutations (I336I/N337N): {(wgs_data['ftsI_mutations'] != '').sum()}")
    else:
        print("Usage: python AMRFinder_ftsI-mutations_ID_Cov.py <project_dir> <griphin_summary_file> <identity_threshold> <coverage_threshold>")
        print("Example: python AMRFinder_ftsI-mutations_ID_Cov.py /path/to/project GRiPHin_Summary.tsv 90 90")
