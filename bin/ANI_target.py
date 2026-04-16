import os
import glob
import pandas as pd
import sys
import argparse
import re

def get_ani_targets(project_dir, ani_cutoff, coverage_cutoff=70.0):
    """
    Scans project directory for ANI files and reports species > cutoff.
    
    Args:
        project_dir: Path to the project directory containing sample folders.
        ani_cutoff: Minimum ANI percentage (e.g., 95.0).
        coverage_cutoff: Minimum coverage percentage (default: 70.0).
        
    Returns:
        DataFrame with 'ID' and 'ANI_Target_Species' columns.
    """
    
    results = []

    # Iterate through all subdirectories in the project folder
    # Assuming each subdirectory is a sample ID
    if not os.path.exists(project_dir):
        print(f"Error: Project directory '{project_dir}' does not exist.")
        return pd.DataFrame()

    print(f"Scanning {project_dir} for ANI files...")
    print(f"ANI cutoff: {ani_cutoff}%, Coverage cutoff: {coverage_cutoff}%")

    # Get list of immediate subdirectories (sample IDs)
    sample_dirs = sorted([d for d in os.listdir(project_dir) if os.path.isdir(os.path.join(project_dir, d))])

    for sample_id in sample_dirs:
        # Construct path to the ANI file
        # Pattern: project_dir/sample_id/ANI/sample_id_REFSEQ_*.ani.txt
        ani_folder = os.path.join(project_dir, sample_id, 'ANI')
        
        # Use glob to find the specific .ani.txt file as the date part might change
        ani_file_pattern = os.path.join(ani_folder, f"{sample_id}_REFSEQ_*.ani.txt")
        ani_files = glob.glob(ani_file_pattern)

        # default empty
        species_str = ""

        if ani_files:
            ani_file = ani_files[0]  # Take the first match
            try:
                # Read as headerless tab-separated (ANI output format)
                # Columns: Query, Reference, ANI, AlignedFragments, TotalFragments
                df = pd.read_csv(ani_file, header=None, sep='\t', engine='python')
                if df.shape[1] >= 5:
                    df.columns = ['Query', 'Reference', 'ANI', 'AlignedFragments', 'TotalFragments']
                elif df.shape[1] >= 3:
                    extra_cols = [f'col{i}' for i in range(4, df.shape[1] + 1)]
                    df.columns = ['Query', 'Reference', 'ANI'] + extra_cols
                else:
                    raise ValueError(f"Unexpected ANI file format ({df.shape[1]} columns): {ani_file}")

                # Ensure ANI is numeric; convert proportion to percent if needed
                df['ANI'] = pd.to_numeric(df['ANI'], errors='coerce')
                if df['ANI'].max(skipna=True) is not None and df['ANI'].max() <= 1.0:
                    df['ANI'] = df['ANI'] * 100.0

                # Calculate coverage from AlignedFragments / TotalFragments * 100
                if 'AlignedFragments' in df.columns and 'TotalFragments' in df.columns:
                    df['AlignedFragments'] = pd.to_numeric(df['AlignedFragments'], errors='coerce')
                    df['TotalFragments'] = pd.to_numeric(df['TotalFragments'], errors='coerce')
                    df['Coverage'] = (df['AlignedFragments'] / df['TotalFragments'] * 100).round(2)
                else:
                    df['Coverage'] = None

                # Extract species from reference path (basename, take first two tokens before GCF/ASM)
                def ref_to_species(refpath):
                    b = os.path.basename(str(refpath))
                    b = re.sub(r'\.fna(\.gz)?$', '', b)
                    parts = b.split('_')
                    stop = next((i for i, p in enumerate(parts) if p.upper().startswith(('GCF', 'ASM'))), len(parts))
                    name_parts = parts[:min(2, stop)] if stop > 0 else parts[:2]
                    return '_'.join(name_parts)

                df['SpeciesName'] = df['Reference'].apply(ref_to_species)

                # Filter out unspecific species (e.g., Acinetobacter_sp., Citrobacter_sp.)
                df = df[~df['SpeciesName'].str.endswith('_sp.')]

                # Filter by ANI cutoff AND coverage cutoff
                ani_mask = df['ANI'] >= ani_cutoff
                if 'Coverage' in df.columns and df['Coverage'].notna().any():
                    cov_mask = df['Coverage'] >= coverage_cutoff
                    report_df = df[ani_mask & cov_mask].copy()
                else:
                    report_df = df[ani_mask].copy()

                if not report_df.empty:
                    # If same species appears multiple times, keep the row with max ANI
                    # Use idxmax to keep corresponding coverage for the max ANI row
                    idx_max = report_df.groupby('SpeciesName')['ANI'].idxmax()
                    max_df = report_df.loc[idx_max].sort_values('ANI', ascending=False)
                    
                    # Exclude the top hit (best match) — it is already shown in
                    # FastANI_Organism / FastANI Match (%) / Bases Aligned columns of taxa table
                    if len(max_df) > 1:
                        additional_df = max_df.iloc[1:]  # Skip first row (best match)
                        entries = []
                        for _, row in additional_df.iterrows():
                            if pd.notna(row.get('Coverage')):
                                entries.append(f"{row['SpeciesName']}({row['ANI']:.2f}%; {row['Coverage']:.2f}%)")
                            else:
                                entries.append(f"{row['SpeciesName']}({row['ANI']:.2f}%)")
                        species_str = "|".join(entries)
                    else:
                        # Only one species found (the top hit) — nothing additional to report
                        species_str = ""
                else:
                    print(f"{sample_id}: No ANI hits >= {ani_cutoff}% ANI and >= {coverage_cutoff}% coverage.")
                    species_str = ""
            except Exception as e:
                print(f"Error processing {sample_id}: {e}")
                species_str = ""
        else:
            print(f"No ANI file found for {sample_id}. Expected pattern: {ani_file_pattern}")
            continue

        results.append({
            'ID': sample_id,
            'ANI_Target_Species': species_str
        })

    # Create DataFrame
    results_df = pd.DataFrame(results)
    return results_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract species from ANI files above a specific cutoff.")
    parser.add_argument("-d", "--dir", required=True, help="Path to the project directory containing sample folders")
    parser.add_argument("-c", "--cutoff", type=float, default=95.0, help="ANI cutoff percentage (default: 95.0)")
    parser.add_argument("-v", "--coverage", type=float, default=70.0, help="Coverage cutoff percentage (default: 70.0)")
    parser.add_argument("-o", "--output", help="Output CSV filename (optional)")

    args = parser.parse_args()

    df = get_ani_targets(args.dir, args.cutoff, args.coverage)
    
    if not df.empty:
        # Determine output filename
        if args.output:
            outfile = args.output
        else:
            outfile = f"ANI_targets_{args.cutoff}.csv"
            
        df.to_csv(outfile, index=False)
        print(f"Results saved to {outfile}")
        
        # summary stats
        print(f"Total samples processed: {len(df)}")
        print(f"Samples with targets >= {args.cutoff}% ANI and >= {args.coverage}% coverage: {len(df[df['ANI_Target_Species'] != ''])}")
    else:
        print("No data processed.")