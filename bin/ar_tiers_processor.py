import pandas as pd
import sys
import os
import glob

def get_ar_tiers(project_dir, griphin_summary_file):
    """
    Extract AR gene data from abritamr.txt files and categorize into tiers
    """
    # Read the GRiPHin summary to get WGS_ID and FastANI_Organism
    griphin_df = pd.read_csv(griphin_summary_file, sep='\t')
    wgs_data = griphin_df[['WGS_ID', 'FastANI_Organism']].copy()
    # Find the index of the first row where all values are NA
    na_row_index = wgs_data.index[wgs_data.isna().all(axis=1)][0]
    # Keep only rows BEFORE that row
    wgs_data = wgs_data.loc[:na_row_index - 1]
    wgs_data = wgs_data.rename(columns={'WGS_ID': 'ID', 'FastANI_Organism': 'Species'})
    # Hardcoded exception genes - add your exception genes here
    ex_gene = [
        'vanH-A',
        'blaTEM-1',
        'blaVEB-6'
    ]
    
    # Define tier columns based on abritamr.txt structure
    tier_1_columns = ['Carbapenemase (MBL)', 'Carbapenemase', 'Methicillin', 'Vancomycin']
    tier_2_columns = ['AmpC', 'ESBL']
    tier_3_columns = ['Carbapenemase (OXA-51 family)', 'Beta-lactam', 'Cephalosporin', 
                     'Aminoglycosides (Ribosomal methyltransferase)', 
                     'Azithromycin/Erythromycin/Spiramycin/Telithromycin',
                     'Chloramphenicol/Clindamycin/Florfenicol/Linezolid/Streptogramin b/Tiamulin/Virginiamycin',
                     'Florfenicol/Oxazolidinone']
    
    # Initialize result columns
    wgs_data['Tier 1'] = ''
    wgs_data['Tier 2'] = ''
    wgs_data['Tier 3'] = ''
    wgs_data['other'] = ''
    
    # Process each WGS_ID
    for idx, row in wgs_data.iterrows():
        wgs_id = str(row['ID'])
        abritamr_file = os.path.join(project_dir, wgs_id, 'AMRFinder', f'{wgs_id}.abritamr.txt')
        if os.path.exists(abritamr_file):
            try:
                # Read abritamr file
                abritamr_df = pd.read_csv(abritamr_file, sep='\t')
                
                # Exclude 'Isolate' column from processing
                if 'Isolate' in abritamr_df.columns:
                    abritamr_df = abritamr_df.drop(columns=['Isolate'])
                
                # Get the data row (assuming first data row after header)
                if len(abritamr_df) > 0:
                    data_row = abritamr_df.iloc[0]
                    
                    def process_tier_genes(tier_columns):
                        tier_genes = []
                        for col in tier_columns:
                            if col in data_row and pd.notna(data_row[col]) and str(data_row[col]).strip():
                                genes = [g.strip() for g in str(data_row[col]).split(',') if g.strip()]
                                tier_genes.extend(genes)
                        return tier_genes
                    
                    # Process each tier
                    tier_1_genes = process_tier_genes(tier_1_columns)
                    tier_2_genes = process_tier_genes(tier_2_columns)
                    tier_3_genes = process_tier_genes(tier_3_columns)
                    
                    # Process other genes (columns not in tier 1, 2, or 3)
                    other_genes = []
                    all_tier_columns = tier_1_columns + tier_2_columns + tier_3_columns
                    for col in data_row.index:
                        if col not in all_tier_columns and pd.notna(data_row[col]) and str(data_row[col]).strip():
                            genes = [g.strip() for g in str(data_row[col]).split(',') if g.strip()]
                            other_genes.extend(genes)
                    
                    # Apply exception gene formatting and join with |
                    def format_genes(gene_list):
                        if not gene_list:
                            return ''
                        # Remove special characters like *, ^, etc. and format
                        cleaned_genes = []
                        for g in gene_list:
                            clean_g = g.replace('*', '').replace('^', '').strip()
                            if clean_g:
                                formatted_g = f"{{{clean_g}}}" if clean_g in ex_gene else clean_g
                                cleaned_genes.append(formatted_g)
                        return '|'.join(cleaned_genes)
                    
                    wgs_data.at[idx, 'Tier 1'] = format_genes(tier_1_genes)
                    wgs_data.at[idx, 'Tier 2'] = format_genes(tier_2_genes)
                    wgs_data.at[idx, 'Tier 3'] = format_genes(tier_3_genes)
                    wgs_data.at[idx, 'other'] = format_genes(other_genes)
                    
            except Exception as e:
                print(f"Error processing {abritamr_file}: {e}")
        else:
            print(f"File not found: {abritamr_file}")
    return wgs_data

if __name__ == "__main__":
    if len(sys.argv) == 4:
        project_dir = sys.argv[1]
        griphin_summary_file = sys.argv[2]
        output_csv_file = sys.argv[3]
        
        # Run the get_ar_tiers function and save wgs_data to CSV
        wgs_data = get_ar_tiers(project_dir, griphin_summary_file)
        wgs_data.to_csv(output_csv_file, index=False)
        print(f"Saved to: {output_csv_file}")
    else:
        print("Usage: python csc_category_v2.py <project_dir> <griphin_summary_file> <output_csv_file>")