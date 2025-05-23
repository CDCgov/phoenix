#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re

##################################  Centar functions   #########################################
def transform_value(value):
    # Split the value into components
    if (isinstance(value, float) and np.isnan(value)):
        return ""
    elif value != "NA|NA|NA" and value != "[NA|NA|NA]":
        parts = value.split('|')
        if len(parts) == 3:
            # Reorder and format the components into the desired format
            nuc_identity = parts[0][:-2].replace("[","")   # Extract the number from '98NT'
            aa_identity = parts[1][:-2]   # Extract the number from '98AA'
            coverage = parts[2].replace("COV]","")           # Extract the coverage number
            # Format the new string
            return f'[{nuc_identity}NT/{aa_identity}AA/{coverage}]G'
        elif len(parts) == 2:
            # Reorder and format the components into the desired format
            nuc_identity = parts[0][:-2].replace("[","")   # Extract the number from '98NT'
            coverage = parts[1].replace("COV[]","")           # Extract the coverage number
            # Format the new string
            return f'[{nuc_identity}NT/{coverage}]G'
    else:
        return ""
    return value

# ******Note: this function is called in GRiPHin.py any changes to it need to be tested in -entry CDC_PHOENIX and PHOENIX******
def clean_and_format_centar_dfs(centar_df):
    '''If Centar was run get info to add to the dataframe.'''
    cols_to_transform = [x for x in centar_df.columns if '[%Nuc_Identity' in x ]
    for col in cols_to_transform:
        centar_df[col] = centar_df[col].apply(transform_value)
    #drop presence/absence columns
    columns_to_drop = [col for col in centar_df.columns if 'presence'  in col ]
    clean_centar_df = centar_df.drop(columns=columns_to_drop)
    # Remove the substring from all column headers
    clean_centar_df.rename(columns=lambda x: re.sub(r'\[%Nuc_Identity \| %AA_Identity \| %Coverage\]', '', x).strip(), inplace=True)
    clean_centar_df.rename(columns=lambda x: re.sub(r'\[%Nuc_Identity \| %Coverage\]', '', x).strip(), inplace=True)
    clean_centar_df.rename(columns=lambda x: re.sub(r'Diffbase_', '', x).strip(), inplace=True)
    #Replace empty strings with NaN and drop columns that are completely blank. ## .infer_objects(copy=False) ensures pandas infers correct types without creating a copy (retains legacy behavior)
    clean_centar_df = clean_centar_df.replace('', np.nan).infer_objects(copy=False).dropna(axis=1, how='all')
    #separate dataframes
    RB_type = [ "CEMB RT Crosswalk", "Inferred RT", "Probability", "ML Note", "Plasmid Info" ]
    RB_type_col = [col for col in clean_centar_df.columns if any(substring in col for substring in RB_type) ]
    RB_type_len = len(RB_type_col)
    A_B_Tox = [ "Toxinotype", "Toxin-A_sub-type", "tcdA", "Toxin-B_sub-type", "tcdB"]
    A_B_Tox_col = [col for col in clean_centar_df.columns if any(substring in col for substring in A_B_Tox) ]
    A_B_Tox_len = len(A_B_Tox_col)
    other_Tox = [ "tcdC", "tcdR", "tcdE", "cdtA", "cdtB", "cdtR", "cdtAB1", "cdtAB2", "non-tox", "PaLoc" ]
    other_Tox_col = [col for col in clean_centar_df.columns if any(substring in col for substring in other_Tox) ]
    other_Tox_len = len(other_Tox_col)
    mutants = [ 'gyr','dac','feo','fur','gdp','gly','hem','hsm','Isc','mur', 'mur','nifJ','PNim','rpo','sda','thi','Van','mutations' ]
    mutations_col = [col for col in clean_centar_df.columns if any(substring in col for substring in mutants) ]
    # List of mutation names to remove
    mutations_to_remove = ['tcdC other mutations', 'cdtR other mutations', 'PaLoc_NonTox other mutations']
    # Remove each mutation name if it exists in mutations_col
    mutations_col = [mutation for mutation in mutations_col if mutation not in mutations_to_remove]
    mutant_len = len(mutations_col)
    #if "MLST Clade" in clean_centar_df.columns:
    existing_columns_in_order = ["MLST Clade"] + A_B_Tox_col + other_Tox_col + mutations_col + RB_type_col
    #else:
    #   existing_columns_in_order = A_B_Tox_col + other_Tox_col + mutations_col + RB_type_col 
    #if clean_centar_df.empty or clean_centar_df.columns.tolist() == ["WGS_ID"]: #for cases where centar wasn't run for that sample - not c. diff or a qc failure sample
    #    clean_centar_df = pd.DataFrame(columns = existing_columns_in_order) # Assign the headers to the DataFrame
    ordered_centar_df = clean_centar_df[existing_columns_in_order]
    return ordered_centar_df, A_B_Tox_len, other_Tox_len, mutant_len, RB_type_len

def create_centar_combined_df(directory, sample_name):
    '''If Centar was run get info to add to the dataframe.'''
    # Filtering rows where Final_Taxa_ID contains 'shigella' or 'e.coli'
    # if there is a trailing / remove it
    directory = directory.rstrip('/')
    # create file names
    centar_summary = directory + "/CENTAR/" + sample_name + "_centar_output.tsv"
    #clean up the dataframe
    try: # handling for samples that failed and didn't get centar files created
        centar_df = pd.read_csv(centar_summary, sep='\t', header=0)
        centar_df["WGS_ID"] = sample_name
    except FileNotFoundError:
        try: #retry looking at a different location
            centar_summary = "./" + sample_name + "_centar_output.tsv"
            centar_df = pd.read_csv(centar_summary, sep='\t', header=0)
            centar_df["WGS_ID"] = sample_name
        except FileNotFoundError: 
            print("Warning: " + sample_name + "_centar_output.tsv file not found")
            # Add a row to the DataFrame with the WGS_ID column set to sample_name
            centar_df = pd.DataFrame({"WGS_ID": [sample_name]})
    # make NOT_FOUND and - to blank to keep inline with the AR/PF/HV calls.
    centar_df.replace("NOT_FOUND", "", inplace=True)
    centar_df.replace("-", "", inplace=True)
    return centar_df

######################################## ShigaPass functions ##############################################
def create_shiga_df(directory, sample_name, shiga_df, taxa):
    '''If Shigapass was run get info to add to the dataframe.'''
    # if there is a trailing / remove it
    directory = directory.rstrip('/')
    if "Escherichia" in taxa or "Shigella" in taxa:
        # create file names
        shiga_summary = directory + "/ANI/" + sample_name + "_ShigaPass_summary.csv"
        # Create a dictionary to store row information
        row_data = { "WGS_ID": sample_name, "ShigaPass_Organism": ""}
        row_data["WGS_ID"] = sample_name
        with open(shiga_summary) as shiga_file:
            for line in shiga_file.readlines()[1:]:  # Skip the first line
                if line.split(";")[9] == 'Not Shigella/EIEC\n':
                    row_data["ShigaPass_Organism"] = "Not Shigella/EIEC"
                else:
                    row_data["ShigaPass_Organism"] = line.split(";")[7]
                # Convert the row data into a DataFrame and concatenate with the main DataFrame
                shiga_df = pd.concat([shiga_df, pd.DataFrame([row_data])], ignore_index=True)
        # Define the mapping of short strings to longer strings
        mapping_dict = {'SB': 'Shigella boydii',
                        'SD': 'Shigella dysenteriae',
                        'SS': 'Shigella sonnei',
                        'SF1-5': 'Shigella flexneri'}
        # Apply the mapping using map()
        shiga_df['ShigaPass_Organism'] = shiga_df['ShigaPass_Organism'].replace(mapping_dict)
    else:
        # Add a row to the DataFrame with the WGS_ID column set to sample_name
        shiga_df = pd.DataFrame({"WGS_ID": [sample_name]})
    return shiga_df

def double_check_taxa_id(shiga_df, phx_df):
    # Merge the DataFrames on 'WGS_ID'
    print(shiga_df.columns)
    merged_df = pd.merge(phx_df, shiga_df, on='WGS_ID', how='left')
    # Identify the position of the insertion point
    insert_position = merged_df.columns.get_loc("FastANI_Organism")
    # Reorder the columns: place the new columns at the desired position
    columns = list(merged_df.columns)
    # Reorder columns to insert the new columns between 'Column_A' and 'Column_B'
    new_columns = ['ShigaPass_Organism']
    # Reorder columns: place the new columns between 'Column_A' and 'Column_B'
    columns_reordered = (
        columns[:insert_position] +  # Columns before the insertion point
        new_columns +                # New columns to be inserted
        [col for col in columns if col not in new_columns and col not in columns[:insert_position]]) # Remaining columns
    # Reorder the merged DataFrame columns
    merged_df = merged_df[columns_reordered]
    # Apply the custom function to fill the Taxa_ID column
    merged_df['Final_Taxa_ID'] = merged_df.apply(fill_taxa_id, axis=1)
    return merged_df

# Define the custom function to update the Taxa_ID based on conditions
def fill_taxa_id(row):
    if row['Taxa_Source'] == 'ANI_REFSEQ':
        return row['FastANI_Organism']
    elif row['Taxa_Source'] == 'kraken2_wtasmbld':
        genus = row['Kraken_ID_WtAssembly_%'].split(" ")[0]
        species = row['Kraken_ID_WtAssembly_%'].split(" ")[2]
        return genus + " " + species
    elif row['Taxa_Source'] == 'kraken2_trimmed':
        genus = row['Kraken_ID_Raw_Reads_%'].split(" ")[0]
        species = row['Kraken_ID_Raw_Reads_%'].split(" ")[2]
        return genus + " " + species
    elif row['Taxa_Source'] == 'ShigaPass':
        return row['ShigaPass_Organism']
    else:
        return 'Unknown'  # Default case if no condition matches

#def main():
#    directory = "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Jill_DIR/PHX_v2/v2.2.0-dev/centar/cdc_centar_newer"
#    sample_names = [ "2022GL-00907", "2022GL-00947", "2022GL-01162" ]
#    centar_dfs = []
#    for sample_name in sample_names:
#        centar_df = create_centar_combined_df(directory, sample_name)
#        centar_dfs.append(centar_df)
#    full_centar_df = pd.concat(centar_dfs, ignore_index=True)
#    ordered_centar_df = clean_and_format_centar_dfs(full_centar_df)
