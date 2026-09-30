#!/usr/bin/env python3

import pandas as pd
import numpy as np
import re
import os

def try_paths(path1, path2):
    """Function to try the first path, then fall back to second path if file doesn't exist"""
    if os.path.exists(path1):
        return path1
    elif os.path.exists(path2):
        return path2
    else:
        return path1  # Return the first path even if it doesn't exist, for consistent error handling

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

# Custom function to handle nulls better
def safe_str_convert(x):
    return str(x) if pd.notna(x) and x != '' else x

def print_df(df_toprint, label, all):
    pd.set_option('display.max_rows', 200)
    pd.set_option('display.max_columns', 500)
    pd.set_option('display.width', 1000)
    print(label+"  ------------------------------------------------------")
    for col in df_toprint.columns.tolist():
        print(col)
#    print("Columns:", df_toprint.columns)
    if all == True:
        with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000, 'display.colheader_justify', 'center', 'display.precision', 2, 'display.max_colwidth', 100):  # more options can be specified also
            print(df_toprint)

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
    clean_centar_df['CEMB RT Crosswalk'] = clean_centar_df['CEMB RT Crosswalk'].astype(str).replace('NaN', '') # Ensure leading zeroes are kept before it makes it to writing excel file
#    print_df(clean_centar_df, "Cleaned Centar DF", True)
#    print(clean_centar_df['CEMB RT Crosswalk'].values)
    #Replace empty strings with NaN and drop columns that are completely blank. ## .infer_objects(copy=False) ensures pandas infers correct types without creating a copy (retains legacy behavior)
    clean_centar_df = clean_centar_df.replace('', np.nan).infer_objects(copy=False).dropna(axis=1, how='all')
    #separate dataframes
    RB_type = [ "CEMB RT Crosswalk", "Inferred RT", "Probability", "ML Note", "Plasmid Info" ]
    RB_type_col = [col for col in clean_centar_df.columns if any(substring in col for substring in RB_type) ]
    RB_type_len = len(RB_type_col)
    A_B_Tox = [ "Toxinotype", "Toxin-A_sub-type", "tcdA", "Toxin-B_sub-type", "tcdB"]
    A_B_Tox_col = [col for col in clean_centar_df.columns if any(substring in col for substring in A_B_Tox) ]
    A_B_Tox_len = len(A_B_Tox_col)
    other_Tox = [ "tcdC", "tcdR", "tcdE", "cdtA", "cdtB", "cdtR", "cdtAB1", "cdtAB2", "cdt_NonTox", "PaLoc" ]
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

def create_centar_combined_df(directory1, sample_name, directory2):
    '''If Centar was run get info to add to the dataframe.'''
    # Filtering rows where Final_Taxa_ID contains 'shigella' or 'e.coli'
    # if there is a trailing / remove it
    directory1 = directory1.rstrip('/')
    # create file names
    print(directory1 + "/" + sample_name + "_centar_output.tsv")
    print(directory2 + "/" + sample_name + "/CENTAR/" + sample_name + "_centar_output.tsv")
    centar_summary = try_paths( directory1 + "/" + sample_name + "_centar_output.tsv", directory2 + "/" + sample_name + "/CENTAR/" + sample_name + "_centar_output.tsv" )
    print(centar_summary)
    #clean up the dataframe
    try: # handling for samples that failed and didn't get centar files created
#        centar_df = pd.read_csv(centar_summary, sep='\t', header=0)
        centar_df = pd.read_csv(centar_summary, sep='\t', header=0,
                converters={'CEMB RT Crosswalk': safe_str_convert})
        centar_df["WGS_ID"] = sample_name
#        print_df(centar_df, "Centar DF for " + sample_name, True)
    except FileNotFoundError:
        try: #retry looking at a different location
            centar_summary = "./" + sample_name + "_centar_output.tsv"
#            centar_df = pd.read_csv(centar_summary, sep='\t', header=0)
            centar_df = pd.read_csv(centar_summary, sep='\t', header=0,
                converters={'CEMB RT Crosswalk': safe_str_convert})
            centar_df["WGS_ID"] = sample_name
        except FileNotFoundError: 
            print("Warning: " + sample_name + "_centar_output.tsv file not found")
            # Add a row to the DataFrame with the WGS_ID column set to sample_name
            centar_df = pd.DataFrame({"WGS_ID": [sample_name]})
    # make NOT_FOUND and - to blank to keep inline with the AR/PF/HV calls.
    #if 'ML Note' in centar_df.columns:
    #    centar_df.replace("NOT_FOUND", "", inplace=True)
        #centar_df['column_name'] = centar_df['column_name'].str.replace("NOT_FOUND", "")
        #print(centar_df.columns)
        #print(centar_df)
        # Want to keep the dashes in the ML_Notes section so, temporarily replace them with double dashes
    #    centar_df.loc[centar_df['ML Note'] == '-', 'ML Note'] = '--'
    #    centar_df.replace("-", "", inplace=True)
        # Changiong back the double dashes to single dashes in the ML_Notes section
    #    centar_df.loc[centar_df['ML Note'] == '--', 'ML Note'] = '-'
    return centar_df

######################################## ShigaPass functions ##############################################

# def create_shiga_df(directory1, sample_name, shiga_df, taxa, directory2):
#     '''If Shigapass was run get info to add to the dataframe. This preserves ShigaPass's
#     raw call as-read (ShigaPass_Organism) for easy auditing/spot-checking, independent of
#     whatever check_taxa.py ultimately decided the final corrected taxa should be.'''
#     directory1 = directory1.rstrip('/')
#     if "Escherichia" in taxa or "Shigella" in taxa:
#         shiga_summary = try_paths( directory1 + "/" + sample_name + "_ShigaPass_summary.csv", directory2 + "/" + sample_name + "/ANI/" + sample_name + "_ShigaPass_summary.csv" )
#         row_data = { "WGS_ID": sample_name, "ShigaPass_Organism": ""}
#         try:
#             with open(shiga_summary) as shiga_file:
#                 for line in shiga_file.readlines()[1:]:
#                     if line.split(";")[9] == 'Not Shigella/EIEC\n':
#                         row_data["ShigaPass_Organism"] = "Not Shigella/EIEC"
#                     else:
#                         row_data["ShigaPass_Organism"] = line.split(";")[7]
#                     shiga_df = pd.concat([shiga_df, pd.DataFrame([row_data])], ignore_index=True)
#             mapping_dict = {
#                 r'^SB[\w-]*$': 'Shigella boydii',
#                 r'^SD[\w-]*$': 'Shigella dysenteriae',
#                 r'^SS[\w-]*$': 'Shigella sonnei',
#                 r'^SF[\w-]*$': 'Shigella flexneri'
#             }
#             shiga_df['ShigaPass_Organism'] = shiga_df['ShigaPass_Organism'].replace(mapping_dict, regex=True)
#         except FileNotFoundError:
#             print("Warning: ShigaPass file for " + sample_name + " not found")
#             new_row = pd.DataFrame({"WGS_ID": [sample_name], "ShigaPass_Organism": [""]})
#             shiga_df = pd.concat([shiga_df, new_row], ignore_index=True)
#     else:
#         new_row = pd.DataFrame({"WGS_ID": [sample_name], "ShigaPass_Organism": [""]})
#         shiga_df = pd.concat([shiga_df, new_row], ignore_index=True)
#     return shiga_df

def create_shiga_df(directory1, sample_name, shiga_df, taxa, tax_file, directory2):
    directory1 = directory1.rstrip('/')
    if "Escherichia" in taxa or "Shigella" in taxa:
        shiga_summary = try_paths( directory1 + "/" + sample_name + "_ShigaPass_summary.csv", directory2 + "/" + sample_name + "/ANI/" + sample_name + "_ShigaPass_summary.csv" )
        row_data = { "WGS_ID": sample_name, "ShigaPass_Organism": ""}
        try:
            with open(shiga_summary) as shiga_file:
                for line in shiga_file.readlines()[1:]:
                    fields = line.strip().split(";")
                    if any("Not Shigella/EIEC" in f for f in fields):
                        row_data["ShigaPass_Organism"] = "Not Shigella/EIEC"
                    else:
                        row_data["ShigaPass_Organism"] = fields[7] if len(fields) > 7 else ""
                    shiga_df = pd.concat([shiga_df, pd.DataFrame([row_data])], ignore_index=True)
            mapping_dict = {
                r'^SB[\w-]*$': 'Shigella boydii',
                r'^SD[\w-]*$': 'Shigella dysenteriae',
                r'^SS[\w-]*$': 'Shigella sonnei',
                r'^SF[\w-]*$': 'Shigella flexneri'
            }
            shiga_df['ShigaPass_Organism'] = shiga_df['ShigaPass_Organism'].replace(mapping_dict, regex=True)

            # Append ShigaPass-computed %ID;%Coverage from the tax file header, if present
            id_cov = None
            try:
                with open(tax_file, "r") as tf:
                    first_line = tf.readline()
                if first_line.startswith("ShigaPass\t"):
                    parts = first_line.split("\t")
                    if len(parts) > 1 and ";" in parts[1]:
                        id_cov = parts[1]
            except FileNotFoundError:
                pass
            if id_cov:
                mask = shiga_df["WGS_ID"] == sample_name
                shiga_df.loc[mask, "ShigaPass_Organism"] = shiga_df.loc[mask, "ShigaPass_Organism"] + ";" + id_cov

        except FileNotFoundError:
            new_row = pd.DataFrame({"WGS_ID": [sample_name], "ShigaPass_Organism": [""]})
            shiga_df = pd.concat([shiga_df, new_row], ignore_index=True)
    else:
        new_row = pd.DataFrame({"WGS_ID": [sample_name], "ShigaPass_Organism": [""]})
        shiga_df = pd.concat([shiga_df, new_row], ignore_index=True)
    return shiga_df


def get_corrected_taxa_from_tax_file(tax_file):
    '''Reads the already-corrected G:/s: genus and species directly from the .tax file,
    which check_taxa.py has already reconciled against ShigaPass/ANI. Used to populate
    Final_Taxa_ID for ShigaPass-adjudicated samples, separately from the raw
    ShigaPass_Organism audit column (which intentionally preserves ShigaPass's literal
    output, including negative results like "Not Shigella/EIEC", for spot-checking).'''
    genus = None
    species = None
    try:
        with open(tax_file, "r") as f:
            for line in f:
                if line.startswith("G:"):
                    genus = line.split("\t")[1].strip()
                elif line.startswith("s:"):
                    species = line.split("\t")[1].strip()
    except FileNotFoundError:
        return ""
    if genus and species:
        return f"{genus} {species}"
    return genus or ""


def double_check_taxa_id(shiga_df, phx_df, tax_files_by_sample):
    '''tax_files_by_sample: dict mapping WGS_ID -> tax_file path, so Final_Taxa_ID can be
    populated from the corrected tax file rather than ShigaPass_Organism's raw value.'''
    merged_df = pd.merge(phx_df, shiga_df, on='WGS_ID', how='left')
    insert_position = merged_df.columns.get_loc("FastANI_Organism")
    columns = list(merged_df.columns)
    new_columns = ['ShigaPass_Organism']
    columns_reordered = (
        columns[:insert_position] +
        new_columns +
        [col for col in columns if col not in new_columns and col not in columns[:insert_position]])
    merged_df = merged_df[columns_reordered]
    merged_df['Final_Taxa_ID'] = merged_df.apply(
        lambda row: fill_taxa_id(row, tax_files_by_sample.get(row['WGS_ID'], "")), axis=1
    )
    return merged_df


def fill_taxa_id(row, tax_file=None):
    if row['Taxa_Source'] == 'ANI_REFSEQ':
        return row['FastANI_Organism']
    elif row['Taxa_Source'] == 'kraken2_wtasmbld':
        genus = row['Kraken_ID_WtAssembly_%'].split(" ")[0]
        species = row['Kraken_ID_WtAssembly_%'].split(" ")[2]
        return genus + " " + species
    elif row['Taxa_Source'] == 'kraken2_trimmed':
        if 'Kraken_ID_Trimmed_Reads_%' in row and row['Kraken_ID_Trimmed_Reads_%']:
            target_column = row['Kraken_ID_Trimmed_Reads_%']
        elif 'Kraken_ID_Raw_Reads_%' in row and row['Kraken_ID_Raw_Reads_%']:
            target_column = row['Kraken_ID_Raw_Reads_%']
        genus = target_column.split(" ")[0]
        species = target_column.split(" ")[2]
        return genus + " " + species
    elif row['Taxa_Source'].startswith('ShigaPass'):
        # Final_Taxa_ID comes from the corrected tax file, NOT the raw ShigaPass_Organism
        # column -- that column intentionally preserves ShigaPass's literal output
        # (including "Not Shigella/EIEC") for auditing, which isn't a displayable organism name.
        if tax_file:
            corrected = get_corrected_taxa_from_tax_file(tax_file)
            if corrected:
                return corrected
        return row['ShigaPass_Organism']  # fallback if tax_file unavailable
    else:
        return 'Unknown'

#def main():
#    directory = "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Jill_DIR/PHX_v2/v2.2.0-dev/centar/cdc_centar_newer"
#    sample_names = [ "2022GL-00907", "2022GL-00947", "2022GL-01162" ]
#    centar_dfs = []
#    for sample_name in sample_names:
#        centar_df = create_centar_combined_df(directory, sample_name)
#        centar_dfs.append(centar_df)
#    full_centar_df = pd.concat(centar_dfs, ignore_index=True)
#    ordered_centar_df = clean_and_format_centar_dfs(full_centar_df)
