#!/usr/bin/env python3

#disable cache usage in the Python so __pycache__ isn't formed. If you don't do this using 'nextflow run cdcgov/phoenix...' a second time will causes and error
import sys
sys.dont_write_bytecode = True
import glob
import pandas as pd
import numpy as np
import argparse
import xlsxwriter as ws
from xlsxwriter.utility import xl_rowcol_to_cell
from openpyxl import load_workbook
from itertools import chain
from GRiPHin import order_ar_gene_columns, Combine_dfs, big5_check, write_to_excel, convert_excel_to_tsv, sort_qc_through_spec2_dataframe

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python GRiPHin.py -g1 ./samplesheet.csv -a ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output --phoenix --scaffolds
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='''Script to create new griphin excel sheet by combining two griphin summaries. The -g2 is considered the "new" file and thus when 
    samples with data in the -g1 file will have values overwritten to be the values in the -g2 file.''')
    parser.add_argument('-g1', '--old_griphin', default=None, required=False, dest='griphin_old', help='The first griphin excel file to combine.')
    parser.add_argument('-g2', '--new_griphin', required=False, dest='griphin_new', help='The second griphin excel file to combine.')
    parser.add_argument('-o', '--output', required=False, default=None,dest='output', help='Name of output file default is GRiPHin_Summary.xlsx.')
    parser.add_argument('-s', '--samplesheet', required=False, default=None,dest='samplesheet', help='samplesheet with sample,directory columns. Used to doublecheck sample names.')
    parser.add_argument('--griphin_list', required=False, action='store_true',default=False, dest='griphin_list', help='pass instead of -g1/-g2 when you want to combine more than 2 griphins.')
    parser.add_argument('--coverage', default=30, required=False, dest='set_coverage', help='The coverage cut off default is 30x.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def read_excels(file_path1, file_path2, samplesheet):
    #get number of footer lines for each file to skip when reading in as a pd
    footer_lines1 = detect_footer_lines(file_path1)
    footer_lines2 = detect_footer_lines(file_path2)
    # Read the Excel file, skipping the first row and using the second row as the header
    try: #check that this is an excel file
        df_1 = pd.read_excel(file_path1,
            skiprows=1,  # Skip the first header row
            header=0,    # Use the second row as the header
            skipfooter=footer_lines1,engine='openpyxl')
    except Exception as e:
        raise ValueError(f"The input file is not a valid Excel file: {file_path1}")
    try: #check that this is an excel file
        df_2 = pd.read_excel(file_path2,
            skiprows=1,  # Skip the first header row
            header=0,    # Use the second row as the header
            skipfooter=footer_lines2,engine='openpyxl')
    except Exception as e:
        raise ValueError(f"The input file is not a valid Excel file: {file_path2}")
    # because the dir is pulled into the module when making a griphin file some samples that aren't in the samplesheet can end up in the griphin file so we will remove those
    if samplesheet != None:
        samplesheet = pd.read_csv(samplesheet, header=0)
        samples =  samplesheet['sample'].tolist()
        #filters the DataFrame df_1 to include only the rows where the values in the 'WGS_ID' column are present in the samples list.
        df_1 = df_1[df_1['WGS_ID'].isin(samples)]
        df_2 = df_2[df_2['WGS_ID'].isin(samples)]
    # check that we have the files from the same entry
    phoenix, shiga, centar_1, centar_2, all_centar = check_column_presence(df_1, file_path1, df_2, file_path2)
    ordered_centar_df_1, centar_df_lens_1, centar_df_column_names_1 = split_centar_df(centar_1, file_path1)
    if centar_1 == True:
        #parse old griphin df
        df1_qc, df1_gene_hits_with_centar = split_dataframe(df_1,"Toxinotype")
        df1_centar, df1_gene_hits = split_dataframe(df1_gene_hits_with_centar,"AR_Database")
        df1_ar, df1_pf_hv = split_dataframe(df1_gene_hits,"HV_Database")
        df1_hv, df1_pf = split_dataframe(df1_pf_hv,"Plasmid_Replicon_Database")
    else:
        df1_qc, df1_gene_hits = split_dataframe(df_1,"AR_Database")
        df1_ar, df1_pf_hv = split_dataframe(df1_gene_hits,"HV_Database")
        df1_hv, df1_pf = split_dataframe(df1_pf_hv,"Plasmid_Replicon_Database")
    ordered_centar_df_2, centar_df_lens_2, centar_df_column_names_2 = split_centar_df(centar_2, file_path2)
    if centar_2 == True:
        #parse new griphin df
        df2_qc, df2_gene_hits_with_centar = split_dataframe(df_2,"Toxinotype")
        df2_centar, df2_gene_hits = split_dataframe(df2_gene_hits_with_centar,"AR_Database")
        df2_ar, df2_pf_hv = split_dataframe(df2_gene_hits,"HV_Database")
        df2_hv, df2_pf = split_dataframe(df2_pf_hv,"Plasmid_Replicon_Database")
    else:
        df2_qc, df2_gene_hits = split_dataframe(df_2,"AR_Database")
        df2_ar, df2_pf_hv = split_dataframe(df2_gene_hits,"HV_Database")
        df2_hv, df2_pf_tar = split_dataframe(df2_pf_hv,"Plasmid_Replicon_Database")
        df2_pf, df2_blanktar = split_dataframe(df2_pf_tar,"Toxinotype")
    #combine qc columns
    combined_df_qc = combine_qc_dataframes(df1_qc, df2_qc)
    #combine ar columns
    combined_df_ar, samples_to_add = combine_gene_dataframes(df1_ar, df2_ar)
    #make sure the order of the ar genes is correct
    order_combined_ar_df = order_ar_gene_columns(combined_df_ar)
    # combine pf and hv dataframes
    combined_df_pf, samples_to_add = combine_gene_dataframes(df1_pf, df2_pf)
    combined_df_hv, samples_to_add = combine_gene_dataframes(df1_hv, df2_hv)
    #combine centar dataframes
    ordered_centar_df, centar_df_lens, centar_df_column_names = combine_centar(ordered_centar_df_1, centar_df_lens_1, centar_df_column_names_1, ordered_centar_df_2, centar_df_lens_2, centar_df_column_names_2)
    return combined_df_qc, order_combined_ar_df, combined_df_pf, combined_df_hv, phoenix, shiga, all_centar, ordered_centar_df, centar_df_lens, centar_df_column_names

def combine_centar(ordered_centar_df_1, centar_df_lens_1, centar_df_column_names_1, ordered_centar_df_2, centar_df_lens_2, centar_df_column_names_2):
    # Combine the ordered_centar_df DataFrames
    reordered_centar_df = pd.concat([ordered_centar_df_1, ordered_centar_df_2], ignore_index=True, sort=False)
    # Combine the centar_df_lens lists
    centar_df_lens = [max(lens) for lens in zip(centar_df_lens_1, centar_df_lens_2)]
    #separate dataframes
    RB_type = [ "CEMB RT Crosswalk", "Inferred RT", "Probability", "ML Note", "Plasmid Info" ]
    RB_type_col = [col for col in reordered_centar_df.columns if any(substring in col for substring in RB_type) ]
    A_B_Tox = [ "Toxinotype", "Toxin-A_sub-type", "tcdA", "Toxin-B_sub-type", "tcdB"]
    A_B_Tox_col = [col for col in reordered_centar_df.columns if any(substring in col for substring in A_B_Tox) ]
    other_Tox = [ "tcdC", "tcdR", "tcdE", "cdtA", "cdtB", "cdtR", "cdtAB1", "cdtAB2", "non-tox", "PaLoc" ]
    other_Tox_col = [col for col in reordered_centar_df.columns if any(substring in col for substring in other_Tox) ]
    mutants = [ 'gyr','dac','feo','fur','gdp','gly','hem','hsm','Isc','mur', 'mur','nifJ','PNim','rpo','sda','thi','Van','mutations' ]
    mutations_col = [col for col in reordered_centar_df.columns if any(substring in col for substring in mutants) ]
    # List of mutation names to remove
    mutations_to_remove = ['tcdC other mutations', 'cdtR other mutations', 'PaLoc_NonTox other mutations']
    # Remove each mutation name if it exists in mutations_col
    mutations_col = [mutation for mutation in mutations_col if mutation not in mutations_to_remove]
    existing_columns_in_order = A_B_Tox_col + other_Tox_col + mutations_col + RB_type_col
    ordered_centar_df = reordered_centar_df[existing_columns_in_order]
    # Combine the centar_df_column_names lists
    #centar_df_column_names = [list(set(chain(*columns))) for columns in zip(centar_df_column_names_1, centar_df_column_names_2)]
    return ordered_centar_df, centar_df_lens, existing_columns_in_order


def split_centar_df(centar, excel):
    footer_lines1 = detect_footer_lines(excel)
    df = pd.read_excel(excel,
        header=[0, 1],    # Use the second row as the header
        skipfooter=footer_lines1,engine='openpyxl')
    if centar == True:
        df_1 = pd.read_excel(excel,
            skiprows=1,  # Skip the first header row
            header=0,    # Use the second row as the header
            skipfooter=footer_lines1,engine='openpyxl')
        # Get the column names
        columns = df_1.columns.tolist()
        # Find the indices of the specified columns
        ## Might need to fiddle with MLST Clade here eventually...might not if original griphins get made right
        start_index = columns.index('Toxinotype')
        end_index = columns.index('AR_Database')
        # Get the columns between the specified columns
        columns_between = [columns[0]]+columns[start_index:end_index]
        # Subset the DataFrame
        ordered_centar_df = df_1[columns_between]
        #get length of sub dataframes
        centar_df_lens = [ len(df['Toxin A/B Variants'].columns), len(df['Other Toxins'].columns), len(df['C. difficile Specific AR Mutations'].columns), len(df['ML Predicted Ribotype'].columns) ] # have to have it in a list for sum() later
        centar_df_column_names = [ df['Toxin A/B Variants'].columns, df['Other Toxins'].columns, df['C. difficile Specific AR Mutations'].columns, df['ML Predicted Ribotype'].columns ] # have to have it in a list for sum() later
    else:
        centar_df_lens = [0,0,0,0]
        centar_df_column_names = [[],[],[],[]]
        ordered_centar_df = pd.DataFrame()
    return ordered_centar_df, centar_df_lens, centar_df_column_names

def combine_gene_dataframes(old_df, new_df):
    # Ensure the first column is 'WGS_ID'
    if new_df.columns[0] != 'WGS_ID' or old_df.columns[0] != 'WGS_ID':
        raise ValueError("The first column in both dataframes must be 'WGS_ID'")
    df1_gene = old_df.set_index('WGS_ID')
    df2_gene = new_df.set_index('WGS_ID')
    # Identify samples to add and print them
    samples_to_add = df2_gene.index.difference(df1_gene.index)
    # Combine dataframes, prioritizing new_df values and aligning columns
    combined_df = df1_gene.combine_first(df2_gene)
    # Add missing columns from new_df to old_df
    columns_to_add = [col for col in df2_gene.columns if col not in df1_gene.columns and col not in ["AR_Database", "HV_Database","Plasmid_Replicon_Database", "WGS_ID", "No_AR_Genes_Found", "No_HVGs_Found", "No_Plasmid_Markers"]]
    for col in columns_to_add:
        combined_df[col] = combined_df[col].fillna(df2_gene[col])
    # Revert index to column
    combined_df.reset_index(inplace=True)
    return combined_df, samples_to_add

def combine_qc_dataframes(df1_qc, df2_qc):
    # Convert the 'WGS_ID' columns to the index to facilitate updating
    df1 = df1_qc.set_index('WGS_ID')
    df2 = df2_qc.set_index('WGS_ID')
    # Update df1_qc with matching rows from df2_qc
    df1.update(df2)
    # Append non-matching rows from df2_qc to df1_qc using pd.concat
    combined_df = pd.concat([df1, df2[~df2.index.isin(df1.index)]])
    # Reset the index to restore the 'WGS_ID' column
    combined_df.reset_index(inplace=True)
    return combined_df

def add_blank_centar_columns(df, reference_df):
    # Get the list of column names from the reference dataframe
    columns = reference_df.columns.tolist()
    # Find the indices of the "Toxinotype" and "AR_Database" columns
    try:
        start_index = columns.index("Toxinotype")
        end_index = columns.index("AR_Database")
    except ValueError:
        raise ValueError("One of the required columns ('Toxinotype' or 'AR_Database') is missing in the reference dataframe.")

    # Get all the columns between "Toxinotype" (inclusive) and "AR_Database" (exclusive)
    centar_columns = columns[start_index:end_index]
    # Add blank columns for each centar-related column
    for col in centar_columns:
        if col not in df.columns:
            df[col] = pd.NA  # Fill with NaNs or blank
    return df

def check_column_presence(df1, df1_path, df2, df2_path):
    df1_has_column = "BUSCO_Lineage" in df1.columns
    df2_has_column = "BUSCO_Lineage" in df2.columns
    if df1_has_column and df2_has_column: # set that its a CDC PHOENIX run
        phoenix = False
    elif not df1_has_column and not df2_has_column: # set that its a PHOENIX run
        phoenix = True
    else:
        if df1_has_column:
            raise ValueError(f"{CRED}The old griphin file was produced from -entry CDC_PHOENIX and the new griphin summary wasn't. These files aren't compatible.{CEND}")
        else:
            raise ValueError(f"{CRED}The new griphin file was produced from -entry CDC_PHOENIX and the old griphin summary wasn't. These files aren't compatible.{CEND}")
    #check for centar
    if "Toxin-A_sub-type" in df1.columns or "Toxin-A_sub-type" in df2.columns:
        if "Toxin-A_sub-type" not in df1.columns:
            print("Adding centar blank columns to " + df1_path)
            df1 = add_blank_centar_columns(df1, df2)
            centar_1 = False
        else:
            centar_1 = True
        if "Toxin-A_sub-type" not in df2.columns:
            print("Adding centar blank columns to " + df2_path)
            df2 = add_blank_centar_columns(df2, df1)
            centar_2 = False
        else:
            centar_2 = True
    else:
        centar_1 = centar_2 = False
    if centar_1 == True or centar_2 == True:
        all_centar = True
    else:
        all_centar = False
    #check for shigapass
    if "ShigaPass_Organism" in df1.columns or "ShigaPass_Organism" in df2.columns:
        shiga = True
    else:
        shiga = False
    return phoenix, shiga, centar_1, centar_2, all_centar

def split_dataframe(df, split_column):
    # Ensure the first column is 'WGS_ID'
    if df.columns[0] != 'WGS_ID':
        raise ValueError("The first column must be 'WGS_ID'")
    # Find the index of the split_column
    split_index = df.columns.get_loc(split_column)
    # Create two DataFrames based on the split index, including 'WGS_ID' at the beginning
    df_before = df.iloc[:, [0] + list(range(1, split_index))]  # All columns before the split column + 'WGS_ID'
    #df_before = df.iloc[:, list(range(1, split_index))]  # All columns before the split column + 'WGS_ID'
    df_after = df.iloc[:, [0] + list(range(split_index, df.shape[1]))]  # The split column and all columns after it + 'WGS_ID'
    return df_before, df_after

def add_and_concatenate(df1, df2):
    # Step 1: Preserve column order from df1 and append missing columns from df2
    all_columns = list(df1.columns) + [col for col in df2.columns if col not in df1.columns]
    # Step 1: Find the union of columns in both DataFrames
    #all_columns = set(df1.columns).union(set(df2.columns))
    # Step 2: Add missing columns to each DataFrame at once by reindexing
    df1 = df1.reindex(columns=all_columns)
    df2 = df2.reindex(columns=all_columns)
    # Step 3: Concatenate the two DataFrames row-wise
    result = pd.concat([df1, df2], ignore_index=True, sort=False)
    # Step 4: Reset the index and return the final DataFrame
    result.reset_index(drop=True, inplace=True)
    return result

def detect_footer_lines(file_path, sheet_name=0):
    # Load the workbook
    wb = load_workbook(file_path, read_only=True, data_only=True)
    ws = wb[sheet_name if isinstance(sheet_name, str) else wb.sheetnames[sheet_name]]
    footer_lines = 3 #set to 3 because review by row has data in another cell and there are two all blank spacing rows
    is_footer_section = False
    # Count footer lines
    for row in reversed(tuple(ws.rows)):
        first_cell = row[0].value
        other_cells = [cell.value for cell in row[1:]]
        # Check if the first cell has text and all other cells are empty
        if isinstance(first_cell, str) and all(cell is None for cell in other_cells):
            footer_lines += 1
            is_footer_section = True
        elif is_footer_section:
            break
    return footer_lines

def read_excel(file_path, old_phoenix, reference_qc_df, reference_centar_df):
    #get number of footer lines for each file to skip when reading in as a pd
    footer_lines = detect_footer_lines(file_path)
    # Read the Excel file, skipping the first row and using the second row as the header
    try: #check that this is an excel file
        df = pd.read_excel(file_path,
            skiprows=1,  # Skip the first header row
            header=0,    # Use the second row as the header
            skipfooter=footer_lines,engine='openpyxl')

    except Exception as e:
        raise ValueError(f"The input file is not a valid Excel file: {file_path}")
    # check that we have the files from the same entry
    if "BUSCO_Lineage" in df.columns: # set that its a CDC PHOENIX run
        phoenix = False
    else:
        phoenix = False
    if old_phoenix != phoenix: 
        raise ValueError(f"{CRED}The one griphin file in your set was produced from -entry CDC_PHOENIX and some were not wasn't. These files aren't compatible to combine. {CEND}")
    #check for centar
    if "Toxin-A_sub-type" not in df.columns:
        print("Adding centar blank columns to " + file_path)
        for col in reference_centar_df.columns:
            if col not in df.columns:
                df[col] = pd.NA  # Fill with NaNs or blank
        centar = all_centar = False
    else:
        centar = all_centar = True
    #check for shigapass
    if "ShigaPass_Organism" in df.columns:
        shiga = True
    else:
        shiga = False
    #parse old griphin df
    ordered_centar_df, centar_df_lens, centar_df_column_names = split_centar_df(centar, file_path)
    #centar_df_lens[0] += 1
    #centar_df_lens = [x + 1 for x in centar_df_lens]
    if centar == True:
        df_qc, df_gene_hits_with_centar = split_dataframe(df,"Toxinotype")
        df_centar, df_gene_hits = split_dataframe(df_gene_hits_with_centar,"AR_Database")
        df_ar, df_pf_hv = split_dataframe(df_gene_hits,"HV_Database")
        df_hv, df_pf = split_dataframe(df_pf_hv,"Plasmid_Replicon_Database")
    else:
        df_qc, df_gene_hits = split_dataframe(df,"AR_Database")
        df_ar, df_pf_hv = split_dataframe(df_gene_hits,"HV_Database")
        df_hv, df_pf = split_dataframe(df_pf_hv,"Plasmid_Replicon_Database")
        # Ensure the first column is 'WGS_ID'
    if df_qc.columns[0] != 'WGS_ID' or reference_qc_df.columns[0] != 'WGS_ID':
        raise ValueError("The first column in both dataframes must be 'WGS_ID'")
    # Set WGS_ID as index
    df1_qc = reference_qc_df.set_index('WGS_ID')
    ref_qc = df_qc.set_index('WGS_ID')
    # Identify samples to add and print them
    samples_to_add = ref_qc.index.difference(df1_qc.index)
    #make sure the order of the ar genes is correct
    order_ar_df = order_ar_gene_columns(df_ar)
    print("Adding sample(s) to the GRiPHin summary:", samples_to_add.tolist())
    return df_qc, order_ar_df, df_pf, df_hv, phoenix, shiga, centar, ordered_centar_df, centar_df_lens, centar_df_column_names

def update_centar_columns(centar_df_column_names_final, centar_df_column_names):
    # Iterate over each corresponding pair of lists in centar_df_lens_final and centar_df_lens
    for i, (final_list, new_list) in enumerate(zip(centar_df_column_names_final, centar_df_column_names)):
        # Convert both lists to sets to get unique column names
        unique_columns = set(final_list).union(set(new_list))
        # Update centar_df_lens_final with the unique columns
        centar_df_column_names_final[i] = list(unique_columns)
    # After processing all lists, calculate the total number of unique names in each list
    total_unique_columns = [len(columns) for columns in centar_df_column_names_final]
    return found_headlines, total_unique_columns

def main():
    # looping through excel files
    args = parseArgs()
    #figure out the name of the output file 
    if args.output != None:
        output_file = args.output
    else:
        # Derive output file name from input file name
        output_file = args.griphin_new.replace("_GRiPHin_Summary.xlsx", "")
    # checking what the input type is
    if args.griphin_list != False:
        griphin_files = glob.glob("*_GRiPHin_Summary.xlsx")
        if len(griphin_files) < 2:
            raise ValueError(f"{CRED}Need at least two GRiPHin files for combination when using --griphin_list.{CEND}")
        # Initialize the combination with the first file
        base_file = griphin_files.pop(0)
        #combine first two files
        combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, phoenix_final, shiga_final, centar_final, ordered_centar_df, centar_df_lens_final, centar_df_column_names_final = read_excels(base_file, griphin_files[0], args.samplesheet)
        combined_dataframes_1 = [ combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, ordered_centar_df ]
        # Iterate over remaining files, progressively combining them with the base
        for count, next_file in enumerate(griphin_files[1:], start=1):
            # checking
            combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, phoenix, shiga, centar, ordered_centar_df, centar_df_lens_final, centar_df_column_names = read_excel(next_file, phoenix_final, combined_df_qc, ordered_centar_df)
            combined_dataframes = [ combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, ordered_centar_df ]
            # Update flags
            phoenix_final = phoenix_final or phoenix
            shiga_final = shiga_final or shiga
            centar_final = centar_final or centar
            if count == 1:
                # Iterate over both tuples and combine corresponding DataFrames
                combined_dataframes_final = [ add_and_concatenate(df1, df2) for df1, df2 in zip(combined_dataframes_1, combined_dataframes) ]
                if centar_final == True:
                    # Update centar_df_lens1 with the maximum of each corresponding number
                    #centar_df_lens_final = [max(val1, val2) for val1, val2 in zip(centar_df_lens_final, centar_df_lens)]
                    centar_df_column_names_final, centar_df_lens_final = update_centar_columns(centar_df_column_names_final, centar_df_column_names)
            else:
                combined_dataframes_final = [ add_and_concatenate(df1, df2) for df1, df2 in zip(combined_dataframes_final, combined_dataframes) ]
                if centar_final == True:
                    # Update centar_df_lens1 with the maximum of each corresponding number
                    #centar_df_lens_final = [max(val1, val2) for val1, val2 in zip(centar_df_lens_final, centar_df_lens)]
                    centar_df_column_names_final, centar_df_lens_final = update_centar_columns(centar_df_column_names_final, centar_df_column_names)
            # Unpack the tuple into individual DataFrames
            combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, ordered_centar_df = combined_dataframes_final
    else:
        combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, phoenix_final, shiga_final, centar_final, ordered_centar_df, centar_df_lens_final, centar_df_column_names = read_excels(args.griphin_new, args.griphin_old, args.samplesheet)
    # add centar_df to the almost_final_ar_df
    if centar_final == True:
        # Reset the indices of both DataFrames during concat so they are aligned and we don't get NaNs in the final row of the dataframe
        combined_df_qc = pd.concat([combined_df_qc, ordered_centar_df], axis=1)
    # call function from griphin script to combine all dfs
    #combined_df_ar["AR_Database"] = combined_df_ar["AR_Database"].fillna("No AR Genes Found")
    final_df, ar_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db = Combine_dfs(combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, pd.DataFrame(), True)
    #get other information for excel writing
    (qc_max_row, qc_max_col) = combined_df_qc.shape
    # Mst account for WGS_ID removal later, so decreasing size by one here before i forget
    qc_max_row -= 1
    pf_max_col = combined_df_pf.shape[1] - 1 #remove one for the WGS_ID column
    hv_max_col = combined_df_hv.shape[1] - 1 #remove one for the WGS_ID column
    #write excel sheet
    write_to_excel(args.set_coverage, output_file, final_df, qc_max_col, ar_max_col, pf_max_col, hv_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db, phoenix_final, shiga_final, centar_final, centar_df_lens_final)
    #write tsv from excel
    convert_excel_to_tsv(output_file)

if __name__ == '__main__':
    main()