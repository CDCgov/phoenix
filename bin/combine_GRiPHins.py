#!/usr/bin/env python3

#disable cache usage in the Python so __pycache__ isn't formed. If you don't do this using 'nextflow run cdcgov/phoenix...' a second time will causes and error
import sys
sys.dont_write_bytecode = True
import pandas as pd
import argparse
import xlsxwriter as ws
from xlsxwriter.utility import xl_rowcol_to_cell
from openpyxl import load_workbook
from itertools import chain
from GRiPHin import order_ar_gene_columns, Combine_dfs, big5_check, write_to_excel, convert_excel_to_tsv

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
    parser.add_argument('-o', '--output', required=False, dest='output', help='Name of output file default is GRiPHin_Summary.xlsx.')
    parser.add_argument('--coverage', default=30, required=False, dest='set_coverage', help='The coverage cut off default is 30x.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def read_excels(file_path1, file_path2):
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
    # check that we have the files from the same entry
    phoenix = check_column_presence(df_1, df_2)
    #parse old griphin df
    df1_qc, df1_gene_hits = split_dataframe(df_1,"AR_Database")
    df1_ar, df1_pf_hv = split_dataframe(df1_gene_hits,"HV_Database")
    print(df1_ar)
    df1_hv, df1_pf = split_dataframe(df1_pf_hv,"Plasmid_Replicon_Database")
    #parse new griphin df
    df2_qc, df2_gene_hits = split_dataframe(df_2,"AR_Database")
    df2_ar, df2_pf_hv = split_dataframe(df2_gene_hits,"HV_Database")
    df2_hv, df2_pf = split_dataframe(df2_pf_hv,"Plasmid_Replicon_Database")
    #combine qc columns
    combined_df_qc = combine_qc_dataframes(df1_qc, df2_qc)
    #combine ar columns
    combined_df_ar, samples_to_add = combine_gene_dataframes(df1_ar, df2_ar)
    #make sure the order of the ar genes is correct
    order_combined_ar_df = order_ar_gene_columns(combined_df_ar)
    print("Adding sample(s) to the GRiPHin summary:", samples_to_add.tolist())
    # combine pf and hv dataframes
    combined_df_pf, samples_to_add = combine_gene_dataframes(df1_pf, df2_pf)
    combined_df_hv, samples_to_add = combine_gene_dataframes(df1_hv, df2_hv)
    return combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, phoenix

def combine_gene_dataframes(old_df, new_df):
    # Ensure the first column is 'WGS_ID'
    if new_df.columns[0] != 'WGS_ID' or old_df.columns[0] != 'WGS_ID':
        raise ValueError("The first column in both dataframes must be 'WGS_ID'")
    # Set WGS_ID as index
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

def check_column_presence(df1, df2):
    df1_has_column = "BUSCO_Lineage" in df1.columns
    df2_has_column = "BUSCO_Lineage" in df2.columns
    if df1_has_column and df2_has_column: # set that its a CDC PHOENIX run
        phoenix = False
        return phoenix
    elif not df1_has_column and not df2_has_column: # set that its a PHOENIX run
        phoenix = True
        return phoenix
    else:
        if df1_has_column:
            raise ValueError(f"{CRED}The old griphin file was produced from -entry CDC_PHOENIX and the new griphin summary wasn't. These files aren't compatible.{CEND}")
        else:
            raise ValueError(f"{CRED}The new griphin file was produced from -entry CDC_PHOENIX and the old griphin summary wasn't. These files aren't compatible.{CEND}")

def split_dataframe(df, split_column):
    # Ensure the first column is 'WGS_ID'
    if df.columns[0] != 'WGS_ID':
        raise ValueError("The first column must be 'WGS_ID'")
    # Find the index of the split_column
    split_index = df.columns.get_loc(split_column)
    # Create two DataFrames based on the split index, including 'WGS_ID' at the beginning
    df_before = df.iloc[:, [0] + list(range(1, split_index))]  # All columns before the split column + 'WGS_ID'
    df_after = df.iloc[:, [0] + list(range(split_index, df.shape[1]))]  # The split column and all columns after it + 'WGS_ID'
    return df_before, df_after

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


def main():
    args = parseArgs()
    #figure out what we shuld name the output file 
    if args.output:
        output_file = args.output
    else:
        # Derive output file name from input file name
        output_file = args.griphin_new.replace("_GRiPHin_Summary.xlsx", "")
    combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, phoenix = read_excels(args.griphin_new, args.griphin_old)
    # call function from griphin scrript to combine all dfs
    #print(combined_df_ar)
    final_df, ar_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db = Combine_dfs(combined_df_qc, combined_df_ar, combined_df_pf, combined_df_hv, pd.DataFrame(), True)
    #get other information for excel writing
    (qc_max_row, qc_max_col) = combined_df_qc.shape
    pf_max_col = combined_df_pf.shape[1] - 1 #remove one for the WGS_ID column
    hv_max_col = combined_df_hv.shape[1] - 1 #remove one for the WGS_ID column
    #write excel sheet
    write_to_excel(args.set_coverage, output_file, final_df, qc_max_col, ar_max_col, pf_max_col, hv_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db, phoenix)
    #write tsv from excel
    convert_excel_to_tsv(output_file)

if __name__ == '__main__':
    main()