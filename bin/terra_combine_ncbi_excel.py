#!/usr/bin/env python3

# importing the required modules
import glob
import pandas as pd
import argparse
import xlsxwriter as ws
from xlsxwriter.utility import xl_rowcol_to_cell
import re
from re import search
from itertools import chain

##Makes a summary Excel file when given a series of ncbi excel files produced by PHoeNIx using --create_ncbi_sheet
##Usage: >python terra_combine_ncbi_excel.py --sra_output <sra_excel_prefix> --biosample_output <biosample_excel_prefix>
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a combined NCBI summary excel sheet')
    parser.add_argument('--biosample_output', dest='biosample_output', required=False, default="BiosampleAttributes_Microbe.1.0.xlsx", help='prefix for biosample final file.')
    parser.add_argument('--sra_output', dest='sra_output', required=False, default="Sra_Microbe.1.0.xlsx", help='prefix for sra final file.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    parser.add_argument('files', nargs=argparse.REMAINDER)
    return parser.parse_args()
 
def combine_second_rows(excel_files, final_file):
    # Read the header from the first Excel file
    header_df = pd.read_excel(excel_files[0], header=0, nrows=0)
    file_type = excel_files[0].split('_')[0]

    # Initialize an empty DataFrame to store combined second rows
    combined_second_rows = pd.DataFrame()

    # Identify invalid files where all columns in the second row contain NaN values
    invalid_files = [file for file in excel_files if pd.read_excel(file, header=0, nrows=1).isna().all().all()]
    #print results of screening for na columns
    if invalid_files:
        print("{} files with all NaN columns:".format(file_type), invalid_files)
    else:
        print("No {} files found with all NaN columns.".format(file_type))

    # remove invalid files from list
    valid_excel_files = [file for file in excel_files if file not in invalid_files]

    #combine each row to a combined dataframe
    for file in valid_excel_files:
        df = pd.read_excel(file, header=0, nrows=1)
        combined_second_rows = pd.concat([df, combined_second_rows], ignore_index=True, axis=0)
        combined_second_rows = combined_second_rows[list(header_df.columns)]

    if len(valid_excel_files) < 0:  # Check if there are any valid second rows
        print("No valid second rows found in any {} file.".format(file_type))

    # Write the combined second rows along with the header to a new Excel file
    combined_second_rows.to_excel(final_file, index=False,header=list(header_df.columns))

def add_disclaimer(input_excel, input_sheet_name):
    df = pd.read_excel(input_excel, header=0)
    # Load the Excel file using pandas ExcelWriter
    with pd.ExcelWriter(input_excel, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name=input_sheet_name, index=False)
        workbook  = writer.book
        worksheet = writer.sheets[input_sheet_name]
        # getting values to set column widths automatically
        for idx, col in enumerate(df):  # loop through all columns
            series = df[col]
            max_len = max((
            series.astype(str).map(len).max(),  # len of largest item
                len(str(series.name))  # len of column name/header
                )) + 1  # adding a little extra space
            worksheet.set_column(idx, idx, max_len)  # set column width
        # Add text to the cell
        # Change the font color to red
        red_format = workbook.add_format({'color': 'red', 'bold': True, 'text_wrap': True})
        # Change the font color to orange
        orange_format = workbook.add_format({'color': 'orange', 'bold': True, 'text_wrap': True})
        biosample_delete_warning = """Do the following before upload:
1. Delete this row and the rows below!
2. At minimum fill out the following columns: 
    - Host: e.g., Homo sapiens, animal, environmental, other
    - Collection Date: Specimen collection year only
    - Isolate Source: Specimen source. 
        * Homo sapiens or animal: blood, urine, etc
        * environmental: device or surface """
        sra_delete_warning = """Do the following before upload:
1. Delete this row and the rows below!
2. Fill out 'design_description' column with a short description of our library prep info and any other pertinent information. Ex: Sequenced using Nextera XT library prep kit, 2 x 250.
3. Fill out 'instrument_model' column with your illumina model type and number (if it has one). Ex: Illumina HiSeq 1500."""
        disclaimer_text = """As a reminder, please do not submit raw sequencing data to the CDC HAI-Seq BioProject (531911) that is auto populated in this sheet unless you are a state public health laboratory, a CDC partner or have been directed to do so by DHQP. The BioProject accession IDs in this file are specifically designated for domestic HAI bacterial pathogen sequencing data, \
including from the Antimicrobial Resistance Laboratory Network (AR Lab Network), state public health labs, surveillance programs, and outbreaks. For inquiries about the appropriate BioProject location for your data, please contact HAISeq@cdc.gov."""
        # Determine the number of rows already filled
        num_rows = df.shape[0] + 2
        # Add text to the cell
        if input_sheet_name == "SRA_data":
            worksheet.merge_range('A' + str(num_rows+1) + ':J' + str(num_rows+4), sra_delete_warning, orange_format)
            worksheet.merge_range('A' + str(num_rows+5) + ':J' + str(num_rows+8), disclaimer_text, red_format)
        else:
            worksheet.merge_range('A' + str(num_rows+1) + ':J' + str(num_rows+6), biosample_delete_warning, orange_format)
            worksheet.merge_range('A' + str(num_rows+7) + ':J' + str(num_rows+9), disclaimer_text, red_format)

def main():
    args = parseArgs()
    # get files in the path
    biosample_list = glob.glob("BiosampleAttributes_*_Microbe.1.0.xlsx")
    sra_list = glob.glob("Sra_*_Microbe.1.0.xlsx")
    # check that the file_list isn't empty
    if len(biosample_list) != 0:
        print(str(len(biosample_list)) + " BiosampleAttributes_Microbe.1.0.xlsx files were found using BiosampleAttributes_*_Microbe.1.0.xlsx!")
        if args.biosample_output != "BiosampleAttributes_Microbe.1.0.xlsx":
            final_file = args.biosample_output + "_BiosampleAttributes_Microbe.1.0.xlsx"
        else:
            final_file = args.biosample_output
        combine_second_rows(biosample_list, final_file)
        add_disclaimer(final_file, "Sheet1")
    else:
        print("No BiosampleAttributes_Microbe.1.0.xlsx files were found, you need to pass at least one!")
        exit()
    if len(sra_list) != 0:
        print(str(len(sra_list)) + " Sra_Microbe.1.0.xlsx files were found using Sra_*_Microbe.1.0.xlsx!")
        if args.sra_output != "Sra_Microbe.1.0.xlsx":
            final_file = args.sra_output + "_Sra_Microbe.1.0.xlsx"
        else:
            final_file = args.sra_output
        combine_second_rows(sra_list, final_file)
        add_disclaimer(final_file, "SRA_data")
    else:
        print("No Sra_Microbe.1.0.xlsx files were found, you need to pass at least one!")
        exit()


if __name__ == '__main__':
    main()