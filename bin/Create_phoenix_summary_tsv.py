#!/usr/bin/env python3

import sys
import glob
import os
import pandas as pd
import fnmatch
from decimal import *
getcontext().prec = 4
import argparse
import yaml

##Makes a summary Excel file when given a series of output summary line files from PHoeNIx
##Usage: >python Phoenix_Summary_tsv_06-10-22.py -o Output_Report.tsv 
## Written by Rich Stanton (njr5@cdc.gov), updates by Jill Hagey (qpk9@cdc.gov)
# Set display options to show all rows and columns
pd.set_option('display.max_rows', None)  # Show all rows
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.max_colwidth', None)  # Show all columns

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CYELLOW = '\033[93m'
CEND = '\033[0m'

# Function to get the script version
def get_version():
    return "2.2.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet')
    parser.add_argument('-o', '--out', dest='output_file', required=True, help='output file name')
    parser.add_argument('--updater', dest='updater', action='store_true', help='Pass if running updater pipeline.')
    parser.add_argument('-b', '--busco', action='store_true', help='parameter to know if busco was run')
    parser.add_argument('--software_versions', dest='software_versions', help='This will update old files with correct PHX_Version')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    parser.add_argument('files', nargs=argparse.REMAINDER)
    return parser.parse_args()

# Create new column based on Taxa_Source
def get_organism_value(row):
    if row['Taxa_Source'] == 'ANI_REFSEQ':
        return row['FastANI_Organism']
    elif row['Taxa_Source'] == 'Shigapass':
        return row['ShigaPass_Organism']
    else:
        return ''

# Create new column based on Taxa_Source
def get_fastani_org(row):
    print(f"Taxa_Source: {row['Taxa_Source']}")
    if row['Taxa_Source'] == 'ANI_REFSEQ':
        print(f"Final_Taxa_ID: {row['Final_Taxa_ID']}")
        return row['Final_Taxa_ID']
    else:
        return ''

def List_TSV(output_file, input_list, busco, updater, phx_version):
    with open(output_file, 'w') as f:
        if not busco: # if --busco is not passed
            f.write('WGS_ID\tPHX_Version\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tFinal_Taxa_ID\tTaxa_Source\tFastANI_Organism\tFastANI_%ID\tFastANI_%Coverage\tShigaPass_Organism\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason\n')
        if busco:
            f.write('WGS_ID\tPHX_Version\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tBUSCO\tBUSCO_DB\tFinal_Taxa_ID\tTaxa_Source\tFastANI_Organism\tFastANI_%ID\tFastANI_%Coverage\tShigaPass_Organism\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason\n')
        try: #if there are numbers in the name then use that to sort
            print(f"Sorting input files by numeric values in their names")
            if updater == True:
                def extract_sample_id(filepath):
                    filename = os.path.basename(filepath)
                    if filename.endswith('_summaryline.tsv'):
                        return filename.replace("_summaryline.tsv","")  # Remove '_summaryline.tsv' suffix
                    return filename
                input_list_sorted=sorted(input_list, key=extract_sample_id)
            else:
                input_list_sorted=sorted(input_list, key=lambda x: int("".join([i for i in x if i.isdigit()])))
        except: #if no numbers then use only alphabetically
            input_list_sorted=sorted(input_list)
        for file in input_list_sorted:
            print(CRED + file + CEND)
            with open(file, "r") as f2:
                header = next(f2) # skip the first line of the samplesheet
                if "Final_Taxa_ID" not in header: # FOR BACKWARDS COMPATIBILITY
                    #print(f"Older file format detected in {file}, updating")
                    # If you have existing data
                    df = pd.read_csv(file, sep="\t")
                    df["Taxa_Confidence"] = df["Taxa_Confidence"].astype(str).str.split(' ').str[0]
                    # Rename columns
                    df = df.rename(columns={ "ID": "WGS_ID", "Species": "FastANI_Organism", "Taxa_Confidence": "FastANI_%ID", "Taxa_Coverage": "FastANI_%Coverage" })
                    #move Taxa_Source column
                    taxa_source_col = df.pop("Taxa_Source")
                    if busco and not "BUSCO" in df.columns.tolist():
                        add_index = df.columns.get_loc("GC_%") + 6
                    else:
                        add_index = df.columns.get_loc("GC_%") + 4
                    df.insert(add_index, "Taxa_Source", taxa_source_col)
                    # Insert PHX_Version as second column (index 1)
                    df.insert(1, "PHX_Version", phx_version)
                    # Move FastANI_%ID column to after Taxa_Source
                    fastani_id_col = df.pop("FastANI_%ID")
                    #taxa_source_idx = df.columns.get_loc("Taxa_Source")
                    df.insert(df.columns.get_loc("Taxa_Source") + 2, "FastANI_%ID", fastani_id_col)
                    # Insert ShigaPass_Organism column
                    df.insert(df.columns.get_loc("Taxa_Source") + 3, "ShigaPass_Organism", "")
                    # Insert Final_Taxa_ID column
                    df.insert(df.columns.get_loc("Taxa_Source") + 1, "Final_Taxa_ID", "")
                    # fill in Final_Taxa_ID based on Taxa_Source
                    print(df['FastANI_Organism'])
                    df['Final_Taxa_ID'] = df.apply(get_organism_value, axis=1)
                    print(f"Final_Taxa_ID: {df['Final_Taxa_ID']}")
                    # fill in FastANI_Organism based on Taxa_Source
                    df['FastANI_Organism'] = df.apply(get_fastani_org, axis=1)
                    headers = df.columns.tolist()
                    if busco and not "BUSCO" in headers:
                        headers.insert(df.columns.get_loc("GC_%") + 1, "BUSCO")
                        headers.insert(df.columns.get_loc("GC_%") + 2, "BUSCO_DB")
                    # Write the updated data
                    for _, row in df.iterrows():
                        f.write('\t'.join(str(row[col]) if col in row and pd.notna(row[col]) else '' for col in headers) + '\n')
                    continue
                else:
                    df = pd.read_csv(file, sep="\t")
                for line in f2:
                    f.write(line.strip('\n') + '\n')

def collect_files(updater):
    summary_files = glob.glob('*.tsv')
    files_to_remove = ['empty_summaryline.tsv', '*_centar_output.tsv', '*Phoenix_Summary.tsv']
    summary_files = [file for file in summary_files
                        if not any(fnmatch.fnmatch(file, pattern) for pattern in files_to_remove)]
    if updater is True:
        old_summary_files = glob.glob('*/*/*_summaryline.tsv')
        # Filter out items containing any of the substrings
        filtered_list = [file for file in old_summary_files if not any(substring in file for substring in summary_files)]
        summary_files = filtered_list + summary_files
    return summary_files

# Function to get the script version from YAML file
def get_old_version(yaml_file_path):
    """Extract version information from YAML file for cdcgov/phoenix."""
    try:
        with open(yaml_file_path, 'r') as file:
            data = yaml.safe_load(file)
            
        # Navigate through the YAML structure to find cdcgov/phoenix version
        if 'Workflow' in data and 'cdcgov/phoenix' in data['Workflow']:
            return data['Workflow']['cdcgov/phoenix']
        
        # Alternative: search through all keys in case structure varies
        for key, value in data.items():
            if isinstance(value, dict) and 'cdcgov/phoenix' in value:
                return value['cdcgov/phoenix']
    except FileNotFoundError:
        print(f"Warning: YAML file {yaml_file_path} not found, using default version")
    except yaml.YAMLError as e:
        print(f"Warning: Error parsing YAML file: {e}, using default version")
    except Exception as e:
        print(f"Warning: Unexpected error reading version: {e}, using default version")
    # Return default version if file not found or parsing fails
    return ""

def main():
    args = parseArgs()
    summary_files = collect_files(args.updater)
    #if args.updater:
    #    phx_version = get_old_version(args.software_versions)
    List_TSV(args.output_file, summary_files, args.busco, args.updater, args.software_versions)

if __name__ == '__main__':
    main()