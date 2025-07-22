#!/usr/bin/env python3

import sys
import glob
import os
import pandas as pd
import fnmatch
from decimal import *
getcontext().prec = 4
import argparse

##Makes a summary Excel file when given a series of output summary line files from PHoeNIx
##Usage: >python Phoenix_Summary_tsv_06-10-22.py -o Output_Report.tsv 
## Written by Rich Stanton (njr5@cdc.gov), updates by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "2.2.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet')
    parser.add_argument('-o', '--out', dest='output_file', required=True, help='output file name')
    parser.add_argument('--updater', dest='updater', action='store_true', help='Pass if running updater pipeline.')
    parser.add_argument('-b', '--busco', action='store_true', help='parameter to know if busco was run')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    parser.add_argument('files', nargs=argparse.REMAINDER)
    return parser.parse_args()

def List_TSV(output_file, input_list, busco):
    with open(output_file, 'w') as f:
        if not busco: # if --busco is not passed
            f.write('WGS_ID\tPHX_Version\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tFinal_Taxa_ID\tTaxa_Source\tFastANI_Organism\tFastANI_%ID\tFastANI_%Coverage\tShigaPass_Organism\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason\n')
        if busco:
            f.write('WGS_ID\tPHX_Version\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tBUSCO\tBUSCO_DB\tFinal_Taxa_ID\tTaxa_Source\tFastANI_Organism\tFastANI_%ID\tFastANI_%Coverage\tShigaPass_Organism\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason\n')
        try: #if there are numbers in the name then use that to sort
            input_list_sorted=sorted(input_list, key=lambda x: int("".join([i for i in x if i.isdigit()])))
        except: #if no numbers then use only alphabetically
            input_list_sorted=sorted(input_list)
        for entry in input_list_sorted:
            with open(entry, "r") as f2:
                header = next(f2) # skip the first line of the samplesheet
                if "Final_Taxa_ID" not in header: # FOR BACKWARDS COMPATIBILITY
                    print(f"Older file format detected in {entry}, updating")
                    # If you have existing data
                    df = pd.read_csv(entry, sep="\t")
                    # Rename columns
                    df = df.rename(columns={ "ID": "WGS_ID", "Species": "FastANI_Organism", "Taxa_Confidence": "FastANI_%ID", "Taxa_Coverage": "FastANI_%Coverage" })
                    # Insert PHX_Version as second column (index 1)
                    df.insert(1, "PHX_Version", "")
                    headers = df.columns.tolist()
                    taxa_source_idx = headers.index("Taxa_Source")
                    headers.insert(taxa_source_idx + 4, "ShigaPass_Organism")
                    fastani_org_idx = headers.index("FastANI_Organism")
                    headers.insert(fastani_org_idx - 1, "Final_Taxa_ID")
                    if busco and not "BUSCO" in headers:
                        headers.insert(fastani_org_idx + 1, "BUSCO")
                        headers.insert(fastani_org_idx + 2, "BUSCO_DB")
                    # Write the updated data
                    for _, row in df.iterrows():
                        f.write('\t'.join(str(row[col]) if col in row and pd.notna(row[col]) else '' for col in headers) + '\n')
                    continue
                for line in f2:
                    f.write(line.strip('\n') + '\n')

def collect_files(updater):
    summary_files = glob.glob('*.tsv')
    files_to_remove = ['empty_summaryline.tsv', '*_centar_output.tsv', 'Phoenix_Summary.tsv']
    summary_files = [file for file in summary_files
                        if not any(fnmatch.fnmatch(file, pattern) for pattern in files_to_remove)]
    if updater is True:
        old_summary_files = glob.glob('*/*/*_summaryline.tsv')
        # Filter out items containing any of the substrings
        filtered_list = [file for file in old_summary_files if not any(substring in file for substring in summary_files)]
        summary_files = filtered_list + summary_files
    return summary_files

def main():
    args = parseArgs()
    summary_files = collect_files(args.updater)
    List_TSV(args.output_file, summary_files, args.busco)

if __name__ == '__main__':
    main()