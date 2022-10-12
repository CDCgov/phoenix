#!/usr/bin/env python3

import sys
import glob
import os
from decimal import *
getcontext().prec = 4
import argparse

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python Phoenix_Summary_tsv_06-10-22.py -o Output_Report.tsv Line_File_1 Line_File_2 Line_File3
## Written by Rich Stanton (njr5@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet')
    parser.add_argument('-o', '--out', required=True, help='output file name')
    parser.add_argument('-b', '--busco', action='store_true', help='parameter to know if busco was run')
    parser.add_argument('files', nargs=argparse.REMAINDER)
    return parser.parse_args()

def List_TSV(output_file, input_list, busco):
    Out = open(output_file, 'w')
    if not busco: # if --busco is not passed
        Out.write('ID\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tSpecies\tTaxa_Confidence\tTaxa_Source\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason\n')
    if busco:
        Out.write('ID\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tBUSCO\tBUSCO_DB\tSpecies\tTaxa_Confidence\tTaxa_Source\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason\n')
    try: #if there is numbers in the name then use that to sort
        input_list_sorted=sorted(input_list, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    except: #if no numbers then use only alphabetically
        input_list_sorted=sorted(input_list)
    for entry in input_list_sorted:
        f = open(entry, 'r')
        String1 = f.readline()
        Out.write(String1 + '\n')
        f.close()
    Out.close()

args = parseArgs()
List_TSV(args.out, args.files, args.busco)
