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
    parser.add_argument('files', nargs=argparse.REMAINDER)
    return parser.parse_args()

def List_TSV(output_file, input_list):
    Out = open(output_file, 'w')
    Out.write('ID\tQC_Outcome\tCoverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Contigs_>500bp\tSpecies\tTaxa_Confidence\tTaxa_Source\tMLST_Scheme\tMLST\tGC_%\tKraken2_Trimd\tKraken2_Weighted\tBeta_Lactam_Resistance_Genes\tOther_AR_Genes\tHypervirulence_Genes\tAMRFinder_Point_Mutations\tQC_Reason\n')
    for entry in input_list:
        f = open(entry, 'r')
        String1 = f.readline()
        Out.write(String1 + '\n')
        f.close()
    Out.close()

args = parseArgs()
List_TSV(args.out, args.files)