#!/usr/bin/env python3

import sys
import glob
import os
import argparse
import csv

##Script to generate a samplesheet with sample,directory columns
##Usage: >python create_samplesheet.py -d input_directory
## Written by Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a samplesheet with sample,directory columns.')
    parser.add_argument('-d', '--directory', default=None, required=False, dest='directory', help='Will create a samplesheet for all samples in the directory (expects PHoeNIx style directory structure).')
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def create_samplesheet(directory):
    """Function will create a samplesheet from samples in a directory if -d argument passed."""
    directory = os.path.abspath(directory) # make sure we have an absolute path to start with
    with open("GRiPHin_samplesheet_created.csv", "w") as samplesheet:
        samplesheet.write('sample,directory\n')
        dirs = os.listdir(directory)
        skip_list_a = glob.glob(directory + "/*_GRiPHin_Summary.xlsx") # for if griphin is run on a folder that already has a report in it
        skip_list_a2 = glob.glob(directory + "/*_GRiPHin_Summary.tsv") # for if griphin is run on a folder that already has a report in it
        skip_list_a = [ item.split('/')[-1] for item in skip_list_a ]  # just get the excel name not the full path
        skip_list_a2 = [ item.split('/')[-1] for item in skip_list_a2 ]  # just get the excel name not the full path
        skip_list_b = [ "Phoenix_Summary.tsv", "pipeline_info", "GRiPHin_Summary.xlsx","BiosampleAttributes_Microbe.1.0.xlsx", "Sra_Microbe.1.0.xlsx", "multiqc", "samplesheet_converted.csv", "GRiPHin_samplesheet.csv"]
        skip_list = skip_list_a + skip_list_a2 + skip_list_b
        for sample in dirs:
            if sample not in skip_list:
                if directory[-1] != "/": # if directory doesn't have trailing / add one
                    directory = directory + "/"
                samplesheet.write(sample + "," + directory + sample + '\n')
    samplesheet = "GRiPHin_samplesheet_created.csv"
    return samplesheet

def main():
    args = parseArgs()
    # If a directory is given then create a samplesheet from it if not use the samplesheet passed
    if args.directory !=None:
        samplesheet = create_samplesheet(args.directory)
    else:
        sys.exit(CRED + "You MUST pass a stop directory of PHoeNIx output to create a samlesheet.\n" + CEND)

if __name__ == '__main__':
    main()