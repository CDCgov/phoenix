#!/usr/bin/env python3

import sys
import glob
import os
import argparse
import csv

##Given a samplesheet from GRiPHin's wf a new samplesheet
##Usage: >python create_samplesheet.py -s ./samplesheet.csv -a ../PHX/phoenix/assets/databases/ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output
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
        skip_list = [ "Phoenix_Output_Report.tsv", "pipeline_info", "GRiPHin_Report.xlsx", "multiqc", "samplesheet_converted.csv", "GRiPHin_samplesheet.csv"]
        for sample in dirs:
            if sample not in skip_list:
                #with open("GRiPHin_samplesheet_created.csv", "a") as samplesheet:
                    if directory[-1] != "/": # if directory doesn't have trailing / add one
                        directory = directory + "/"
                    #print(sample + "," + directory + sample + '\n')
                    samplesheet.write(sample + "," + directory + sample + '\n')
                    #print(directory)
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