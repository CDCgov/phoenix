#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse
import glob
from pathlib import Path
import pandas as pd

# Function to get the script version
def get_version():
    return "1.0.0"

def parse_args(args=None):
    Description = "Reformat cdcgov/phoenix samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

def check_for_duplicates(file_in):
    df = pd.read_csv(file_in)
    duplicate_mask = df.duplicated()
    if duplicate_mask.any():
        duplicate_rows = df[duplicate_mask]
        for index, row in duplicate_rows.iterrows():
            raise ValueError(f"Duplicate row found: {row}. The pair of sample name and directory must be unique.")

def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,fastq_1,fastq_2
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz
    SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,


    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "directory"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        sample_name_list = [] # used to check if sample name has been used before
        Read_list = []
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check for duplicate sample names
            sample_name = line.split(",")[0]
            if sample_name in sample_name_list:
                print_error(
                    "The sample id {} is used multiple times! IDs need to be unique.".format(sample_name),
                    "Line",
                    line,
                )
            else:
                sample_name_list.append(sample_name)

            files = []
            # Define the file path
            dir = line.split(",")[1]
            sample_folder = line.split(",")[0]
            if str(dir).strip().endswith('/'):
                path = str(dir).strip()[:-1]
                sample_path = path
                project_path = "/".join(sample_path.split("/")[:-1])
            else:
                path = str(dir).strip()
                sample_path = path
                project_path = "/".join(sample_path.split("/")[:-1])
            #files.append(path + "/" + sample_folder + "/file_integrity/" + sample_name + "_scaffolds_summary.txt")
            files.append(sample_path + "/fastp_trimd/" + sample_name + "_1.trim.fastq.gz")
            files.append(sample_path + "/fastp_trimd/" + sample_name + "_2.trim.fastq.gz")
            files.append(sample_path + "/assembly/" + sample_name + ".filtered.scaffolds.fa.gz")
            files.append(sample_path + "/annotation/" + sample_name + ".faa")
            files.append(sample_path + "/annotation/" + sample_name + ".gff")
            files.append(sample_path + "/" + sample_name + ".tax")
            files.append(sample_path + "/" + sample_name + "_summaryline.tsv")
            files.append(sample_path + "/" + sample_name + ".synopsis")
            files.append(project_path + "/" + "Phoenix_Summary.tsv")
            # Handle glob searches with potential errors
            #try:
                # Find the position of the last occurrence of "/"
                #last_slash_index = path.rfind('/')
                # Slice the string up to and including the last slash
                #project_path = path[:last_slash_index + 1]
            try:
                full_path = path + "/file_integrity/"
                fairy_file = glob.glob(path + "/"+ sample_name + "/file_integrity/" + sample_name + "*summary.txt")[0]
                print(f"fairy_file found at " + fairy_file)
                # check that the sample did not fail the file_integrity check
                with open(fairy_file, 'r') as file:
                    for line in file:
                        if 'FAILED' in line:
                            print("The file {} states the file failed integrity checks and should not be included in the analysis. Please remove this sample from the analysis.".format(fairy_file))
                            sys.exit(1)
            except IndexError:
                #raise ValueError(f"No *_summary.tsv file found in {full_path}.")
                print(f"No *_summary.tsv file found in {full_path}.")
            try:
                glob.glob(project_path + "/*_GRiPHin_Summary.tsv")[0]
            except IndexError:
                raise ValueError(f"No *_GRiPHin_Summary.tsv file found in {project_path}.")
            for file_path in files:
                # Check if the file exists
                if Path(file_path).exists():
                    pass
                else:
                    raise ValueError("The file {} does not exist and is required for the pipeline. Please remove this sample from the analysis.".format(file_path))

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, directory = lspl[: len(HEADER)]
             # Define the file path
            if str(directory).strip().endswith('/'):
                directory = str(directory).strip()[:-1]
            else:
                directory = str(directory).strip()
            sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Create sample mapping dictionary = { sample: [ directory ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [directory]
            else:
                if directory in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(directory)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "directory"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "Sample: {}".format(sample))
                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write("{},{}\n".format(sample, val))
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_for_duplicates(args.FILE_IN)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)

if __name__ == "__main__":
    sys.exit(main())
