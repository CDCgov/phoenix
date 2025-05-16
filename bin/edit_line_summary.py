#!/usr/bin/env python3

import pandas as pd
import argparse

##Edits summaryline file
##Usage: >python edit_line_summary.py -i phoenix_summaryline.tsv
## Written by Jill Hagey (qpk9@cdc.gov)


# Function to get the script version
def get_version():
    return "2.0.0"


def parseArgs(args=None):
    parser = argparse.ArgumentParser(
        description="Edits a PhoeNix summary line for samples that fail the scaffold count check."
    )
    parser.add_argument(
        "-i", "--input", dest="file", required=False, help="PhoeNix summary line"
    )
    parser.add_argument(
        "--version", action="version", version=get_version()
    )  # Add an argument to display the version
    return parser.parse_args()


def edit_line(file):
    with open(file, "a+") as f:
        summary_file = pd.read_csv(file, sep="\t")
        summary_file["Auto_QC_Failure_Reason"] = "No scaffolds left after filtering!"
        summary_file["Auto_QC_Outcome"] = "FAIL"
        summary_file.to_csv(file, index=False, header=True, sep="\t")


def main():
    args = parseArgs()
    # if the output file already exists remove it
    edit_line(args.file)


if __name__ == "__main__":
    main()
