#!/usr/bin/env python3

import pandas as pd
import argparse

##Edits summaryline file
##Usage: >python edit_line_summary.py -i phoenix_summaryline.tsv
## Written by Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Edits a PhoeNix summary line for samples that fail the scaffold count check.')
    parser.add_argument('-i', '--input', dest="file", required=False, help='PhoeNix summary line')
    return parser.parse_args()

def edit_line(file):
    with open(file, "a+") as f:
        summary_file = pd.read_csv(file, sep="\t")
        summary_file["Auto_QC_Failure_Reason"] = "No scaffolds left after filtering!"
        summary_file["Auto_QC_Outcome"] = "FAIL"
        summary_file.to_csv(file, index=False, header=True, sep='\t')

def main():
    args = parseArgs()
    # if the output file already exists remove it
    edit_line(args.file)

if __name__ == '__main__':
    main()