#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import csv
import os, sys

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python GRiPHin.py -s ./samplesheet.csv -a ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output --phoenix --scaffolds
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-o', '--output', dest='output', default="", required=False, help='Name of output file default is GRiPHin_Summary.xlsx.')
    parser.add_argument('-f','--format_ani_file', dest='format_ani_file', required=False, help='The coverage cut off default is 30x.')
    parser.add_argument('-a','--ani_file', dest='ani_file', required=False, help='The coverage cut off default is 30x.')
    parser.add_argument('-s','--shigapass_file', dest="shigapass_file", default=False, help='Turn on with --scaffolds to keep samples from failing/warnings/alerts that are based on trimmed data. Default is off.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def main():
    args = parseArgs()
    # Step 1: Open CSV file and check the second line to see if shigella was identified or not
    try:
        with open(args.shigapass_file, "r") as csv_file:
            reader = csv.reader(csv_file)
            next(reader)  # Skip the header row
            second_line = next(reader, None)

            if second_line and "Not Shigella/EIEC" not in second_line[0]:
                print("Taxa Identification was correct. Exiting.")
                try:
                    os.rename(args.shigapass_file, args.output)
                    sys.exit(0)
                except OSError as e:
                    print(f"Error renaming file: {e}")
                    exit(1)

        # Step 2: If the string is present, find the required line
        with open(args.ani_file, "r") as csv_file:
            for line in csv_file:
                if "Escherichia_coli" in line:
                    escherichia_coli_line = line
                    break
            else:
                raise ValueError("No line with 'Escherichia_coli' found.")

        # Parse the line by tabs
        parts = escherichia_coli_line.strip().split("\t")
        if len(parts) < 5:
            raise ValueError("Unexpected format in Escherichia_coli line.")

        genome = parts[1].replace("reference_dir/", "")
        percent_ani_match = float(parts[2])
        fragment_matches = int(parts[3])
        total_fragments = int(parts[4])

        # Calculate best_coverage using bash logic
        best_coverage = round((100 * fragment_matches / total_fragments), 2)

        # Step 3: Update the file using pandas
        df = pd.read_csv(args.format_ani_file, sep="\t")
        df["Source File"] = genome
        df["Organism"] = "Escherichia coli"
        df["% ID"] = round(percent_ani_match, 2)
        df["% Coverage"] = round(best_coverage, 2)

        # Save updated DataFrame to a new file
        df.to_csv(args.output, sep="\t", index=False)
        print("File updated successfully.")

    except FileNotFoundError as e:
        print(f"File not found: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == '__main__':
    main()