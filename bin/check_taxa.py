#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import csv
import os, sys
import fileinput

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python GRiPHin.py -s ./samplesheet.csv -a ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output --phoenix --scaffolds
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-o', '--format_ani_output', dest='format_ani_output', default="", required=False, help='Name of output file for formatted ANI results. Default is empty, which will not save the file.')
    parser.add_argument('-f','--format_ani_file', dest='format_ani_file', required=False, help='The coverage cut off default is 30x.')
    parser.add_argument('-t','--tax_file', dest='tax_file', required=False, help='The .tax file.')
    parser.add_argument('-a','--ani_file', dest='ani_file', required=False, help='The coverage cut off default is 30x.')
    parser.add_argument('-s','--shigapass_file', dest="shigapass_file", default=False, help='Turn on with --scaffolds to keep samples from failing/warnings/alerts that are based on trimmed data. Default is off.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def main(shigapass_file, format_ani_file, ani_file, tax_file):
    args = parseArgs()
    # Step 1: Open CSV file and check the second line to see if shigella was identified or not
    with open(shigapass_file, "r") as csv_file:
        reader = csv.reader(csv_file)
        next(reader)  # Skip the header row
        second_line = next(reader, None)

        if second_line and "Not Shigella/EIEC" not in second_line[0]:
            print("Taxa Identification was correct. Exiting.")
            try:
                os.rename(format_ani_file, args.output)
                sys.exit(0)
            except OSError as e:
                print(f"Error renaming file: {e}")
                exit(1)

    # Step 2: If the string is present, find the required line
    with open(ani_file, "r") as csv_file:
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
    df = pd.read_csv(format_ani_file, sep="\t")
    df["Source File"] = genome
    df["Organism"] = "Escherichia coli"
    df["% ID"] = round(percent_ani_match, 2)
    df["% Coverage"] = round(best_coverage, 2)

    # Save updated DataFrame to a new file
    df.to_csv(args.output, sep="\t", index=False)
    print(f"{args.output} updated successfully.")

    # Step 4: Update the taxonomy file
    # Read the tax file and update Shigella to Escherichia and species to coli
    with open(tax_file, "r") as old_tax_file:
        lines = old_tax_file.readlines()

    # Update the lines - only change G: and s: lines, preserve everything else
    updated_lines = []
    for line in lines:
        if line.startswith("G:") and "Shigella" in line:
            # Change Shigella to Escherichia (G:620 Shigella -> G:561 Escherichia)
            updated_lines.append("G:561\tEscherichia\n")
        elif line.startswith("s:"):
            # Change any species to coli (s:623 flexneri -> s:562 coli)
            updated_lines.append("s:562\tcoli\n")
        else:
            # Keep all other lines exactly as they are
            updated_lines.append(line)

    # Write the updated content back to the tax file
    with open(args.tax_file, "w") as updated_tax_file:
        updated_tax_file.writelines(updated_lines)

    print(f"{args.tax_file} updated successfully.")

def check_tax(shigapass_file, tax_file):
    with open(tax_file, "r") as f:
        for line in f:
            if line.startswith("G:"):
                tax_genus = line.split("\t")[1].strip()
            elif line.startswith("s:"):
                tax_species = line.split("\t")[1].strip()
    with open (shigapass_file, "r") as csv_file:
        reader = csv.reader(csv_file)
        next(reader)  # Skip the header row
        second_line = next(reader, None)
        #in nextflow we checked if tax in .tax was Escherichia and if summary contained "Not Shigella/EIEC" -> so the only option if "Not Shigella/EIEC" is in the file is for the tax to be Shigella
        if tax_genus == "Shigella" and "Not Shigella/EIEC" in second_line[0]:
            shiga_to_ecoli = True
            diff_shiga = ecoli_to_shiga = False
        elif tax_genus == "Escherichia" and "Not Shigella/EIEC" not in second_line[0]:
            shiga_to_ecoli = diff_shiga = False
            ecoli_to_shiga = True
        else:
            diff_shiga = True
            shiga_to_ecoli = ecoli_to_shiga = False
    return diff_shiga, shiga_to_ecoli, ecoli_to_shiga

def convert_ecoli_to_shiga_or_update_shiga(shigapass_file, format_ani_file, ani_file, tax_file):
    percent_id = update_ani_file(format_ani_file, ani_file, "Shigella_")
    #step 1: update taxonomy file
    # Find species by checking file content once
    with open(shigapass_file) as f:
            content = f.read()
    for marker, sp in [("SS", "s:624\tsonnei\n"), ("SF1-5", "s:623\tflexneri\n"), ("SB", "s:621\tboydii\n"), ("SD", "s:622\tdysenteriae\n")]:
        if marker in content:
            species = sp
            break
    # Write taxonomy file
    with open(tax_file, 'w') as f:
        f.write(f"ShigaPass\t{percent_id}\t{shigapass_file}\nK:2\tBacteria\nP:1224\tPseudomonadota\nC:1236\tGammaproteobacteria\nO:91347\tEnterobacterales\nF:543\tEnterobacteriaceae\nG:620\tShigella\n")
        if species:
            f.write(f"{species}\n")


def convert_shiga_to_ecoli(format_ani_file, ani_file, tax_file):
    percent_id = update_ani_file(format_ani_file, ani_file, "Escherichia_coli")
    # Write taxonomy file
    with open(tax_file, 'r') as f:
        # Update the lines - only change G: and s: lines, preserve everything else
        updated_lines = []
        lines = f.readlines()
    with open(tax_file, 'w') as f2:
        for line in lines:
            if line.startswith("G:") and "Shigella" in line:
                # Change Shigella to Escherichia (G:620 Shigella -> G:561 Escherichia)
                updated_lines.append("G:561\tEscherichia\n")
            elif line.startswith("s:"):
                # Change any species to coli (s:623 flexneri -> s:562 coli)
                updated_lines.append("s:562\tcoli\n")
            else:
                # Keep all other lines exactly as they are, but modify the first line if needed
                if len(updated_lines) == 0:  # This is the first line
                    # Modify the first line here
                    new_line = "ANI_REFSEQ\t" + str(percent_id) + "\t" + line.split("\t")[2]
                    modified_line = line.replace(line, new_line)  # Example modification
                    updated_lines.append(modified_line)
                else:
                    # Keep all other lines exactly as they are
                    updated_lines.append(line)
        f2.writelines(updated_lines)

def update_ani_file(format_ani_file, ani_file, taxa_string):
    # Step 2: update the ani_file - find the first line with Escherichia_coli and collect values
    with open(ani_file, "r") as csv_file:
        for line in csv_file:
            if taxa_string in line:
                matching_line = line
                break
        else:
            raise ValueError("No line with " + taxa_string + " found.")

    # Parse the line by tabs
    parts = matching_line.strip().split("\t")
    if len(parts) < 5:
        raise ValueError("Unexpected format in " + taxa_string + " line.")

    genome = parts[1].replace("reference_dir/", "")
    organism = genome.split("_")[0] + " " + genome.split("_")[1] # Extract organism name from the genome
    percent_ani_match = float(parts[2])
    fragment_matches = int(parts[3])
    total_fragments = int(parts[4])

    # Calculate best_coverage using bash logic
    best_coverage = round((100 * fragment_matches / total_fragments), 2)

    # Step 3: Update the file using pandas
    df = pd.read_csv(format_ani_file, sep="\t")
    df["Source File"] = genome
    df["Organism"] = organism
    df["% ID"] = round(percent_ani_match, 2)
    df["% Coverage"] = round(best_coverage, 2)

    # Save updated DataFrame to a new file
    df.to_csv(args.format_ani_output, sep="\t", index=False)
    print(f"{args.format_ani_output} updated successfully.")

    #return the percent_ani_match to add to the .tax file
    return round(percent_ani_match, 2)


if __name__ == '__main__':
    args = parseArgs()
    sample_id = args.shigapass_file.replace("_ShigaPass_summary.csv","")  # Extract sample ID from the file name
    diff_shiga, shiga_to_ecoli, ecoli_to_shiga = check_tax(args.shigapass_file, args.tax_file)
    if ecoli_to_shiga or diff_shiga:
        convert_ecoli_to_shiga_or_update_shiga(args.shigapass_file, args.format_ani_file, args.ani_file, args.tax_file)
    elif shiga_to_ecoli:
        convert_shiga_to_ecoli(args.format_ani_file, args.ani_file, args.tax_file)
    else:
        print("No updates needed for taxa identification.")
        sys.exit(0)