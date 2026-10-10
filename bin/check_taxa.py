#!/usr/bin/env python3

import pandas as pd
import argparse
import sys
import shutil

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python GRiPHin.py -s ./samplesheet.csv -a ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output --phoenix --scaffolds
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
__version__ = "1.1.0"

#                      ShigaPass
#                          │
#               ┌──────────┴──────────┐
#               │                     │
#     contains 'EIEC' key words    NOT contains 'EIEC' key words
#               │                     │
#         target = E. coli      target = Shigella
#               │                     │
#      ┌────────┴───────┐      ┌──────┴──────────┐
#  FastANI          FastANI   FastANI          FastANI
#  Shigella         E. coli   E. coli          Shigella
#     │                │         │                │
#  convert            no       convert          update
#  to E. coli       change    to Shigella       species 

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-o', '--format_ani_output', dest='format_ani_output', default="", required=False, help='Name of output file for formatted ANI results. Default is empty, which will not save the file.')
    parser.add_argument('-f','--format_ani_file', dest='format_ani_file', required=False, help='The coverage cut off default is 30x.')
    parser.add_argument('-t','--tax_file', dest='tax_file', required=False, help='The .tax file.')
    parser.add_argument('-a','--ani_file', dest='ani_file', required=False, help='The coverage cut off default is 30x.')
    parser.add_argument('-s','--shigapass_file', dest="shigapass_file", default=False, help='Turn on with --scaffolds to keep samples from failing/warnings/alerts that are based on trimmed data. Default is off.')
    parser.add_argument("-V", "--version",  action="version", version=f"%(prog)s: {__version__}")# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def convert_ecoli_to_shiga_or_update_shiga(shigapass_file, format_ani_file, ani_file, tax_file):
    percent_id = update_ani_file(format_ani_file, ani_file, "Shigella_")
    #step 1: update taxonomy file
    # Find species by checking file content once
    with open(shigapass_file) as f:
            second_line = f.readlines()[1]
            Predicted_Serotype = second_line.split(";")[7]
    species = None
    for marker, sp in [("SS", "s:624\tsonnei\n"), ("SF", "s:623\tflexneri\n"), ("SB", "s:621\tboydii\n"), ("SD", "s:622\tdysenteriae\n"), ("Shigella spp.", "s:625\tShigella sp.\n")]:
        print(f"Checking for marker '{marker}' in tax file...")
        if marker in Predicted_Serotype:
            species = sp
            break
    if species is None:
        raise ValueError(
            f"Unable to determine Shigella species from ShigaPass result: {Predicted_Serotype}"
        )
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
   
    with open (args.shigapass_file, "r") as file:
        lines = [line.strip() for line in file if line.strip()]
        if len(lines) < 2:
            sys.exit(f"Error: The file '{args.shigapass_file}' is missing header or summary data.")
        elif len(lines) > 2:
            sys.exit(f"Error: The file '{args.shigapass_file}' has more than 2 lines, indicating multiple samples.")

        # ShigaPass is only run for samples classified as Escherichia or Shigella
        # by FastANI, so these are the only two possible genera at this stage.
        # We use the ShigaPass result to reconcile the FastANI taxonomy.
        #
        # ShigaPass produces three relevant result types:
        #   - Shigella type       -> Shigella
        #   - EIEC                -> Escherichia
        #   - Not Shigella/EIEC   -> Escherichia
        #
        # Note that "Not Shigella/EIEC" does not mean "not Escherichia coli".
        # It means the sample is neither Shigella nor enteroinvasive E. coli (EIEC);
        # it may still be a non-EIEC E. coli strain.
        #
        # Therefore, if the ShigaPass result contains "EIEC", the sample is treated
        # as Escherichia; otherwise, it is treated as Shigella.
        shigapass_line = lines[1]
        shigapass_genus = "Escherichia" if "EIEC" in shigapass_line else "Shigella"

    tax_genus = None
    with open(args.tax_file, "r") as f:
        for line in f:
            if line.startswith("G:"):
                tax_genus = line.split("\t")[1].strip()
                break
    
    if tax_genus is None:
        sys.exit(f"Error: Could not find genus line (G:) in {args.tax_file}")

    if tax_genus not in {"Escherichia", "Shigella"}:
        sys.exit(
            f"Error: Unexpected FastANI genus '{tax_genus}'. "
            "ShigaPass reconciliation expects Escherichia or Shigella."
        )

    if shigapass_genus == "Shigella" and tax_genus == "Escherichia":
        print(f"{CRED}Taxa Identification changed from Escherichia coli to Shigella for sample {sample_id}.{CEND}")
        convert_ecoli_to_shiga_or_update_shiga(args.shigapass_file, args.format_ani_file, args.ani_file, args.tax_file)
    elif shigapass_genus == "Shigella" and tax_genus == "Shigella":
        print(f"{CRED}Taxa Identification updated for Shigella species for sample {sample_id}.{CEND}")
        convert_ecoli_to_shiga_or_update_shiga(args.shigapass_file, args.format_ani_file, args.ani_file, args.tax_file)
    elif shigapass_genus == "Escherichia" and tax_genus == "Shigella":
        print(f"{CRED}Taxa Identification changed from Shigella to Escherichia coli for sample {sample_id}.{CEND}")
        convert_shiga_to_ecoli(args.format_ani_file, args.ani_file, args.tax_file)
    else:
        print("No updates needed for taxa identification.")
        if args.format_ani_output and args.format_ani_file != args.format_ani_output:
            shutil.copyfile(args.format_ani_file, args.format_ani_output)
            print(f"{args.format_ani_output} copied successfully.")
        sys.exit(0)