#!/usr/bin/env python3

import pandas as pd
import numpy as np
import argparse
import sys

__version__ = "1.1.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--format_ani_output', dest='format_ani_output', default="", required=False)
    parser.add_argument('-f','--format_ani_file', dest='format_ani_file', required=False)
    parser.add_argument('-t','--tax_file', dest='tax_file', required=False)
    parser.add_argument('-a','--ani_file', dest='ani_file', required=False)
    parser.add_argument('-s','--shigapass_file', dest="shigapass_file", default=False)
    parser.add_argument("-V", "--version", action="version", version=f"%(prog)s: {__version__}")
    return parser.parse_args()

CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def check_tax(shigapass_file, tax_file):
    tax_genus = None
    tax_species = None
    with open(tax_file, "r") as f:
        for line in f:
            if line.startswith("G:"):
                tax_genus = line.split("\t")[1].strip()
            elif line.startswith("s:"):
                tax_species = line.split("\t")[1].strip()

    with open(shigapass_file, "r") as f:
        lines = f.readlines()
        second_line = lines[1].strip().split(";")
#        print(f"DEBUG second_line = {second_line}")

        is_ecoli_per_shigapass = any("EIEC" in field or "Not Shigella/EIEC" in field for field in second_line)
#        print(f"DEBUG is_ecoli_per_shigapass = {is_ecoli_per_shigapass}")

        predicted_flex_serotype = second_line[8] if len(second_line) > 8 else ""
        shigapass_species_map = {"SS": "sonnei", "SF": "flexneri", "SB": "boydii", "SD": "dysenteriae"}
        shigapass_species = next((sp for marker, sp in shigapass_species_map.items() if marker in predicted_flex_serotype), None)

#    print(f"DEBUG tax_genus = '{tax_genus}', tax_species = '{tax_species}'")

    if tax_genus == "Shigella" and is_ecoli_per_shigapass:
        return "shiga_to_ecoli"
    elif tax_genus == "Escherichia" and not is_ecoli_per_shigapass:
        return "ecoli_to_shiga"
    elif tax_genus == "Shigella" and not is_ecoli_per_shigapass and shigapass_species and tax_species != shigapass_species:
        return "diff_shiga"
    elif tax_genus == "Escherichia" and is_ecoli_per_shigapass:
        return "ani_only"
    else:
        return "none"


# def convert_ecoli_to_shiga_or_update_shiga(shigapass_file, format_ani_file, ani_file, tax_file):
#     percent_id = update_ani_file(format_ani_file, ani_file, "Shigella_")
#     percent_id_str = str(percent_id) if percent_id is not None else "N/A (ShigaPass override; no matching ANI hit found)"
#     with open(shigapass_file) as f:
#         second_line = f.readlines()[1]
#         fields = second_line.strip().split(";")
#         Predicted_FlexSerotype = fields[8] if len(fields) > 8 else ""
#     species = None
#     for marker, sp in [("SS", "s:624\tsonnei\n"), ("SF", "s:623\tflexneri\n"), ("SB", "s:621\tboydii\n"), ("SD", "s:622\tdysenteriae\n"), ("Shigella spp.", "s:625\tShigella sp.\n")]:
#         if marker in Predicted_FlexSerotype:
#             species = sp
#             break
#     with open(tax_file, 'w') as f:
#         f.write(f"ShigaPass\t{percent_id_str}\t{shigapass_file}\nK:2\tBacteria\nP:1224\tPseudomonadota\nC:1236\tGammaproteobacteria\nO:91347\tEnterobacterales\nF:543\tEnterobacteriaceae\nG:620\tShigella\n")
#         if species:
#             f.write(f"{species}\n")


# def convert_shiga_to_ecoli(format_ani_file, ani_file, tax_file):
#     percent_id = update_ani_file(format_ani_file, ani_file, "Escherichia_coli")
#     percent_id_str = str(percent_id) if percent_id is not None else "N/A (ShigaPass override; no matching ANI hit found)"
#     with open(tax_file, 'r') as f:
#         lines = f.readlines()
#     updated_lines = []
#     species_line_written = False
#     for line in lines:
#         if line.startswith("G:") and "Shigella" in line:
#             updated_lines.append("G:561\tEscherichia\n")
#         elif line.startswith("s:"):
#             updated_lines.append("s:562\tcoli\n")
#             species_line_written = True
#         else:
#             if len(updated_lines) == 0:
#                 new_line = "ShigaPass\t" + percent_id_str + "\t" + line.split("\t")[2]
#                 updated_lines.append(new_line)
#             else:
#                 updated_lines.append(line)
#     if not species_line_written:
#         updated_lines.append("s:562\tcoli\n")
#     with open(tax_file, 'w') as f2:
#         f2.writelines(updated_lines)

def convert_ecoli_to_shiga_or_update_shiga(shigapass_file, format_ani_file, ani_file, tax_file):
    percent_id, coverage = update_ani_file(format_ani_file, ani_file, "Shigella_")
    if percent_id is not None:
        percent_id_str = f"{percent_id};{coverage}"
    else:
        percent_id_str = "N/A (ShigaPass override; no matching ANI hit found)"
    with open(shigapass_file) as f:
        second_line = f.readlines()[1]
        fields = second_line.strip().split(";")
        Predicted_FlexSerotype = fields[8] if len(fields) > 8 else ""
    species = None
    for marker, sp in [("SS", "s:624\tsonnei\n"), ("SF", "s:623\tflexneri\n"), ("SB", "s:621\tboydii\n"), ("SD", "s:622\tdysenteriae\n"), ("Shigella spp.", "s:625\tShigella sp.\n")]:
        if marker in Predicted_FlexSerotype:
            species = sp
            break
    with open(tax_file, 'w') as f:
        f.write(f"ShigaPass\t{percent_id_str}\t{shigapass_file}\nK:2\tBacteria\nP:1224\tPseudomonadota\nC:1236\tGammaproteobacteria\nO:91347\tEnterobacterales\nF:543\tEnterobacteriaceae\nG:620\tShigella\n")
        if species:
            f.write(f"{species}\n")


def convert_shiga_to_ecoli(format_ani_file, ani_file, tax_file):
    percent_id, coverage = update_ani_file(format_ani_file, ani_file, "Escherichia_coli")
    if percent_id is not None:
        percent_id_str = f"{percent_id};{coverage}"
    else:
        percent_id_str = "N/A (ShigaPass override; no matching ANI hit found)"
    with open(tax_file, 'r') as f:
        lines = f.readlines()
    updated_lines = []
    species_line_written = False
    for line in lines:
        if line.startswith("G:") and "Shigella" in line:
            updated_lines.append("G:561\tEscherichia\n")
        elif line.startswith("s:"):
            updated_lines.append("s:562\tcoli\n")
            species_line_written = True
        else:
            if len(updated_lines) == 0:
                new_line = "ShigaPass\t" + percent_id_str + "\t" + line.split("\t")[2]
                updated_lines.append(new_line)
            else:
                updated_lines.append(line)
    if not species_line_written:
        updated_lines.append("s:562\tcoli\n")
    with open(tax_file, 'w') as f2:
        f2.writelines(updated_lines)


# def correct_ani_only(format_ani_file, ani_file, tax_file):
#     tax_genus = None
#     with open(tax_file, "r") as f:
#         for line in f:
#             if line.startswith("G:"):
#                 tax_genus = line.split("\t")[1].strip()
#                 break
#     taxa_string = "Escherichia_coli" if tax_genus == "Escherichia" else "Shigella_"
#     percent_id = update_ani_file(format_ani_file, ani_file, taxa_string)
#     percent_id_str = str(percent_id) if percent_id is not None else "N/A (ShigaPass override; no matching ANI hit found)"

#     with open(tax_file, 'r') as f:
#         lines = f.readlines()
#     if lines:
#         first_line_parts = lines[0].split("\t")
#         remainder = "\t".join(first_line_parts[2:]) if len(first_line_parts) > 2 else ""
#         lines[0] = f"ShigaPass\t{percent_id_str}\t{remainder}" if remainder else f"ShigaPass\t{percent_id_str}\n"
#     with open(tax_file, 'w') as f2:
#         f2.writelines(lines)

def correct_ani_only(format_ani_file, ani_file, tax_file):
    tax_genus = None
    with open(tax_file, "r") as f:
        for line in f:
            if line.startswith("G:"):
                tax_genus = line.split("\t")[1].strip()
                break
    taxa_string = "Escherichia_coli" if tax_genus == "Escherichia" else "Shigella_"
    update_ani_file(format_ani_file, ani_file, taxa_string)


# def update_ani_file(format_ani_file, ani_file, taxa_string):
#     matching_line = None
#     with open(ani_file, "r") as csv_file:
#         for line in csv_file:
#             if taxa_string in line:
#                 matching_line = line
#                 break

#     if matching_line is None:
#         df = pd.read_csv(format_ani_file, sep="\t")
#         df["Source File"] = "N/A - ShigaPass override"
#         df["Organism"] = taxa_string.rstrip("_").replace("_", " ") + " (ShigaPass override; no matching ANI hit in top results)"
#         df["% ID"] = np.nan
#         df["% Coverage"] = np.nan
#         df.to_csv(args.format_ani_output, sep="\t", index=False)
#         return None

#     parts = matching_line.strip().split("\t")
#     if len(parts) < 5:
#         raise ValueError("Unexpected format in " + taxa_string + " line.")

#     genome = parts[1].replace("reference_dir/", "")
#     organism = genome.split("_")[0] + " " + genome.split("_")[1]
#     percent_ani_match = float(parts[2])
#     fragment_matches = int(parts[3])
#     total_fragments = int(parts[4])
#     best_coverage = round((100 * fragment_matches / total_fragments), 2)

#     df = pd.read_csv(format_ani_file, sep="\t")
#     df["Source File"] = genome
#     df["Organism"] = organism
#     df["% ID"] = round(percent_ani_match, 2)
#     df["% Coverage"] = round(best_coverage, 2)
#     df.to_csv(args.format_ani_output, sep="\t", index=False)

#     return round(percent_ani_match, 2)

def update_ani_file(format_ani_file, ani_file, taxa_string):
    matching_line = None
    with open(ani_file, "r") as csv_file:
        for line in csv_file:
            if taxa_string in line:
                matching_line = line
                break

    if matching_line is None:
        df = pd.read_csv(format_ani_file, sep="\t")
        df["Source File"] = "N/A - ShigaPass override"
        df["Organism"] = taxa_string.rstrip("_").replace("_", " ") + " (ShigaPass override; no matching ANI hit in top results)"
        df["% ID"] = np.nan
        df["% Coverage"] = np.nan
        df.to_csv(args.format_ani_output, sep="\t", index=False)
        return None, None

    parts = matching_line.strip().split("\t")
    if len(parts) < 5:
        raise ValueError("Unexpected format in " + taxa_string + " line.")

    genome = parts[1].replace("reference_dir/", "")
    organism = genome.split("_")[0] + " " + genome.split("_")[1]
    percent_ani_match = float(parts[2])
    fragment_matches = int(parts[3])
    total_fragments = int(parts[4])
    best_coverage = round((100 * fragment_matches / total_fragments), 2)

    df = pd.read_csv(format_ani_file, sep="\t")
    df["Source File"] = genome
    df["Organism"] = organism
    df["% ID"] = round(percent_ani_match, 2)
    df["% Coverage"] = best_coverage
    df.to_csv(args.format_ani_output, sep="\t", index=False)

    return round(percent_ani_match, 2), best_coverage


if __name__ == '__main__':
    args = parseArgs()
    sample_id = args.shigapass_file.replace("_ShigaPass_summary.csv","")
    outcome = check_tax(args.shigapass_file, args.tax_file)
    print(f"DEBUG outcome = {outcome}")

    if outcome == "ecoli_to_shiga":
        convert_ecoli_to_shiga_or_update_shiga(args.shigapass_file, args.format_ani_file, args.ani_file, args.tax_file)
    elif outcome == "diff_shiga":
        convert_ecoli_to_shiga_or_update_shiga(args.shigapass_file, args.format_ani_file, args.ani_file, args.tax_file)
    elif outcome == "shiga_to_ecoli":
        convert_shiga_to_ecoli(args.format_ani_file, args.ani_file, args.tax_file)
    elif outcome == "ani_only":
        correct_ani_only(args.format_ani_file, args.ani_file, args.tax_file)
    else:
        print("No updates needed for taxa identification.")
        sys.exit(0)