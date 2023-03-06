#!/usr/bin/env python3

import glob
import shutil
import argparse
import re

## Makes a summary Excel file when given a run folder from PhoeNiX
## Usage: >python get_taxa_for_amrfinder.py -t tax_file -m mutation_tab_file -o out_File
## Written by Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    """Takes in a taxa file and creates a file with the taxa found"""
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary line')
    parser.add_argument('-t', '--taxa', dest="taxa_file", required=True, help='.tax file that comes from determine_taxID.sh')
    parser.add_argument('-o', '--output', dest="output", required=True, help='name of output file')
    return parser.parse_args()

def get_taxa(taxa_file):
    """Open Taxonomy file and get genus, species and a combo variable"""
    with open(taxa_file, "r") as f:
        for line in f:
            if line.startswith("G:"):
                genus_match = line.split(":")[1]
                genus_match = re.sub( "\d+|\n|\s|\.", '', genus_match)
            if line.startswith("s:"):
                species_match = line.split(":")[1]
                species_match = re.sub( "\d+|\n|\.|\s", '', species_match)
        gen_sp = genus_match + "_" + species_match
    return genus_match, species_match, gen_sp

def taxa_check(genus_match, species_match, gen_sp):
    # species and genus lists
    #species = ['Enterococcus_faecalis','Enterococcus_faecium','Staphylococcus_aureus','Staphylococcus_pseudintermedius','Streptococcus_agalactiae', 'Streptococcus_pyogenes', \
    #'Clostridioides_difficile', 'Pseudomonas_aeruginosa', 'Vibrio_cholerae','Neisseria_gonorrhoeae', 'Neisseria_meningitidis'] # for update
    species = ['Enterococcus_faecalis','Enterococcus_faecium','Staphylococcus_aureus','Staphylococcus_pseudintermedius','Streptococcus_agalactiae', 'Streptococcus_pyogenes', \
    'Clostridioides_difficile', 'Pseudomonas_aeruginosa', 'Vibrio_cholerae']
    genus = ['Salmonella']
    #complex_genera = ['Acinetobacter', 'Neisseria', 'Escherichia', 'Klebsiella','Campylobacter', 'Shigella', 'Streptococcus'] # for update
    complex_genera = ['Acinetobacter', 'Neisseria', 'Escherichia', 'Klebsiella','Campylobacter', 'Neisseria', 'Shigella', 'Streptococcus']
    if gen_sp in species:
        taxa=gen_sp
    elif genus_match in genus:
        taxa=genus_match
    elif genus_match in complex_genera:
        taxa=complicated_genus(gen_sp, genus_match)
    else:
        taxa=str("No Match Found")
    return taxa

def complicated_genus(gen_sp, genus_match):
    # Genus lists
    Escherichia = ["Escherichia", "Shigella"]
    # Genus_species lists
    Klebsiella = ["Klebsiella_pneumoniae", "Klebsiella_aerogenes", "Klebsiella_oxytoca"]
    Neisseria = ['Neisseria_gonorrhea', 'Neisseria_meningitidis']  # remove for update
    Campylobacter = ['Campylobacter_coli','Campylobacter_jejuni']
    Streptococcus_pneumoniae = ['Streptococcus_pneumoniae', 'Streptococcus_mitis']
    Acinetobacter = ['Acinetobacter_baumannii','Acinetobacter_calcoaceticus', 'Acinetobacter_lactucae', 'Acinetobacter_nosocomialis', 'Acinetobacter_pittii', 'Acinetobacter_seifertii']
    # create Acinetobacter sp. names
    AB_complex_sp = []
    AB_sp = [ '0000-0051-1906','0000-0051-1909','0000-0051-8448','0000-0052-5682','0000-0052-7200','0000-0053-2562','0000-0053-7518','0000-0060-4277','0000-0082-5590','0000-0091-0263',
    '1000160','1130196','1179249','1239920','1245249','1245593','1264765','1281984','1289694','1294243','1294596','1295259','1396970','1424608','1461402','1475718','1542444','1564232',
    '1566109','1578804','1592897','216872','21871','225588','230853','259052','25977_1','25977_10','25977_2','25977_3','25977_4','25977_6','25977_7','25977_8','263903-1','263903-2',
    '272263','478810','479375','694762','716(2011)','723929','72431','735(2011)','742879','766875','796380-1375','809848','826659','869535','883425','88816','907131','983759','ABNIH27',
    'AR_0276','CB15','FDAARGOS_131','LUH 07045','LUH 10726','LUH 14530','OIFC021','RUH 14531','TG27347','WC-141','WC-323']
    for sp in AB_sp:
        name="Acinetobacter_sp._" + sp
        AB_complex_sp.append(name)
    # Look for correct match
    if genus_match in Escherichia:
        taxa=str("Escherichia")
    elif gen_sp in Klebsiella:
        taxa=str("Klebsiella")
#        if gen_sp == "Klebsiella_pneumoniae": # for update
#            taxa=str("Klebsiella_pneumoniae")
#        elif gen_sp == "Klebsiella_aerogenes":
#            taxa=str("Klebsiella_oxytoca")
#        elif gen_sp == "Klebsiella_oxytoca":
#            taxa=str("Klebsiella_oxytoca")
    elif gen_sp in Streptococcus_pneumoniae:
        taxa=str("Streptococcus_pneumoniae")
    elif gen_sp in Campylobacter:
        taxa=str("Campylobacter")
    elif gen_sp in Acinetobacter:
        taxa=str("Acinetobacter_baumannii")
    elif gen_sp in AB_complex_sp:
        taxa=str("Acinetobacter_baumannii")
    elif gen_sp in Neisseria: #remove for update
        taxa=str("Neisseria") #remove for update
    else:
        taxa="No Match Found"
    return taxa

def write_file(taxa, output_file):
    print(taxa)
    with open(output_file, 'w') as f:
        f.write(taxa)

def main():
    args = parseArgs()
    genus_match, species_match, gen_sp = get_taxa(args.taxa_file)
    taxa = taxa_check(genus_match, species_match, gen_sp)
    write_file(taxa, args.output)

if __name__ == '__main__':
    main()