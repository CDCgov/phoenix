#!/usr/bin/env python3

from datetime import date
import os
import os.path
import re
from decimal import *
getcontext().prec = 4
import argparse
import yaml
import pandas as pd

## updates a sample's readme file when --pipeline update_phoenix is run
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary line')
    parser.add_argument('-d', '--sample_directory', dest="sample_directory",required=True, help="Directory for the sample's readme")
    parser.add_argument('-m', '--mlst_db', dest="mlst_db",required=True, help='New MLST database')
    parser.add_argument('-a', '--amrfinder_db', dest="amrfinder_db", required=True, help='New AMRFinder database')
    parser.add_argument('-g', '--ar_db', dest="ar_db",required=True, help='New combined AR db')
    parser.add_argument('--old_gamma', dest="old_gamma",required=True, help='Old <sample_id>_ResGANNCBI_<date>_srst2.gamma file.')
    parser.add_argument('--new_gamma', dest="new_gamma",required=True, help='New <sample_id>_ResGANNCBI_<date>_srst2.gamma file.')
    parser.add_argument('--old_ncbi', dest="old_ncbi",required=True, help='Old <sample_id>_<date>_all_genes.tsv file.')
    parser.add_argument('--new_ncbi', dest="new_ncbi",required=True, help='New <sample_id>_<date>_all_genes.tsv file.')
    parser.add_argument('--old_tax', dest="old_tax",required=False, help='Old <sample_id>.tax file.')
    parser.add_argument('--new_tax', dest="new_tax",required=False, help='New <sample_id>_updater_log.tax file.')
    parser.add_argument('-p', '--pipeline_info', dest="pipeline_info", default=True, help='software_versions.yml file') # Need this for when you call -entry CDC_PHOENIX or CDC_SCAFFOLDS, but spades fails
    parser.add_argument('-o', '--out', dest="output", required=True, help='output file name')
    parser.add_argument('-v', dest="current_phx_version", required=True, help='current phx version.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

def get_old_database_IDs(software_versions):
    with open(software_versions) as yml_file:
        yml_data = yaml.safe_load(yml_file)
        if 'GAMMA_AR' in yml_data:
            gamma_ar_data = yml_data['GAMMA_AR']
            # Extract values into variables
            old_gamma_ar_db = gamma_ar_data.get('Database', '')
            #gamma_ver = gamma_ar_data.get('gamma', '')
            #gamma_container = gamma_ar_data.get('gamma_container', '')
        if 'MLST' in yml_data:
            mlst_data = yml_data['MLST']
            # Extract values into variables
            #mlst_ver = mlst_data.get('mlst', '')
            old_mlst_db = mlst_data.get('mlst_db', '')
            #mlst_container = mlst_data.get('mlst_container', '')
        if 'AMRFINDERPLUS_RUN' in yml_data:
            amrfinder_data = yml_data['AMRFINDERPLUS_RUN']
            # Extract values into variables
            #amrfinderplus_ver = amrfinder_data.get('amrfinderplus', '')
            old_amrfinderplus_db = amrfinder_data.get('amrfinderplus_db_version', '')
            #amrfinderplus_container = amrfinder_data.get('amrfinderplus_container', '')
        if 'Workflow' in yml_data:
            phoenix_data = yml_data['Workflow']
            # Extract values into variables
            phoenix_ver = phoenix_data.get('cdcgov/phoenix', '')
    return old_gamma_ar_db, old_mlst_db, old_amrfinderplus_db, phoenix_ver

def get_new_database_IDs(mlst_db, amrfinder_db):
    # Use regular expression to extract the date part from the filename
    amrfinderplus_match = re.search(r'(\d{4})(\d{2})(\d{2})\.(\d{1})', amrfinder_db)
    new_amrfinderplus_db = amrfinderplus_match.group(1) + "-" + amrfinderplus_match.group(2) + "-" + amrfinderplus_match.group(3) +"." + amrfinderplus_match.group(4) 
    # doing the same for mlst
    mlst_match = re.search(r'(\d{4})(\d{2})(\d{2})', mlst_db)
    new_mlst_db = mlst_match.group(1) + "-" + mlst_match.group(2) + "-" + mlst_match.group(3) 
    return new_mlst_db, new_amrfinderplus_db

def extract_gamma_gene_name(gene_string):
    """
    Extract the gene name from the Gene column. We want the gene_name portion (second part after splitting by __)
    Format: ID__gene_name__additional_info__category__source
    """
    parts = gene_string.split('__')
    if len(parts) >= 3:
        # Return the gene name portion (second element) and the accession (third element)
        #return f"{parts[1]}__{parts[2]}"
        allele_name = parts[2].split('_')[0]
        return f"{allele_name}"
    return gene_string

def compare_ar_files(old_file, new_file, file_type):
    """Compare two AR files and identify added and dropped AR genes."""
    # Read the files
    old_df = pd.read_csv(old_file, sep='\t')
    new_df = pd.read_csv(new_file, sep='\t')
    # Determine gene column and extraction function based on file type
    if file_type.lower() == "gamma":
        gene_column = "Gene"
        extract_gene_name = extract_gamma_gene_name
        file_description = "gamma"
    elif file_type.lower() == "ncbi":
        old_df.rename(columns={'Gene_symbol': 'Element_symbol', 'Protein_identifier': 'Protein_id', 'Element_type': 'Type', 'Sequence_name': 'Element_name', 'HMM_id': 'HMM_accession', '%_Coverage_of_reference_sequence': '%_Coverage_of_reference',
                                '%_Identity_to_reference_sequence': '%_Identity_to_reference','Element_subtype': 'Subtype', 'Accession_of_closest_sequence': 'Closest_reference_accession', 'Name_of_closest_sequence': 'Closest_reference_name'}, inplace=True) 
        gene_column = "Element_symbol"
        extract_gene_name = lambda x: x  # NCBI files use the column value directly
        file_description = "NCBI"
    else:
        raise ValueError(f"Unsupported file_type: {file_type}. Must be 'gamma' or 'ncbi'")
    added_ar_genes = []
    dropped_ar_genes = []
    
    # Process old file genes
    print(f"Processing old {file_description} file...")
    for idx, old_row in old_df.iterrows():
        old_gene_name = extract_gene_name(old_row[gene_column])
        #print(f"Old gene: {old_gene_name}")
        
        # Find matching gene in new file
        match_found = False
        for new_idx, new_row in new_df.iterrows():
            new_gene_name = extract_gene_name(new_row[gene_column])
            
            # Check if gene names match (substring matching)
            if old_gene_name in new_gene_name or new_gene_name in old_gene_name:
                match_found = True
                # for debugging purposes, we can print that a match was found
                #print(f"  Found match: {new_gene_name}")
                
                # Compare entire rows for differences (excluding gene column for exact comparison)
                # We'll compare all columns except the gene column since we already know they might differ slightly
                old_row_data = old_row.drop(gene_column).to_dict()
                new_row_data = new_row.drop(gene_column).to_dict()
                
                if old_row_data != new_row_data:
                    print(f"  Row differences found for {old_gene_name}")
                    # Show the differences
                    for col in old_row_data:
                        if old_row_data[col] != new_row_data[col]:
                            print(f"    {col}: {old_row_data[col]} -> {new_row_data[col]}")
                    # Add to added_ar_genes if there are differences
                    added_ar_genes.append(old_gene_name)
                #else:
                    #for debugging purposes, we can print that no differences were found
                    #print(f"  No differences found for {old_gene_name}")
                break
        
        if not match_found:
            print(f"  No match found for {old_gene_name} - adding to dropped_ar_genes")
            dropped_ar_genes.append(old_gene_name)
    
    # Process new file genes to find truly new genes
    print(f"\nProcessing new {file_description} file for new genes...")
    for idx, new_row in new_df.iterrows():
        new_gene_name = extract_gene_name(new_row[gene_column])
        # Check if this gene exists in old file
        match_found = False
        for old_idx, old_row in old_df.iterrows():
            old_gene_name = extract_gene_name(old_row[gene_column])
            # Check if gene names match (substring matching)
            if new_gene_name in old_gene_name or old_gene_name in new_gene_name:
                match_found = True
                break
        if not match_found:
            print(f"New gene found: {new_gene_name} - adding to added_ar_genes")
            added_ar_genes.append(new_gene_name)
    
    print(f"\nAdded AR Genes ({len(added_ar_genes)}):")
    for gene in added_ar_genes:
        print(f"  + {gene}")

    print(f"\nDropped AR Genes ({len(dropped_ar_genes)}):")
    for gene in dropped_ar_genes:
        print(f"  - {gene}")

    return added_ar_genes, dropped_ar_genes

def write_readme(old_gamma, old_mlst, old_amrfinder, new_ar_db, new_mlst_db, new_amrfinderplus_db, output, phoenix_ver, sample_directory, current_phx_version, added_gamma_ar_genes, dropped_gamma_ar_genes, added_ncbi_ar_genes, dropped_ncbi_ar_genes, new_tax, old_tax):
    Gamma_db_updated = old_gamma + " --> " + new_ar_db.strip()
    amrfinder_db_updated = old_amrfinder + " --> " + new_amrfinderplus_db.strip()
    MLST_db_updated = old_mlst + " --> " + new_mlst_db.strip()
    phx_versions = phoenix_ver + " --> " + current_phx_version.strip()
    # If both old_tax and new_tax are empty strings, set tax_updated to empty string
    if old_tax == "" and new_tax == "":
        tax_updated = ""
    else:
        tax_updated = old_tax + " --> " + new_tax
    # Convert date to string format
    date_string = date.today().strftime('%Y-%m-%d')
    # Check if the file exists
    if os.path.exists(sample_directory + "/"+ output):
        # File exists, append a new line
        with open(output, 'a') as f:
            f.write(date_string + "\t" + phx_versions + "\t" + tax_updated + "\t" + MLST_db_updated + "\t" + Gamma_db_updated + "\t" + amrfinder_db_updated + '\t' + ','.join(added_gamma_ar_genes) + '\t' + ','.join(dropped_gamma_ar_genes) + '\t' + ','.join(added_ncbi_ar_genes) + '\t' + ','.join(dropped_ncbi_ar_genes) + '\n')
    else:
        # File doesn't exist, write header and new line
        with open(output, 'w') as f:
            f.write('Date\tPhoenix_Version\tTaxa\tMLST_DB\tGAMMA_DB\tAMRFinderPlus_DB\tAdded_GAMMA_AR_genes\tDropped/Changed_GAMMA_AR_genes\tAdded_NCBI_AR_genes\tDropped/Changed_NCBI_AR_genes\n')
            f.write(date_string + "\t" + phx_versions + "\t" + tax_updated + "\t" + MLST_db_updated + "\t" + Gamma_db_updated + "\t" + amrfinder_db_updated + '\t' + ','.join(added_gamma_ar_genes) + '\t' + ','.join(dropped_gamma_ar_genes) + '\t' + ','.join(added_ncbi_ar_genes) + '\t' + ','.join(dropped_ncbi_ar_genes) + '\n')

def compare_tax_files(old_tax_file, new_tax_file):
    # If tax files aren't provided, return empty strings
    if not old_tax_file or not new_tax_file:
        return "", ""
    # if the files are passed then get the taxa info
    with open(new_tax_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("G:"):
                new_genus = line.split("\t")[1]
            elif line.startswith("s:"):
                new_species = line.split("\t")[1]
    with open(old_tax_file, 'r') as f2:
        lines = f2.readlines()
        for line in lines:
            if line.startswith("G:"):
                old_genus = line.split("\t")[1]
            elif line.startswith("s:"):
                old_species = line.split("\t")[1]
    old_tax = old_genus.strip() + " " + old_species.strip()
    new_tax = new_genus.strip() + " " + new_species.strip()
    return new_tax, old_tax

def main():
    args = parseArgs()
    old_gamma, old_mlst, old_amrfinder, phoenix_ver = get_old_database_IDs(args.pipeline_info)
    new_mlst_db, new_amrfinderplus_db = get_new_database_IDs(args.mlst_db, args.amrfinder_db)
    added_gamma_ar_genes, dropped_gamma_ar_genes = compare_ar_files(args.old_gamma, args.new_gamma, file_type="gamma")
    added_ncbi_ar_genes, dropped_ncbi_ar_genes = compare_ar_files(args.old_ncbi, args.new_ncbi, file_type="ncbi")
    new_tax, old_tax = compare_tax_files(args.old_tax, args.new_tax)
    write_readme(old_gamma, old_mlst, old_amrfinder, args.ar_db, new_mlst_db, new_amrfinderplus_db, args.output, phoenix_ver, args.sample_directory, args.current_phx_version, added_gamma_ar_genes, dropped_gamma_ar_genes, added_ncbi_ar_genes, dropped_ncbi_ar_genes, new_tax, old_tax)


if __name__ == '__main__':
    main()