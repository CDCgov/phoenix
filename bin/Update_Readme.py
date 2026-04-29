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
    parser.add_argument('--pf_db', dest="pf_db",required=True, help='New PF-Replicons database')
    parser.add_argument('--old_gamma', dest="old_gamma",required=True, help='Old <sample_id>_ResGANNCBI_<date>_srst2.gamma file.')
    parser.add_argument('--new_gamma', dest="new_gamma",required=True, help='New <sample_id>_ResGANNCBI_<date>_srst2.gamma file.')
    parser.add_argument('--old_pf', dest="old_pf",required=True, help='Old <sample_id>_PF-Replicons_<date>.gamma file.')
    parser.add_argument('--new_pf', dest="new_pf",required=True, help='New <sample_id>_PF-Replicons_<date>.gamma file.')
    parser.add_argument('--old_ncbi', dest="old_ncbi",required=True, help='Old <sample_id>_<date>_all_genes.tsv file.')
    parser.add_argument('--new_ncbi', dest="new_ncbi",required=True, help='New <sample_id>_<date>_all_genes.tsv file.')
    parser.add_argument('--old_tax', dest="old_tax",required=False, help='Old <sample_id>.tax file.')
    parser.add_argument('--new_tax', dest="new_tax",required=False, help='New <sample_id>_updater_log.tax file.')
    parser.add_argument('-p', '--pipeline_info', dest="pipeline_info", default=True, help='software_versions.yml file') # Need this for when you call -entry CDC_PHOENIX or CDC_SCAFFOLDS, but spades fails
    parser.add_argument('-o', '--out', dest="output", required=True, help='output file name')
    parser.add_argument('-v', dest="current_phx_version", required=True, help='current phx version.')
    parser.add_argument('--old_software_version_file', dest="old_software_version_file", required=False, help='Intermediate software versions file to show what the version of updater run was')
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
        else:
            old_gamma_ar_db = 'Missing'
        if 'GAMMA_PF' in yml_data:
            gamma_pf_data = yml_data['GAMMA_PF']
            # Extract values into variables
            old_gamma_pf_db = gamma_pf_data.get('Database', '')
            #gamma_ver = gamma_ar_data.get('gamma', '')
            #gamma_container = gamma_ar_data.get('gamma_container', '')
        else:
            old_gamma_pf_db = 'Missing'
        if 'MLST' in yml_data:
            mlst_data = yml_data['MLST']
            # Extract values into variables
            #mlst_ver = mlst_data.get('mlst', '')
            old_mlst_db = mlst_data.get('mlst_db', '')
            #mlst_container = mlst_data.get('mlst_container', '')
        else:
            old_mlst_db = 'Missing'
        if 'AMRFINDERPLUS_RUN' in yml_data:
            amrfinder_data = yml_data['AMRFINDERPLUS_RUN']
            # Extract values into variables
            #amrfinderplus_ver = amrfinder_data.get('amrfinderplus', '')
            old_amrfinderplus_db = amrfinder_data.get('amrfinderplus_db_version', '')
            #amrfinderplus_container = amrfinder_data.get('amrfinderplus_container', '')
        else:
            old_amrfinderplus_db = 'Missing'
        if 'Workflow' in yml_data:
            phoenix_data = yml_data['Workflow']
            # Extract values into variables
            phoenix_ver = phoenix_data.get('cdcgov/phoenix', '')
        else:
            phoenix_ver = 'Missing'
    return old_gamma_ar_db, old_mlst_db, old_amrfinderplus_db, phoenix_ver, old_gamma_pf_db

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
    old_df = pd.read_csv(old_file, sep='\t')
    new_df = pd.read_csv(new_file, sep='\t')

    rename_map = {
        'Gene_symbol': 'Element_symbol', 'Protein_identifier': 'Protein_id',
        'Element_type': 'Type', 'Sequence_name': 'Element_name',
        'HMM_id': 'HMM_accession',
        '%_Coverage_of_reference_sequence': '%_Coverage_of_reference',
        '%_Identity_to_reference_sequence': '%_Identity_to_reference',
        'Element_subtype': 'Subtype',
        'Accession_of_closest_sequence': 'Closest_reference_accession',
        'Name_of_closest_sequence': 'Closest_reference_name'
    }

    if file_type.lower() == "gamma":
        gene_column = "Gene"
        extract_gene_name = extract_gamma_gene_name
    elif file_type.lower() == "ncbi":
        old_df.rename(columns=rename_map, inplace=True)
        new_df.rename(columns=rename_map, inplace=True)  # Fix: rename both
        gene_column = "Element_symbol"
        extract_gene_name = lambda x: x
    else:
        raise ValueError(f"Unsupported file_type: {file_type}. Must be 'gamma' or 'ncbi'")

    # Build dicts keyed by gene name for O(1) lookup
    old_genes = {extract_gene_name(row[gene_column]): row.drop(gene_column).to_dict()
                 for _, row in old_df.iterrows()}
    new_genes = {extract_gene_name(row[gene_column]): row.drop(gene_column).to_dict()
                 for _, row in new_df.iterrows()}

    added_or_changed = []
    dropped = []

    for gene, old_data in old_genes.items():
        if gene not in new_genes:
            print(f"Dropped: {gene}")
            dropped.append(gene)
        elif old_data != new_genes[gene]:
            print(f"Changed: {gene}")
            for col in old_data:
                if old_data[col] != new_genes[gene].get(col):
                    print(f"  {col}: {old_data[col]} -> {new_genes[gene][col]}")
            added_or_changed.append(gene)

    for gene in new_genes:
        if gene not in old_genes:
            print(f"Added: {gene}")
            added_or_changed.append(gene)

    print(f"\nAdded/Changed ({len(added_or_changed)}): {added_or_changed}")
    print(f"Dropped ({len(dropped)}): {dropped}")
    return added_or_changed, dropped

def write_readme(old_gamma, old_mlst, old_amrfinder, new_ar_db, new_mlst_db, new_amrfinderplus_db, output, phoenix_ver, sample_directory, current_phx_version, added_gamma_ar_genes, dropped_gamma_ar_genes, added_ncbi_ar_genes, dropped_ncbi_ar_genes, new_tax, old_tax, old_pf, added_pf_genes, dropped_pf_genes, pf_db):
    Gamma_db_updated = old_gamma + " --> " + new_ar_db.strip()
    gamma_pf_db_updated = old_pf + " --> " + pf_db.strip()
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

    header = 'Date\tPhoenix_Version\tTaxa\tMLST_DB\tGAMMA_DB\tAMRFinderPlus_DB\tGAMMA_PF_DB\tAdded/Changed_GAMMA_AR_genes\tDropped_GAMMA_AR_genes\tAdded/Changed_NCBI_AR_genes\tDropped_NCBI_AR_genes\tAdded/Changed_PF_genes\tDropped_PF_genes\n'
    canonical_cols = header.strip().split('\t')

    new_line = date_string + "\t" + phx_versions + "\t" + tax_updated + "\t" + MLST_db_updated + "\t" + Gamma_db_updated + "\t" + amrfinder_db_updated + '\t' + gamma_pf_db_updated + '\t' + ','.join(added_gamma_ar_genes) + '\t' + ','.join(dropped_gamma_ar_genes) + '\t' + ','.join(added_ncbi_ar_genes) + '\t' + ','.join(dropped_ncbi_ar_genes) + '\t' + ','.join(added_pf_genes) + '\t' + ','.join(dropped_pf_genes) + '\n'

    # Build a single-row DataFrame for the new line explicitly by column name
    new_row = pd.DataFrame([{
        'Date':                          date_string,
        'Phoenix_Version':               phx_versions,
        'Taxa':                          tax_updated,
        'MLST_DB':                       MLST_db_updated,
        'GAMMA_DB':                      Gamma_db_updated,
        'AMRFinderPlus_DB':              amrfinder_db_updated,
        'GAMMA_PF_DB':                   gamma_pf_db_updated,
        'Added/Changed_GAMMA_AR_genes':          ','.join(added_gamma_ar_genes),
        'Dropped_GAMMA_AR_genes':','.join(dropped_gamma_ar_genes),
        'Added/Changed_NCBI_AR_genes':           ','.join(added_ncbi_ar_genes),
        'Dropped_NCBI_AR_genes': ','.join(dropped_ncbi_ar_genes),
        'Added/Changed_PF_genes':                ','.join(added_pf_genes),
        'Dropped_PF_genes':      ','.join(dropped_pf_genes),
    }], dtype=str)

    if os.path.exists(sample_directory + "/" + output):
        # File exists — read, normalize columns, reindex, append new row, write back
        df = pd.read_csv(output, sep='\t', dtype=str, on_bad_lines='warn', engine='python')

        # Normalize column names to handle any case differences
        col_map = {col: match for col in df.columns for match in canonical_cols if col.lower() == match.lower()}
        df = df.rename(columns=col_map)

        # Reindex to canonical columns, adding any missing ones as empty
        df = df.reindex(columns=canonical_cols)

        # Append new row and write back
        df = pd.concat([df, new_row], ignore_index=True)
        df.to_csv(output, sep='\t', index=False)
    else:
        # File doesn't exist — write header and new line directly
        with open(output, 'w') as f:
            f.write(header)
            f.write(new_line)

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
    if args.old_software_version_file:
        old_gamma, old_mlst, old_amrfinder, phoenix_ver, old_pf = get_old_database_IDs(args.old_software_version_file)
    elif args.pipeline_info:
        old_gamma, old_mlst, old_amrfinder, phoenix_ver, old_pf = get_old_database_IDs(args.pipeline_info)
    else:
        print("No software_versions.yml file provided, cannot extract old database versions. Please provide a software_versions.yml file using the -p argument or an intermediate software version file using --old_software_version_file.")
    new_mlst_db, new_amrfinderplus_db = get_new_database_IDs(args.mlst_db, args.amrfinder_db)
    added_gamma_ar_genes, dropped_gamma_ar_genes = compare_ar_files(args.old_gamma, args.new_gamma, file_type="gamma")
    added_ncbi_ar_genes, dropped_ncbi_ar_genes = compare_ar_files(args.old_ncbi, args.new_ncbi, file_type="ncbi")
    added_pf_genes, dropped_pf_genes = compare_ar_files(args.old_pf, args.new_pf, file_type="gamma")
    new_tax, old_tax = compare_tax_files(args.old_tax, args.new_tax)
    write_readme(old_gamma, old_mlst, old_amrfinder, args.ar_db, new_mlst_db, new_amrfinderplus_db, args.output, phoenix_ver, args.sample_directory, args.current_phx_version, added_gamma_ar_genes, dropped_gamma_ar_genes, added_ncbi_ar_genes, dropped_ncbi_ar_genes, new_tax, old_tax, old_pf, added_pf_genes, dropped_pf_genes, args.pf_db)


if __name__ == '__main__':
    main()
