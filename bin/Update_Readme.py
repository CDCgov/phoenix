#!/usr/bin/env python3

from datetime import date
import os
import os.path
import re
from decimal import *
getcontext().prec = 4
import argparse
import yaml

## updates a sample's readme file when -entry UPDATE is run
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

def write_readme(old_gamma, old_mlst, old_amrfinder, new_ar_db, new_mlst_db, new_amrfinderplus_db, output, phoenix_ver, sample_directory, current_phx_version):
    Gamma_db_updated = old_gamma + " --> " + new_ar_db
    amrfinder_db_updated = old_amrfinder + " --> " + new_amrfinderplus_db
    MLST_db_updated = old_mlst + " --> " + new_mlst_db
    phx_versions = phoenix_ver + " --> " + current_phx_version
    # Convert date to string format
    date_string = date.today().strftime('%Y-%m-%d')
    # Check if the file exists
    if os.path.exists(sample_directory + "/"+ output):
        # File exists, append a new line
        with open(output, 'a') as f:
            f.write(date_string + "\t" + phoenix_ver + "\t" + MLST_db_updated + "\t" + Gamma_db_updated + "\t" + amrfinder_db_updated + '\n')
    else:
        # File doesn't exist, write header and new line
        with open(output, 'w') as f:
            f.write('Date\tPhoenix_Version\tMLST_DB\tGAMMA_DB\tAMRFINDERPLUS_DB\n')
            f.write(date_string + "\t" + phx_versions + "\t" + MLST_db_updated + "\t" + Gamma_db_updated + "\t" + amrfinder_db_updated + '\n')

def main():
    args = parseArgs()
    old_gamma, old_mlst, old_amrfinder, phoenix_ver = get_old_database_IDs(args.pipeline_info)
    new_mlst_db, new_amrfinderplus_db = get_new_database_IDs(args.mlst_db, args.amrfinder_db)
    write_readme(old_gamma, old_mlst, old_amrfinder, args.ar_db, new_mlst_db, new_amrfinderplus_db, args.output, phoenix_ver, args.sample_directory, args.current_phx_version)


if __name__ == '__main__':
    main()