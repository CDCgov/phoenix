#!/usr/bin/env python3
import argparse
import csv
import os
import subprocess

rscript_path="/opt/conda/bin/Rscript"

# Function to get the script version
def get_version():
    return "0.1.0"
    
def run_visualization_script(query_fasta, ref_fasta, visualization_file, rscript_location):
    command = f'{rscript_path} {rscript_location} {query_fasta} {ref_fasta} {visualization_file}'
    subprocess.run(command, shell=True)

def main(tsv_file, contig_folder, accession_folder, visual_files_folder, rscript_location):
    # Read the TSV file
    with open(tsv_file, 'r') as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter='\t')
        for row in reader:
            accession_number = row['Accession Number']
            contig_id = row['Contig']
            
            # Construct paths to the query and reference FASpdf( paste(fastANI_visual_file,".pdf",sep="") )TA files
            query_fasta = os.path.join(contig_folder, f'{contig_id}.fasta')
            ref_fasta = os.path.join(accession_folder, f'{accession_number}.fasta')
            
            # Construct path to the visualization file
            visualization_file = os.path.join(visual_files_folder, f'{accession_number}_vs_{contig_id}.visual')
            
            # Check if the visualization file exists
            if os.path.exists(visualization_file):
                # Run the visualization script
                run_visualization_script(query_fasta, ref_fasta, visualization_file, rscript_location)
            else:
                print(f"Visualization file not found for {accession_number} and {contig_id}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Script to parse a TSV file, find sample pairs, and run an R script for visualization.')
    parser.add_argument('tsv_file', help='Path to the TSV file')
    parser.add_argument('contig_folder', help='Path to the contig folder')
    parser.add_argument('accession_folder', help='Path to the accession folder')
    parser.add_argument('visual_files_folder', help='Path to the .visual files location')
    parser.add_argument('rscript_location', help='Path to the R script location')
    parser.add_argument('--version', action='version', version=get_version()) # Add an argument to display the version

    args = parser.parse_args()

    main(args.tsv_file, args.contig_folder, args.accession_folder, args.visual_files_folder, args.rscript_location)

