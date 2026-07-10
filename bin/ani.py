#!/usr/bin/env python3
import csv
import os
import subprocess
from Bio import Entrez, SeqIO
import argparse
import re
import tempfile

# Replace 'Your_Email@example.com' with your email address
Entrez.email = 'Your_Email@example.com'

# Function to get the script version
def get_version():
    return "0.1.0"

# Function to fetch the FASTA sequence given an accession number or file path
def fetch_fasta_sequence(accession_number_or_path):
    try:
        # Check if the input is a local file path
        if os.path.isfile(accession_number_or_path):
            with open(accession_number_or_path, 'r') as file:
                record = SeqIO.read(file, 'fasta')
            return record
        else:
            handle = Entrez.efetch(db="nucleotide", id=accession_number_or_path, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return record
    except Exception as e:
        print(f"Error fetching sequence for {accession_number_or_path}: {e}")
        return None

# Function to extract contig sequences from an assembly.fasta file
def extract_contig_sequences(assembly_file, contig_ids, output_directory):
    for record in SeqIO.parse(assembly_file, "fasta"):
        contig_id = record.id
        if contig_id in contig_ids:
            output_file_path = os.path.join(output_directory, f'{contig_id}.fasta')
            with open(output_file_path, 'w') as output_file:
                SeqIO.write(record, output_file, 'fasta')
            print(f"Contig {contig_id} saved to {output_file_path}")

# Function to calculate fastANI between two sequences
def calculate_fastani(reference_sequence, query_sequence, output_name=None):
    ref_file_name = None
    query_file_name = None
    try:
        # Use the original file names for reference and query sequences
        ref_file_name = reference_sequence.id + '_reference.fasta'
        query_file_name = query_sequence.id + '_query.fasta'

        with open(ref_file_name, 'w') as ref_file:
            ref_file.write(f'>{reference_sequence.id}\n{str(reference_sequence.seq)}\n')

        with open(query_file_name, 'w') as query_file:
            query_file.write(f'>{query_sequence.id}\n{str(query_sequence.seq)}\n')

        # Use subprocess to call the fastANI tool
        output_option = f'-o {output_name}' if output_name else ''
        output_txt_name = f'{output_name}.txt'  # Add .txt extension to the output file
        command = f'fastANI -q {query_file_name} -r {ref_file_name} --visualize {output_option}'

        result = subprocess.run(command, shell=True, capture_output=True, text=True)

        # Debugging: Print the fastANI output and errors
        #print("fastANI Output:")
        #print(result.stdout)
        #print("fastANI Errors:")
        #print(result.stderr)

        # Check for errors in stderr
        #if result.returncode != 0:
            #print(f"Error calculating fastANI: {result.stderr}")
            #return None
        #else:
            # Extract the fastANI values from the output using regular expressions
            #ani_match = re.search(r'ANI\s*:\s*([\d.]+)%', result.stdout)

            #if ani_match:
                #ani_value = float(ani_match.group(1))
                #return ani_value
            #else:
                #print(f"No ANI value found in the output")
                #return None
    except Exception as e:
        print(f"Error calculating fastANI: {e}")
        return None
    finally:
        # Remove temporary files
        os.remove(ref_file_name)
        os.remove(query_file_name)
        if output_name:
            os.rename(output_name, output_txt_name)  # Rename output file with .txt extension


def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description='Combined script for fetching sequences, extracting contigs, and calculating fastANI.')
    parser.add_argument('-i', '--input', help='Path to the TSV file', required=True)
    parser.add_argument('-a', '--assembly', help='Path to the assembly.fasta file', required=True)
    parser.add_argument('--version', action='version', version=get_version()) # Add an argument to display the version
    

    # Parse command-line arguments
    args = parser.parse_args()

    # Set the path to the TSV file
    tsv_file_path = args.input

    # Set the working directory as the base path for both accession and contig folders
    working_directory = os.getcwd()
    accession_folder = os.path.join(working_directory, 'accession')
    contig_folder = os.path.join(working_directory, 'contig')

    # Create the output directories if they don't exist
    os.makedirs(accession_folder, exist_ok=True)
    os.makedirs(contig_folder, exist_ok=True)

    # Read the TSV file
    with open(tsv_file_path, 'r', newline='') as file:
        reader = csv.DictReader(file, delimiter='\t')

        # Iterate through each row in the TSV file
        for row in reader:
            accession_number = row['Accession Number']
            contig_id = row['Contig']

            # Fetch sequences for each accession number
            accession_sequence = fetch_fasta_sequence(accession_number)
            if accession_sequence:
                output_file_path = os.path.join(accession_folder, f'{accession_number}.fasta')
                with open(output_file_path, 'w') as output_file:
                    SeqIO.write(accession_sequence, output_file, 'fasta')
                print(f"Sequence for {accession_number} saved to {output_file_path}")

            # Extract contig sequences from the assembly file
            contig_ids = [contig_id]
            extract_contig_sequences(args.assembly, contig_ids, contig_folder)

            # Construct the paths to the accession and contig FASTA files
            accession_fasta_file = os.path.join(accession_folder, f'{accession_number}.fasta')
            contig_fasta_file = os.path.join(contig_folder, f'{contig_id}.fasta')

            # Fetch sequences from the FASTA files
            accession_sequence = fetch_fasta_sequence(accession_fasta_file)
            contig_sequence = fetch_fasta_sequence(contig_fasta_file)

            if accession_sequence and contig_sequence:
                # Generate the output name based on accession number and contig ID
                output_name = f'{accession_number}_vs_{contig_id}'

                # Calculate fastANI
                result = calculate_fastani(accession_sequence, contig_sequence, output_name)

                #if result is not None:
                    #ani_value, visualization_file_name = result
                    #print(f"fastANI between {accession_number} and {contig_id}: {ani_value:.4f}%")
                    #if visualization_file_name:
                        #print(f"Visualization saved to: {visualization_file_name}")
                #else:
                    #print(f"Error calculating fastANI for {accession_number} and {contig_id}")

if __name__ == "__main__":
    main()
