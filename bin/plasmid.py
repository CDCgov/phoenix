#!/usr/bin/env python3
import argparse
import subprocess
import os
from Bio import SeqIO

# Function to get the script version
def get_version():
    return "0.1.0"

def read_config_data(filename):
    """Reads configuration data from a tab-delimited file"""
    with open(filename) as config_file:
        data = []
        for line in config_file:
            fields = line.strip().split("\t")
            if fields:
                data.append(
                    {
                        "accession_number": fields[0],
                        "plasmid_name": fields[1],
                        "genome_accession": fields[2],
                    }
                )
    return data


def parse_config(config_file):
    """Parses the configuration file"""
    data = read_config_data(config_file)
    return data


def parse_blast_output(blast_output_file, config_data):
    """Parses BLAST output and returns a dictionary of contigs with their top 3 hits based on length and identity"""
    top_3_hits = {}
    with open(blast_output_file, "r") as f:
        for line in f:
            data = line.strip().split("\t")
            contig = data[0]
            accession_number = data[1]
            identity = float(data[2])
            length = int(data[3])

            # Find plasmid name from config data
            plasmid_name = None
            for entry in config_data:
                if entry["accession_number"] == accession_number:
                    plasmid_name = entry["plasmid_name"]
                    break

            # Append data to top_3_hits dictionary
            if plasmid_name:
                if contig not in top_3_hits:
                    top_3_hits[contig] = []
                top_3_hits[contig].append(
                    (accession_number, plasmid_name, identity, length)
                )

        # Sort hits by descending order of length and identity
        for contig, hits in top_3_hits.items():
            sorted_hits = sorted(hits, key=lambda x: (x[3], x[2]), reverse=True)[:3]
            top_3_hits[contig] = sorted_hits

        return top_3_hits


def get_contig_sizes(assembly_file):
    """Returns a dictionary mapping contig IDs to their sizes"""
    contig_sizes = {}
    for record in SeqIO.parse(assembly_file, "fasta"):
        contig_sizes[record.id] = len(record.seq)
    return contig_sizes


def run_blast(
    genome_file, database_file, reference_seq=None, output_file="blast_output.txt"
):
    """Runs BLAST command and retrieves matching plasmid name and highest identity"""
    with open(output_file, "w") as blast_output_file:
        blast_process = subprocess.Popen(
            ["blastn", "-query", genome_file, "-db", database_file, "-outfmt", "6"],
            stdout=blast_output_file,
        )
        blast_process.communicate()  # Wait for process to finish

    return output_file

# def run_blast(
#     genome_file, database_file, reference_seq=None, output_file="blast_output.txt"
# ):
#     """Runs BLAST command and retrieves matching plasmid name and highest identity"""
#     with open(output_file, "w") as blast_output_file:
#         if genome_file.endswith(".gz"):
#             # gunzip -c genome_file | blastn -query - -db database_file -outfmt 6
#             gunzip_proc =subprocess.Popen(
#                 ["gunzip", "-c", genome_file],
#                 stdout=subprocess.PIPE,
#             )
#             blast_proc = subprocess.Popen(
#                 ["blastn", "-query", "-", "-db", database_file, "-outfmt", "6"],
#                 stdin=gunzip_proc.stdout,
#                 stdout=blast_output_file,
#             )
#             gunzip_proc.stdout.close()
#             blast_proc.communicate()
#             gunzip_proc.wait()

#             if gunzip_proc.returncode != 0 or blast_proc.returncode != 0:
#                 raise subprocess.CalledProcessError(
#                     blast_proc.returncode or gunzip_proc.returncode,
#                     "gunzip|blastn pipeline",
#                 )
#         else:
#             blast_proc = subprocess.Popen(
#                 ["blastn", "-query", genome_file, "-db", database_file, "-outfmt", "6"],
#                 stdout=blast_output_file,
#             )
#             blast_proc.communicate()
#             if blast_proc.returncode != 0:
#                 raise subprocess.CalledProcessError(blast_proc.returncode, "blastn")

#     return output_file

def extract_contig_sequences(assembly_file, contig_ids, output_directory):
    """Extracts specified contig sequences from an assembly file"""
    os.makedirs(output_directory, exist_ok=True)

    for record in SeqIO.parse(assembly_file, "fasta"):
        if record.id in contig_ids:
            output_file_path = os.path.join(output_directory, record.id + ".fasta")
            with open(output_file_path, "w") as output_file:
                SeqIO.write(record, output_file, "fasta")
            print(f"Contig {record.id} saved to {output_file_path}")


def main():
    """Main script execution"""
    parser = argparse.ArgumentParser(
        description="Get highest identity for each unique contig from BLAST output"
    )
    parser.add_argument("-i", "--input", help="Input genome file path", required=True)
    parser.add_argument(
        "-db", "--database", help="Database FASTA file path", required=True
    )
    parser.add_argument("-c", "--config", help="Configuration file path", required=True)
    parser.add_argument("-r", "--reference", help="Reference sequence (optional)")
    parser.add_argument(
        "-o",
        "--output-file",
        help="BLAST output file name (default: blast_output.txt)",
        default="blast_output.txt",
    )
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    
    args = parser.parse_args()

    parsed_config = parse_config(args.config)

    # Run BLAST and retrieve the output file
    output_file = run_blast(
        args.input, args.database, args.reference, output_file=args.output_file
    )

    # Parse the BLAST output to get contigs and their hits
    parsed_data = parse_blast_output(output_file, parsed_config)

    # Get contig sizes from the genome assembly file
    contig_sizes = get_contig_sizes(args.input)

    # Extract contigs from assembly file if provided
    contigs_to_extract = list(parsed_data.keys())
    extract_contig_sequences(args.input, contigs_to_extract, "contigs")

    # Construct the formatted output with classification for top 3 hits
    output_string = (
        "Contig\tContig Size\tBest Hit Plasmid Size\tAccession Number\tPlasmid Name\t"
        "Identity\tLength\tContig_Plasmid_Ratio\tConfidence\tChromosome/Plasmid\n"
    )

    for contig, hits in parsed_data.items():
        contig_size = contig_sizes.get(contig, "N/A")
        # Classify the contig based on size
        classification = "chromosomal" if contig_size > 600000 else "plasmid"

        for hit in hits:  # Only iterating over the top 3 hits (already sorted)
            accession_number, plasmid_name, identity, length = hit
            plasmid_size = length  # Length of the plasmid from the best hit
            # Calculate the ratio of contig length to plasmid length
            ratio = contig_size / plasmid_size if plasmid_size else "N/A"
            confidence = "high" if 0.8 <= ratio <= 1.2 else "low"

            # Append all the information including classification
            output_string += (
                f"{contig}\t{contig_size}\t{plasmid_size}\t{accession_number}\t{plasmid_name}\t"
                f"{identity}\t{length}\t{ratio:.2f}\t{confidence}\t{classification}\n"
            )

    # Write the output to a file
    with open(args.output_file, "w") as formatted_output_file:
        formatted_output_file.write(output_string)


if __name__ == "__main__":
    main()
