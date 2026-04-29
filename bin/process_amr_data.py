#!/usr/bin/env python3

import csv
import re
import os
import glob
import argparse

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

def clean_gene_name(gene):
    """Remove special characters from gene names"""
    # Remove asterisks, carets, and other special characters
    return re.sub(r'[*^#@!&$%+\-=<>{}[\]|\\/:;.,?~`]', '', gene)

def process_abritamr_file(abritamr_file):
    """Process abritamr.txt file to extract clean gene names"""
    with open(abritamr_file, 'r') as f:
        lines = f.readlines()
    
    # Skip header line and get the results line
    results_line = lines[1].strip()
    
    # Split by tab and get the gene columns (skip the first 'Isolate' column)
    columns = results_line.split('\t')[1:]
    
    # Extract all genes from all columns
    all_genes = []
    for column in columns:
        if column.strip():  # Skip empty columns
            genes = column.split(',')
            for gene in genes:
                gene = gene.strip()
                if gene:
                    clean_gene = clean_gene_name(gene)
                    if clean_gene:
                        all_genes.append(clean_gene)
    
    return list(set(all_genes))  # Remove duplicates

def process_amrfinder_file(amrfinder_file):
    """Process amrfinder.out file to extract coverage and identity data"""
    gene_data = {}
    
    with open(amrfinder_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            gene_symbol = row['Gene symbol']
            sequence_name = row['Name of closest sequence']
            coverage = row['% Coverage of reference sequence']
            identity = row['% Identity to reference sequence']
            
            # Clean the gene symbol
            clean_gene = clean_gene_name(gene_symbol)
            
            # Store data with the clean gene symbol as key
            gene_data[clean_gene] = {
                'Coverage': coverage,
                'Identity': identity
            }
            
            # Also try to extract gene name from sequence description for better matching
            if sequence_name and sequence_name != 'NA':
                # Look for patterns like "AAC(6')-Ib4" in the sequence name
                import re
                gene_matches = re.findall(r'[A-Za-z]+\([^)]*\)[-\w]*', sequence_name)
                for match in gene_matches:
                    clean_match = clean_gene_name(match)
                    if clean_match and clean_match not in gene_data:
                        gene_data[clean_match] = {
                            'Coverage': coverage,
                            'Identity': identity
                        }
    
    return gene_data

def create_summary(abritamr_file, amrfinder_file, output_file, sample_name):
    """Create summary file with gene, coverage, and identity data"""
    try:
        # Get clean gene names from abritamr file
        abritamr_genes = process_abritamr_file(abritamr_file)
        
        # Get coverage and identity data from amrfinder file
        amrfinder_data = process_amrfinder_file(amrfinder_file)
        
        # Create summary data
        summary_data = []
        
        for gene in abritamr_genes:
            found = False
            
            # First try exact match
            if gene in amrfinder_data:
                summary_data.append([
                    sample_name,
                    gene,
                    amrfinder_data[gene]['Coverage'],
                    amrfinder_data[gene]['Identity']
                ])
                found = True
            else:
                # Try fuzzy matching - look for partial matches
                for amr_gene in amrfinder_data.keys():
                    # Check if the gene names are similar (allowing for minor differences)
                    if (gene.lower().replace('(', '').replace(')', '').replace("'", '') == 
                        amr_gene.lower().replace('(', '').replace(')', '').replace("'", '')):
                        summary_data.append([
                            sample_name,
                            gene,
                            amrfinder_data[amr_gene]['Coverage'],
                            amrfinder_data[amr_gene]['Identity']
                        ])
                        found = True
                        break
            
            if not found:
                # Gene not found in amrfinder data
                summary_data.append([
                    sample_name,
                    gene,
                    'N/A',
                    'N/A'
                ])
        
        # Write to CSV file
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            # Write header
            writer.writerow(['Sample', 'Gene', 'Coverage', 'Identity'])
            # Write data
            writer.writerows(summary_data)
        
        print(f"  ✓ Summary file created: {output_file}")
        print(f"  ✓ Total genes processed: {len(summary_data)}")
        
        return summary_data
        
    except Exception as e:
        print(f"  ✗ Error processing {sample_name}: {str(e)}")
        return None

def find_sample_files():
    """Find sample files in AMRFinder folder"""
    samples = {}
    
    # Look for .abritamr.txt files
    abritamr_files = glob.glob("*.abritamr.txt")
    
    for abritamr_file in abritamr_files:
        # Extract sample name (everything before .abritamr.txt)
        sample_name = os.path.basename(abritamr_file).replace('.abritamr.txt', '')
        
        # Look for corresponding .amrfinder.out file
        amrfinder_file = f"{sample_name}.amrfinder.out"
        
        if os.path.exists(amrfinder_file):
            samples[sample_name] = {
                'abritamr': abritamr_file,
                'amrfinder': amrfinder_file
            }
    
    return samples

def process_multiple_folders():
    """Process all AMRFinder folders in the given"""
    total_processed = 0

    # Find sample files in this folder
    samples = find_sample_files()
        
    if not samples:
        print("  No matching sample files found.")
    print(f"  Found {len(samples)} sample(s): {', '.join(samples.keys())}")
        
    # Process each sample
    for sample_name, files in samples.items():
        print(f"\n  Processing sample: {sample_name}")
        # Create output filename
        output_file = os.path.join(f"{sample_name}_AR_details.csv")

        # Process the sample
        result = create_summary(
            files['abritamr'], 
            files['amrfinder'], 
            output_file,
            sample_name
        )
            
        if result is not None:
            total_processed += 1
    
    print(f"\n" + "="*50)
    print(f"Processing complete!")
    print(f"Total samples processed: {total_processed}")
    print(f"Output files created with pattern: [sample]_AR_details.csv")

def main():
    args = parseArgs()
    # Process all AMRFinder folders in the given path
    process_multiple_folders()

if __name__ == "__main__":
    main()

