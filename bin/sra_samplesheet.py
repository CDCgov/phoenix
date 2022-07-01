#!/usr/bin/env python
# v1.0 (06/28/2022)
#
# Adds filepaths for SRA FastQs automatically downloaded by Phoenix
# Auto mapping of filepath for user 
# ##Usage: >python sra_samplesheet.py <FILENAME>
# Written by Maria Diaz (lex0@cdc.gov)

import os 
import argparse
import pandas as pd

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a complete PhoeNix paired reads samplesheet')
    parser.add_argument('file', nargs=argparse.REMAINDER)
    return parser.parse_args()

cwd = os.getcwd()
fastqLoc = "results/fastq_files/"
suffA = "_R1_001.fastq.gz"
suffB = "_R2_001.fastq.gz"

def formatFilesSamplesheet(partialCsv):
    
    df = pd.read_csv(partialCsv[0], names= ['sample','fastq_1','fastq_2'])
    
    df['fastq_1'] = fastqLoc + df['sample'] + suffA
    
    df['fastq_2'] = fastqLoc + df['sample'] + suffB

    df = df.to_csv(cwd + "/samplesheet.csv", index=False, header=True)

args = parseArgs()
formatFilesSamplesheet(args.file)