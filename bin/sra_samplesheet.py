#!/usr/bin/env python
# v1.0 (06/28/2022)
#
# Adds filepaths for SRA FastQs automatically downloaded by Phoenix
# Auto mapping of filepath for user 
# ##Usage: >python sra_samplesheet.py -d <directory where fastq files are found>
# Written by Maria Diaz (lex0@cdc.gov) and Jill Hagey (qpk9@cdc.gov)

import argparse
import pandas as pd
import glob

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a complete PhoeNix paired reads samplesheet')
    parser.add_argument('-d', '--directory', required=True, dest='directory', help='Will create a samplesheet for all samples in the directory.')
    return parser.parse_args()

def write_samplesheet(directory, metadata_df):
    suff_R1 = "_R1_001.fastq.gz"
    suff_R2 = "_R2_001.fastq.gz"
    if directory[-1] != "/": # if directory doesn't have trailing / add one
        directory = directory + "/"
    with open("samplesheet.csv", "w") as samplesheet:
        samplesheet.write('sample' + ',' + 'fastq_1' + ',' + 'fastq_2' + ',' + 'sra_number\n')
        for index, row in metadata_df.iterrows(): #loop over each row
            sample = row['SampleName']
            sra_number = row['Run']
            samplesheet.write(sample + "," + directory + sample + "/raw_fastqs/" + sample + suff_R1 + "," + directory + sample + "/raw_fastqs/" + sample + suff_R2 + "," + sra_number + '\n')

def get_metadata():
    metadata_df = pd.DataFrame()
    for filename in glob.glob('*_sra_metadata.csv'):
        df = pd.read_csv(filename, header=0, dtype='str')
        df_less = df[['Run','SampleName']]
        metadata_df = pd.concat([metadata_df, df_less], axis=0, ignore_index=True)
    return metadata_df

def main():
    args = parseArgs()
    metadata_df = get_metadata()
    write_samplesheet(args.directory, metadata_df)

if __name__ == '__main__':
    main()
