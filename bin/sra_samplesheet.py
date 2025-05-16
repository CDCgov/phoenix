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


# Function to get the script version
def get_version():
    return "2.0.0"


def parseArgs(args=None):
    parser = argparse.ArgumentParser(
        description="Script to generate a complete PhoeNix paired reads samplesheet"
    )
    parser.add_argument(
        "-d",
        "--directory",
        required=True,
        dest="directory",
        help="Will create a samplesheet for all samples in the directory.",
    )
    parser.add_argument(
        "-s",
        "--use_srr",
        dest="use_srr",
        default=False,
        action="store_true",
        help="Will create a samplesheet for all samples in the directory.",
    )
    parser.add_argument(
        "--version", action="version", version=get_version()
    )  # Add an argument to display the version
    return parser.parse_args()


def write_samplesheet(directory, metadata_df, duplicates, use_srr):
    suff_R1 = "_R1_001.fastq.gz"
    suff_R2 = "_R2_001.fastq.gz"
    if directory[-1] != "/":  # if directory doesn't have trailing / add one
        directory = directory + "/"
    with open("sra_samplesheet.csv", "w") as samplesheet:
        samplesheet.write(
            "sample" + "," + "fastq_1" + "," + "fastq_2" + "," + "sra_number\n"
        )
        for index, row in metadata_df.iterrows():  # loop over each row
            sample = row["SampleName"]
            sample = sample.replace(" ", "_").replace("/", "_")
            sra_number = row["Run"]
            if duplicates == True or use_srr == True:
                samplesheet.write(
                    sra_number
                    + ","
                    + directory
                    + sra_number
                    + "/raw_fastqs/"
                    + sra_number
                    + suff_R1
                    + ","
                    + directory
                    + sra_number
                    + "/raw_fastqs/"
                    + sra_number
                    + suff_R2
                    + ","
                    + sample
                    + "\n"
                )
            else:
                samplesheet.write(
                    sample
                    + ","
                    + directory
                    + sample
                    + "/raw_fastqs/"
                    + sample
                    + suff_R1
                    + ","
                    + directory
                    + sample
                    + "/raw_fastqs/"
                    + sample
                    + suff_R2
                    + ","
                    + sra_number
                    + "\n"
                )


def get_metadata():
    metadata_df = pd.DataFrame()  # make empty dataframe to add too
    for filename in glob.glob("*_sra_metadata.csv"):
        # get sra number aka the run
        sra_number = filename.replace("_sra_metadata.csv", "")
        df = pd.read_csv(filename, header=0, dtype="str")
        # when there are duplicate samples names this can mean there is multiple runs that show up in the metadata, but we only want the one row
        df = df.loc[df["Run"] == sra_number]
        df_less = df[
            ["Run", "SampleName"]
        ]  # reduce dataframe to only the information we need
        metadata_df = pd.concat([metadata_df, df_less], axis=0, ignore_index=True)
    duplicates = metadata_df["SampleName"].duplicated().any()
    return metadata_df, duplicates


def main():
    args = parseArgs()
    metadata_df, duplicates = get_metadata()
    write_samplesheet(args.directory, metadata_df, duplicates, args.use_srr)


if __name__ == "__main__":
    main()
