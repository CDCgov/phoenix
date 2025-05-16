#!/usr/bin/env python
# v1.0 (11/25/2022)
#
# Adds filepaths for scaffolds stored in PHoeNIx scaffolds folder
# Auto mapping of filepath for user
# ##Usage: >python scaffolds_samplesheet.py <FILENAME>
# Written by Maria Diaz (lex0@cdc.gov)

import os
import argparse
import pandas as pd


def parseArgs(args=None):
    parser = argparse.ArgumentParser(
        description="Script to generate a complete PhoeNix scaffolds samplesheet"
    )
    parser.add_argument("file", nargs=argparse.REMAINDER)
    return parser.parse_args()


cwd = os.getcwd()
scaffLoc = "scaffolds_files/"
suffA = ".scaffolds.fa.gz"


def formatFilesSamplesheet(partialCsv):

    df = pd.read_csv(partialCsv[0], names=["sample", "scaffolds_file"])

    df["scaffolds_file"] = scaffLoc + df["sample"] + suffA

    df = df.to_csv(cwd + "/samplesheet.csv", index=False, header=True)


args = parseArgs()
formatFilesSamplesheet(args.file)
