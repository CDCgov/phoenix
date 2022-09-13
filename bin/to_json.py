#!/usr/bin/env python
# v1.0 (09/13/2022)
#
# Creates JSON objects of CSV or TSV files created during Phoenix analysis
# ##Usage: >python to_json.py <FILENAME>
# Written by Maria Diaz (lex0@cdc.gov)
#

import argparse
import pandas as pd

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to convert Phoenix Excel output to a JSON object.')
    parser.add_argument('file')
    return parser.parse_args()

def jsonConverter(filePath):
    
    if filePath.endswith(".csv") or filePath.endswith(".txt"):#convert CSV files
        df = pd.read_csv(filePath, header=0)
        df.to_json(filePath[:-4] + ".json")
    else: #convert tab delim files
        df = pd.read_csv(filePath, sep="\t")
        df.to_json(filePath[:-4]  + ".json")

args = parseArgs()
jsonConverter(args.file)
