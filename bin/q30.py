#!/usr/bin/env python3

## orginal script from https://github.com/dayedepps/q30/blob/master/fastq.py
## modified to be python3 not 2 and added read count.
## by Jill Hagey qpk9@cdc.gov 4/3/2023

import os, sys

# disable cache usage in the Python so __pycache__ isn't formed. If you don't do this using 'nextflow run cdcgov/phoenix...' a second time will causes and error
sys.dont_write_bytecode = True  # needs to be before the import fastq
import fastq
import time
import argparse


# Function to get the script version
def get_version():
    return "2.0.0"


def parseArgs(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--version", action="version", version=get_version()
    )  # Add an argument to display the version
    parser.add_argument(
        "-i", "--input", dest="input", required=False, help="input fasta filename"
    )
    parser.add_argument("files", nargs=argparse.REMAINDER)
    return parser.parse_args()


def qual_stat(qstr):
    q20 = 0
    q30 = 0
    for q in qstr:
        qual = ord(chr(q)) - 33
        # qual = ord(q) - 33 #python2 version
        if qual >= 30:
            q30 += 1
            q20 += 1
        elif qual >= 20:
            q20 += 1
    return q20, q30


def stat(filename):
    reader = fastq.Reader(filename)
    total_read_count = 0
    total_base_count = 0
    q20_count = 0
    q30_count = 0
    while True:
        read = reader.nextRead()
        if read == None:
            break
        total_read_count = total_read_count + 1
        total_base_count += len(read[3])
        q20, q30 = qual_stat(read[3])
        q20_count += q20
        q30_count += q30

    print("total reads:", total_read_count)
    print("total bases:", total_base_count)
    print("q20 bases:", q20_count)
    print("q30 bases:", q30_count)
    print("q20 percents:", 100 * float(q20_count) / float(total_base_count))
    print("q30 percents:", 100 * float(q30_count) / float(total_base_count))


def main():
    args = parseArgs()
    stat(args.input)


if __name__ == "__main__":
    time1 = time.time()
    main()
    time2 = time.time()
    print("Time used: " + str(time2 - time1))
