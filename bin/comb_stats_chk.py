#!/usr/bin/env python3
import pandas as pd
import argparse

## Output check for messages indicating read pairs
## that do not match
## Written by Maria Diaz
## v.1.0.0
def parseArgs(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--raw_read', required=True)
    return parser.parse_args()

## Compare the GET_RAW_STATS module combined reads output
def reads_compare(file):
    approved = "PASS"
    failure = "FAILED"
    aggr_read_stats = pd.read_csv(file, sep="\t")
    
    # Confirm number of R1 reads are the same as R2 reads
    if aggr_read_stats["R1[reads]"].equals(aggr_read_stats["R2[reads]"]) is True:
        outcome = approved
    else:
        outcome = failure

    prefix = file.split("_raw")[0]

    #write to .txt file
    filename = prefix + "_result.txt"
    with open(filename, "w") as tmp:
        tmp.write(outcome)
        tmp.close()
    
def main():
    args = parseArgs()
    reads_compare(args.raw_read)

if __name__ == '__main__':
    main() 

    
