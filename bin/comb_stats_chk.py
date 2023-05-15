
import pandas as pd
import argparse

## Output check for messages indicating read pairs
## that do not match
## Written by Maria Diaz
## v.1.0.0

def parseArgs(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('path', nargs=argparse.REMAINDER)
    return parser.parse_args(args)

## Compare the GET_RAW_STATS module combined reads output
def reads_compare(file):
    approved = "PASS"
    failure = "FAILED"
    aggr_read_stats = pd.read_csv(file, sep="\t")
    
    # Confirm number of R1 reads are the same as R2 reads
    if aggr_read_stats["R1[reads]"].equals(aggr_read_stats["R2[reads]"]) is True:
        return approved
    else:
        return failure

args = parseArgs()
reads_compare(args.path)
    
