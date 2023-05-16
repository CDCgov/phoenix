#!/usr/bin/env python3
import pandas as pd

## Output check for messages indicating read pairs
## that do not match
## Written by Maria Diaz
## v.1.0.0

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
    filename = prefix + "_results.txt"
    tmp = open(filename, "w")
    tmp.write(outcome)
    tmp.close()
    


    
