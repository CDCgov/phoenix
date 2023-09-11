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
    failure = "FAILED: Read pairs are NOT the same!" #essentially should never show up
    aggr_read_stats = pd.read_csv(file, sep="\t")
    prefix = file.split("_raw")[0]
    
    # Confirm number of R1 reads are the same as R2 reads
    if aggr_read_stats["R1[reads]"].equals(aggr_read_stats["R2[reads]"]) is True:
        outcome = approved
    else:
        outcome = failure
        data = {'ID': [prefix],'Auto_QC_Outcome': ['FAIL'],'Warning_Count': ['Unknown'],'Estimated_Coverage': ['Unknown'],'Genome_Length': ['Unknown'],
        'Assembly_Ratio_(STDev)': ['Unknown'],'#_of_Scaffolds_>500bp': ['Unknown'],'GC_%,Species,Taxa_Confidence': ['Unknown'],
        'Taxa_Coverage': ['Unknown'],'Taxa_Source': ['Unknown'],'Kraken2_Trimd': ['Unknown'],'Kraken2_Weighted': ['Unknown'],'MLST_Scheme_1': ['Unknown'],'MLST_1': ['Unknown'],
        'MLST_Scheme_2': ['Unknown'],'MLST_2': ['Unknown'],'GAMMA_Beta_Lactam_Resistance_Genes': ['Unknown'],'GAMMA_Other_AR_Genes': ['Unknown'],
        'AMRFinder_Point_Mutations': ['Unknown'],'Hypervirulence_Genes': ['Unknown'],'Plasmid_Incompatibility_Replicons': ['Unknown'],
        'Auto_QC_Failure_Reason': [failure]}
        df = pd.DataFrame(data)
        df.to_csv(prefix + '_summaryline_failure.tsv', sep="\t",index=False)

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
