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
        column_names = ['ID','Auto_QC_Outcome','Warning_Count','Estimated_Coverage','Genome_Length',
        'Assembly_Ratio_(STDev)','#_of_Scaffolds_>500bp','GC_%,Species,Taxa_Confidence',
        'Taxa_Coverage','Taxa_Source','Kraken2_Trimd','Kraken2_Weighted','MLST_Scheme_1','MLST_1',
        'MLST_Scheme_2','MLST_2','GAMMA_Beta_Lactam_Resistance_Genes','GAMMA_Other_AR_Genes',
        'AMRFinder_Point_Mutations','Hypervirulence_Genes','Plasmid_Incompatibility_Replicons',
        'Auto_QC_Failure_Reason']
        data = [[prefix,'FAIL','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown',
        'Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown',
        'Unknown','Unknown','Unknown',failure]]
        df = pd.DataFrame(data, columns=column_names)
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
