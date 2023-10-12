#!/usr/bin/env python3
import pandas as pd
import argparse
import os
from datetime import date

## Output check for messages indicating read pairs
## that do not match
## Written by Maria Diaz
## v.1.0.0
def parseArgs(args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--raw_read', dest="raw_read", required=True)
    parser.add_argument('-f', '--results_file', dest="results_file",required=True)
    parser.add_argument('-b', '--busco', dest="busco", default=False, action='store_true', help='Pass to make true for -entry cdc pipelines') # Need this for when you call -entry CDC_PHOENIX or CDC_SCAFFOLDS, but spades fails
    return parser.parse_args()

## Compare the GET_RAW_STATS module combined reads output
def reads_compare(file, filename, busco):
    approved = "PASSED: Read pairs are equal."
    failure = "FAILED: The number of reads in R1/R2 are NOT the same!" #essentially should never show up
    aggr_read_stats = pd.read_csv(file, sep="\t")
    prefix = file.split("_raw")[0]

    # Confirm number of R1 reads are the same as R2 reads
    if aggr_read_stats["R1[reads]"].equals(aggr_read_stats["R2[reads]"]) is True:
        outcome = approved
    else:
        outcome = failure
        raw_length_R1, raw_length_R2, raw_reads, raw_pairs, raw_Q30_R1_rounded, raw_Q30_R2_rounded = get_read_stats(aggr_read_stats)
        warning_count = write_synopsis(prefix, busco, raw_length_R1, raw_length_R2, raw_reads, raw_pairs, raw_Q30_R1_rounded, raw_Q30_R2_rounded)
        write_summary_line(prefix, busco, warning_count)
    #write to end of *_results.txt file
    filename = prefix + "_results_old.txt"
    with open(filename, "a") as tmp:
        tmp.write(outcome)
    os.rename(filename, prefix + "_results.txt")

def get_read_stats(aggr_read_stats):
    raw_length_R1 = str(aggr_read_stats["R1[bp]"].values[0])
    raw_length_R2 = str(aggr_read_stats["R2[bp]"].values[0])
    raw_reads = int(aggr_read_stats["Total_Sequenced_[reads]"].values[0])
    raw_pairs = str(raw_reads/2)
    raw_Q30_R1_rounded = round((float(aggr_read_stats["Q30_R1_[%]"].values[0])*100), 2)
    raw_Q30_R2_rounded = round((float(aggr_read_stats["Q30_R2_[%]"].values[0])*100), 2)
    return raw_length_R1, raw_length_R2, raw_reads, raw_pairs, raw_Q30_R1_rounded, raw_Q30_R2_rounded


def write_synopsis(sample_name, busco, raw_length_R1, raw_length_R2, raw_reads, raw_pairs, raw_Q30_R1_rounded, raw_Q30_R2_rounded):
    status="FAILED"
    warning_count=0
    #create synopsis file
    with open(sample_name +".synopsis", "a") as f:
        # Creates and prints header info for the sample being processed
        today = date.today()
        f.write("---------- Checking " + sample_name + " for successful completion on ----------\n")
        f.write("Summarized                    : SUCCESS : " + str(today) + "\n")
        #Write out QC counts as failures
        f.write("FASTQs                        : SUCCESS :  R1: " + raw_length_R1 + "bps R2: " + raw_length_R2 + "bps\n")
        if raw_reads <= 1000000 and raw_reads >= 1:
            f.write("RAW_READ_COUNTS               : WARNING : Low individual read count before trimming: " + str(raw_reads) + "(" + raw_pairs + " paired reads)\n")
            warning_count = warning_count + 1
        elif raw_reads <=0:
            f.write("RAW_READ_COUNTS               : FAILED  : No individual read count before trimming: " + str(raw_reads) + "(" + raw_pairs + " paired reads)\n")
        else:
            f.write("RAW_READ_COUNTS               : SUCCESS : " + str(raw_reads) + " individual reads found in sample (" + raw_pairs + " paired reads)\n")
        if raw_Q30_R1_rounded < 90:
            f.write("RAW_Q30_R1%                   : WARNING : Q30_R1% at " + str(raw_Q30_R1_rounded) + "% (Threshold is 90%)\n")
            warning_count = warning_count + 1
        else:
            f.write("RAW_Q30_R1%                   : SUCCESS : Q30_R1% at " + str(raw_Q30_R1_rounded) + "% (Threshold is 90%)\n")
        if raw_Q30_R2_rounded < 70:
            f.write("RAW_Q30_R1%                   : WARNING : Q30_R1% at " + str(raw_Q30_R2_rounded) + "% (Threshold is 70%)\n")
            warning_count = warning_count + 1
        else:
            f.write("RAW_Q30_R2%                   : SUCCESS : Q30_R2% at " + str(raw_Q30_R2_rounded) + "% (Threshold is 70%)\n")
        f.write("TRIMMED_R1                    : FAILED  : " + sample_name + "_trimmed_read_counts.txt not found. No trimming performed -- Unequal number of reads in R1/R2!\n")
        f.write("TRIMMED_R2                    : FAILED  : " + sample_name + "_trimmed_read_counts.txt not found. No trimming performed -- Unequal number of reads in R1/R2!\n")
        f.write("TRIMMED_FASTQs                : FAILED  : " + sample_name + "_trimmed_read_counts.txt not found. No trimming performed -- Unequal number of reads in R1/R2!\n")
        f.write("TRIMMED_READ_COUNTS           : FAILED  : " + sample_name + "_trimmed_read_counts.txt not found. No trimming performed -- Unequal number of reads in R1/R2!\n")
        f.write("TRIMMED_Q30_R1%               : FAILED  : " + sample_name + "_trimmed_read_counts.txt not found. No trimming performed -- Unequal number of reads in R1/R2!\n")
        f.write("TRIMMED_Q30_R2%               : FAILED  : " + sample_name + "_trimmed_read_counts.txt not found. No trimming performed -- Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_READS                 : FAILED  : " + sample_name + ".kraken2_trimd.report.txt not found -- Unequal number of reads in R1/R2!\n")
        f.write("KRONA_READS                   : FAILED  : kraken2 reads do not exist -- Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_CLASSIFY_READS        : FAILED  : There are no classified reads -- Unequal number of reads in R1/R2!\n")
        f.write("QC_COUNTS                     : FAILED  : Unequal number of reads in R1/R2!\n")
        f.write("Q30_STATS                     : FAILED  : Unequal number of reads in R1/R2!\n")
        f.write("BBDUK                         : FAILED  : Unequal number of reads in R1/R2!\n")
        f.write("TRIMMING                      : FAILED  : Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_READS                 : FAILED  : Unequal number of reads in R1/R2!\n")
        f.write("KRONA_READS                   : FAILED  : kraken2 reads do not exist -- Unequal number of reads in R1/R2!")
        f.write("ASSEMBLY                      : FAILED  : " + sample_name + ".scaffolds.fa.gz not found -- Unequal number of reads in R1/R2!\n")
        if busco == True:
            f.write("SRST2                         : FAILED  : " + sample_name + "__fullgenes__*.txt file does not exist. SRST2 was not run. -- Unequal number of reads in R1/R2!\n")
        f.write("SCAFFOLD_TRIM                 : FAILED  : " + sample_name + ".filtered.scaffolds.fa.gz not found -- Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_ASMBLD                : FAILED  : " + sample_name + ".kraken2_asmbld.report.txt not found -- Unequal number of reads in R1/R2!\n")
        f.write("KRONA_ASMBLD                  : FAILED  : kraken2 unweighted not performed -- Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_CLASSIFY_ASMBLD       : FAILED  : kraken2 assembly not performed -- Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_ASMBLD_CONTAM         : FAILED  : kraken2 assembly not performed -- Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_WEIGHTED              : FAILED  : kraken2 weighted not performed -- Unequal number of reads in R1/R2!\n")
        f.write("KRONA_WEIGHTED                : FAILED  : kraken2 weighted not performed -- Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_CLASSIFY_WEIGHTED     : FAILED  : kraken2 weighted not performed -- Unequal number of reads in R1/R2!\n")
        f.write("KRAKEN2_WEIGHTED_CONTAM       : FAILED  : kraken2 weighted not performed -- Unequal number of reads in R1/R2!\n")
        f.write("QUAST                         : FAILED  : " + sample_name + "_report.tsv does not exist. QUAST was not run, no assembly -- Unequal number of reads in R1/R2!\n")
        f.write("TAXA                          : FAILED  : None of the classifiers completed successfully-- Unequal number of reads in R1/R2!\n")
        f.write("ASSEMBLY_RATIO(SD)            : FAILED  : No Ratio File exists -- Unequal number of reads in R1/R2!\n")
        f.write("COVERAGE                      : FAILED  : No trimmed reads to review for coverage -- Unequal number of reads in R1/R2!\n")
        if busco == True:
            f.write("BUSCO                         : FAILED  : BUSCO was not performed -- Unequal number of reads in R1/R2!\n")
        f.write("FASTANI_REFSEQ                : FAILED  : FASTANI was not performed -- Unequal number of reads in R1/R2!\n")
        f.write("MLST                          : FAILED  : MLST was not performed -- Unequal number of reads in R1/R2!\n")
        f.write("GAMMA_AR                      : FAILED  : GAMMA_AR was not performed -- Unequal number of reads in R1/R2!\n")
        f.write("PLASMID_REPLICONS             : FAILED  : GAMMA was not performed -- Unequal number of reads in R1/R2!\n")
        f.write("HYPERVIRULENCE                : FAILED  : GAMMA was not performed -- Unequal number of reads in R1/R2!\n")
        f.write("Auto PASS/FAIL                : FAILED  : FASTQs couldn't be unzipped!\n")
        f.write("---------- " + sample_name + " completed as " + status + " ----------\n")
        f.write("Unequal number of reads in R1/R2. Please fix and rerun PHoeNIx.\n")
        f.write("PHoeNIx will only analyze FASTQ files that have an equal number of reads.\n")
        f.write("WARNINGS: out of line with what is expected and MAY cause problems downstream.\n")
        f.write("ALERT: something to note, does not mean it is a poor-quality assembly.")
    return warning_count

def write_summary_line(prefix, busco, warning_count):
    if busco == True:
        column_names = ['ID','Auto_QC_Outcome','Warning_Count','Estimated_Coverage','Genome_Length',
        'Assembly_Ratio_(STDev)','#_of_Scaffolds_>500bp','GC_%', 'BUSCO', 'BUSCO_DB', 'Species','Taxa_Confidence',
        'Taxa_Coverage','Taxa_Source','Kraken2_Trimd','Kraken2_Weighted','MLST_Scheme_1','MLST_1',
        'MLST_Scheme_2','MLST_2','GAMMA_Beta_Lactam_Resistance_Genes','GAMMA_Other_AR_Genes',
        'AMRFinder_Point_Mutations','Hypervirulence_Genes','Plasmid_Incompatibility_Replicons',
        'Auto_QC_Failure_Reason']
        data = [[prefix,'FAIL',warning_count,'Unknown','Unknown','Unknown','Unknown','Unknown','Unknown', 'Unknown','Unknown','Unknown',
        'Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','The number of reads in R1/R2 are NOT equal!']]
    else:
        column_names = ['ID','Auto_QC_Outcome','Warning_Count','Estimated_Coverage','Genome_Length',
        'Assembly_Ratio_(STDev)','#_of_Scaffolds_>500bp','GC_%', 'Species','Taxa_Confidence',
        'Taxa_Coverage','Taxa_Source','Kraken2_Trimd','Kraken2_Weighted','MLST_Scheme_1','MLST_1',
        'MLST_Scheme_2','MLST_2','GAMMA_Beta_Lactam_Resistance_Genes','GAMMA_Other_AR_Genes',
        'AMRFinder_Point_Mutations','Hypervirulence_Genes','Plasmid_Incompatibility_Replicons',
        'Auto_QC_Failure_Reason']
        data = [[prefix,'FAIL',warning_count,'Unknown','Unknown','Unknown','Unknown','Unknown','Unknown', 'Unknown'
        'Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','Unknown','The number of reads in R1/R2 are NOT equal!']]
    df = pd.DataFrame(data, columns=column_names)
    df.to_csv(prefix + '_summaryline_failure.tsv', sep="\t",index=False)

def main():
    args = parseArgs()
    reads_compare(args.raw_read, args.results_file, args.busco)

if __name__ == '__main__':
    main() 
