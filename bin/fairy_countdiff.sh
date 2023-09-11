#!/bin/bash -l
#
# Description: script to create a custom .synopsis file for mismatched read counts
#
# Usage: ./fairy_countdiff.sh raw_read_counts_file 
#
# v.1.0.0 (08/03/2023)
#
# Created by Maria Diaz (lex0@cdc.gov)
#

sfx="_raw_read_counts.txt"
fname="${1}"
prefix=${fname%"$sfx"}
prefix=${prefix%%_*}
#create synopsis file
sample_name=${prefix}
# Creates and prints header info for the sample being processed
today=$(date)
echo "----------Checking ${sample_name} for successful completion on ----------"  > "${sample_name}.synopsis"

printf "%-30s: %-8s : %s\\n" "Summarized" "SUCCESS" "${today}"  >> "${sample_name}.synopsis"

#Write out QC counts as failures
raw_length_R1=$(tail -n1  "${fname}" | cut -d$'\t' -f3)
raw_length_R2=$(tail -n1  "${fname}" | cut -d$'\t' -f5)
printf "%-30s: %-8s : %s\\n" "FASTQs" "FAILED" "R1: ${raw_length_R1}bps R2: ${raw_length_R2}bps -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "RAW_READ_COUNTS" "FAILED" "Read pairs are NOT the same!"  >> "${sample_name}.synopsis"

printf "%-30s: %-8s : %s\\n" "RAW_Q30_R1%" "FAILED" "Cannot use these FASTQ files for analysis -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "RAW_Q30_R2%" "FAILED" "Cannot use these FASTQ files for analysis  -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
	
printf "%-30s: %-8s : %s\\n" "TRIMMED_R1" "FAILED" "No trimming performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "TRIMMED_R2" "FAILED" "No trimming performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"

printf "%-30s: %-8s : %s\\n" "TRIMMED_FASTQs" "FAILED" "Trimmed reads QC file does not exist -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "TRIMMED_READ_COUNTS" "FAILED" "Trimmed reads QC file does not exist -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
	
printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R1%" "FAILED" "No trimming performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R2%" "FAILED" "No trimming performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"

printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS" "FAILED" "${sample_name}.kraken2_trimd.report.txt not found -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"

printf "%-30s: %-8s : %s\\n" "KRONA_READS" "FAILED" "kraken2 reads do not exist -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
status="FAILED"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "FAILED" "There are no classified reads -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"

printf "%-30s: %-8s : %s\\n" "QC_COUNTS" "FAILED" "Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "Q30_STATS" "FAILED" "Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "BBDUK" "FAILED" "Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "TRIMMING" "FAILED" "Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS" "FAILED" "Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KROFAILED_READS" "FAILED" "Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
	
printf "%-30s: %-8s : %s\\n" "ASSEMBLY" "FAILED" "${sample_name}.scaffolds.fa.gz not found -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"

printf "%-30s: %-8s : %s\\n" "SRST2" "FAILED" "SPAdes not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "SCAFFOLD_TRIM" "FAILED" "${sample_name}.filtered.scaffolds.fa.gz not found -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD" "FAILED" "${sample_name}.kraken2_asmbld.report.txt not found -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRONA_ASMBLD" "FAILED" "kraken2 unweighted not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "FAILED" "kraken2 assembly not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD_CONTAM" "FAILED" "kraken2 assembly not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED" "FAILED" "kraken2 weighted not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRONA_WEIGHTED" "FAILED" "kraken2 weighted not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "FAILED" "kraken2 weighted not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED_CONTAM" "FAILED" "kraken2 weighted not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "QUAST" "FAILED" "QUAST not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "TAXA-${tax_source}" "FAILED" "No Taxa File found -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "ASSEMBLY_RATIO(SD)" "FAILED" "No Ratio File exists -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "COVERAGE" "FAILED" "No trimmed reads to review for coverage -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "BUSCO" "FAILED" "BUSCO was not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "FAILED" "FASTANI was not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "MLST" "FAILED" "MLST was not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "GAMMA_AR" "FAILED" "GAMMA_AR was not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "PLASMID_REPLICONS" "FAILED" "GAMMA was not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
printf "%-30s: %-8s : %s\\n" "HYPERVIRULENCE" "FAILED" "GAMMA was not performed -- Read pairs are NOT the same!"  >> "${sample_name}.synopsis"

printf "%-30s: %-8s : %s\\n" "Auto Pass/FAIL" "FAIL" "Read pairs are NOT the same!"  >> "${sample_name}.synopsis"

echo "---------- ${sample_name} completed as ${status} ----------"  >> "${sample_name}.synopsis"
echo "WARNINGS: Read pairs are NOT the same!"  >> "${sample_name}.synopsis"
echo "ALERT: PHoeNIx will only analyze FASTQ files that have equivalent read counts for R1 and R2."  >> "${sample_name}.synopsis"

