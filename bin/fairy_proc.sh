#!/bin/bash -l
set +e
#
# Description: script to check for file integrity and log errors 
#
# Usage: ./fairy_proc.sh reads_file software_version
#
# v.1.0.0 (07/13/2023)
#
# Created by Maria Diaz (lex0@cdc.gov)
#

sfx=".fastq.gz"
fname="${1}"
prefix=${fname%"$sfx"}
gzip -t $1 2>> ${prefix}.txt


if grep -q -e "error" -e "unexpected" ${prefix}.txt ]; then
	prefix=${prefix%%_*}
	echo "FAILED CORRUPTION CHECK! CANNOT UNZIP FASTQ FILE. CHECK FASTQ FILE(S) FOR CORRUPTION!" >> ${prefix}_results.txt
	
	#error warning for line_summary channel
	echo "ID	Auto_QC_Outcome	Warning_Count	Estimated_Coverage	Genome_Length	Assembly_Ratio_(STDev)	#_of_Scaffolds_>500bp	GC_%	Species	Taxa_Confidence	Taxa_Coverage	Taxa_Source	Kraken2_Trimd	Kraken2_Weighted	MLST_Scheme_1	MLST_1	MLST_Scheme_2	MLST_2	GAMMA_Beta_Lactam_Resistance_Genes	GAMMA_Other_AR_Genes	AMRFinder_Point_Mutations	Hypervirulence_Genes	Plasmid_Incompatibility_Replicons	Auto_QC_Failure_Reason" > ${prefix}_summaryline_failure.tsv
	#file contents
	echo "${prefix}	FAIL	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	kraken2_trimmed	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown		CANNOT UNZIP FASTQ FILE. CHECK FASTQ FILE(S) FOR CORRUPTION!" | tr -d '\n' >> ${prefix}_summaryline_failure.tsv
	#create synopsis file
	sample_name=${prefix}
	# Creates and prints header info for the sample being processed
	today=$(date)
	echo "----------Checking ${sample_name} for successful completion on ----------"  > "${sample_name}.synopsis"

	printf "%-30s: %-8s : %s\\n" "Summarized" "SUCCESS" "${today}"  >> "${sample_name}.synopsis"

	#Write out QC counts as failures
    printf "%-30s: %-8s : %s\\n" "FASTQs" "FAILED" "${sample_name}_raw_read_counts reads QC file does not exist -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "RAW_READ_COUNTS" "FAILED" "${sample_name}_raw_read_counts reads QC file does not exist -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"

	printf "%-30s: %-8s : %s\\n" "RAW_Q30_R1%" "FAILED" "${sample_name}_raw_read_counts.txt not found -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "RAW_Q30_R2%" "FAILED" "${sample_name}_raw_read_counts.txt not found -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	
	printf "%-30s: %-8s : %s\\n" "TRIMMED_R1" "FAILED" "No trimming performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "TRIMMED_R2" "FAILED" "No trimming performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"

    printf "%-30s: %-8s : %s\\n" "TRIMMED_FASTQs" "FAILED" "Trimmed reads QC file does not exist -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "TRIMMED_READ_COUNTS" "FAILED" "Trimmed reads QC file does not exist -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	
	printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R1%" "FAILED" "No trimming performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
    printf "%-30s: %-8s : %s\\n" "TRIMMED_Q30_R2%" "FAILED" "No trimming performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"

	printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS" "FAILED" "${sample_name}.kraken2_trimd.report.txt not found -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"

	printf "%-30s: %-8s : %s\\n" "KRONA_READS" "FAILED" "kraken2 reads do not exist -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	status="FAILED"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_READS" "FAILED" "There are no classified reads -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"

	printf "%-30s: %-8s : %s\\n" "QC_COUNTS" "FAILED" "FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "Q30_STATS" "FAILED" "FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "BBDUK" "FAILED" "FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "TRIMMING" "FAILED" "FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_READS" "FAILED" "FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KROFAILED_READS" "FAILED" "FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	
	printf "%-30s: %-8s : %s\\n" "ASSEMBLY" "FAILED" "${sample_name}.scaffolds.fa.gz not found -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"

	printf "%-30s: %-8s : %s\\n" "SRST2" "FAILED" "SPAdes not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "SCAFFOLD_TRIM" "FAILED" "${sample_name}.filtered.scaffolds.fa.gz not found -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD" "FAILED" "${sample_name}.kraken2_asmbld.report.txt not found -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRONA_ASMBLD" "FAILED" "kraken2 unweighted not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_ASMBLD" "FAILED" "kraken2 assembly not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_ASMBLD_CONTAM" "FAILED" "kraken2 assembly not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED" "FAILED" "kraken2 weighted not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRONA_WEIGHTED" "FAILED" "kraken2 weighted not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_CLASSIFY_WEIGHTED" "FAILED" "kraken2 weighted not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "KRAKEN2_WEIGHTED_CONTAM" "FAILED" "kraken2 weighted not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "QUAST" "FAILED" "QUAST not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "TAXA-${tax_source}" "FAILED" "No Taxa File found -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "ASSEMBLY_RATIO(SD)" "FAILED" "No Ratio File exists -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "COVERAGE" "FAILED" "No trimmed reads to review for coverage -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "BUSCO" "FAILED" "BUSCO was not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "FASTANI_REFSEQ" "FAILED" "FASTANI was not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "MLST" "FAILED" "MLST was not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "GAMMA_AR" "FAILED" "GAMMA_AR was not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "PLASMID_REPLICONS" "FAILED" "GAMMA was not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	printf "%-30s: %-8s : %s\\n" "HYPERVIRULENCE" "FAILED" "GAMMA was not performed -- FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"

	printf "%-30s: %-8s : %s\\n" "Auto Pass/FAIL" "FAIL" "FASTQs couldn't be unzipped!"  >> "${sample_name}.synopsis"
	
	echo "---------- ${sample_name} completed as ${status} ----------"  >> "${sample_name}.synopsis"
	echo "WARNINGS: File corruption detected in one or more isolate FASTQs. Please check the FAIry result folder for details."  >> "${sample_name}.synopsis"
	echo "ALERT: PHoeNIx will only analyze FASTQ files that can be unzipped without error."  >> "${sample_name}.synopsis"
else
	echo "PASS" >> ${prefix%%_*}_results.txt
fi
