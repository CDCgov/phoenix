#!/bin/bash -l
set +e
#
# Description: script to check for file integrity and log errors 
#
# Usage: ./fairy_proc.sh reads_file 
#
# v.1.0.0 (07/13/2023)
#
# Created by Maria Diaz (lex0@cdc.gov)
#

sfx=".fastq.gz"
fname="${1}"
prefix=${fname%"$sfx"}
gzip -t $1 2>> ${prefix}.txt

if grep "error" ${prefix}.txt || grep "error" ${prefix}.txt; then
	echo "FAILED CORRUPTION CHECK! CANNOT UNZIP FASTQ FILE. CHECK FASTQ FILE(S) FOR CORRUPTION!" > ${prefix}_results.txt
	#error warning for line_summary channel
	echo "ID	Auto_QC_Outcome	Warning_Count	Estimated_Coverage	Genome_Length	Assembly_Ratio_(STDev)	#_of_Scaffolds_>500bp	GC_%	Species	Taxa_Confidence	Taxa_Coverage	Taxa_Source	Kraken2_Trimd	Kraken2_Weighted	MLST_Scheme_1	MLST_1	MLST_Scheme_2	MLST_2	GAMMA_Beta_Lactam_Resistance_Genes	GAMMA_Other_AR_Genes	AMRFinder_Point_Mutations	Hypervirulence_Genes	Plasmid_Incompatibility_Replicons	Auto_QC_Failure_Reason" > ${prefix}_summaryline_failure.tsv
	#file contents
	echo "${prefix}	FAIL	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	kraken2_trimmed	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown	Unknown		CANNOT UNZIP FASTQ FILE. CHECK FASTQ FILE(S) FOR CORRUPTION!" | tr -d '\n' >> ${prefix}_summaryline_failure.tsv
else
	echo "PASS" > ${prefix}_results.txt
fi

#proceed to cumulative read counts if files aren't corrupt
if grep -Fx "PASS" ${prefix}_results.txt; then
	q30.py $1 > ${prefix}_stats.txt
fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
	python: \$(python --version | sed 's/Python //g')
END_VERSIONS