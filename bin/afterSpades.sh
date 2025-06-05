#!/bin/bash

#
# Description: Script to clean up the output of SPAdes
#
# Usage:
#
# Output location:
#
# Modules required: None
#
# V1.0 (07/03/2022)
#
# Created by Jill Hagey (qpk9@cdc.gov)

version=2.0 # (11/15/2023) Changed to signify adoption of CLIA minded versioning. This version is equivalent to previous version 1.0 (07/03/2022)

# Parse command line options
options_found=0
while getopts ":Vh" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
			show_help
			exit 0
			;;
		V)
			show_version="True";;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

if [[ "${show_version}" = "True" ]]; then
    echo "afterspades.sh: ${version}"
    exit
fi

log=$(find *.spades.log)
prefix=$(basename $log .spades.log)
if [ -f scaffolds.fasta ]; then
    mv scaffolds.fasta ${prefix}.scaffolds.fa
    gzip -n ${prefix}.scaffolds.fa
    spades_complete=scaffolds_created
    echo ,$spades_complete | tr -d "\n" >> ${prefix}_spades_outcome.csv
else
    spades_complete=no_scaffolds
    echo ,$spades_complete | tr -d "\n" >> ${prefix}_spades_outcome.csv
    mv ${prefix}_summary_old_3.txt ${prefix}_trimstats_summary.txt
fi
if [ -f contigs.fasta ]; then
    mv contigs.fasta ${prefix}.contigs.fa
    gzip -n ${prefix}.contigs.fa
    spades_complete=contigs_created
    echo ,$spades_complete | tr -d "\n" >> ${prefix}_spades_outcome.csv
else
    spades_complete=no_contigs
    echo ,$spades_complete | tr -d "\n" >> ${prefix}_spades_outcome.csv
    mv ${prefix}_summary_old_3.txt ${prefix}_trimstats_summary.txt
fi
if [ -f assembly_graph_with_scaffolds.gfa ]; then
    mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
    gzip -n ${prefix}.assembly.gfa
fi