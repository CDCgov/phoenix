#!/bin/bash -l

#
# Description: Script to clean up the output of SPAdes
#
# Usage:
#
# Output location:
#
# Modules required: None
#
# V1.0 (04/18/2023)
#
# Created by Jill Hagey (qpk9@cdc.gov)
#
#  Function to print out help blurb
show_help () {
	echo "Usage is ./beforeSpades.sh -d path_to_output -n report.tsv"
	echo "required: -d = path to specific sorted database with statistics related to entries from NCBI"
	echo "required: -n = sample name of file"
	echo "required: -k = kraken trimmed best hit summary"
	echo "required: -s = synopsis file"
}

# Parse command line options
options_found=0
while getopts ":h?n" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

#literally just creating an empty json
echo '{\n	\"summary\": {\n		\"after_filtering\": {\n			\"total_reads\":0,\n			\"total_bases\":0,\n			\"q20_bases\":0,\n			\"q30_bases\":0,\n			\"q20_rate\":0,\n			\"q30_rate\":0,\n			\"read1_mean_length\":0,\n			\"gc_content\":0\n		}\n	}\n}' > ${sample_name}_singles.fastp.json