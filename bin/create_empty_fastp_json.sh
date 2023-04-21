#!/bin/bash -l

#
# Description: Script to create empty json for single reads info to be used downstream
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
	echo "Usage is ./crate_empty_fastp_json.sh  -n sample_name"
	echo "required: -n = sample name of isolate"
	}

# Parse command line options
options_found=0
while getopts ":h?n:" option; do
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
echo -e '{\n\t"summary": {\n\t\t"after_filtering": {\n\t\t\t"total_reads":0,\n\t\t\t"total_bases":0,\n\t\t\t"q20_bases":0,\n\t\t\t"q30_bases":0,\n\t\t\t"q20_rate":0,\n\t\t\t"q30_rate":0,\n\t\t\t"read1_mean_length":0,\n\t\t\t"gc_content":0\n\t\t}\n\t}\n}' > ${sample_name}_singles.fastp.json