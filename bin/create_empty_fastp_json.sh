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
# Created by Jill Hagey (qpk9@cdc.gov)
#
#  Function to print out help blurb

version=2.0 # (11/15/2023) Changed to signify adoption of CLIA minded versioning. This version is equivalent to previous version 1.0 (04/18/2023)

show_help () {
	echo "Usage is ./crate_empty_fastp_json.sh  -n sample_name [-V version]"
	echo "required: -n = sample name of isolate"
	echo "optional: -V = show version and exit"
	echo ""
	echo "version: ${version}"
}

# Parse command line options
options_found=0
while getopts ":h?n:V" option; do
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
	echo "create_empty_fastp_json.sh: ${version}"
	exit
fi

#literally just creating an empty json
echo -e '{\n\t"summary": {\n\t\t"after_filtering": {\n\t\t\t"total_reads":0,\n\t\t\t"total_bases":0,\n\t\t\t"q20_bases":0,\n\t\t\t"q30_bases":0,\n\t\t\t"q20_rate":0,\n\t\t\t"q30_rate":0,\n\t\t\t"read1_mean_length":0,\n\t\t\t"gc_content":0\n\t\t}\n\t}\n}' > ${sample_name}_singles.fastp.json