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
# V1.0 (07/03/2022)
#
# Created by Jill Hagey (qpk9@cdc.gov)
#
#  Function to print out help blurb
show_help () {
	echo "Usage is ./beforeSpades.sh -d path_to_output -n report.tsv"
	echo "required: -d = path to specific sorted database with statistics related to entries from NCBI"
	echo "required: -n = sample name of file"
	echo "required: -k = kraken trimmed best hit summary"
}

# Parse command line options
options_found=0
while getopts ":h?d:n:k:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			output_path=${OPTARG};;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
		k)
			echo "Option -k triggered, argument = ${OPTARG}"
			k2_bh_summary=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

genus=$(grep -R 'G:' $k2_bh_summary | cut -d ' ' -f3 | tr -d '\n')
gpercent=$(grep -R 'G:' $k2_bh_summary | cut -d ' ' -f2 | tr -d '\n')
species=$(grep -R 's:' $k2_bh_summary | cut -d ' ' -f3 | tr -d '\n')
spercent=$(grep -R 's:' $k2_bh_summary | cut -d ' ' -f2 | tr -d '\n')
name=$(echo "${genus}(${gpercent}%) ${species}(${spercent}%)")

echo "${sample_name}\tFAIL\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\t${name}\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tSPAdes_Failure" | tr -d '\n' > ${sample_name}_summaryline_failure.tsv
cp ${sample_name}_summaryline_failure.tsv ${output_path}/${sample_name}/
# copy the synopsis file
cp ${sample_name}.synopsis ${output_path}/${sample_name}