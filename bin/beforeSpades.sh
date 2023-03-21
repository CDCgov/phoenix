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
	echo "required: -s = synopsis file"
}

# Parse command line options
options_found=0
while getopts ":h?d:n:k:s:" option; do
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
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			synopsis=${OPTARG};;
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
species_col=$(echo "${genus} ${species}")


#get the number of warnings in the synopsis file
warning_count=$(grep ": WARNING  :" $synopsis | wc -l)

#header
echo "ID\tAuto_QC_Outcome\tWarning_Count\tEstimated_Coverage\tGenome_Length\tAssembly_Ratio_(STDev)\t#_of_Scaffolds_>500bp\tGC_%\tSpecies\tTaxa_Confidence\tTaxa_Source\tKraken2_Trimd\tKraken2_Weighted\tMLST_Scheme_1\tMLST_1\tMLST_Scheme_2\tMLST_2\tGAMMA_Beta_Lactam_Resistance_Genes\tGAMMA_Other_AR_Genes\tAMRFinder_Point_Mutations\tHypervirulence_Genes\tPlasmid_Incompatibility_Replicons\tAuto_QC_Failure_Reason" > ${sample_name}_summaryline_failure.tsv
#file contents
echo "${sample_name}\tFAIL\t${warning_count}\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\t${species_col}\tUnknown\tkraken2_trimmed\t${name}\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tUnknown\tSPAdes_Failure" | tr -d '\n' >> ${sample_name}_summaryline_failure.tsv
cp ${sample_name}_summaryline_failure.tsv ${output_path}/${sample_name}/
# copy the synopsis file
cp ${sample_name}.synopsis ${output_path}/${sample_name}