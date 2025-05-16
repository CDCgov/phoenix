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
# Created by Jill Hagey (qpk9@cdc.gov)
#
#  Function to print out help blurb

version=2.0 # (11/15/2023) Changed to signify adoption of CLIA minded versioning. This version is equivalent to previous version 1.0 (07/03/2022)

show_help () {
    echo "Usage is ./beforeSpades.sh -d path_to_output -n report.tsv"
    echo "required: -d = path to specific sorted database with statistics related to entries from NCBI"
    echo "required: -n = sample name of file"
    echo "required: -k = kraken trimmed best hit summary"
    echo "required: -s = synopsis file"
    echo "optional: -V = show version"
}

# Parse command line options
options_found=0
while getopts ":h?d:n:k:s:cV" option; do
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
        c)
            echo "Option -c triggered"
            cdc_extended_qc="true";;
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
    echo "beforespades.sh: ${version}"
    exit
fi

genus=$(grep -R 'G:' $k2_bh_summary | cut -d ' ' -f3 | tr -d '\n')
gpercent=$(grep -R 'G:' $k2_bh_summary | cut -d ' ' -f2 | tr -d '\n')
species=$(grep -R 's:' $k2_bh_summary | cut -d ' ' -f3 | tr -d '\n')
spercent=$(grep -R 's:' $k2_bh_summary | cut -d ' ' -f2 | tr -d '\n')
name=$(echo "${genus}(${gpercent}%) ${species}(${spercent}%)")
species_col=$(echo "${genus} ${species}")


#get the number of warnings in the synopsis file
warning_count=$(grep ": WARNING  :" $synopsis | wc -l)

if [[ "${cdc_extended_qc}" == "true" ]]; then
    #for cdc_phoenix or cdc_scaffolds entry
    echo "ID    Auto_QC_Outcome Warning_Count   Estimated_Coverage  Genome_Length   Assembly_Ratio_(STDev)  #_of_Scaffolds_>500bp   GC_%    Species Taxa_Confidence Taxa_Coverage   Taxa_Source Kraken2_Trimd   Kraken2_Weighted    MLST_Scheme_1   MLST_1  MLST_Scheme_2   MLST_2  GAMMA_Beta_Lactam_Resistance_Genes  GAMMA_Other_AR_Genes    AMRFinder_Point_Mutations   Hypervirulence_Genes    Plasmid_Incompatibility_Replicons   Auto_QC_Failure_Reason" > ${sample_name}_summaryline_failure.tsv
    #file contents
    echo "${sample_name}    FAIL    ${warning_count}    Unknown Unknown Unknown Unknown Unknown Unknown Unknown ${species_col}  ${spercent}% Reads_assigned Unknown kraken2_trimmed ${name} Unknown Unknown Unknown Unknown Unknown Unknown Unknown Unknown Unknown Unknown SPAdes_Failure" | tr -d '\n' >> ${sample_name}_summaryline_failure.tsv
else
    #for phoenix or scaffolds entry
    #header
    echo "ID    Auto_QC_Outcome Warning_Count   Estimated_Coverage  Genome_Length   Assembly_Ratio_(STDev)  #_of_Scaffolds_>500bp   GC_%    Species Taxa_Confidence Taxa_Coverage   Taxa_Source Kraken2_Trimd   Kraken2_Weighted    MLST_Scheme_1   MLST_1  MLST_Scheme_2   MLST_2  GAMMA_Beta_Lactam_Resistance_Genes  GAMMA_Other_AR_Genes    AMRFinder_Point_Mutations   Hypervirulence_Genes    Plasmid_Incompatibility_Replicons   Auto_QC_Failure_Reason" > ${sample_name}_summaryline_failure.tsv
    #file contents
    echo "${sample_name}    FAIL    ${warning_count}    Unknown Unknown Unknown Unknown Unknown ${species_col}  ${spercent}% Reads_assigned Unknown kraken2_trimmed ${name} Unknown Unknown Unknown Unknown Unknown Unknown Unknown Unknown Unknown Unknown SPAdes_Failure" | tr -d '\n' >> ${sample_name}_summaryline_failure.tsv
fi

cp ${sample_name}_summaryline_failure.tsv ${output_path}/${sample_name}/
# copy the synopsis file
cp ${sample_name}.synopsis ${output_path}/${sample_name}
