#!/bin/bash

#$ -o 	CENTAR_CONSOLIDATE.out
#$ -e 	CENTAR_CONSOLIDATE.err
#$ -N 	CENTAR_CONSOLIDATE
#$ -cwd
#$ -q short.q

#
# Description: 'More' formal way to extract the specific GOIs for Cdiff samples from a PHoeNIx source isolate
#
# Usage: ./Centar_consolidater.sh -i path_to_phx_isolate -output_file
#
# Output location: Varies on contents
#
# v1.0.0 (07/15/2024)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ../Centar_Consolidater.sh -t path_to_phx_isolate_GAMMA_tox_file -c path_to_clade_file -o output_file -y toxinotype_file -r Ribotype_file -p plasmid_file -s sample_name -a gamma_ar_file -n gammas_nt_file -x crosswalk_STRT_file"
}

version="1.0"

tox_input=""
clade_input=""
# Parse command line options
options_found=0
while getopts ":h?t:c:o:y:a:r:p:Vs:n:x:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
            show_help
            exit 0
            ;;
		t)
			echo "Option -t triggered, argument = ${OPTARG}"
			tox_input=${OPTARG};;
		c)
			echo "Option -c triggered, argument = ${OPTARG}"
			clade_input=${OPTARG};;
        o)
			echo "Option -o triggered, argument = ${OPTARG}"
			output=${OPTARG};;
        y)
			echo "Option -y triggered, argument = ${OPTARG}"
			ttype_file=${OPTARG};;
        a)
			echo "Option -a triggered, argument = ${OPTARG}"
			aa_mut_file=${OPTARG};;
        n)
			echo "Option -n triggered, argument = ${OPTARG}"
			nt_mut_file=${OPTARG};;
        r)
			echo "Option -r triggered, argument = ${OPTARG}"
			rt_file=${OPTARG};;
        s)
			echo "Option -s triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
        p)
			echo "Option -p triggered, argument = ${OPTARG}"
			plasmid_file=${OPTARG};;
        x)
			echo "Option -x triggered, argument = ${OPTARG}"
			xwalkRT_file=${OPTARG};;
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

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi


if [[ "${show_version}" = "True" ]]; then
	echo "${version}"
	exit
fi

# Little func to pull out elements of a GAMA output file. $1 is the entry you want and $2 is the GAMA line to parse
function extract_gama_elements {
    re='^[0-9]+$'
    if ! [[ $1 =~ $re ]] ; then
        echo "error: Input is not a number"; exit 1
    fi
    formatted_string=$(echo $2 | sed 's/__/!/g')
    new_string=$(echo "${formatted_string}" | cut -d'!' -f$1)
    echo "${new_string}"
}

# Little func to compare gamma output to what the mutation type is, e.g. AA vs NT and how well they match. $1 is the gamma line you want to check
function confirm_mutation_type {
    #re='^[0-9]+$'
    #if ! [[ $1 =~ $re ]] ; then
    #    echo "error: Input is not a number"; exit 1
    #fi
    formatted_string=$(echo $2 | sed 's/__/!/g')
    new_string=$(echo "${formatted_string}" | cut -d'!' -f$1)
    echo "${new_string}"
}

# Little function to clean up mutation strings...mostly of the trailing comma
function clean_mut {
    description=$(echo "${1}" | cut -d$'\t' -f6)
    if [[ ${description: -1} = "," ]]; then
        new_string="${description::-1}"
    else
        new_string="${description}"
    fi
    echo "${new_string}"
}


if [[ -f "${rt_file}" ]]; then
    rt="unset"
    prob="unset"
    type="unset"
    source="unset"
    note=""
    clusterx="unset"
    clustery="unset"

    ## Parse the RT file once ready and figure out how to make it expandable
     line=$(tail -n1 ${rt_file})
     ML_RT=""
    # rt=$(echo "${line}" | cut -d$'\t' -f2)
    # IFS=',' read -a rts <<< "${rt}"
    # rt_count=${#rts[@]}
    # if [ ${rt_count} -gt 0 ]; then
    #     new_rt=""
    #     for rib in "${rts[@]}"; do
    #         rt_length=${#rib}
    #         if [ $rt_length -eq 2 ] && [[ "${rib}" != "NA" ]]; then
    #             new_rt="0${rib}"
    #         elif [ $rt_length -eq 1 ]; then
    #             new_rt="00${rib}"
    #         else
    #             new_rt="${rib}"
    #         fi
    #     done
    #     if [[ ${new_rt:0:1} == "," ]] ; then
    #         new_rt="${new_rt:1}"
    #     fi
    # else
    #     new_rt="unset"
    # fi
    raw_prob=$(echo "${line}" | cut -d$'\t' -f5)
    temp_prob=$(echo "$raw_prob * 100" | bc)
    prob=$(printf "%2.2f" "${temp_prob}")
    echo "PROB-$raw_prob-$temp_prob-$prob"}
    clusterx=$(echo "${line}" | cut -d$'\t' -f2)
    clustery=$(echo "${line}" | cut -d$'\t' -f3)
    source=$(echo "${line}" | cut -d$'\t' -f6)
    rt=$(echo "${line}" | cut -d$'\t' -f4)
    note=$(echo "${line}" | cut -d$'\t' -f7)
    if [[ "${note}" = "" ]]; then
        ML_RT="${rt}\t${prob}\t${source}"
    else
        ML_RT="${rt}\t${prob}\t${source}\t${note}"
    fi
else
    echo "No Ribotype file"
    ML_RT="NO_RT_FILE"
fi

#// Removing until we add this feature in, no more printing placeholder
#if [[ -f "${plasmid_file}" ]]; then
#    ## Parse the plasmid file once ready and figure out how to make it expandable
#    echo "No plasmid info file yet"
#    plasmids="NO_PLASMID_FILE_YET"
#else
#    echo "No plasmid info file"
#    plasmids="NO_PLASMID_FILE_YET"
#fi

if [[ -f "${ttype_file}" ]]; then
    line=0
    
    while IFS= read -r var; do
        if [[ "${line}" -eq 0 ]]; then
            line=$(( line + 1 ))
            continue
        else
            if [[ "${var}" = "Toxinotype:"* ]]; then
                toxinotype=$(echo "${var}" | cut -d$'\t' -f2)
            else
                current_type=$(echo "${var}" | cut -d$'\t' -f2)
                current_subtype=$(echo "${var}" | cut -d$'\t' -f3)
               # echo "${current_type}, ${current_subtype}"
                if [[ "${current_type}" = "Toxin-A" ]]; then
                    if [[ "${subtype_A}" = "" ]]; then
                        subtype_A="${current_subtype}"
                    else
                        subtype_A="${subtype_A},${current_subtype}"
                    fi
                elif [[ "${current_type}" = "Toxin-B" ]]; then
                    if [[ "${subtype_B}" = "" ]]; then
                        subtype_B="${current_subtype}"
                    else
                        subtype_B="${subtype_B},${current_subtype}"
                    fi
                else
                    echo "Unknown type experienced in Centar_consolidator.sh when looking up toxinotypes"
                fi
            fi
        fi
        line=$(( line + 1 ))
    done < "${ttype_file}"

    if [[ "${subtype_A}" = "" ]]; then
        subtype_A="Unknown"
    fi
    if [[ "${subtype_B}" = "" ]]; then
        subtype_B="Unknown"
    fi
else
    echo "No toxinotype input file"
    subtype_A="NO_ttype_file"
    subtype_B="NO_ttype_file"
    ttype="NO_ttype_file"
fi

if [[ -f "${tox_input}" ]]; then
    # Check and format for tcdA count
    tcdA_count=$(grep -c 'tcdA' "${tox_input}")
    if [[ "${tcdA_count}" -eq 0 ]]; then
        tcdA=0
        tcdA_set="0\tNOT_FOUND\tNOT_FOUND"
    elif [[ "${tcdA_count}" -eq 1 ]]; then
        tcdA_line=$(grep 'tcdA' "${tox_input}")
        tcdA_matchtype=$(echo "${tcdA_line}" | cut -d$'\t' -f5)
        tcdA_allele=$(extract_gama_elements 3 "${tcdA_line}")
        tcdA_raw_IDA=$(echo "${tcdA_line}" | cut -d$'\t' -f10)
        tcdA_raw_IDN=$(echo "${tcdA_line}" | cut -d$'\t' -f11)
        tcdA_raw_length=$(echo "${tcdA_line}" | cut -d$'\t' -f12)
        tcdA_IDA=$(echo "$tcdA_raw_IDA * 100" | bc | cut -d'.' -f1)
        tcdA_IDN=$(echo "$tcdA_raw_IDN * 100" | bc | cut -d'.' -f1)
        tcdA_length=$(echo "$tcdA_raw_length * 100" | bc | cut -d'.' -f1)
        if [[ "${tcdA_matchtype}" = *"Trunc"* ]]; then
            tcdA="TRUNC:CODON[${tcdA_IDA}AA]"
        elif (( $(echo "$tcdA_length*100 < 90" | bc -l) )); then
            tcdA="TRUNC:LENGTH<90[${tcdA_length}COV]]"
        else
            tcdA=1
        fi
        tcdA_set="${tcdA}\t${subtype_A}\t[${tcdA_IDN}NT|${tcdA_IDA}AA|${tcdA_length}COV]"
    else
        tcdA=1
        tcdA_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"tcdA"* ]]; then
                tcdA_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        tcdA_count_array=${#tcdA_arr[@]}

        if [[ ${tcdA_count} != ${tcdA_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdA_stats=""
        tcdA_allele=""
        while [[ ${count} -lt ${tcdA_count_array} ]]; do
            tcdA_line="${tcdA_arr[count]}"
            tcdA_allele_current=$(extract_gama_elements 3 "${tcdA_line}")
            tcdA_raw_IDA=$(echo "${tcdA_line}" | cut -d$'\t' -f10)
            tcdA_raw_IDN=$(echo "${tcdA_line}" | cut -d$'\t' -f11)
            tcdA_raw_length=$(echo "${tcdA_line}" | cut -d$'\t' -f12)
            tcdA_IDA=$(echo "$tcdA_raw_IDA * 100" | bc | cut -d'.' -f1)
            tcdA_IDN=$(echo "$tcdA_raw_IDN * 100" | bc | cut -d'.' -f1)
            tcdA_length=$(echo "$tcdA_raw_length * 100" | bc | cut -d'.' -f1)
            tcdA_matchtype=$(echo "${tcdA_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${tcdA} = "" ]]; then
                if [[ "${tcdA_matchtype}" = *"Trunc"* ]]; then
                    tcdA="TRUNC:CODON[${tcdA_IDA}AA]"
                elif (( $(echo "$tcdA_length*100 < 90" | bc -l) )); then
                    tcdA="TRUNC:LENGTH<90[${tcdA_length}COV]"
                else
                    tcdA=1
                fi
            else
                if [[ "${tcdA_matchtype}" = *"Trunc"* ]]; then
                    tcdA="${tcdA}-TRUNC:CODON[${tcdA_IDA}AA]"
                elif (( $(echo "$tcdA_length*100 < 90" | bc -l) )); then
                    tcdA="${tcdA}-TRUNC:LENGTH<90[${tcdB_length}COV]]"
                else
                    tcdA="${tcdA}-1"
                fi
            fi
            if [[ "${tcdA_stats}" == "" ]]; then
                tcdA_stats="[${tcdA_IDN}NT|${tcdA_IDA}AA|${tcdA_length}COV]"
            else
                tcdA_stats="${tcdA_stats}-[${tcdA_IDN}NT|${tcdA_IDA}AA|${tcdA_length}COV]"
            fi
            count=$(( count + 1 ))
        done
        tcdA_set="${tcdA}\t${subtype_A}\t${tcdA_stats}"
    fi

    # Check and format for tcdB count
    tcdB_count=$(grep -c 'tcdB' "${tox_input}")
    if [[ "${tcdB_count}" -eq 0 ]]; then
        tcdB=0
        tcdB_set="0\tNOT_FOUND\tNOT_FOUND"
    elif [[ "${tcdB_count}" -eq 1 ]]; then
        tcdB_line=$(grep 'tcdB' "${tox_input}")
        tcdB_matchtype=$(echo "${tcdB_line}" | cut -d$'\t' -f5)
        tcdB_allele=$(extract_gama_elements 3 "${tcdB_line}")
        tcdB_raw_IDA=$(echo "${tcdB_line}" | cut -d$'\t' -f10)
        tcdB_raw_IDN=$(echo "${tcdB_line}" | cut -d$'\t' -f11)
        tcdB_raw_length=$(echo "${tcdB_line}" | cut -d$'\t' -f12)
        tcdB_IDA=$(echo "$tcdB_raw_IDA * 100" | bc | cut -d'.' -f1)
        tcdB_IDN=$(echo "$tcdB_raw_IDN * 100" | bc | cut -d'.' -f1)
        tcdB_length=$(echo "$tcdB_raw_length * 100" | bc | cut -d'.' -f1)
        if [[ "${tcdB_matchtype}" = *"Trunc"* ]]; then
            tcdB="TRUNC:CODON[${tcdB_IDA}AA]"
        elif (( $(echo "$tcdB_length*100 < 90" | bc -l) )); then
            tcdB="TRUNC:LENGTH<90[${tcdB_length}COV]]"
        else
            tcdB=1
        fi
        tcdB_set="${tcdB}\t${subtype_B}\t[${tcdB_IDN}NT|${tcdB_IDA}AA|${tcdB_length}COV]"
    else
        tcdB=1
        tcdB_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"tcdB"* ]]; then
                tcdB_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        tcdB_count_array=${#tcdB_arr[@]}

        if [[ ${tcdB_count} != ${tcdB_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdB_stats=""
        tcdB_allele=""
        while [[ ${count} -lt ${tcdB_count_array} ]]; do
            tcdB_line="${tcdB_arr[count]}"
            tcdB_allele_current=$(extract_gama_elements 3 "${tcdB_line}")
            tcdB_raw_IDA=$(echo "${tcdB_line}" | cut -d$'\t' -f10)
            tcdB_raw_IDN=$(echo "${tcdB_line}" | cut -d$'\t' -f11)
            tcdB_raw_length=$(echo "${tcdB_line}" | cut -d$'\t' -f12)
            tcdB_IDA=$(echo "$tcdB_raw_IDA * 100" | bc | cut -d'.' -f1)
            tcdB_IDN=$(echo "$tcdB_raw_IDN * 100" | bc | cut -d'.' -f1)
            tcdB_length=$(echo "$tcdB_raw_length * 100" | bc | cut -d'.' -f1)
            tcdB_matchtype=$(echo "${tcdB_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${tcdB} = "" ]]; then
                if [[ "${tcdB_matchtype}" = *"Trunc"* ]]; then
                    tcdB="TRUNC:CODON[${tcdB_IDA}AA]"
                elif (( $(echo "$tcdB_length*100 < 90" | bc -l) )); then
                    tcdB="TRUNC:LENGTH<90[${tcdB_length}COV]"
                else
                    tcdB=1
                fi
            else
                if [[ "${tcdB_matchtype}" = *"Trunc"* ]]; then
                    tcdB="${tcdB}-TRUNC:CODON[${tcdB_IDA}AA]"
                elif (( $(echo "$tcdB_length*100 < 90" | bc -l) )); then
                    tcdB="${tcdB}-TRUNC:LENGTH<90[${tcdB_length}COV]"
                else
                    tcdB="${tcdB}-1"
                fi
            fi
            if [[ "${tcdB_stats}" == "" ]]; then
                tcdB_stats="[${tcdB_IDN}NT|${tcdB_IDA}AA|${tcdB_length}COV]"
            else
                echo "extra"
                tcdB_stats="${tcdB_stats}-[${tcdB_IDN}NT|${tcdB_IDA}AA|${tcdB_length}COV]"
            fi
            #echo "${count}-${tcdB}\t${tcdB_count_array}\t${tcdB_stats}"
            count=$(( count + 1 ))
        done
        tcdB_set="${tcdB}\t${subtype_B}\t${tcdB_stats}"
    fi

    # Check and format for tcdC count
    tcdC_count=$(grep -c 'tcdC' "${tox_input}")
    if [[ "${tcdC_count}" -eq 0 ]]; then
        tcdC=0
        tcdC_set="0\tNOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    elif [[ "${tcdC_count}" -eq 1 ]]; then
        tcdC_line=$(grep 'tcdC' "${tox_input}")
        tcdC_matchtype=$(echo "${tcdC_line}" | cut -d$'\t' -f5)
        tcdC_allele=$(extract_gama_elements 3 "${tcdC_line}")
        tcdC_muts=$(clean_mut "${tcdC_line}")
        tcdC_raw_IDA=$(echo "${tcdC_line}" | cut -d$'\t' -f10)
        tcdC_raw_IDN=$(echo "${tcdC_line}" | cut -d$'\t' -f11)
        tcdC_raw_length=$(echo "${tcdC_line}" | cut -d$'\t' -f12)
        tcdC_IDA=$(echo "$tcdC_raw_IDA * 100" | bc | cut -d'.' -f1)
        tcdC_IDN=$(echo "$tcdC_raw_IDN * 100" | bc | cut -d'.' -f1)
        tcdC_length=$(echo "$tcdC_raw_length * 100" | bc | cut -d'.' -f1)
        if [[ "${tcdC_matchtype}" = *"Trunc"* ]]; then
            if [[ "${tcdC_IDA}" = "100" ]]; then
                tcdC=1
            else
                tcdC="TRUNC:CODON[${tcdC_IDA}AA]"
            fi
        elif (( $(echo "$tcdC_length*100 < 90" | bc -l) )); then
            tcdC="TRUNC:LENGTH<90[${tcdC_length}COV]"
        else
            tcdC=1
        fi
        if [[ "${tcdC_muts}" = "No coding mutations" ]]; then
            tcdC_muts='-'
        fi
        
        tcdC_set="${tcdC}\t${tcdC_allele}\t${tcdC_muts}\t[${tcdC_IDN}NT|${tcdC_IDA}AA|${tcdC_length}COV]"
    else
        tcdC=""
        tcdC_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"tcdC"* ]]; then
                tcdC_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        tcdC_count_array=${#tcdC_arr[@]}

        if [[ ${tcdC_count} != ${tcdC_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdC_stats=""
        tcdC_allele=""
        tcdC_muts=""
        while [[ ${count} -lt ${tcdC_count_array} ]]; do
            tcdC_line="${tcdC_arr[count]}"
            tcdC_allele_current=$(extract_gama_elements 3 "${tcdC_line}")
            tcdC_raw_IDA=$(echo "${tcdC_line}" | cut -d$'\t' -f10)
            tcdC_raw_IDN=$(echo "${tcdC_line}" | cut -d$'\t' -f11)
            tcdC_raw_length=$(echo "${tcdC_line}" | cut -d$'\t' -f12)
            tcdC_muts_current=$(clean_mut "${tcdC_line}")
            tcdC_IDA=$(echo "$tcdC_raw_IDA * 100" | bc | cut -d'.' -f1)
            tcdC_IDN=$(echo "$tcdC_raw_IDN * 100" | bc | cut -d'.' -f1)
            tcdC_length=$(echo "$tcdC_raw_length * 100" | bc | cut -d'.' -f1)
            tcdC_matchtype=$(echo "${tcdC_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${tcdC} = "" ]]; then
                if [[ "${tcdC_matchtype}" = *"Trunc"* ]]; then
                    if [[ "${tcdC_IDA}" = "100" ]]; then
                        tcdC=1
                    else
                        tcdC="TRUNC:CODON[${tcdC_IDA}AA]"
                    fi
                elif (( $(echo "$tcdC_length*100 < 90" | bc -l) )); then
                    tcdC="TRUNC:LENGTH<90[${tcdC_length}COV]"
                else
                    tcdC=1
                fi
            else
                if [[ "${tcdC_matchtype}" = *"Trunc"* ]]; then
                    if [[ "${tcdC_IDA}" = "100" ]]; then
                        tcdC="${tcdC}-1"
                    else
                        tcdC="${tcdC}-TRUNC:CODON[${tcdC_IDA}AA]"
                    fi
                elif (( $(echo "$tcdC_length*100 < 90" | bc -l) )); then
                    tcdC="${tcdC}-TRUNC:LENGTH<90[${tcdC_length}COV]"
                else
                    tcdC="${tcdC}-1"
                fi
            fi
            if [[ "${tcdC_allele}" = "" ]]; then
                tcdC_allele=${tcdC_allele_current}
            else
                tcdC_allele="${tcdC_allele}-${tcdC_allele_current}"
            fi
            if [[ "${tcdC_stats}" = "" ]]; then
                tcdC_stats="[${tcdC_IDN}NT|${tcdC_IDA}AA|${tcdC_length}COV]"
            else
                tcdC_stats="${tcdC_stats}-[${tcdC_IDN}NT|${tcdC_IDA}AA|${tcdC_length}COV]"
            fi
            if [[ "${tcdC_muts_current}" = "No coding mutations" ]]; then
                tcdC_muts_current='-'
            fi
            if [[ "${tcdC_muts}" = "" ]]; then
                tcdC_muts="${tcdC_muts_current}"
            else
                tcdC_muts="${tcdC_muts}-${tcdC_muts_current}"
            fi
            count=$(( count + 1 ))
        done
        tcdC_set="${tcdC}\t${tcdC_allele}\t${tcdC_muts}\t${tcdC_stats}"
    fi

    # Check and format for tcdD count
    tcdD_count=$(grep -c 'tcdD' "${tox_input}")
    if [[ "${tcdD_count}" -eq 0 ]]; then
        tcdD=0
        tcdD_set="0\tNOT_FOUND"
    elif [[ "${tcdD_count}" -eq 1 ]]; then
        tcdD_line=$(grep 'tcdD' "${tox_input}")
        tcdD_matchtype=$(echo "${tcdD_line}" | cut -d$'\t' -f5)
        tcdD_raw_IDA=$(echo "${tcdD_line}" | cut -d$'\t' -f10)
        tcdD_raw_IDN=$(echo "${tcdD_line}" | cut -d$'\t' -f11)
        tcdD_raw_length=$(echo "${tcdD_line}" | cut -d$'\t' -f12)
        tcdD_IDA=$(echo "$tcdD_raw_IDA * 100" | bc | cut -d'.' -f1)
        tcdD_IDN=$(echo "$tcdD_raw_IDN * 100" | bc | cut -d'.' -f1)
        tcdD_length=$(echo "$tcdD_raw_length * 100" | bc | cut -d'.' -f1)
        if [[ "${tcdD_matchtype}" = *"Trunc"* ]]; then
            tcdD="TRUNC:CODON[${tcdD_IDA}AA]"
        elif (( $(echo "$tcdD_length*100 < 90" | bc -l) )); then
            tcdD="TRUNC:LENGTH<90[${tcdD_length}COV]"
        else
            tcdD=1
        fi
        tcdD_set="${tcdD}\t[${tcdD_IDN}NT|${tcdD_IDA}AA|${tcdD_length}COV]"
    else
        tcdD=1
        tcdD_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"tcdD"* ]]; then
                tcdD_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        tcdD_count_array=${#tcdD_arr[@]}

        if [[ ${tcdD_count} != ${tcdD_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdD_stats=""
        tcdD_allele=""
        while [[ ${count} -lt ${tcdD_count_array} ]]; do
            tcdD_line="${tcdD_arr[count]}"
            tcdD_allele_current=$(extract_gama_elements 3 "${tcdD_line}")
            tcdD_raw_IDA=$(echo "${tcdD_line}" | cut -d$'\t' -f10)
            tcdD_raw_IDN=$(echo "${tcdD_line}" | cut -d$'\t' -f11)
            tcdD_raw_length=$(echo "${tcdD_line}" | cut -d$'\t' -f12)
            tcdD_IDA=$(echo "$tcdD_raw_IDA * 100" | bc | cut -d'.' -f1)
            tcdD_IDN=$(echo "$tcdD_raw_IDN * 100" | bc | cut -d'.' -f1)
            tcdD_length=$(echo "$tcdD_raw_length * 100" | bc | cut -d'.' -f1)
            tcdD_matchtype=$(echo "${tcdD_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${tcdD} = "" ]]; then
                if [[ "${tcdD_matchtype}" = *"Trunc"* ]]; then
                    tcdD="TRUNC:CODON[${tcdD_IDA}AA]"
                elif (( $(echo "$tcdD_length*100 < 90" | bc -l) )); then
                    tcdD="TRUNC:LENGTH<90[${tcdD_length}COV]"
                else
                    tcdD=1
                fi
            else
                if [[ "${tcdD_matchtype}" = *"Trunc"* ]]; then
                    tcdD="${tcdD}-TRUNC:CODON[${tcdD_IDA}AA]"
                elif (( $(echo "$tcdD_length*100 < 90" | bc -l) )); then
                    tcdD="${tcdD}-TRUNC:LENGTH<90[${tcdD_length}COV]"
                else
                    tcdD="${tcdD}-1"
                fi
            fi
            if [[ "${tcdD_stats}" == "" ]]; then
                tcdD_stats="[${tcdD_IDN}NT|${tcdD_IDA}AA|${tcdD_length}COV]"
            else
                tcdD_stats="${tcdD_stats}-[${tcdD_IDN}NT|${tcdD_IDA}AA|${tcdD_length}COV]"
            fi
            count=$(( count + 1 ))
        done
        tcdD_set="${tcdD}\t${tcdD_stats}"
    fi

    # Check and format for tcdE count
    tcdE_count=$(grep -c 'tcdE' "${tox_input}")
    if [[ "${tcdE_count}" -eq 0 ]]; then
        tcdE=0
        tcdE_set="0\tNOT_FOUND"
    elif [[ "${tcdE_count}" -eq 1 ]]; then
        tcdE_line=$(grep 'tcdE' "${tox_input}")
        tcdE_matchtype=$(echo "${tcdE_line}" | cut -d$'\t' -f5)
        tcdE_raw_IDA=$(echo "${tcdE_line}" | cut -d$'\t' -f10)
        tcdE_raw_IDN=$(echo "${tcdE_line}" | cut -d$'\t' -f11)
        tcdE_raw_length=$(echo "${tcdE_line}" | cut -d$'\t' -f12)
        tcdE_IDA=$(echo "$tcdE_raw_IDA * 100" | bc | cut -d'.' -f1)
        tcdE_IDN=$(echo "$tcdE_raw_IDN * 100" | bc | cut -d'.' -f1)
        tcdE_length=$(echo "$tcdE_raw_length * 100" | bc | cut -d'.' -f1)
        if [[ "${tcdE_matchtype}" = *"Trunc"* ]]; then
            tcdE="TRUNC:CODON[${tcdE_IDA}AA]"
        elif (( $(echo "$tcdE_length*100 < 90" | bc -l) )); then
            tcdE="TRUNC:LENGTH<90[${tcdE_length}COV]"
        else
            tcdE=1
        fi
        tcdE_set="${tcdE}\t[${tcdE_IDN}NT|${tcdE_IDA}AA|${tcdE_length}COV]"
    else
        tcdE=1
        tcdE_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"tcdE"* ]]; then
                tcdE_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        tcdE_count_array=${#tcdE_arr[@]}

        if [[ ${tcdE_count} != ${tcdE_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdE_stats=""
        tcdE_allele=""
        while [[ ${count} -lt ${tcdE_count_array} ]]; do
            tcdE_line="${tcdE_arr[count]}"
            tcdE_allele_current=$(extract_gama_elements 3 "${cdtE_line}")
            tcdD_raw_IDA=$(echo "${tcdD_line}" | cut -d$'\t' -f10)
            tcdD_raw_IDN=$(echo "${tcdD_line}" | cut -d$'\t' -f11)
            tcdD_raw_length=$(echo "${tcdD_line}" | cut -d$'\t' -f12)
            tcdD_IDA=$(echo "$tcdD_raw_IDA * 100" | bc | cut -d'.' -f1)
            tcdD_IDN=$(echo "$tcdD_raw_IDN * 100" | bc | cut -d'.' -f1)
            tcdD_length=$(echo "$tcdE_raw_length * 100" | bc | cut -d'.' -f1)
            tcdE_matchtype=$(echo "${tcdE_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${tcdE} = "" ]]; then
                if [[ "${tcdE_matchtype}" = *"Trunc"* ]]; then
                    tcdE="TRUNC:CODON[${tcdE_IDA}AA]"
                elif (( $(echo "$tcdE_length*100 < 90" | bc -l) )); then
                    tcdE="TRUNC:LENGTH<90[${tcdE_length}COV]"
                else
                    tcd=1
                fi
            else
                if [[ "${tcdE_matchtype}" = *"Trunc"* ]]; then
                    tcdE="${tcdE}-TRUNC:CODON[${tcdE_IDA}AA]"
                elif (( $(echo "$tcdE_length*100 < 90" | bc -l) )); then
                    tcdE="${tcdE}-TRUNC:LENGTH<90[${tcdE_length}COV]"
                else
                    tcdE="${tcdE}-1"
                fi
            fi
            if [[ "${tcdE_stats}" == "" ]]; then
                tcdE_stats="[${tcdE_IDN}NT|${tcdE_IDA}AA|${tcdE_length}COV]"
            else
                tcdE_stats="${tcdE_stats}-[${tcdE_IDN}NT|${tcdE_IDA}AA|${tcdE_length}COV]"
            fi
            count=$(( count + 1 ))
        done
        tcdE_set="${tcdE}\t${tcdE_stats}"
    fi
        # Check and format for PaLoc_NonTox count
    PaLoc_NonTox_count=$(grep -c 'PaLoc_NonTox' "${tox_input}")
    if [[ "${PaLoc_NonTox_count}" -eq 0 ]]; then
        PaLoc_NonTox=0
        PaLoc_NonTox_set="0\tNOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    elif [[ "${PaLoc_NonTox_count}" -eq 1 ]]; then
        PaLoc_NonTox_line=$(grep 'PaLoc_NonTox' "${tox_input}")
        PaLoc_NonTox_matchtype=$(echo "${PaLoc_NonTox_line}" | cut -d$'\t' -f5)
        PaLoc_NonTox_allele=$(extract_gama_elements 3 "${PaLoc_NonTox_line}")
        PaLoc_NonTox_muts=$(clean_mut "${PaLoc_NonTox_line}")
        PaLoc_NonTox_raw_IDA=$(echo "${PaLoc_NonTox_line}" | cut -d$'\t' -f10)
        PaLoc_NonTox_raw_IDN=$(echo "${PaLoc_NonTox_line}" | cut -d$'\t' -f11)
        PaLoc_NonTox_raw_length=$(echo "${PaLoc_NonTox_line}" | cut -d$'\t' -f12)
        PaLoc_NonTox_IDA=$(echo "$PaLoc_NonTox_raw_IDA * 100" | bc | cut -d'.' -f1)
        PaLoc_NonTox_IDN=$(echo "$PaLoc_NonTox_raw_IDN * 100" | bc | cut -d'.' -f1)
        PaLoc_NonTox_length=$(echo "$PaLoc_NonTox_raw_length * 100" | bc | cut -d'.' -f1)
        #Dont check for codon trunc because it has one in the reference
        if [[ "${PaLoc_NonTox_matchtype}" = *"Trunc"* ]]; then
            if [[ "${PaLoc_NonTox_IDA}" = "100" ]]; then
                PaLoc_NonTox=1
            else
                PaLoc_NonTox="TRUNC:CODON[${PaLoc_NonTox_IDA}AA]"
            fi
        elif (( $(echo "$PaLoc_NonTox_length*100 < 90" | bc -l) )); then
            PaLoc_NonTox="TRUNC:LENGTH<90[${PaLoc_NonTox_length}COV]"
        else
            PaLoc_NonTox=1
        fi
        if [[ "${PaLoc_NonTox_muts}" = "No coding mutations" ]]; then
            PaLoc_NonTox_muts='-'
        fi
        PaLoc_NonTox_set="${PaLoc_NonTox}\t${PaLoc_NonTox_allele}\t${PaLoc_NonTox_muts}\t[${PaLoc_NonTox_IDN}NT|${PaLoc_NonTox_IDA}AA|${PaLoc_NonTox_length}COV]"
    else
        PaLoc_NonTox=""
        PaLoc_NonTox_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"PaLoc_NonTox"* ]]; then
                PaLoc_NonTox_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        PaLoc_NonTox_count_array=${#PaLoc_NonTox_arr[@]}

        if [[ ${PaLoc_NonTox_count} != ${PaLoc_NonTox_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        PaLoc_NonTox_stats=""
        PaLoc_NonTox_allele=""
        PaLoc_NonTox_muts=""
        while [[ ${count} -lt ${PaLoc_NonTox_count_array} ]]; do
            PaLoc_NonTox_line="${PaLoc_NonTox_arr[count]}"
            PaLoc_NonTox_allele_current=$(extract_gama_elements 3 "${PaLoc_NonTox_line}")
            PaLoc_NonTox_raw_IDA=$(echo "${PaLoc_NonTox_line}" | cut -d$'\t' -f10)
            PaLoc_NonTox_raw_IDN=$(echo "${PaLoc_NonTox_line}" | cut -d$'\t' -f11)
            PaLoc_NonTox_raw_length=$(echo "${PaLoc_NonTox_line}" | cut -d$'\t' -f12)
            PaLoc_NonTox_muts_current=$(clean_mut "${PaLoc_NonTox_line}")
            PaLoc_NonTox_IDA=$(echo "$PaLoc_NonTox_raw_IDA * 100" | bc | cut -d'.' -f1)
            PaLoc_NonTox_IDN=$(echo "$PaLoc_NonTox_raw_IDN * 100" | bc | cut -d'.' -f1)
            PaLoc_NonTox_length=$(echo "$PaLoc_NonTox_raw_length * 100" | bc | cut -d'.' -f1)
            PaLoc_NonTox_matchtype=$(echo "${PaLoc_NonTox_arr[count]}" | cut -d$'\t' -f5)
            #Dont check for codon trunc since it has one in the reference
            if [[ ${PaLoc_NonTox} = "" ]]; then
                if [[ "${PaLoc_NonTox_matchtype}" = *"Trunc"* ]]; then
                    if [[ "${PaLoc_NonTox_IDA}" = "100" ]]; then
                        PaLoc_NonTox=1
                    else
                        PaLoc_NonTox="TRUNC:CODON[${PaLoc_NonTox_IDA}AA]"
                    fi
                elif (( $(echo "$PaLoc_NonTox_length*100 < 90" | bc -l) )); then
                    PaLoc_NonTox="TRUNC:LENGTH<90[${PaLoc_NonTox_length}COV]"
                else
                    PaLoc_NonTox=1
                fi
            else
                if [[ "${PaLoc_NonTox_matchtype}" = *"Trunc"* ]]; then
                    if [[ "${PaLoc_NonTox_IDA}" = "100" ]]; then
                        PaLoc_NonTox="${PaLoc_NonTox}"-1
                    else
                        PaLoc_NonTox="${PaLoc_NonTox}-TRUNC:CODON[${PaLoc_NonTox_IDA}AA]"
                    fi
                elif (( $(echo "$PaLoc_NonTox_length*100 < 90" | bc -l) )); then
                    PaLoc_NonTox="${PaLoc_NonTox}-TRUNC:LENGTH<90[${PaLoc_NonTox_length}COV]"
                else
                    PaLoc_NonTox="${PaLoc_NonTox}-1"
                fi
            fi
            if [[ "${PaLoc_NonTox_allele}" = "" ]]; then
                PaLoc_NonTox_allele=${PaLoc_NonTox_allele_current}
            else
                PaLoc_NonTox_allele="${PaLoc_NonTox_allele}-${PaLoc_NonTox_allele_current}"
            fi
            if [[ "${PaLoc_NonTox_stats}" = "" ]]; then
                PaLoc_NonTox_stats="[${PaLoc_NonTox_IDN}NT|${PaLoc_NonTox_IDA}AA|${PaLoc_NonTox_length}COV]"
            else
                PaLoc_NonTox_stats="${PaLoc_NonTox_stats}-[${PaLoc_NonTox_IDN}NT|${PaLoc_NonTox_IDA}AA|${PaLoc_NonTox_length}COV]"
            fi
            if [[ "${PaLoc_NonTox_muts_current}" = "No coding mutations" ]]; then
                PaLoc_NonTox_muts_current='-'
            fi
            if [[ "${PaLoc_NonTox_muts}" = "" ]]; then
                PaLoc_NonTox_muts="${PaLoc_NonTox_muts_current}"
            else
                PaLoc_NonTox_muts="${PaLoc_NonTox_muts}-${PaLoc_NonTox_muts_current}"
            fi
            count=$(( count + 1 ))
        done
        PaLoc_NonTox_set="${PaLoc_NonTox}\t${PaLoc_NonTox_allele}\t${PaLoc_NonTox_muts}\t${PaLoc_NonTox_stats}"
    fi

    # Check and format for cdtA count
    cdtA_count=$(grep -c 'cdtA_' "${tox_input}")
    if [[ "${cdtA_count}" -eq 0 ]]; then
        cdtA=0
        cdtA_set="0\tNOT_FOUND"
    elif [[ "${cdtA_count}" -eq 1 ]]; then
        cdtA_line=$(grep 'cdtA_' "${tox_input}")
        cdtA_matchtype=$(echo "${cdtA_line}" | cut -d$'\t' -f5)
        cdtA_raw_IDA=$(echo "${cdtA_line}" | cut -d$'\t' -f10)
        cdtA_raw_IDN=$(echo "${cdtA_line}" | cut -d$'\t' -f11)
        cdtA_raw_length=$(echo "${cdtA_line}" | cut -d$'\t' -f12)
        cdtA_IDA=$(echo "$cdtA_raw_IDA * 100" | bc | cut -d'.' -f1)
        cdtA_IDN=$(echo "$cdtA_raw_IDN * 100" | bc | cut -d'.' -f1)
        cdtA_length=$(echo "$cdtA_raw_length * 100" | bc | cut -d'.' -f1)
        if [[ "${cdtA_matchtype}" = *"Trunc"* ]]; then
            cdtA="TRUNC:CODON[${cdtA_IDA}AA]"
        elif (( $(echo "$cdtA_length*100 < 90" | bc -l) )); then
            cdtA="TRUNC:LENGTH<90[${cdtA_length}COV]"
        else
            cdtA=1
        fi
        cdtA_set="${cdtA}\t[${cdtA_IDN}NT|${cdtA_IDA}AA|${cdtA_length}COV]"
    else
        cdtA=1
        cdtA_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"cdtA"* ]]; then
                cdtA_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        cdtA_count_array=${#cdtA_arr[@]}

        if [[ ${cdtA_count} != ${cdtA_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtA_stats=""
        cdtA_allele=""
        while [[ ${count} -lt ${cdtA_count_array} ]]; do
            cdtA_line="${cdtA_arr[count]}"
            cdtA_allele_current=$(extract_gama_elements 3 "${cdtA_line}")
            cdtA_raw_IDA=$(echo "${cdtA_line}" | cut -d$'\t' -f10)
            cdtA_raw_IDN=$(echo "${cdtA_line}" | cut -d$'\t' -f11)
            cdtA_raw_length=$(echo "${cdtA_line}" | cut -d$'\t' -f12)
            cdtA_IDA=$(echo "$cdtA_raw_IDA * 100" | bc | cut -d'.' -f1)
            cdtA_IDN=$(echo "$cdtA_raw_IDN * 100" | bc | cut -d'.' -f1)
            cdtA_length=$(echo "$cdtA_raw_length * 100" | bc | cut -d'.' -f1)
            cdtA_matchtype=$(echo "${cdtA_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${cdtA} = "" ]]; then
                if [[ "${cdtA_matchtype}" = *"Trunc"* ]]; then
                    cdtA="TRUNC:CODON[${cdtA_IDA}AA]"
                elif (( $(echo "$cdtA_length*100 < 90" | bc -l) )); then
                    cdtA="TRUNC:LENGTH<90[${cdtA_length}COV]"
                else
                    cdtA=1
                fi
            else
                if [[ "${cdtA_matchtype}" = *"Trunc"* ]]; then
                    cdtA="${cdtA}-TRUNC:CODON[${cdtA_IDA}AA]"
                elif (( $(echo "$cdtA_length*100 < 90" | bc -l) )); then
                    cdtA="${cdtA}-TRUNC:LENGTH<90[${cdtA_length}COV]"
                else
                    cdtA="${cdtA}-1"
                fi
            fi
            if [[ "${cdtA_stats}" == "" ]]; then
                cdtA_stats="[${cdtA_IDN}NT|${cdtA_IDA}AA|${cdtA_length}COV]"
            else
                cdtA_stats="${cdtA_stats}-[${cdtA_IDN}NT|${cdtA_IDA}AA|${cdtA_length}COV]"
            fi
            count=$(( count + 1 ))
        done
        cdtA_set="${cdtA}\t${cdtA_stats}"
    fi

    # Check and format for cdtB count
    cdtB_count=$(grep -c 'cdtB_' "${tox_input}")
    if [[ "${cdtB_count}" -eq 0 ]]; then
        cdtB=0
        cdtB_set="0\tNOT_FOUND"
    elif [[ "${cdtB_count}" -eq 1 ]]; then
        cdtB_line=$(grep 'cdtB_' "${tox_input}")
        cdtB_matchtype=$(echo "${cdtB_line}" | cut -d$'\t' -f5)
        cdtB_raw_IDA=$(echo "${cdtB_line}" | cut -d$'\t' -f10)
        cdtB_raw_IDN=$(echo "${cdtB_line}" | cut -d$'\t' -f11)
        cdtB_raw_length=$(echo "${cdtB_line}" | cut -d$'\t' -f12)
        cdtB_IDA=$(echo "$cdtB_raw_IDA * 100" | bc | cut -d'.' -f1)
        cdtB_IDN=$(echo "$cdtB_raw_IDN * 100" | bc | cut -d'.' -f1)
        cdtB_length=$(echo "$cdtB_raw_length * 100" | bc | cut -d'.' -f1)
        if [[ "${cdtB_matchtype}" = *"Trunc"* ]]; then
            cdtB="TRUNC:CODON[${cdtB_IDA}AA]"
        elif (( $(echo "$cdtB_length*100 < 90" | bc -l) )); then
            cdtB="TRUNC:LENGTH<90[${cdtB_length}COV]"
        else
            cdtB=1
        fi
        cdtB_set="${cdtB}\t[${cdtB_IDN}NT|${cdtB_IDA}AA|${cdtB_length}COV]"
    else
        cdtA=1
        cdtB_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"cdtB"* ]]; then
                cdtB_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        cdtB_count_array=${#cdtB_arr[@]}

        if [[ ${cdtB_count} != ${cdtB_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtB_stats=""
        cdtB_allele=""
        while [[ ${count} -lt ${cdtB_count_array} ]]; do
            cdtB_line="${cdtB_arr[count]}"
            cdtB_allele_current=$(extract_gama_elements 3 "${cdtB_line}")
            cdtB_raw_IDA=$(echo "${cdtB_line}" | cut -d$'\t' -f10)
            cdtB_raw_IDN=$(echo "${cdtB_line}" | cut -d$'\t' -f11)
            cdtB_raw_length=$(echo "${cdtB_line}" | cut -d$'\t' -f12)
            cdtB_IDA=$(echo "$cdtB_raw_IDA * 100" | bc | cut -d'.' -f1)
            cdtB_IDN=$(echo "$cdtB_raw_IDN * 100" | bc | cut -d'.' -f1)
            cdtB_length=$(echo "$cdtB_raw_length * 100" | bc | cut -d'.' -f1)
            cdtB_matchtype=$(echo "${cdtB_arr[count]}" | cut -d$'\t' -f5)
           if [[ ${cdtB} = "" ]]; then
                if [[ "${cdtB_matchtype}" = *"Trunc"* ]]; then
                    cdtB="TRUNC:CODON[${cdtB_IDA}AA]"
                elif (( $(echo "$cdtB_length*100 < 90" | bc -l) )); then
                    cdtB="TRUNC:LENGTH<90[${cdtB_length}COV]"
                else
                    cdtB=1
                fi
            else
                if [[ "${cdtB_matchtype}" = *"Trunc"* ]]; then
                    cdtB="${cdtB}-TRUNC:CODON[${cdtB_IDA}AA]"
                elif (( $(echo "$cdtB_length*100 < 90" | bc -l) )); then
                    cdtB="${cdtB}-TRUNC:LENGTH<90[${cdtB_length}COV]"
                else
                    cdtB="${cdtB}-1"
                fi
            fi
            if [[ "${cdtB_stats}" == "" ]]; then
                cdtB_stats="[${cdtB_IDN}NT|${cdtB_IDA}AA|${cdtB_length}COV]"
            else
                cdtB_stats="${cdtB_stats}-[${cdtB_IDN}NT|${cdtB_IDA}AA|${cdtB_length}COV]"
            fi
            count=$(( count + 1 ))
        done
        cdtB_set="${cdtB}\t${cdtB_stats}"
    fi

    # Check and format for cdtR count
    cdtR_count=$(grep -c 'cdtR' "${tox_input}")
    if [[ "${cdtR_count}" -eq 0 ]]; then
        cdtA=0
        cdtR_set="0\tNOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    elif [[ "${cdtR_count}" -eq 1 ]]; then
        cdtR_line=$(grep 'cdtR' "${tox_input}")
        cdtR_matchtype=$(echo "${cdtR_line}" | cut -d$'\t' -f5)
        cdtR_allele=$(extract_gama_elements 3 "${cdtR_line}")
        cdtR_muts=$(clean_mut "${cdtR_line}")
        cdtR_raw_IDA=$(echo "${cdtR_line}" | cut -d$'\t' -f10)
        cdtR_raw_IDN=$(echo "${cdtR_line}" | cut -d$'\t' -f11)
        cdtR_raw_length=$(echo "${cdtR_line}" | cut -d$'\t' -f12)
        cdtR_IDA=$(echo "$cdtR_raw_IDA * 100" | bc | cut -d'.' -f1)
        cdtR_IDN=$(echo "$cdtR_raw_IDN * 100" | bc | cut -d'.' -f1)
        cdtR_length=$(echo "$cdtR_raw_length * 100" | bc | cut -d'.' -f1)
        if [[ "${cdtR_matchtype}" = *"Trunc"* ]]; then
            if [[ "${cdtR_IDA}" = "100" ]]; then
                cdtR=1
            else
                cdtR="TRUNC:CODON[${cdtR_IDA}AA]"
            fi
        elif (( $(echo "$cdtR_length*100 < 90" | bc -l) )); then
            cdtR="TRUNC:LENGTH<90[${cdtR_length}COV]"
        else
            cdtR=1
        fi
        if [[ "${cdtR_muts}" = "No coding mutations" ]]; then
            cdtR_muts='-'
        fi 
        cdtR_set="${cdtR}\t${cdtR_allele}\t${cdtR_muts}\t[${cdtR_IDN}NT|${cdtR_IDA}AA|${cdtR_length}COV]"
    else
        cdtR=1
        cdtR_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"cdtR"* ]]; then
                cdtR_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        cdtR_count_array=${#cdtR_arr[@]}

        if [[ ${cdtR_count} != ${cdtR_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtR_stats=""
        cdtR_allele=""
        cdtR_muts=""
        while [[ ${count} -lt ${cdtR_count_array} ]]; do
            cdtR_line="${cdtR_arr[count]}"
            cdtR_allele_current=$(extract_gama_elements 3 "${cdtR_line}")
            cdtR_muts_current=$(clean_mut "${cdtR_line}")
            cdtR_raw_IDA=$(echo "${cdtR_line}" | cut -d$'\t' -f10)
            cdtR_raw_IDN=$(echo "${cdtR_line}" | cut -d$'\t' -f11)
            cdtR_raw_length=$(echo "${cdtR_line}" | cut -d$'\t' -f12)
            cdtR_IDA=$(echo "$cdtR_raw_IDA * 100" | bc | cut -d'.' -f1)
            cdtR_IDN=$(echo "$cdtR_raw_IDN * 100" | bc | cut -d'.' -f1)
            cdtR_length=$(echo "$cdtR_raw_length * 100" | bc | cut -d'.' -f1)
            cdtR_matchtype=$(echo "${cdtR_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${cdtR} = "" ]]; then
                if [[ "${cdtR_matchtype}" = *"Trunc"* ]]; then
                    if [[ "${cdtR_IDA}" = "100" ]]; then
                        cdtR=1
                    else
                        cdtR="TRUNC:CODON[${cdtR_IDA}AA]"
                    fi
                elif (( $(echo "$cdtR_length*100 < 90" | bc -l) )); then
                    cdtR="TRUNC:LENGTH<90[${cdtR_length}COV]"
                else
                    cdtR=1
                fi
            else
                if [[ "${cdtR_matchtype}" = *"Trunc"* ]]; then
                    if [[ "${cdtR_IDA}" = "100" ]]; then
                        cdtR="${cdtR}-1"
                    else
                        cdtR="${cdtR}-TRUNC:CODON[${cdtR_IDA}AA]"
                    fi
                elif (( $(echo "$cdtR_length*100 < 90" | bc -l) )); then
                    cdtR="${cdtR}-TRUNC:LENGTH<90[${cdtR_length}COV]"
                else
                    cdtR="${cdtR}-1"
                fi
            fi
          
            if [[ "${cdtR_allele}" = "" ]]; then
                cdtR_allele=${cdtR_allele_current}
            else
                cdtR_allele="${cdtR_allele}-${cdtR_allele_current}"
            fi
            if [[ "${cdtR_muts_current}" = "No coding mutations" ]]; then
                cdtR_muts_current='-'
            fi
            if [[ "${cdtR_muts}" = "" ]]; then
                cdtR_muts="${cdtR_muts_current}"
            else
                cdtR_muts="${cdtR_muts}-${cdtR_muts_current}"
            fi
            if [[ "${cdtR_stats}" == "" ]]; then
                cdtR_stats="[${cdtR_IDN}NT|${cdtR_IDA}AA|${cdtR_length}COV]"
            else
                cdtR_stats="${cdtR_stats}-[${cdtR_IDN}NT|${cdtR_IDA}AA|${cdtR_length}COV]"
            fi

            count=$(( count + 1 ))
        done
        cdtR_set="${cdtR}\t${cdtR_allele}\t${cdtR_muts}\t${cdtR_stats}"
    fi

    # Check and format for cdtAB1 count
    cdtAB1_count=$(grep -c 'cdtAB1' "${tox_input}")
    if [[ "${cdtAB1_count}" -eq 0 ]]; then
        cdtAB1=0
        cdtAB1_set="0\tNOT_FOUND"
    elif [[ "${cdtAB1_count}" -eq 1 ]]; then
        cdtAB1_line=$(grep 'cdtAB1' "${tox_input}")
        cdtAB1_matchtype=$(echo "${cdtAB1_line}" | cut -d$'\t' -f5)
        cdtAB1_raw_IDA=$(echo "${cdtAB1_line}" | cut -d$'\t' -f10)
        cdtAB1_raw_IDN=$(echo "${cdtAB1_line}" | cut -d$'\t' -f11)
        cdtAB1_raw_length=$(echo "${cdtAB1_line}" | cut -d$'\t' -f12)
        cdtAB1_IDA=$(echo "$cdtAB1_raw_IDA * 100" | bc | cut -d'.' -f1)
        cdtAB1_IDN=$(echo "$cdtAB1_raw_IDN * 100" | bc | cut -d'.' -f1)
        cdtAB1_length=$(echo "$cdtAB1_raw_length * 100" | bc | cut -d'.' -f1)
        # Dont need normal trunc check since we cdtAB1 is a codon truncation already
        if (( $(echo "$cdtAB1_length*100 < 90" | bc -l) )); then
            cdtAB1="TRUNC:LENGTH<90[${cdtAB1_length}COV]"
        else
            cdtAB1=1
        fi
        cdtAB1_set="${cdtAB1}\t[${cdtAB1_IDN}NT|${cdtAB1_IDA}AA|${cdtAB1_length}COV]"
    else
        cdtAB1=""
        cdtAB1_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"cdtAB1"* ]]; then
                cdtAB1_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        cdtAB1_count_array=${#cdtAB1_arr[@]}

        if [[ ${cdtAB1_count} != ${cdtAB1_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtAB1_stats=""
        cdtAB1_allele=""
        while [[ ${count} -lt ${cdtAB1_count_array} ]]; do
            cdtAB1_line="${cdtAB1_arr[count]}"
            cdtAB1_allele_current=$(extract_gama_elements 3 "${cdtAB1_line}")
            cdtAB1_raw_IDA=$(echo "${cdtAB1_line}" | cut -d$'\t' -f10)
            cdtAB1_raw_IDN=$(echo "${cdtAB1_line}" | cut -d$'\t' -f11)
            cdtAB1_raw_length=$(echo "${cdtAB1_line}" | cut -d$'\t' -f12)
            cdtAB1_IDA=$(echo "$cdtAB1_raw_IDA * 100" | bc | cut -d'.' -f1)
            cdtAB1_IDN=$(echo "$cdtAB1_raw_IDN * 100" | bc | cut -d'.' -f1)
            cdtAB1_length=$(echo "$cdtAB1_raw_length * 100" | bc | cut -d'.' -f1)
            cdtAB1_matchtype=$(echo "${cdtAB1_arr[count]}" | cut -d$'\t' -f5)
            # Dont need trunc codon check since we know it is already part of the sequence
            if [[ ${cdtAB1} = "" ]]; then
                if (( $(echo "$cdtAB1_length*100 < 90" | bc -l) )); then
                    cdtAB1="TRUNC:LENGTH<90[${cdtAB1_length}COV]"
                else
                    cdtAB1=1
                fi
            else
                if (( $(echo "$cdtAB1_length*100 < 90" | bc -l) )); then
                    cdtAB1="${cdtAB1}-TRUNC:LENGTH<90[${cdtAB1_length}COV]"
                else
                    cdtAB1="${cdtAB1}-1"
                fi
            fi
            if [[ "${cdtAB1_stats}" == "" ]]; then
                cdtAB1_stats="[${cdtAB1_IDN}NT|${cdtAB1_IDA}AA|${cdtAB1_length}COV]"
            else
                cdtAB1_stats="${cdtAB1_stats}-[${cdtAB1_IDN}NT|${cdtAB1_IDA}AA|${cdtAB1_length}COV]"
            fi
            count=$(( count + 1 ))
        done
        cdtAB1_set="${cdtAB1}\t${cdtAB1_stats}"
    fi

    # Check and format for cdtAB2 count
    cdtAB2_count=$(grep -c 'cdtAB2' "${tox_input}")
    if [[ "${cdtAB2_count}" -eq 0 ]]; then
        cdtAB2=0
        cdtAB2_set="0\tNOT_FOUND"
    elif [[ "${cdtAB2_count}" -eq 1 ]]; then
        cdtAB2_line=$(grep 'cdtAB2' "${tox_input}")
        cdtAB2_matchtype=$(echo "${cdtAB2_line}" | cut -d$'\t' -f5)
        cdtAB2_raw_IDA=$(echo "${cdtAB2_line}" | cut -d$'\t' -f10)
        cdtAB2_raw_IDN=$(echo "${cdtAB2_line}" | cut -d$'\t' -f11)
        cdtAB2_raw_length=$(echo "${cdtAB2_line}" | cut -d$'\t' -f12)
        cdtAB2_IDA=$(echo "$cdtAB2_raw_IDA * 100" | bc | cut -d'.' -f1)
        cdtAB2_IDN=$(echo "$cdtAB2_raw_IDN * 100" | bc | cut -d'.' -f1)
        cdtAB2_length=$(echo "$cdtAB2_raw_length * 100" | bc | cut -d'.' -f1)
        # Dont need normal trunc check since we cdtAB1 is a codon truncation already
        if (( $(echo "$cdtAB2_length*100 < 90" | bc -l) )); then
            cdtAB2="TRUNC:LENGTH<90[${cdtAB2_length}COV]"
        else
            cdtAB2=1
        fi
        cdtAB2_set="${cdtAB2}\t[${cdtAB2_IDN}NT|${cdtAB2_IDA}AA|${cdtAB2_length}COV]"
    else
        cdtAB2=""
        cdtAB2_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"cdtAB2"* ]]; then
                cdtAB2_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        cdtAB2_count_array=${#cdtAB2_arr[@]}

        if [[ ${cdtAB2_count} != ${cdtAB2_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtAB2_stats=""
        cdtAB2_allele=""
        while [[ ${count} -lt ${cdtAB2_count_array} ]]; do
            cdtAB2_line="${cdtAB2_arr[count]}"
            cdtAB2_allele_current=$(extract_gama_elements 3 "${cdtAB2_line}")
            cdtAB2_raw_IDA=$(echo "${cdtAB2_line}" | cut -d$'\t' -f10)
            cdtAB2_raw_IDN=$(echo "${cdtAB2_line}" | cut -d$'\t' -f11)
            cdtAB2_raw_length=$(echo "${cdtAB2_line}" | cut -d$'\t' -f12)
            cdtAB2_IDA=$(echo "$cdtAB2_raw_IDA * 100" | bc | cut -d'.' -f1)
            cdtAB2_IDN=$(echo "$cdtAB2_raw_IDN * 100" | bc | cut -d'.' -f1)
            cdtAB2_length=$(echo "$cdtAB2_raw_length * 100" | bc | cut -d'.' -f1)
            cdtAB2_matchtype=$(echo "${cdtAB2_arr[count]}" | cut -d$'\t' -f5)
            # Dont need trunc codon check since we know it is already part of the sequence
            if [[ ${cdtAB2} = "" ]]; then
                if (( $(echo "$cdtAB2_length*100 < 90" | bc -l) )); then
                    cdtAB2="TRUNC:LENGTH<90[${cdtAB2_length}COV]"
                else
                    cdtAB2=1
                fi
            else
                if (( $(echo "$cdtAB2_length*100 < 90" | bc -l) )); then
                    cdtAB2="${cdtAB2}-TRUNC:LENGTH<90[${cdtAB2_length}COV]"
                else
                    cdtAB2="${cdtAB2}-1"
                fi
            fi
            if [[ "${cdtAB2_stats}" == "" ]]; then
                cdtAB2_stats="[${cdtAB2_IDN}NT|${cdtAB2_IDA}AA|${cdtAB2_length}COV]"
            else
                cdtAB2_stats="${cdtAB2_stats}-[${cdtAB2_IDN}NT|${cdtAB2_IDA}AA|${cdtAB2_length}COV]"
            fi
            count=$(( count + 1 ))
        done
        cdtAB2_set="${cdtAB2}\t${cdtAB2_stats}"
    fi

        # Check and format for nontox count
    nontox_count=$(grep -c 'cdtNonTox' "${tox_input}")
    if [[ "${nontox_count}" -eq 0 ]]; then
        nontox=0
        nontox_set="0\tNOT_FOUND"
    elif [[ "${nontox_count}" -eq 1 ]]; then
        nontox_line=$(grep 'cdtNonTox' "${tox_input}")
        nontox_matchtype=$(echo "${nontox_line}" | cut -d$'\t' -f5)
        nontox_raw_IDA=$(echo "${nontox_line}" | cut -d$'\t' -f10)
        nontox_raw_IDN=$(echo "${nontox_line}" | cut -d$'\t' -f11)
        nontox_raw_length=$(echo "${nontox_line}" | cut -d$'\t' -f12)
        nontox_IDA=$(echo "$nontox_raw_IDA * 100" | bc | cut -d'.' -f1)
        nontox_IDN=$(echo "$nontox_raw_IDN * 100" | bc | cut -d'.' -f1)
        nontox_length=$(echo "$nontox_raw_length * 100" | bc | cut -d'.' -f1)
        # Don't need a codon check as nontox is a known truncation
        if [[ "${nontox_matchtype}" = *"Trunc"* ]]; then
            if [[ ${nontox_IDA} = "100" ]]; then
                nontox=1
            else
                nontox="TRUNC:CODON[${nontox_IDA}AA]"
            fi
        elif (( $(echo "$nontox_length*100 < 90" | bc -l) )); then
            nontox="TRUNC:LENGTH<90[${nontox_length}COV]"
        else
            nontox=1
        fi
        nontox_set="${nontox}\t[${nontox_IDN}NT|${nontox_IDA}AA|${nontox_length}COV]"
    else
        nontox=1
        nontox_array=()
        tc=0
        while IFS= read -r var; do
            hit=$(echo "${var}" | cut -d$'\t' -f1)
            if [[ "${hit}" == *"nontox"* ]]; then
                nontox_arr[$tc]="${var}"
                tc=$(( tc + 1))
            fi
        done < "${tox_input}"
        nontox_count_array=${#nontox_arr[@]}

        if [[ ${nontox_count} != ${nontox_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        nontox_stats=""
        nontox_allele=""
        while [[ ${count} -lt ${nontox_count_array} ]]; do
            nontox_line="${nontox_arr[count]}"
            nontox_allele_current=$(extract_gama_elements 3 "${nontox_line}")
            nontox_raw_IDA=$(echo "${nontox_line}" | cut -d$'\t' -f10)
            nontox_raw_IDN=$(echo "${nontox_line}" | cut -d$'\t' -f11)
            nontox_raw_length=$(echo "${nontox_line}" | cut -d$'\t' -f12)
            nontox_IDA=$(echo "$nontox_raw_IDA * 100" | bc | cut -d'.' -f1)
            nontox_IDN=$(echo "$nontox_raw_IDN * 100" | bc | cut -d'.' -f1)
            nontox_length=$(echo "$nontox_raw_length * 100" | bc | cut -d'.' -f1)
            nontox_matchtype=$(echo "${nontox_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${nontox} = "" ]]; then
                if [[ "${nontox_matchtype}" = *"Trunc"* ]]; then
                    if [[ ${nontox_IDA} = "100" ]]; then
                        nontox=1
                    else
                        nontox="TRUNC:CODON[${nontox_IDA}AA]"
                    fi
                elif (( $(echo "$nontox_length*100 < 90" | bc -l) )); then
                        nontox="TRUNC:LENGTH<90[${nontox_length}COV]"
                    else
                        nontox=1
                    fi
            else
                if [[ "${nontox_matchtype}" = *"Trunc"* ]]; then
                    if [[ ${nontox_IDA} = "100" ]]; then
                        nontox="${nontox}-1"
                    else
                        nontox="${nontox}-TRUNC:CODON[${nontox_IDA}AA]"
                    fi
                elif (( $(echo "$nontox_length*100 < 90" | bc -l) )); then
                    nontox="${nontox}-TRUNC:LENGTH<90[${nontox_length}COV]"
                else
                    nontox="${nontox}-1"
                fi
            fi
            if [[ "${nontox_stats}" == "" ]]; then
                nontox_stats="[${nontox_IDN}NT|${nontox_IDA}AA|${nontox_length}COV]"
            else
                nontox_stats="${nontox_stats}-[${nontox_IDN}NT|${nontox_IDA}AA|${nontox_length}COV]"
            fi
            count=$(( count + 1 ))
        done
        nontox_set="${nontox}\t${nontox_stats}"
    fi
else
    echo "No tox file found"
    tcdA_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file"
    tcdB_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file"
    tcdC_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file"
    tcdD_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file"
    tcdE_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file"
    cdtA_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file"
    cdtB_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file"
    cdtAB1_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file"
    cdtAB2_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file"
    cdtR_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file"
    nontox_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file"
fi

if [[ -f "${clade_input}" ]]; then
    clade_lines=$(cat "${clade_input}" | wc -l)
    if [[ "${clade_lines}" -eq 2 ]]; then
        clade=$(tail -n1 "${clade_input}" | cut -d$'\t' -f2)
        mlst=$(tail -n1 "${clade_input}" | cut -d$'\t' -f3)
        if [[ -f "${xwalkRT_file}" ]]; then
            # Do lookup
            xrt="unset"
            while IFS= read -r var; do
                lst=$(echo "${var}" | cut -d$'\t' -f1)
                echo "X:${lst}-M:${mlst}"
                if [[ "${lst}" = "${mlst}" ]]; then
                    xrt=$(echo "${var}" | cut -d$'\t' -f2)
                    IFS=',' read -a rts <<< "${xrt}"
                    rt_count=${#rts[@]}
                    if [ ${rt_count} -gt 0 ]; then
                        new_xrt=""
                        for rib in "${rts[@]}"; do
                            rt_length=${#rib}
                            if [ $rt_length -eq 2 ]; then
                                new_xrt="${new_rt},0${rib}"
                            elif [ $rt_length -eq 1 ]; then
                                new_xrt="${new_xrt},00${rib}"
                            else
                                new_xrt="${new_xrt},${rib}"
                            fi
                        done
                        if [[ ${new_xrt:0:1} == "," ]] ; then
                            new_xrt="${new_xrt:1}"
                        fi
                    else
                        new_xrt="unset"
                    fi
                    xrt="${new_xrt}"
                    break
                fi
            done < "${xwalkRT_file}"
            #xrt="something"
        else
            xrt="No_lookup_file"
        fi
        if [[ "${xrt}" = "unset" ]]; then
            xrt="${mlst} has no crosswalk match"
        fi
    else
        clade="Clade_file_incorrect"
        xrt="Clade_file_incorrect"
    fi
else
    echo "No clade/mlst file, carry one"
    clade="No_clade/MLST_file"
fi

#declare -A AA_POINTS=( [dacS]="E238D V183A" [fur]="E41K" [gdpP]="E328Stop truncation__at__codon__328__(of__665__codons)" [glyC]="A229T" [lscR]="V76A" [murG]="P109L" [nifJ]="G423E" [rpoB]="V1143D V1143F V1143G V1143L Q1074R Q1074H Q1074K" [rpoC]="D245Y D1127E D237Y Q781R R89G" [thiH]="S328F" [vanR]="T115A" [vanS]="G319D R314L R314H S313F T349I" [gyrA]="A117S A118S A118T A384D A92E D103N D71G D71V D81N E123K L345I P116A R90K T82A T82I T82V V43D" [gyrB]="D426N D426V E466K E466V I139R L444F Q434K R377G R447K R447L S366V S464T V130I")
#declare -A NT_POINTS=( [feoB]="117DelA 1__bp__Deletion__at__119" [hemN]="Y214Stop 1__bp__Deletion__at__641" [hsmA]="372DelA 1__bp__Deletion__at__371" [lscR]="153DelA 1__bp__Deletion__at__152" [marR]="349DelT 1__bp__Deletion__at__355" [PNimB]="T115G" [sdaB]="883DelGCA 879DelACG 3__bp__Deletion__at__886" )
declare -A AA_POINTS=( [dacS]="E238D V183A" [fur]="E41K" [gdpP]="E328Stop truncation__at__codon__328__(of__665__codons)" [glyC]="A229T" [murG]="P109L" [nifJ]="Q803R G423E" [rpoB]="V1143D V1143F V1143G V1143L Q1074R Q1074H Q1074K" [rpoC]="D245Y D1127E D237Y Q781R R89G" [thiH]="S328F" [vanR]="T115A" [vanS]="G319D R314L R314H S313F T349I" [gyrA]="A117S A118S A118T A384D A92E D103N D71G D71V D81N E123K L345I P116A R90K T82A T82I T82V V43D" [gyrB]="D426N D426V E466K E466V I139R L444F Q434K R377G R447K R447L S366V S464T V130I" [feoB]="117DelA 1__bp__Deletion__at__120" [hemN]="Y214Stop 1__bp__Deletion__at__642" [hsmA]="372DelA 1__bp__Deletion__at__372" [lscR]="V76A 153DelA 1__bp__Deletion__at__153" [marR]="349DelT 1__bp__Deletion__at__356" [sdaB]="883DelGCA 879DelACG 3__bp__Deletion__at__887")
declare -A NT_POINTS=(  [PNimB]="T115G" )


#printf "%s\n" "${!AA_POINTS[@]}"
#printf "%s\n" "${!NT_POINTS[@]}"

if [[ -f "${aa_mut_file}" ]]; then
    ## Parse the AR file once ready and figure out how to make it expandable
    gyrA_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    gyrB_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    dacS_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    fur_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    gdpP_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    glyC_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    lscR_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    murG_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    rpoB_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    rpoC_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    thiH_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    vanR_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    vanS_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    nifJ_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    feoB_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    hemN_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    hsmA_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    #lscRNT_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    marR_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    PNimB_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    sdaB_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"

    for i in "${!AA_POINTS[@]}"
        do
        #echo "Doing $i"
        while IFS= read -r var; do
            gene_section=$(echo ${var} | awk '{print $1}')
            #echo "Checking $i against ${gene_section}"
            if [[ "${gene_section}" == *"${i}"* ]]; then
                aa_mut_line="${var}"
                all_muts=$(clean_mut "${aa_mut_line}")
                if [[ "${all_muts}" = "No coding mutations" ]]; then
                    all_muts=""
                fi
                line_raw_IDA=$(echo "${aa_mut_line}" | cut -d$'\t' -f10)
                line_raw_IDN=$(echo "${aa_mut_line}" | cut -d$'\t' -f11)
                line_raw_length=$(echo "${aa_mut_line}" | cut -d$'\t' -f12)
                if [[ "${line_raw_IDA}" = *"."* ]]; then
                    line_IDAA=$(echo "$line_raw_IDA * 100" | bc | cut -d'.' -f1)
                else
                    line_IDAA=$(echo "$line_raw_IDA * 100" | bc)
                fi
                if [[ "${line_raw_IDN}" = *"."* ]]; then
                    line_IDNT=$(echo "$line_raw_IDN * 100" | bc | cut -d'.' -f1)
                else
                    line_IDNT=$(echo "$line_raw_IDN * 100" | bc)
                fi
                if [[ "${line_length}" = *"."* ]]; then
                    line_length=$(echo "$line_raw_length * 100" | bc | cut -d'.' -f1)
                else
                    line_length=$(echo "$line_raw_length * 100" | bc)
                fi
                matchtype=$(echo "${aa_mut_line}" | cut -d$'\t' -f5)
                echo "$i,$line_raw_IDA,$line_IDAA:$line_raw_IDN,$line_IDNT:${line_raw_length},$line_length"
                #echo "AM:${all_muts}"
                IFS=',' read -r -a mut_array <<< "$all_muts"
                known_muts=()
                other_muts=()
                if [[ "${matchtype}" = "Truncation" ]] || [[ "${matchtype}" = "Indel Truncation" ]]; then
                    other_muts[0]="TRUNC:CODON[${line_IDAA}AA]"
                elif [[ "${line_length}" -le 90 ]]; then
                    other_muts[0]="TRUNC:LENGTH<90[${line_length}COV]"
                fi
                for j in "${mut_array[@]}"; do
                    added="False"
                    for k in ${AA_POINTS[$i]}; do
                        trim_j=$(echo -e "${j//__/ }" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
                        trim_k=$(echo -e "${k//__/ }" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
                        #echo "${j}/${trim_j} VS ${k}/${trim_k}"    
                        if [[ "${trim_j//__/ }" = "${trim_k//__/ }" ]]; then
                            #echo "FOUND - ${trim_j//__/ }"
                            added="True"
                            known_muts+=("${trim_j//__/ }")
                        fi
                    done
                    if [[ "${added}" = "False" ]]; then
                        other_muts+=("${trim_j//__/ }")
                    fi
                done
                muts_string=""
                #echo "KM: ${known_muts[@]}"
                for item in "${known_muts[@]}"; do
                    if [[ "${muts_string}" = "" ]]; then
                        muts_string="${item}"
                    else
                        muts_string="${muts_string},${item}"
                    fi
                done
                other_muts_string=""
                other_mut_count=0
                #echo "OM: ${other_muts[@]}"
                for item in "${other_muts[@]}"; do
                    #echo "OMS-pre:${other_muts_string},${item}"
                    if [[ "${item}" = "No mutations" ]]; then
                        # Dont do anything if it is just saying there is nothing there
                        :
                    else
                        if [[ "${other_muts_string}" = "" ]]; then
                            other_muts_string="${item}"
                        else
                            other_muts_string="${other_muts_string},${item}"
                        fi
                        other_mut_count=$(( other_mut_count + 1 ))
                        #echo "OMS-post:${other_muts_string}"
                    fi
                done
                if [[ "${other_mut_count}" -gt 10 ]]; then
                    other_muts_string="${other_mut_count} other mutations"
                fi
                if [[ "${mut_string}" = "" ]]; then
                    mut_string="-"
                fi
                if [[ "${other_mut_string}" = "" ]]; then
                    mut_string="-"
                fi
            # echo "0-${i}"
            # echo "1-${aa_mut_line}"
            # echo "2-${all_muts}"
            # echo "3-${line_IDAA}"
            # echo "4-${line_IDNT}"
            # echo "5-${line_length}"
            # echo "6-${known_muts}"
            # echo "7-${AA_POINTS[$i]}"

                case "$i" in
                    gyrA)
                        gyrA_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    gyrB)
                        gyrB_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    dacS)
                        dacS_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    fur)
                        fur_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    gdpP)
                        gdpP_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    glyC)
                        glyC_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    lscR)
                        lscR_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    murG)
                        murG_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    rpoB)
                        rpoB_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    rpoC)
                        rpoC_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    thiH)
                        thiH_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    vanR)
                        vanR_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    vanS)
                        vanS_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    nifJ)
                        nifJ_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    feoB)
                        feoB_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    hemN)
                        hemN_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    hsmA)
                        hsmA_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    marR)
                        marR_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    sdaB)
                        sdaB_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_IDAA}AA|${line_length}COV]";;
                    *)
                        echo "Unknown gene: $i" ;;
                esac
                continue 1
            else
                # echo "$i not found this line"
                :
            fi
        done < "${aa_mut_file}"
    done
else
    echo "NO_AA_mut_file"
    gyrA_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    gyrB_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    dacS_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    fur_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    gdpP_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    glyC_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    lscRAA_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    murG_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    rpoB_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    rpoC_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    thiH_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    vanR_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    vanS_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
    nifJ_set="NO_AA_mut_file\tNO_AA_mut_file\tNO_AA_mut_file"
fi

if [[ -f "${nt_mut_file}" ]]; then
    ## Parse the AR file once ready and figure out how to make it expandable
    #feoB_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    #hemN_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    #hsmA_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    #lscRNT_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    #marR_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    PNimB_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    #sdaB_set="NOT_FOUND\tNOT_FOUND\tNOT_FOUND"
    #cat "${nt_mut_file}"
    for i in "${!NT_POINTS[@]}"
    do
        while IFS= read -r var; do
            gene_section=$(echo ${var} | awk '{print $1}')
            if [[ "${gene_section}" == *"${i}"* ]]; then
                #echo "found ${i}"
                nt_mut_line=${var}
                all_muts=$(clean_mut "${nt_mut_line}")
                if [[ ${all_muts: -1} = "," ]]; then
                    all_muts="${all_muts::-1}"
                elif [[ "${all_muts}" = "Exact match" ]]; then
                    all_muts=""
                fi
                line_raw_IDA=$(echo "${nt_mut_line}" | cut -d$'\t' -f13)
                line_raw_IDN=$(echo "${nt_mut_line}" | cut -d$'\t' -f14)
                line_raw_length=$(echo "${nt_mut_line}" | cut -d$'\t' -f15)
                line_IDNT=$(echo "$line_raw_IDN * 100" | bc | cut -d'.' -f1)
                line_IDNT=$(echo "$line_raw_IDA * 100" | bc | cut -d'.' -f1)
                line_length=$(echo "$line_raw_length * 100" | bc | cut -d'.' -f1)
                IFS=',' read -r -a mut_array <<< "$all_muts"
                known_muts=()
                other_muts=()
                if [[ "${line_length}" -le 90 ]]; then
                    other_muts[0]="TRUNC:LENGTH<90[${line_length}COV]"
                fi
                for j in "${mut_array[@]}"; do
                    added="False"
                    for k in ${NT_POINTS[$i]}; do
                        trim_j=$(echo -e "${j//__/ }" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
                        trim_k=$(echo -e "${k//__/ }" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
                        #echo "${j}/${trim_j} VS ${k}/${trim_k}"    
                        if [[ "${trim_j//__/ }" = "${trim_k//__/ }" ]]; then
                            #echo "FOUND - ${trim_j//__/ }"
                            added="True"
                            known_muts+=("${trim_j//__/ }")
                        fi
                    done
                    if [[ "${added}" = "False" ]]; then
                        other_muts+=("${trim_j//__/ }")
                    fi
                done
                muts_string=""
                #echo "KM: ${known_muts[@]}"
                for item in "${known_muts[@]}"; do
                    if [[ "${muts_string}" = "" ]]; then
                        muts_string="${item}"
                    else
                        muts_string="${muts_string},${item}"
                    fi
                done
                other_muts_string=""
                other_mut_count=0
                #echo "OM: ${other_muts[@]}"
                for item in "${other_muts[@]}"; do
                    if [[ "${item}" = "No mutations" ]]; then
                        # Dont do anything if it is just saying there is nothing there
                        :
                    else
                        #echo "OMS-pre:${other_muts_string},${item}"
                        if [[ "${other_muts_string}" = "" ]]; then
                            other_muts_string="${item}"
                        else
                            other_muts_string="${other_muts_string},${item}"
                        fi
                        other_mut_count=$(( other_mut_count + 1 ))
                        #echo "OMS-post:${other_muts_string}"
                    fi
                done
                if [[ "${other_mut_count}" -gt 10 ]]; then
                    other_muts_string="${other_mut_count} other mutations"
                fi
                if [[ "${mut_string}" = "" ]]; then
                    mut_string="-"
                fi
                if [[ "${other_mut_string}" = "" ]]; then
                    mut_string="-"
                fi
                #echo "0-${i}"
                #echo "1-${nt_mut_line}"
                #echo "2-${all_muts}"
                #echo "4-${line_IDNT}"
                #echo "5-${line_length}"
                #echo "6-${known_muts}"
                #echo "7-${NT_POINTS[$i]}"
                case "$i" in
                    #feoB)
                    #    feoB_set="${muts_string}\t${other_muts_string}\t${line_IDNT}NT|${line_IDNA}AA|${line_length}COV";;
                    #hemN)
                    #    hemN_set="${muts_string}\t${other_muts_string}\t${line_IDNT}NT|${line_IDNA}AA|${line_length}COV";;
                    #hsmA)
                    #    hsmA_set="${muts_string}\t${other_muts_string}\t${line_IDNT}NT|${line_IDNA}AA|${line_length}COV";;
                    #lscR)
                    #    lscRNT_set="${muts_string}\t${other_muts_string}\t${line_IDNT}NT|${line_IDNA}AA|${line_length}COV";;
                    #marR)
                    #    marR_set="${muts_string}\t${other_muts_string}\t${line_IDNT}NT|${line_IDNA}AA|${line_length}COV";;
                    PNimB)
                        PNimB_set="${muts_string}\t${other_muts_string}\t[${line_IDNT}NT|${line_length}COV]";;
                    #sdaB)
                    #    sdaB_set="${muts_string}\t${other_muts_string}\t${line_IDNT}NT|${line_IDNA}AA|${line_length}COV";;
                    *)
                        echo "Unknown gene: $i" ;;
                esac
                continue 1
            else
                #echo "$i not found this line"
                :
            fi
        done < "${nt_mut_file}"
    done
else
    echo "NO_NT_mut_file"
    #feoB_set="NO_NT_mut_file\tNO_NT_mut_file\tNO_NT_mut_file"
    #hemN_set="NO_NT_mut_file\tNO_NT_mut_file\tNO_NT_mut_file"
    #hsmA_set="NO_NT_mut_file\tNO_NT_mut_file\tNO_NT_mut_file"
    #lscRNT_set="NO_NT_mut_file\tNO_NT_mut_file\tNO_NT_mut_file"
    #marR_set="NO_NT_mut_file\tNO_NT_mut_file\tNO_NT_mut_file"
    PNimB_set="NO_NT_mut_file\tNO_NT_mut_file\tNO_NT_mut_file"
    #sdaB_set="NO_NT_mut_file\tNO_NT_mut_file\tNO_NT_mut_file"
fi

echo -e "clade:${clade}
toxinoptype:${toxinotype}
tcdA:${tcdA_set}
tcdB:${tcdB_set}
tcdc:${tcdC_set}
tcdD:${tcdD_set}
tcdE:${tcdE_set}
cdtA:${cdtA_set}
cdtB:${cdtB_set}
cdtR:${cdtR_set}
cdtAB1:${cdtAB1_set}
cdtAB2:${cdtAB2_set}
Non:${nontox_set}
gyrA:${gyrA_set}
gyrB:${gyrB_set}
dacS:${dacS_set}
feoB:${feoB_set}
fur:${fur_set}
gdpP:${gdpP_set}
glyC:${glyC_set}
hemN:${hemN_set}
hsmA:${hsmA_set}
lscR:${lscR_set}
marR:${marR_set}
murG:${murG_set}
nifJ:${nifJ_set}
PNimG:${PNimB_set}
rpoB:${rpoB_set}
rpoC:${rpoC_set}
sdaB:${sdaB_set}
thiH:${thiH_set}
vanR:${vanR_set}
vanS:${vanS_set}
xwalk:${xrt}
ML_RT:${ML_RT}
plasmids:${plasmids}
"

# Loop through the genes in the list and format to match desired output style
if [[ ! -f "${output}" ]]; then
    echo -e "WGS_ID\tMLST Clade\tDiffbase_Toxinotype\ttcdA_presence\tDiffbase_Toxin-A_sub-type\ttcdA [%Nuc_Identity | %AA_Identity | %Coverage]\ttcdB_presence\tDiffbase_Toxin-B_sub-type\ttcdB [%Nuc_Identity | %AA_Identity | %Coverage]\ttcdC_presence\ttcdC_Variant\ttcdC other mutations\ttcdC [%Nuc_Identity | %AA_Identity | %Coverage]\ttcdR_presence\ttcdR [%Nuc_Identity | %AA_Identity | %Coverage]\ttcdE_presence\ttcdE [%Nuc_Identity | %AA_Identity | %Coverage]\tPaLoc_NonTox_presence\tPaLoc_NonTox_Variant\tPaLoc_NonTox other mutations\tPaLoc_NonTox [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtA_presence\tcdtA [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtB_presence\tcdtB [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtR_presence\tcdtR_Variant\tcdtR other mutations\tcdtR [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtAB1_presence\tcdtAB1 [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtAB2_presence\tcdtAB2 [%Nuc_Identity | %AA_Identity | %Coverage]\tcdt_NonTox_presence\tcdt_NonTox[%Nuc_Identity | %AA_Identity | %Coverage]\tgyrA known mutations\tgyrA other mutations\tgyrA [%Nuc_Identity | %AA_Identity | %Coverage]\tgyrB known mutations\tgyrB other mutations\tgyrB [%Nuc_Identity | %AA_Identity | %Coverage]\tdacS known mutations\tdacS other mutations\tdacS [%Nuc_Identity | %AA_Identity | %Coverage]\tfeoB known mutations\tfeoB other mutations\tfeoB [%Nuc_Identity | %AA_Identity | %Coverage]\tfur known mutations\tfur other mutations\tfur [%Nuc_Identity | %AA_Identity | %Coverage]\tgdpP known mutations\tgdpP other mutations\tgdpP [%Nuc_Identity | %AA_Identity | %Coverage]\tglyC known mutations\tglyC other mutations\tglyC [%Nuc_Identity | %AA_Identity | %Coverage]\themN known mutations\themN other mutations\themN [%Nuc_Identity | %AA_Identity | %Coverage]\thsmA known mutations\thsmA other mutations\thsmA [%Nuc_Identity | %AA_Identity | %Coverage]\tlscR known mutations\tlscR other mutations\tlscR [%Nuc_Identity | %AA_Identity | %Coverage]\tmarR known mutations\tmarR other mutations\tmarR [%Nuc_Identity | %AA_Identity | %Coverage]\tmurG known mutations\tmurG other mutations\tmurG [%Nuc_Identity | %AA_Identity | %Coverage]\tnifJ known mutations\tnifJ other mutations\tnifJ [%Nuc_Identity | %AA_Identity | %Coverage]\tPNimB known mutations\tPNimB other mutations\tPNimB [%Nuc_Identity | %Coverage]\trpoB known mutations\trpoB other mutations\trpoB [%Nuc_Identity | %AA_Identity | %Coverage]\trpoC known mutations\trpoC other mutations\trpoC [%Nuc_Identity | %AA_Identity | %Coverage]\tsdaB known mutations\tsdaB other mutations\tsdaB [%Nuc_Identity | %AA_Identity | %Coverage]\tthiH known mutations\tthiH other mutations\tthiH [%Nuc_Identity | %AA_Identity | %Coverage]\tvanR known mutations\tvanR other mutations\tvanR [%Nuc_Identity | %AA_Identity | %Coverage]\tvanS known mutations\tvanS other mutations\tvanS [%Nuc_Identity | %AA_Identity | %Coverage]\tCEMB RT Crosswalk\tInferred RT\tProbability\tML Method\tML Note" > "${output}"
fi
echo -e "${sample_name}\t${clade}\t${toxinotype}\t${tcdA_set}\t${tcdB_set}\t${tcdC_set}\t${tcdD_set}\t${tcdE_set}\t${PaLoc_NonTox_set}\t${cdtA_set}\t${cdtB_set}\t${cdtR_set}\t${cdtAB1_set}\t${cdtAB2_set}\t${nontox_set}\t${gyrA_set}\t${gyrB_set}\t${dacS_set}\t${feoB_set}\t${fur_set}\t${gdpP_set}\t${glyC_set}\t${hemN_set}\t${hsmA_set}\t${lscR_set}\t${marR_set}\t${murG_set}\t${nifJ_set}\t${PNimB_set}\t${rpoB_set}\t${rpoC_set}\t${sdaB_set}\t${thiH_set}\t${vanR_set}\t${vanS_set}\t${xrt}\t${ML_RT}"  >> "${output}"

# Keep original output for when we add in plasmids again
## Loop through the genes in the list and format to match desired output style
#if [[ ! -f "${output}" ]]; then
#    echo -e "WGS_ID\tMLST Clade\tDiffbase_Toxinotype\ttcdA_presence\tDiffbase_Toxin-A_sub-type\ttcdA [%Nuc_Identity | %AA_Identity | %Coverage]\ttcdB_presence\tDiffbase_Toxin-B_sub-type\ttcdB [%Nuc_Identity | %AA_Identity | %Coverage]\ttcdC_presence\ttcdC_Variant\ttcdC other mutations\ttcdC [%Nuc_Identity | %AA_Identity | %Coverage]\ttcdR_presence\ttcdR [%Nuc_Identity | %AA_Identity | %Coverage]\ttcdE_presence\ttcdE [%Nuc_Identity | %AA_Identity | %Coverage]\tPaLoc_NonTox_presence\tPaLoc_NonTox_Variant\tPaLoc_NonTox other mutations\tPaLoc_NonTox [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtA_presence\tcdtA [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtB_presence\tcdtB [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtR_presence\tcdtR_Variant\tcdtR other mutations\tcdtR [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtAB1_presence\tcdtAB1 [%Nuc_Identity | %AA_Identity | %Coverage]\tcdtAB2_presence\tcdtAB2 [%Nuc_Identity | %AA_Identity | %Coverage]\tcdt_NonTox_presence\tcdt_NonTox[%Nuc_Identity | %AA_Identity | %Coverage]\tgyrA known mutations\tgyrA other mutations\tgyrA [%Nuc_Identity | %AA_Identity | %Coverage]\tgyrB known mutations\tgyrB other mutations\tgyrB [%Nuc_Identity | %AA_Identity | %Coverage]\tdacS known mutations\tdacS other mutations\tdacS [%Nuc_Identity | %AA_Identity | %Coverage]\tfeoB known mutations\tfeoB other mutations\tfeoB [%Nuc_Identity | %AA_Identity | %Coverage]\tfur known mutations\tfur other mutations\tfur [%Nuc_Identity | %AA_Identity | %Coverage]\tgdpP known mutations\tgdpP other mutations\tgdpP [%Nuc_Identity | %AA_Identity | %Coverage]\tglyC known mutations\tglyC other mutations\tglyC [%Nuc_Identity | %AA_Identity | %Coverage]\themN known mutations\themN other mutations\themN [%Nuc_Identity | %AA_Identity | %Coverage]\thsmA known mutations\thsmA other mutations\thsmA [%Nuc_Identity | %AA_Identity | %Coverage]\tlscR known mutations\tlscR other mutations\tlscR [%Nuc_Identity | %AA_Identity | %Coverage]\tmarR known mutations\tmarR other mutations\tmarR [%Nuc_Identity | %AA_Identity | %Coverage]\tmurG known mutations\tmurG other mutations\tmurG [%Nuc_Identity | %AA_Identity | %Coverage]\tnifJ known mutations\tnifJ other mutations\tnifJ [%Nuc_Identity | %AA_Identity | %Coverage]\tPNimB known mutations\tPNimB other mutations\tPNimB [%Nuc_Identity | %Coverage]\trpoB known mutations\trpoB other mutations\trpoB [%Nuc_Identity | %AA_Identity | %Coverage]\trpoC known mutations\trpoC other mutations\trpoC [%Nuc_Identity | %AA_Identity | %Coverage]\tsdaB known mutations\tsdaB other mutations\tsdaB [%Nuc_Identity | %AA_Identity | %Coverage]\tthiH known mutations\tthiH other mutations\tthiH [%Nuc_Identity | %AA_Identity | %Coverage]\tvanR known mutations\tvanR other mutations\tvanR [%Nuc_Identity | %AA_Identity | %Coverage]\tvanS known mutations\tvanS other mutations\tvanS [%Nuc_Identity | %AA_Identity | %Coverage]\tCEMB RT Crosswalk\tInferred RT\tProbability\tML Method\tML Note\tPlasmid Info\t${other_AR_header}" > "${output}"
#fi
#echo -e "${sample_name}\t${clade}\t${toxinotype}\t${tcdA_set}\t${tcdB_set}\t${tcdC_set}\t${tcdD_set}\t${tcdE_set}\t${PaLoc_NonTox_set}\t${cdtA_set}\t${cdtB_set}\t${cdtR_set}\t${cdtAB1_set}\t${cdtAB2_set}\t${nontox_set}\t${gyrA_set}\t${gyrB_set}\t${dacS_set}\t${feoB_set}\t${fur_set}\t${gdpP_set}\t${glyC_set}\t${hemN_set}\t${hsmA_set}\t${lscR_set}\t${marR_set}\t${murG_set}\t${nifJ_set}\t${PNimB_set}\t${rpoB_set}\t${rpoC_set}\t${sdaB_set}\t${thiH_set}\t${vanR_set}\t${vanS_set}\t${xrt}\t${ML_RT}\t${plasmids}"  >> "${output}"