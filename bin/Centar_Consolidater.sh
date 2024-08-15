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
	echo "Usage is ../Centar_Consolidater.sh -t path_to_phx_isolate_GAMMA_tox_file -c path_to_clade_file -o output_file -y toxinotype_file -r Ribotype_file -p plasmid_file -s sample_name"
}

version="1.0"

tox_input=""
clade_input=""
# Parse command line options
options_found=0
while getopts ":h?t:c:o:y:a:r:p:Vs:" option; do
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
			ar_file=${OPTARG};;
        r)
			echo "Option -r triggered, argument = ${OPTARG}"
			rt_file=${OPTARG};;
        s)
			echo "Option -s triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
        p)
			echo "Option -p triggered, argument = ${OPTARG}"
			plasmid_file=${OPTARG};;
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
	echo "Centar_Consolidater.sh: ${version}"
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


if [[ -f "${ar_file}" ]]; then
    ## Parse the AR file once ready and figure out how to make it expandable
    echo "NO_OTHER_AR_FILE_YET"
    other_AR="NO_AR_FILE_YET"
    other_AR_count=0
    other_AR_header="Other_Cdiff_AR"
else
    echo "NO_OTHER_AR_FILE"
    other_AR="NO_AR_FILE"
    other_AR_count=0
    other_AR_header="Other_Cdiff_AR"
fi

if [[ -f "${rt_file}" ]]; then
    ML_RT="unset"
    ## Parse the RT file once ready and figure out how to make it expandable
    echo "No Ribotype file yet"
    ML_RT="NO_RT_FILE_YET\tNO_RT_FILE_YET"
else
    echo "No Ribotype file"
    ML_RT="NO_RT_FILE_YET\tNO_RT_FILE"
fi

if [[ -f "${rt_file}" ]]; then
    ## Parse the RT file once ready and figure out how to make it expandable
    echo "No plasmid info file yet"
    plasmids="NO_PLASMID_FILE_YET"
else
    echo "No plasmid info file"
    plasmids="NO_PLASMID_FILE_YET"
fi

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
                echo "${current_type}, ${current_subtype}"
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
        tcdA_set="0\t0\tNA\tNA\tNA|NA|NA"
    elif [[ "${tcdA_count}" -eq 1 ]]; then
        tcdA_line=$(grep 'tcdA' "${tox_input}")
        tcdA_matchtype=$(echo "${tcdA_line}" | cut -d$'\t' -f5)
        if [[ "${tcdA_matchtype}" = *"Trunc"* ]]; then
            tcdA="Trunc"
        else
            tcdA=1
        fi
        tcdA_allele=$(extract_gama_elements 3 "${tcdA_line}")
        tcdA_IDA=$(echo "${tcdA_line}" | cut -d$'\t' -f10)
        tcdA_IDN=$(echo "${tcdA_line}" | cut -d$'\t' -f11)
        tcdA_length=$(echo "${tcdA_line}" | cut -d$'\t' -f12)
        tcdA_set="${tcdA}\t${tcdA_count}\t${tcdA_allele}\t${subtype_A}\t${tcdA_IDN}|${tcdA_IDA}|${tcdA_length}"
    else
        tcdA=1
        readarray -t tcdA_arr < <(grep "tcdA" "${tox_input}")
        tcdA_count_array=${#tcdA_arr[@]}
        if [[ ${tcdA_count} != ${tcdA_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdA_stats=""
        tcdA_allele=""
        while [[ ${count} -lt ${tcdA_count_array} ]]; do
            tcdA_allele=$(extract_gama_elements 3 "${tcdA_arr[count]}")
            tcdA_IDA=$(echo "${tcdA_arr[count]}" | cut -d$'\t' -f10)
            tcdA_IDN=$(echo "${tcdA_arr[count]}" | cut -d$'\t' -f11)
            tcdA_length=$(echo "${tcdA_arr[count]}" | cut -d$'\t' -f12)
            tcdA_matchtype=$(echo "${tcdA_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${tcdA_matchtype}" = *"Trunc"* ]] || [[ ${tcdA} = "Trunc" ]]; then
                tcdA="Trunc"
            else
                tcdA=1
            fi
            if [[ "${tcdA_allele}" = "" ]]; then
                tcdA_allele=${tcdA_allele_current}
            else
                tcdA_allele="${tcdA_allele}-${tcdA_allele_current}"
            fi
            if [[ "${tcdA_stats}" == "" ]]; then
                tcdA_stats="${tcdA_IDN}|${tcdA_IDA}|${tcdA_length}"
            else
                tcdA_stats="${tcdA_stats}-${tcdA_IDN}|${tcdA_IDA}|${tcdA_length}"
            fi
            count=$(( count + 1 ))
        done
        tcdA_set="${tcdA}\t${tcdA_count_array}\t${tcdA_allele}\t${subtype_A}\t${tcdA_stats}"
    fi

    # Check and format for tcdB count
    tcdB_count=$(grep -c 'tcdB' "${tox_input}")
    echo "${tcdB_count}"
    if [[ "${tcdB_count}" -eq 0 ]]; then
        tcdB=0
        tcdB_set="0\t0\tNA\tNA\tNA|NA|NA"
    elif [[ "${tcdB_count}" -eq 1 ]]; then
        tcdB_line=$(grep 'tcdB' "${tox_input}")
        tcdB_matchtype=$(echo "${tcdB_line}" | cut -d$'\t' -f5)
        if [[ "${tcdB_matchtype}" = *"Trunc"* ]]; then
            tcdB="Trunc"
        else
            tcdB=1
        fi
        tcdB_allele=$(extract_gama_elements 3 "${tcdB_line}")
        tcdB_IDA=$(echo "${tcdB_line}" | cut -d$'\t' -f10)
        tcdB_IDN=$(echo "${tcdB_line}" | cut -d$'\t' -f11)
        tcdB_length=$(echo "${tcdB_line}" | cut -d$'\t' -f12)
        tcdB_set="${tcdB}\t${tcdB_count}\t${tcdB_allele}\t${subtype_B}\t${tcdB_IDN}|${tcdB_IDA}|${tcdB_length}"
    else
        tcdB=1
        readarray -t tcdB_arr < <(grep "tcdB" "${tox_input}")
        tcdB_count_array=${#tcdB_arr[@]}
        if [[ ${tcdB_count} != ${tcdB_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdB_stats=""
        tcdB_allele=""
        while [[ ${count} -lt ${tcdB_count_array} ]]; do
            echo "${count}"
            tcdB_allele=$(extract_gama_elements 3 "${tcdB_arr[count]}")
            tcdB_IDA=$(echo "${tcdB_arr[count]}" | cut -d$'\t' -f10)
            tcdB_IDN=$(echo "${tcdB_arr[count]}" | cut -d$'\t' -f11)
            tcdB_length=$(echo "${tcdB_arr[count]}" | cut -d$'\t' -f12)
            tcdB_matchtype=$(echo "${tcdB_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${tcdB_matchtype}" = *"Trunc"* ]] || [[ ${tcdB} = "Trunc" ]]; then
                tcdB="Trunc"
            else
                tcdB=1
            fi
            if [[ "${tcdB_allele}" = "" ]]; then
                tcdB_allele=${tcdB_allele_current}
            else
                tcdB_allele="${tcdB_allele}-${tcdB_allele_current}"
            fi
            if [[ "${tcdB_stats}" == "" ]]; then
                tcdB_stats="${tcdB_IDN}|${tcdB_IDA}|${tcdB_length}"
            else
                echo "extra"
                tcdB_stats="${tcdB_stats}-${tcdB_IDN}|${tcdB_IDA}|${tcdB_length}"
            fi
            echo "${count}-${tcdB}\t${tcdB_count_array}\t${tcdB_stats}"
            count=$(( count + 1 ))
        done
        tcdB_set="${tcdB}\t${tcdB_count_array}\t${tcdB_allele}\t${subtype_B}\t${tcdB_stats}"
    fi

    # Check and format for tcdC count
    tcdC_count=$(grep -c 'tcdC' "${tox_input}")
    if [[ "${tcdC_count}" -eq 0 ]]; then
        tcdC=0
        tcdC_set="0\t0\tNA\tNA|NA|NA"
    elif [[ "${tcdC_count}" -eq 1 ]]; then
        tcdC_line=$(grep 'tcdC' "${tox_input}")
        tcdC_matchtype=$(echo "${tcdC_line}" | cut -d$'\t' -f5)
        if [[ "${tcdC_matchtype}" = *"Trunc"* ]]; then
            tcdC="Trunc"
        else
            tcdC=1
        fi
        tcdC_allele=$(extract_gama_elements 3 "${tcdC_line}")
        tcdC_IDA=$(echo "${tcdC_line}" | cut -d$'\t' -f10)
        tcdC_IDN=$(echo "${tcdC_line}" | cut -d$'\t' -f11)
        tcdC_length=$(echo "${tcdC_line}" | cut -d$'\t' -f12)
        tcdC_set="${tcdC}\t${tcdC_count}\t${tcdC_allele}\t${tcdC_IDN}|${tcdC_IDA}|${tcdC_length}"
    else
        readarray -t tcdC_arr < <(grep "tcdC" "${tox_input}")
        tcdC_count_array=${#tcdC_arr[@]}
        if [[ ${tcdC_count} != ${tcdC_count_array} ]]; then 
            echo "wut?! ${tcdC_count}!=${tcdC_count_array}"
            echo -e "${tcdC_arr[0]}\n${tcdC_arr[1]}"
        fi
        count=0

        tcdC_allele=""
        tcdC_stats=""
        tcdC=""
        while [[ ${count} -lt ${tcdC_count_array} ]]; do
            tcdC_allele_current=$(extract_gama_elements 3 "${tcdC_arr[count]}")
            tcdC_IDA=$(echo "${tcdC_arr[count]}" | cut -d$'\t' -f10)
            tcdC_IDN=$(echo "${tcdC_arr[count]}" | cut -d$'\t' -f11)
            tcdC_length=$(echo "${tcdC_arr[count]}" | cut -d$'\t' -f12)
            tcdC_matchtype=$(echo "${tcdC_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${tcdC} = "" ]]; then
                if [[ "${tcdC_matchtype_temp}" = *"Trunc"* ]]; then
                    tcdC="Trunc"
                else
                    tcdC=1
                fi
            else
                if [[ "${tcdC_matchtype}" = *"Trunc"* ]]; then
                    tcdC="${tcdC}-Trunc"
                else
                    tcdC="${tcdC}-1"
                fi
            fi
            if [[ "${tcdC_allele}" = "" ]]; then
                tcdC_allele=${tcdC_allele_current}
            else
                tcdC_allele="${tcdC_allele}-${tcdC_allele_current}"
            fi
            if [[ "${tcdC_stats}" == "" ]]; then
                tcdC_stats="${tcdC_IDN}|${tcdC_IDA}|${tcdC_length}"
            else
                tcdC_stats="${tcdC_stats}-${tcdC_IDN}|${tcdC_IDA}|${tcdC_length}"
            fi
            count=$(( count + 1 ))
        done
        tcdC_set="${tcdC}\t${tcdC_count_array}\t${tcdC_allele}\t${tcdC_stats}"
    fi

    # Check and format for tcdD count
    tcdD_count=$(grep -c 'tcdD' "${tox_input}")
    if [[ "${tcdD_count}" -eq 0 ]]; then
        tcdD=0
        tcdD_set="0\t0\tNA|NA|NA"
    elif [[ "${tcdD_count}" -eq 1 ]]; then
        tcdD_line=$(grep 'tcdD' "${tox_input}")
        tcdD_matchtype=$(echo "${tcdD_line}" | cut -d$'\t' -f5)
        if [[ "${tcdD_matchtype}" = *"Trunc"* ]]; then
            tcdD="Trunc"
        else
            tcdD=1
        fi
        tcdD_IDA=$(echo "${tcdD_line}" | cut -d$'\t' -f10)
        tcdD_IDN=$(echo "${tcdD_line}" | cut -d$'\t' -f11)
        tcdD_length=$(echo "${tcdD_line}" | cut -d$'\t' -f12)
        tcdD_set="${tcdD}\t${tcdD_count}\t${tcdD_IDN}|${tcdD_IDA}|${tcdD_length}"
    else
        tcdD=1
        readarray -t tcdD_arr < <(grep "tcdD" "${tox_input}")
        tcdD_count_array=${#tcdD_arr[@]}
        if [[ ${tcdD_count} != ${tcdD_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdD_stats=""
        while [[ ${count} -lt ${tcdD_count_array} ]]; do
            tcdD_IDA=$(echo "${tcdD_arr[count]}" | cut -d$'\t' -f10)
            tcdD_IDN=$(echo "${tcdD_arr[count]}" | cut -d$'\t' -f11)
            tcdD_length=$(echo "${tcdD_arr[count]}" | cut -d$'\t' -f12)
            tcdD_matchtype=$(echo "${tcdD_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${tcdD_matchtype}" = *"Trunc"* ]] || [[ ${tcdD} = "Trunc" ]]; then
                tcdD="Trunc"
            else
                tcdD=1
            fi
            if [[ "${tcdD_stats}" == "" ]]; then
                tcdD_stats="${tcdD_IDN}|${tcdD_IDA}|${tcdD_length}"
            else
                tcdD_stats="${tcdD_stats}-${tcdD_IDN}|${tcdD_IDA}|${tcdD_length}"
            fi
            count=$(( count + 1 ))
        done
        tcdD_set="${tcdD}\t${tcdD_count_array}\t${tcdD_stats}"
    fi

    # Check and format for tcdE count
    tcdE_count=$(grep -c 'tcdE' "${tox_input}")
    if [[ "${tcdE_count}" -eq 0 ]]; then
        tcdE=0
        tcdE_set="0\t0\tNA|NA|NA"
    elif [[ "${tcdE_count}" -eq 1 ]]; then
        tcdE_line=$(grep 'tcdE' "${tox_input}")
        tcdE_matchtype=$(echo "${tcdE_line}" | cut -d$'\t' -f5)
        if [[ "${tcdE_matchtype}" = *"Trunc"* ]]; then
            tcdE="Trunc"
        else
            tcdE=1
        fi
        tcdE_IDA=$(echo "${tcdE_line}" | cut -d$'\t' -f10)
        tcdE_IDN=$(echo "${tcdE_line}" | cut -d$'\t' -f11)
        tcdE_length=$(echo "${tcdE_line}" | cut -d$'\t' -f12)
        tcdE_set="${tcdE}\t${tcdE_count}\t${tcdE_IDN}|${tcdE_IDA}|${tcdE_length}"
    else
        tcdE=1
        readarray -t tcdE_arr < <(grep "tcdE" "${tox_input}")
        tcdE_count_array=${#tcdE_arr[@]}
        if [[ ${tcdE_count} != ${tcdE_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        tcdE_stats=""
        while [[ ${count} -lt ${tcdE_count_array} ]]; do
            tcdE_IDA=$(echo "${tcdE_arr[count]}" | cut -d$'\t' -f10)
            tcdE_IDN=$(echo "${tcdE_arr[count]}" | cut -d$'\t' -f11)
            tcdE_length=$(echo "${tcdE_arr[count]}" | cut -d$'\t' -f12)
            tcdE_matchtype=$(echo "${tcdE_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${tcdE_matchtype}" = *"Trunc"* ]] || [[ ${tcdE} = "Trunc" ]]; then
                tcdE="Trunc"
            else
                tcdE=1
            fi
            if [[ "${tcdE_stats}" == "" ]]; then
                tcdE_stats="${tcdE_IDN}|${tcdE_IDA}|${tcdE_length}"
            else
                tcdE_stats="${tcdE_stats}-${tcdE_IDN}|${tcdE_IDA}|${tcdE_length}"
            fi
            count=$(( count + 1 ))
        done
        tcdE_set="${tcdE}\t${tcdE_count_array}\t${tcdE_stats}"
    fi

    # Check and format for cdtA count
    cdtA_count=$(grep -c 'cdtA_' "${tox_input}")
    if [[ "${cdtA_count}" -eq 0 ]]; then
        cdtA=0
        cdtA_set="0\t0\tNA|NA|NA"
    elif [[ "${cdtA_count}" -eq 1 ]]; then
        cdtA_line=$(grep 'cdtA_' "${tox_input}")
        cdtA_matchtype=$(echo "${cdtA_line}" | cut -d$'\t' -f5)
        if [[ "${cdtA_matchtype}" = *"Trunc"* ]]; then
            cdtA="Trunc"
        else
            cdtA=1
        fi
        cdtA_IDA=$(echo "${cdtA_line}" | cut -d$'\t' -f10)
        cdtA_IDN=$(echo "${cdtA_line}" | cut -d$'\t' -f11)
        cdtA_length=$(echo "${cdtA_line}" | cut -d$'\t' -f12)
        cdtA_set="${cdtA}\t${cdtA_count}\t${cdtA_IDN}|${cdtA_IDA}|${cdtA_length}"
    else
        readarray -t cdtA_arr < <(grep "cdtA" "${tox_input}")
        cdtA_count_array=${#cdtA_arr[@]}
        if [[ ${cdtA_count} != ${cdtA_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtA_stats=""
        cdtA=""
        while [[ ${count} -lt ${cdtA_count_array} ]]; do
            cdtA_IDA=$(echo "${cdtA_arr[count]}" | cut -d$'\t' -f10)
            cdtA_IDN=$(echo "${cdtA_arr[count]}" | cut -d$'\t' -f11)
            cdtA_length=$(echo "${cdtA_arr[count]}" | cut -d$'\t' -f12)
            cdtA_matchtype=$(echo "${cdtA_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${cdtA} = "" ]]; then
                if [[ "${cdtA_matchtype_temp}" = *"Trunc"* ]]; then
                    cdtA="Trunc"
                else
                    cdtA=1
                fi
            else
                if [[ "${cdtA_matchtype}" = *"Trunc"* ]]; then
                    cdtA="${cdtA}-Trunc"
                else
                    cdtA="${cdtA}-1"
                fi
            fi
            if [[ "${cdtA_stats}" == "" ]]; then
                cdtA_stats="${cdtA_IDN}|${cdtA_IDA}|${cdtA_length}"
            else
                cdtA_stats="${cdtA_stats}-${cdtA_IDN}|${cdtA_IDA}|${cdtA_length}"
            fi
            count=$(( count + 1 ))
        done
        cdtA_set="${cdtA}\t${cdtA_count_array}\t${cdtA_stats}"
    fi

    # Check and format for cdtB count
    cdtB_count=$(grep -c 'cdtB_' "${tox_input}")
    if [[ "${cdtB_count}" -eq 0 ]]; then
        cdtB=0
        cdtB_set="0\t0\tNA|NA|NA"
    elif [[ "${cdtB_count}" -eq 1 ]]; then
        cdtB_line=$(grep 'cdtB_' "${tox_input}")
        cdtB_matchtype=$(echo "${cdtB_line}" | cut -d$'\t' -f5)
        if [[ "${cdtB_matchtype}" = *"Trunc"* ]]; then
            cdtB="Trunc"
        else
            cdtB=1
        fi
        cdtB_IDA=$(echo "${cdtB_line}" | cut -d$'\t' -f10)
        cdtB_IDN=$(echo "${cdtB_line}" | cut -d$'\t' -f11)
        cdtB_length=$(echo "${cdtB_line}" | cut -d$'\t' -f12)
        cdtB_set="${cdtB}\t${cdtB_count}\t${cdtB_IDN}|${cdtB_IDA}|${cdtB_length}"
    else
        readarray -t cdtB_arr < <(grep "cdtB" "${tox_input}")
        cdtB_count_array=${#cdtB_arr[@]}
        if [[ ${cdtB_count} != ${cdtB_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtB_stats=""
        cdtB=""
        while [[ ${count} -lt ${cdtB_count_array} ]]; do
            cdtB_IDA=$(echo "${cdtB_arr[count]}" | cut -d$'\t' -f10)
            cdtB_IDN=$(echo "${cdtB_arr[count]}" | cut -d$'\t' -f11)
            cdtB_length=$(echo "${cdtB_arr[count]}" | cut -d$'\t' -f12)
            cdtB_matchtype=$(echo "${cdtB_arr[count]}" | cut -d$'\t' -f5)
            if [[ ${cdtB} = "" ]]; then
                if [[ "${cdtB_matchtype_temp}" = *"Trunc"* ]]; then
                    cdtB="Trunc"
                else
                    cdtB=1
                fi
            else
                if [[ "${cdtB_matchtype}" = *"Trunc"* ]]; then
                    cdtB="${cdtB}-Trunc"
                else
                     cdtB="${cdtB}-1"
                fi
            fi
            if [[ "${cdtB_stats}" == "" ]]; then
                cdtB_stats="${cdtB_IDN}|${cdtB_IDA}|${cdtB_length}"
            else
                cdtB_stats="${cdtB_stats}-${cdtB_IDN}|${cdtB_IDA}|${cdtB_length}"
            fi
            count=$(( count + 1 ))
        done
        cdtB_set="${cdtB}\t${cdtB_count_array}\t${cdtB_stats}"
    fi

    # Check and format for cdtR count
    cdtR_count=$(grep -c 'cdtR' "${tox_input}")
    if [[ "${cdtR_count}" -eq 0 ]]; then
        cdtA=0
        cdtR_set="0\t0\tNA\tNA|NA|NA"
    elif [[ "${cdtR_count}" -eq 1 ]]; then
        cdtR_line=$(grep 'cdtR' "${tox_input}")
        cdtR_matchtype=$(echo "${cdtR_line}" | cut -d$'\t' -f5)
        if [[ "${cdtR_matchtype}" = *"Trunc"* ]]; then
            cdtA="Trunc"
        else
            cdtA=1
        fi
        cdtR_allele=$(extract_gama_elements 3 "${cdtR_line}")
        cdtR_IDA=$(echo "${cdtR_line}" | cut -d$'\t' -f10)
        cdtR_IDN=$(echo "${cdtR_line}" | cut -d$'\t' -f11)
        cdtR_length=$(echo "${cdtR_line}" | cut -d$'\t' -f12)
        cdtR_set="${cdtA}\t${cdtR_count}\t${cdtR_allele}\t${cdtR_IDN}|${cdtR_IDA}|${cdtR_length}"
    else
        cdtA=1
        readarray -t cdtR_arr < <(grep "cdtR" "${tox_input}")
        cdtR_count_array=${#cdtR_arr[@]}
        if [[ ${cdtR_count} != ${cdtR_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtR_stats=""
        cdtR_allele=""
        while [[ ${count} -lt ${cdtR_count_array} ]]; do
            cdtR_allele=$(extract_gama_elements 3 "${cdtR_arr[count]}")
            cdtR_IDA=$(echo "${cdtR_arr[count]}" | cut -d$'\t' -f10)
            cdtR_IDN=$(echo "${cdtR_arr[count]}" | cut -d$'\t' -f11)
            cdtR_length=$(echo "${cdtR_arr[count]}" | cut -d$'\t' -f12)
            cdtR_matchtype=$(echo "${cdtR_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${cdtR_matchtype}" = *"Trunc"* ]] || [[ ${cdtR} = "Trunc" ]]; then
                cdtA="Trunc"
            else
                cdtA=1
            fi
            if [[ "${cdtR_allele}" = "" ]]; then
                cdtR_allele=${cdtR_allele_current}
            else
                cdtR_allele="${cdtR_allele}-${cdtR_allele_current}"
            fi
            if [[ "${cdtR_stats}" == "" ]]; then
                cdtR_stats="${cdtR_IDN}|${cdtR_IDA}|${cdtR_length}"
            else
                cdtR_stats="${cdtR_stats}-${cdtR_IDN}|${cdtR_IDA}|${cdtR_length}"
            fi
            count=$(( count + 1 ))
        done
        cdtR_set="${cdtA}\t${cdtR_count_array}\t${cdtR_stats}"
    fi

    # Check and format for cdtAB1 count
    cdtAB1_count=$(grep -c 'cdtAB1' "${tox_input}")
    if [[ "${cdtAB1_count}" -eq 0 ]]; then
        cdtAB1=0
        cdtAB1_set="0\t0\tNA|NA|NA"
    elif [[ "${cdtAB1_count}" -eq 1 ]]; then
        cdtAB1_line=$(grep 'cdtAB1' "${tox_input}")
        cdtAB1_matchtype=$(echo "${cdtAB1_line}" | cut -d$'\t' -f5)
        if [[ "${cdtAB1_matchtype}" = *"Trunc"* ]]; then
            cdtAB1="Trunc"
        else
            cdtAB1=1
        fi
        cdtAB1_IDA=$(echo "${cdtAB1_line}" | cut -d$'\t' -f10)
        cdtAB1_IDN=$(echo "${cdtAB1_line}" | cut -d$'\t' -f11)
        cdtAB1_length=$(echo "${cdtAB1_line}" | cut -d$'\t' -f12)
        cdtAB1_set="${cdtAB1}\t${cdtAB1_count}\t${cdtAB1_IDN}|${cdtAB1_IDA}|${cdtAB1_length}"
    else
        cdtAB1=1
        readarray -t cdtAB1_arr < <(grep "cdtAB1" "${tox_input}")
        cdtAB1_count_array=${#cdtAB1_arr[@]}
        if [[ ${cdtAB1_count} != ${cdtAB1_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtAB1_stats=""
        while [[ ${count} -lt ${cdtAB1_count_array} ]]; do
            cdtAB1_IDA=$(echo "${cdtAB1_arr[count]}" | cut -d$'\t' -f10)
            cdtAB1_IDN=$(echo "${cdtAB1_arr[count]}" | cut -d$'\t' -f11)
            cdtAB1_length=$(echo "${cdtAB1_arr[count]}" | cut -d$'\t' -f12)
            cdtAB1_matchtype=$(echo "${cdtAB1_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${cdtAB1_matchtype}" = *"Trunc"* ]] || [[ ${cdtAB1} = "Trunc" ]]; then
                cdtAB1="Trunc"
            else
                cdtAB1=1
            fi
            if [[ "${cdtAB1_stats}" == "" ]]; then
                cdtAB1_stats="${cdtAB1_IDN}|${cdtAB1_IDA}|${cdtAB1_length}"
            else
                cdtAB1_stats="${cdtAB1_stats}-${cdtAB1_IDN}|${cdtAB1_IDA}|${cdtAB1_length}"
            fi
            count=$(( count + 1 ))
        done
        cdtAB1_set="${cdtAB1}\t${cdtAB1_count_array}\t${cdtAB1_stats}"
    fi

    # Check and format for cdtAB2 count
    cdtAB2_count=$(grep -c 'cdtAB2' "${tox_input}")
    if [[ "${cdtAB2_count}" -eq 0 ]]; then
        cdtAB2=0
        cdtAB2_set="0\t0\tNA|NA|NA"
    elif [[ "${cdtAB2_count}" -eq 1 ]]; then
        cdtAB2_line=$(grep 'cdtAB2' "${tox_input}")
        cdtAB2_matchtype=$(echo "${cdtAB2_line}" | cut -d$'\t' -f5)
        if [[ "${cdtAB2_matchtype}" = *"Trunc"* ]]; then
            cdtAB2="Trunc"
        else
            cdtAB2=1
        fi
        cdtAB2_IDA=$(echo "${cdtAB2_line}" | cut -d$'\t' -f10)
        cdtAB2_IDN=$(echo "${cdtAB2_line}" | cut -d$'\t' -f11)
        cdtAB2_length=$(echo "${cdtAB2_line}" | cut -d$'\t' -f12)
        cdtAB2_set="${cdtAB2}\t${cdtAB2_count}\t${cdtAB2_IDN}|${cdtAB2_IDA}|${cdtAB2_length}"
    else
        cdtAB2=1
        readarray -t cdtAB2_arr < <(grep "cdtAB2" "${tox_input}")
        cdtAB2_count_array=${#cdtAB2_arr[@]}
        if [[ ${cdtAB2_count} != ${cdtAB2_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        cdtAB2_stats=""
        while [[ ${count} -lt ${cdtAB2_count_array} ]]; do
            cdtAB2_IDA=$(echo "${cdtAB2_arr[count]}" | cut -d$'\t' -f10)
            cdtAB2_IDN=$(echo "${cdtAB2_arr[count]}" | cut -d$'\t' -f11)
            cdtAB2_length=$(echo "${cdtAB2_arr[count]}" | cut -d$'\t' -f12)
            cdtAB2_matchtype=$(echo "${cdtAB2_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${cdtAB2_matchtype}" = *"Trunc"* ]] || [[ ${cdtAB2} = "Trunc" ]]; then
                cdtAB2="Trunc"
            else
                cdtAB2=1
            fi
            if [[ "${cdtAB2_stats}" == "" ]]; then
                cdtAB2_stats="${cdtAB2_IDN}|${cdtAB2_IDA}|${cdtAB2_length}"
            else
                cdtAB2_stats="${cdtAB2_stats}-${cdtAB2_IDN}|${cdtAB2_IDA}|${cdtAB2_length}"
            fi
            count=$(( count + 1 ))
        done
        cdtAB2_set="${cdtAB2}\t${cdtAB2_count_array}\t${cdtAB2_stats}"
    fi

        # Check and format for nontox count
    nontox_count=$(grep -c 'cdtNonTox' "${tox_input}")
    if [[ "${nontox_count}" -eq 0 ]]; then
        nontox=0
        nontox_set="0\t0\tNA|NA|NA"
    elif [[ "${nontox_count}" -eq 1 ]]; then
        nontox_line=$(grep 'cdtNonTox' "${tox_input}")
        nontox_matchtype=$(echo "${nontox_line}" | cut -d$'\t' -f5)
        if [[ "${nontox_matchtype}" = *"Trunc"* ]]; then
            nontox="Trunc"
        else
            nontox=1
        fi
        nontox_IDA=$(echo "${nontox_line}" | cut -d$'\t' -f10)
        nontox_IDN=$(echo "${nontox_line}" | cut -d$'\t' -f11)
        nontox_length=$(echo "${nontox_line}" | cut -d$'\t' -f12)
        nontox_set="${nontox}\t${nontox_count}\t${nontox_IDN}|${nontox_IDA}|${nontox_length}"
    else
        nontox=1
        readarray -t nontox_arr < <(grep "cdtnontox" "${tox_input}")
        nontox_count_array=${#nontox_arr[@]}
        if [[ ${nontox_count} != ${nontox_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        nontox_stats=""
        while [[ ${count} -lt ${nontox_count_array} ]]; do
            nontox_IDA=$(echo "${nontox_arr[count]}" | cut -d$'\t' -f10)
            nontox_IDN=$(echo "${nontox_arr[count]}" | cut -d$'\t' -f11)
            nontox_length=$(echo "${nontox_arr[count]}" | cut -d$'\t' -f12)
            nontox_matchtype=$(echo "${nontox_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${nontox_matchtype}" = *"Trunc"* ]] || [[ ${nontox} = "Trunc" ]]; then
                nontox="Trunc"
            else
                nontox=1
            fi
            if [[ "${nontox_stats}" == "" ]]; then
                nontox_stats="${nontox_IDN}|${nontox_IDA}|${nontox_length}"
            else
                nontox_stats="${nontox_stats}-${nontox_IDN}|${nontox_IDA}|${nontox_length}"
            fi
            count=$(( count + 1 ))
        done
        nontox_set="${nontox}\t${nontox_count_array}\t${nontox_stats}"
    fi

 # Check and format for gyrA count, these are ALL point mutations and must be handled differently
    gyrA_count=$(grep -c 'gyrA' "${tox_input}")
    if [[ "${gyrA_count}" -eq 0 ]]; then
        gyrA=0
        gyrA_set="0\t0\tNA\tN/A\tNA|NA|NA"
    elif [[ "${gyrA_count}" -eq 1 ]]; then
        gyrA_line=$(grep 'gyrA' "${tox_input}")
        gyrA_matchtype=$(echo "${gyrA_line}" | cut -d$'\t' -f5)
        gyrA_AA_codon_ID=$(echo "${gyrA_line}" | cut -d$'\t' -f10)
        gyrA_AA_codon_length=$(echo "${gyrA_line}" | cut -d$'\t' -f12)
        gyrA_allele=$(extract_gama_elements 3 "${gyrA_line}")
        if [[ "${gyrA_matchtype}" = "Native" ]]; then
            if [[ "${gyrA_AA_codon_ID}" -eq 1 ]] && [[ "${gyrA_AA_codon_length}" -eq 1 ]]; then
                gyrA_additional_info=""
            # What to do if it is not a full length match
            else
                gyrA_allele="${gyrA_allele}*"
                gyrA_additional_info="Not Full Length"
            fi
        elif [[ "${gyrA_matchtype}" = *"Trunc"* ]]; then
                gyrA_allele_current="${gyrA_allele_current}*"
                gyrA_additional_info_current=$(echo "${gyrA_line}" | cut -d$'\t' -f6)
                gyrA="Trunc"
        elif [[ "${gyrA_matchtype}" = "Mutant" ]] || [[ "${gyrA_matchtype}" = "Indel" ]]; then
            gyrA_allele="${gyrA_allele}*"
            gyrA_additional_info=$(echo "${gyrA_line}" | cut -d$'\t' -f6)
        # Hit is some other kind of match, i.e. Mutatnt, Indel, Trunc
        else
            echo "UNKNOWN MATCHTYPE, PLEASE ADD TO CHECK (CENTAR CONSOLIDATOR.sh)"
        fi
        if [[ "${gyrA_additional_info: -1}" = "," ]]; then
                gyrA_additional_info="${gyrA_additional_info:0:-1}"
        fi
        
        gyrA_IDA=$(echo "${gyrA_line}" | cut -d$'\t' -f10)
        gyrA_IDN=$(echo "${gyrA_line}" | cut -d$'\t' -f11)
        gyrA_length=$(echo "${gyrA_line}" | cut -d$'\t' -f12)
        gyrA_set="1\t1\t${gyrA_allele}\t${gyrA_additional_info}\t${gyrA_IDN}|${gyrA_IDA}|${gyrA_length}"
    else
        readarray -t gyrA_arr < <(grep "gyrA" "${tox_input}")
        gyrA_count_array=${#gyrA_arr[@]}
        if [[ ${gyrA_count} != ${gyrA_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        gyrA_allele="NA"
        gyrA_additional_info=""
        gyrA_stats=""
        gyrA=""
        while [[ ${count} -lt ${gyrA_count_array} ]]; do
            gyrA_allele_current=$(extract_gama_elements 3 "${gyrA_arr[count]}")
            gyrA_IDA=$(echo "${gyrA_arr[count]}" | cut -d$'\t' -f10)
            gyrA_IDN=$(echo "${gyrA_arr[count]}" | cut -d$'\t' -f11)
            gyrA_length=$(echo "${gyrA_arr[count]}" | cut -d$'\t' -f12)
            gyrA_matchtype=$(echo "${gyrA_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${gyrA_matchtype}" = "Native" ]]; then
                if [[ "${gyrA_AA_codon_ID}" -eq 1 ]] && [[ "${gyrA_AA_codon_length}" -eq 1 ]]; then
                    gyrA_additional_info_current=""
                # What to do if it is not a full length match
                else
                    gyrA_allele_current="${gyrA_allele_current}*"
                    gyrA_additional_info_current="Not Full Length"
                fi
            elif [[ "${gyrA_matchtype}" = *"Trunc"* ]]; then
                gyrA_allele_current="${gyrA_allele_current}*"
                gyrA_additional_info_current=$(echo "${gyrA_line}" | cut -d$'\t' -f6)
            elif [[ "${gyrA_matchtype}" = "Mutant" ]] || [[ "${gyrA_matchtype}" = "Indel" ]]; then
                gyrA_allele_current="${gyrA_allele_current}*"
                gyrA_additional_info_current=$(echo "${gyrA_line}" | cut -d$'\t' -f6)
            # Hit is some other kind of match, i.e. Mutatnt, Indel, Trunc
            else
                echo "UNKNOWN MATCHTYPE, PLEASE ADD TO CHECK (CENTAR CONSOLIDATOR.sh)"
            fi

            if [[ "${gyrA_additional_info_current: -1}" = "," ]]; then
                gyrA_additional_info_current="${gyrA_additional_info_current:0:-1}"
            fi


            if [[ ${gyrA} = "" ]]; then
                if [[ "${gyrA_matchtype_temp}" = *"Trunc"* ]]; then
                    gyrA="Trunc"
                else
                    gyrA=1
                fi
            else
                if [[ "${gyrA_matchtype}" = *"Trunc"* ]]; then
                    gyrA="${gyrA}-Trunc"
                else
                    gyrA="${gyrA}-1"
                fi
            fi
            if [[ "${gyrA_allele}" = "" ]]; then
                gyrA_allele=${gyrA_allele_current}
            else
                gyrA_allele="${gyrA_allele}-${gyrA_allele_current}"
            fi
            if [[ "${gyrA_stats}" == "" ]]; then
                gyrA_stats="${gyrA_IDN}|${gyrA_IDA}|${gyrA_length}"
            else
                gyrA_stats="${gyrA_stats}-${gyrA_IDN}|${gyrA_IDA}|${gyrA_length}"
            fi
            count=$(( count + 1 ))
        done
        gyrA_set="${gyrA}\t${gyrA_count_array}\t${gyrA_allele}\t${gyrA_additional_info}\t${gyrA_stats}"
    fi

    # Check and format for gyrB count, these are ALL point mutations and must be handled differently
    gyrB_count=$(grep -c 'gyrB' "${tox_input}")
    if [[ "${gyrB_count}" -eq 0 ]]; then
        gyrB=0
        gyrB_set="0\t0\tNA\tN/A\tNA|NA|NA"
    elif [[ "${gyrB_count}" -eq 1 ]]; then
        gyrB_line=$(grep 'gyrB' "${tox_input}")
        gyrB_matchtype=$(echo "${gyrB_line}" | cut -d$'\t' -f5)
        gyrB_AA_codon_ID=$(echo "${gyrB_line}" | cut -d$'\t' -f10)
        gyrB_AA_codon_length=$(echo "${gyrB_line}" | cut -d$'\t' -f12)
        gyrB_allele=$(extract_gama_elements 3 "${gyrB_line}")
        if [[ "${gyrB_matchtype}" = "Native" ]]; then
            if [[ "${gyrB_AA_codon_ID}" -eq 1 ]] && [[ "${gyrB_AA_codon_length}" -eq 1 ]]; then
                gyrB_additional_info=""
            # What to do if it is not a full length match
            else
                gyrB_allele="${gyrB_allele}*"
                gyrB_additional_info="Not Full Length"
            fi
        elif [[ "${gyrB_matchtype}" = *"Trunc"* ]]; then
                gyrB_allele_current="${gyrB_allele_current}*"
                gyrB_additional_info_current=$(echo "${gyrB_line}" | cut -d$'\t' -f6)
                gyrB="Trunc"
        elif [[ "${gyrB_matchtype}" = "Mutant" ]] || [[ "${gyrB_matchtype}" = "Indel" ]]; then
            gyrB_allele="${gyrB_allele}*"
            gyrB_additional_info=$(echo "${gyrB_line}" | cut -d$'\t' -f6)
        # Hit is some other kind of match, i.e. Mutatnt, Indel, Trunc
        else
            echo "UNKNOWN MATCHTYPE, PLEASE ADD TO CHECK (CENTAR CONSOLIDATOR.sh)"
        fi
        
        if [[ "${gyrB_additional_info: -1}" = "," ]]; then
                gyrB_additional_info="${gyrB_additional_info:0:-1}"
        fi

        gyrB_IDA=$(echo "${gyrB_line}" | cut -d$'\t' -f10)
        gyrB_IDN=$(echo "${gyrB_line}" | cut -d$'\t' -f11)
        gyrB_length=$(echo "${gyrB_line}" | cut -d$'\t' -f12)
        gyrB_set="1\t1\t${gyrB_allele}\t${gyrB_additional_info}\t${gyrB_IDN}|${gyrB_IDA}|${gyrB_length}"
    else
        readarray -t gyrB_arr < <(grep "gyrB" "${tox_input}")
        gyrB_count_array=${#gyrB_arr[@]}
        if [[ ${gyrB_count} != ${gyrB_count_array} ]]; then 
            echo "wut?!"
        fi
        count=0

        gyrB_allele="NA"
        gyrB_additional_info=""
        gyrB_stats=""
        gyrB=""
        while [[ ${count} -lt ${gyrB_count_array} ]]; do
            gyrB_allele_current=$(extract_gama_elements 3 "${gyrB_arr[count]}")
            gyrB_IDA=$(echo "${gyrB_arr[count]}" | cut -d$'\t' -f10)
            gyrB_IDN=$(echo "${gyrB_arr[count]}" | cut -d$'\t' -f11)
            gyrB_length=$(echo "${gyrB_arr[count]}" | cut -d$'\t' -f12)
            gyrB_matchtype=$(echo "${gyrB_arr[count]}" | cut -d$'\t' -f5)
            if [[ "${gyrB_matchtype}" = "Native" ]]; then
                if [[ "${gyrB_AA_codon_ID}" -eq 1 ]] && [[ "${gyrB_AA_codon_length}" -eq 1 ]]; then
                    gyrB_additional_info_current=""
                # What to do if it is not a full length match
                else
                    gyrB_allele_current="${gyrB_allele_current}*"
                    gyrB_additional_info_current="Not Full Length"
                fi
            elif [[ "${gyrB_matchtype}" = *"Trunc"* ]]; then
                gyrB_allele_current="${gyrB_allele_current}*"
                gyrB_additional_info_current=$(echo "${gyrB_line}" | cut -d$'\t' -f6)
            elif [[ "${gyrB_matchtype}" = "Mutant" ]] || [[ "${gyrB_matchtype}" = "Indel" ]]; then
                gyrB_allele_current="${gyrB_allele_current}*"
                gyrB_additional_info_current=$(echo "${gyrB_line}" | cut -d$'\t' -f6)
            # Hit is some other kind of match, i.e. Mutatnt, Indel, Trunc
            else
                echo "UNKNOWN MATCHTYPE, PLEASE ADD TO CHECK (CENTAR CONSOLIDATOR.sh)"
            fi

            if [[ "${gyrB_additional_info_current: -1}" = "," ]]; then
                gyrB_additional_info_current="${gyrB_additional_info_current:0:-1}"
            fi

            if [[ ${gyrB} = "" ]]; then
                if [[ "${gyrB_matchtype_temp}" = *"Trunc"* ]]; then
                    gyrB="Trunc"
                else
                    gyrB=1
                fi
            else
                if [[ "${gyrB_matchtype}" = *"Trunc"* ]]; then
                    gyrB="${gyrB}-Trunc"
                else
                    gyrB="${gyrB}-1"
                fi
            fi
            if [[ "${gyrB_allele}" = "" ]]; then
                gyrB_allele=${gyrB_allele_current}
            else
                gyrB_allele="${gyrB_allele}-${gyrB_allele_current}"
            fi
            if [[ "${gyrB_stats}" == "" ]]; then
                gyrB_stats="${gyrB_IDN}|${gyrB_IDA}|${gyrB_length}"
            else
                gyrB_stats="${gyrB_stats}-${gyrB_IDN}|${gyrB_IDA}|${gyrB_length}"
            fi
            count=$(( count + 1 ))
        done
        gyrB_set="${gyrB}\t${gyrB_count_array}\t${gyrB_allele}\t${gyrB_additional_info}\t${gyrB_stats}"
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
    gyrA_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file"
    gyrB_set="No_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file\tNo_Tox_file"
fi

if [[ -f "${clade_input}" ]]; then
    clade_lines=$(cat "${clade_input}" | wc -l)
    if [[ "${clade_lines}" -eq 2 ]]; then
        clade=$(tail -n1 "${clade_input}" | cut -d$'\t' -f2)
    else
        clade="Clade_file_incorrect"
    fi
else
    echo "No clade/mlst file, carry one"
    clade="No_clade/MLST_file"
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
plasmids:${plasmids}
"

# Loop through the genes in the list and format to match desired output style
if [[ ! -f "${output}" ]]; then
    echo -e "isolate_ID\tMLST Clade\tDiffbase_Toxinotype\ttcdA_presence\ttcdA_occurences\ttcdA_Variant\tDiffbase_Toxin-A_sub-type\ttcdA Confidence (Coverage_NT|Coverage_AA|Length)\ttcdB_prsence\ttcdB_occurences\ttcdB_Variant\tDiffbase_Toxin-B_sub-type\ttcdB Confidence (Coverage_NT|Coverage_AA|Length)\ttcdC_presence\ttcdC_occurences\ttcdC_Variant\ttcdC Confidence (Coverage_NT|Coverage_AA|Length)\ttcdR_presence\ttcdR_occurences\ttcdR Confidence (Coverage_NT|Coverage_AA|Length)\ttcdE_presence\ttcdE_occurences\ttcdE Confidence (Coverage_NT|Coverage_AA|Length)\tcdtA_presence\tcdtA_occurences\tcdtA Confidence (Coverage_NT|Coverage_AA|Length)\tcdtB_presence\tcdtB_occurences\tcdtB Confidence (Coverage_NT|Coverage_AA|Length)\tcdtR_presence\tcdtR_occurences\tcdtR_Variant\tcdtR Confidence (Coverage_NT|Coverage_AA|Length)\tcdtAB1_presence\tcdtAB1_occurences\tcdtAB1 Confidence (Coverage_NT|Coverage_AA|Length)\tcdtAB2_presence\tcdtAB2_occurences\tcdtAB2 Confidence (Coverage_NT|Coverage_AA|Length)\tnon-tox_presence\tnon-tox_occurences\tnon-tox Confidence (Coverage_NT|Coverage_AA|Length)\tgyrA_presence\tgyrA_occurences\tgyrA_Variant\tgyrA_additional_info\tgyrA Confidence (Coverage_NT|Coverage_AA|Length)\tgyrB_presence\tgyrB_occurences\tgyrB_Variant\tgyrB_additional_info\tgyrB Confidence (Coverage_NT|Coverage_AA|Length)\tInferred RT\tProbability\tPlasmid Info\t${other_AR_header}" > "${output}"
fi
echo -e "${sample_name}\t${clade}\t${toxinotype}\t${tcdA_set}\t${tcdB_set}\t${tcdC_set}\t${tcdD_set}\t${tcdE_set}\t${cdtA_set}\t${cdtB_set}\t${cdtR_set}\t${cdtAB1_set}\t${cdtAB2_set}\t${nontox_set}\t${gyrA_set}\t${gyrB_set}\t${ML_RT}\t${plasmids}\t${other_AR}"  >> "${output}"
