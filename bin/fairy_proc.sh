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
