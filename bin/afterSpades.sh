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
# Created by  (@cdc.gov)
#

fasta=$(find *.spades.log)
prefix=$(basename $fasta .spades.log)
if [ -f scaffolds.fasta ]; then
    mv scaffolds.fasta ${prefix}.scaffolds.fa
    gzip -n ${prefix}.scaffolds.fa
    spades_complete=scaffolds_created
    echo ,$spades_complete | tr -d "\n" >> ${prefix}_spades_outcome.csv
else
    spades_complete=no_scaffolds
    echo ,$spades_complete | tr -d "\n" >> ${prefix}_spades_outcome.csv
fi
if [ -f contigs.fasta ]; then
    mv contigs.fasta ${prefix}.contigs.fa
    gzip -n ${prefix}.contigs.fa
    spades_complete=contigs_created
    echo ,$spades_complete | tr -d "\n" >> ${prefix}_spades_outcome.csv
else
    spades_complete=no_contigs
    echo ,$spades_complete | tr -d "\n" >> ${prefix}_spades_outcome.csv
fi
if [ -f assembly_graph_with_scaffolds.gfa ]; then
    mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
    gzip -n ${prefix}.assembly.gfa
fi