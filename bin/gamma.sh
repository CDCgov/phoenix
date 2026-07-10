#!/bin/bash
#1 db
#2 output
#for i in contig_*.fasta 
#do
GAMMA.py $i *renamed.scaffolds.fa $1 $2_$i
##done