#!/bin/bash

awk 'BEGIN {OFS=","; print "query,reference, ANI%,total,aligned"} {print $1, $2, $3, $4,$5}' *.txt > plasmidANI.csv
