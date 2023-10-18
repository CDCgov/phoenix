#!/usr/bin/env python3

import sys
import os
import json
import argparse
from decimal import *
getcontext().prec = 4
from argparse import ArgumentParser

##Makes pre- and post-filtering QC outputs
##Usage: >python3 FastP_QC.py paired_fastp.json single_fastp_json Isolate_Name
##Written by Rich Stanton (njr5@cdc.gov) and Nick Vlachos (nvx4@cdc.gov)

def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="FastP_QC.py", description="""Script collects stats from fastp jsons.""")
	parser.add_argument("-t", "--trimmed_json", dest="trimmed_json", action="store", required=True, help="Json from fastp output on trimmed reads")
	parser.add_argument("-s", "--single_json", dest="single_json", action="store", required=True, help="Json from fastp output on single reads.")
	parser.add_argument("-n", "--name", dest="name", action="store", required=True, help="Sample name")
	args = parser.parse_args()
	return args

def FastP_QC_before(input_json, output_file, name):
    """Makes a QC output file from an input FastP json"""
    Out = open(output_file, 'w')
    Out.write('Name\tR1[reads]\tR1[bp]\tR2[reads]\tR2[bp]\tQ20_Total_[bp]\tQ30_Total_[bp]\tQ20_R1_[bp]\tQ20_R2_[bp]\tQ20_R1_[%]\tQ20_R2_[%]\tQ30_R1_[bp]\tQ30_R2_[bp]\tQ30_R1_[%]\tQ30_R2_[%]\tTotal_Sequenced_[bp]\tTotal_Sequenced_[reads]\n')
    f = open(input_json)
    data = json.load(f)
    f.close()
    paired_total_info = data['summary']['before_filtering']
    Q20_Total = str(paired_total_info['q20_bases'])
    Q30_Total = str(paired_total_info['q30_bases'])
    Total_Sequenced_bp = str(paired_total_info['total_bases'])
    Total_Sequenced_reads = str(paired_total_info['total_reads'])
    Info1 = data['read1_before_filtering']
    raw_R1_reads = str(Info1['total_reads'])
    raw_R1_bases = str(Info1['total_bases'])
    Q20_R1_bp = str(Info1['q20_bases'])
    Q20_R1_bp = str(Info1['q20_bases'])
    Q30_R1_bp = str(Info1['q30_bases'])
    Info2 = data['read2_before_filtering']
        #check if there are reads otherwise dividing by zero throw
    if Decimal(Info1['total_bases']) == 0:
        Q20_R1_percent, Q20_R2_percent, Q30_R1_percent, Q30_R2_percent = (str(0) for i in range(4))
    else:
        Q20_R1_percent = str(Decimal(Info1['q20_bases']) / Decimal(Info1['total_bases']))
        Q30_R1_percent = str(Decimal(Info1['q30_bases']) / Decimal(Info1['total_bases']))
        Q20_R2_percent = str(Decimal(Info2['q20_bases']) / Decimal(Info2['total_bases']))
        Q30_R2_percent = str(Decimal(Info2['q30_bases']) / Decimal(Info2['total_bases']))
    raw_R2_reads = str(Info2['total_reads'])
    raw_R2_bases = str(Info2['total_bases'])
    Q20_R2_bp = str(Info2['q20_bases'])
    Q30_R2_bp = str(Info2['q30_bases'])
    Line = name + '\t' + raw_R1_reads + '\t' + raw_R1_bases + '\t' + raw_R2_reads + '\t' + raw_R2_bases + '\t' + Q20_Total + '\t' + Q30_Total + '\t' + Q20_R1_bp + '\t' + Q20_R2_bp + '\t' + Q20_R1_percent + '\t' + Q20_R2_percent + '\t' + Q30_R1_bp + '\t' + Q30_R2_bp + '\t' + Q30_R1_percent + '\t' + Q30_R2_percent + '\t' + Total_Sequenced_bp + '\t' + Total_Sequenced_reads
    Out.write(Line)
    Out.close()

def FastP_QC_after(input_trimmed_json, input_singles_json, output_file, name):
    """Makes a QC output file from an input FastP json for orphaned reads"""
    Out = open(output_file, 'w')
    Out.write('Name\tR1[reads]\tR1[bp]\tR2[reads]\tR2[bp]\tUnpaired[reads]\tUnpaired[bps]\tQ20_Total_[bp]\tQ30_Total_[bp]\tQ20_R1_[bp]\tQ20_R2_[bp]\tQ20_unpaired[bp]\tQ20_R1_[%]\tQ20_R2_[%]\tQ20_unpaired[%]\tQ30_R1_[bp]\tQ30_R2_[bp]\tQ30_unpaired[bp]\tQ30_R1_[%]\tQ30_R2_[%]\tQ30_unpaired[%]\tTotal_Sequenced_[bp]\tPaired_Sequenced_[reads]\tTotal_Sequenced_[reads]\n')
    f = open(input_trimmed_json)
    data = json.load(f)
    f.close()
    paired_total_info = data['summary']['after_filtering']
    Q20_Total_trimmed = str(paired_total_info['q20_bases'])
    Q30_Total_trimmed = str(paired_total_info['q30_bases'])
    Trimmed_Sequenced_bp = str(paired_total_info['total_bases'])
    Trimmed_Sequenced_reads = str(paired_total_info['total_reads'])
    Info1 = data['read1_after_filtering']
    raw_R1_reads = str(Info1['total_reads'])
    raw_R1_bases = str(Info1['total_bases'])
    Q20_R1_bp = str(Info1['q20_bases'])
    Q30_R1_bp = str(Info1['q30_bases'])
    Info2 = data['read2_after_filtering']
    #check if there are reads otherwise dividing by zero throw
    if Decimal(Info1['total_bases']) == 0:
        Q20_R1_percent, Q20_R2_percent, Q30_R1_percent, Q30_R2_percent = (str(0) for i in range(4))
    else:
        Q20_R1_percent = str(Decimal(Info1['q20_bases']) / Decimal(Info1['total_bases']))
        Q30_R1_percent = str(Decimal(Info1['q30_bases']) / Decimal(Info1['total_bases']))
        Q20_R2_percent = str(Decimal(Info2['q20_bases']) / Decimal(Info2['total_bases']))
        Q30_R2_percent = str(Decimal(Info2['q30_bases']) / Decimal(Info2['total_bases']))
    raw_R2_reads = str(Info2['total_reads'])
    raw_R2_bases = str(Info2['total_bases'])
    Q20_R2_bp = str(Info2['q20_bases'])
    Q30_R2_bp = str(Info2['q30_bases'])
    f = open(input_singles_json)
    data = json.load(f)
    f.close()
    singles_total_info = data['summary']['after_filtering']
    Q20_Total_singles = str(singles_total_info['q20_bases'])
    Q30_Total_singles = str(singles_total_info['q30_bases'])
    Singles_Sequenced_bp = str(singles_total_info['total_bases'])
    Singles_Sequenced_reads = str(singles_total_info['total_reads'])
    unpaired_reads = str(singles_total_info['total_reads'])
    unpaired_bases = str(singles_total_info['total_bases'])
    ## Found these are already calculated in the FASTP output
    #Q20_unpaired_percent = str(Decimal(singles_total_info['q20_bases']) / Decimal(singles_total_info['total_bases']))
    #Q30_unpaired_percent = str(Decimal(singles_total_info['q30_bases']) / Decimal(singles_total_info['total_bases']))
    Q20_unpaired_percent=str(round(singles_total_info['q20_rate'],4))
    Q30_unpaired_percent=str(round(singles_total_info['q30_rate'],4))


    Total_Sequenced_bp = str(int(Trimmed_Sequenced_bp) + int(Singles_Sequenced_bp))
    Total_Sequenced_reads = str(int(Trimmed_Sequenced_reads) + int(Singles_Sequenced_reads))


    Q20_Total = str(int(Q20_Total_trimmed) + int(Q20_Total_singles))
    Q30_Total = str(int(Q30_Total_trimmed) + int(Q30_Total_singles))

    Line = name + '\t' + raw_R1_reads + '\t' + raw_R1_bases + '\t' + raw_R2_reads + '\t' + raw_R2_bases + '\t' + unpaired_reads + '\t' + unpaired_bases + '\t' + Q20_Total + '\t' + Q30_Total + '\t' + Q20_R1_bp + '\t' + Q20_R2_bp + '\t' + Q20_Total_singles + '\t' + Q20_R1_percent + '\t' + Q20_R2_percent + '\t' + Q20_unpaired_percent + '\t' + Q30_R1_bp + '\t' + Q30_R2_bp + '\t' + Q30_Total_singles + '\t' + Q30_R1_percent + '\t' + Q30_R2_percent + '\t' + Q30_unpaired_percent + '\t' + Total_Sequenced_bp + '\t' + Trimmed_Sequenced_reads + '\t' + Total_Sequenced_reads
    Out.write(Line)
    Out.close()

def FastP_QC_All(input_paired_json, input_singles_json, name):
    """Makes a pre- and post-filtering output of trimmed info"""
    trimmed_output= name + "_trimmed_read_counts.txt"
    FastP_QC_after(input_paired_json, input_singles_json, trimmed_output, name)

def main():
    args = parse_cmdline()
    FastP_QC_All(args.trimmed_json, args.single_json, args.name)

if __name__ == '__main__':
    main()