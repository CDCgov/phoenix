#!/usr/bin/env python3

import argparse
from argparse import ArgumentParser
from decimal import *

##Makes pre- and post-filtering QC outputs
##Usage: >python3 FastP_QC.py paired_fastp.json single_fastp_json Isolate_Name
##Written by Rich Stanton (njr5@cdc.gov) and Nick Vlachos (nvx4@cdc.gov)

def parse_cmdline():
	"""Parse command-line arguments for script."""
	parser = ArgumentParser(prog="FastP_QC.py", description="""Script collects stats from fastp jsons.""")
	parser.add_argument("-r1", "--r1_stats", dest="r1_stats", action="store", required=True, help="Text file with r1 stats, from q30.py script.")
	parser.add_argument("-r2", "--r2_stats", dest="r2_stats", action="store", required=True, help="Text file with r2 stats, from q30.py script.")
	parser.add_argument("-n", "--name", dest="name", action="store", required=True, help="Sample name")
	args = parser.parse_args()
	return args

def write_raw_stats(raw_R1_reads, raw_R1_bases, Q20_R1_bp, Q20_R1_percent, Q30_R1_bp, Q30_R1_percent, raw_R2_reads, raw_R2_bases, Q20_R2_bp, Q20_R2_percent, Q30_R2_bp, Q30_R2_percent, output_file, name):
    """Makes a QC output file from parsed output of q30.py files."""
    with open(output_file, 'w') as f:
        f.write('Name\tR1[reads]\tR1[bp]\tR2[reads]\tR2[bp]\tQ20_Total_[bp]\tQ30_Total_[bp]\tQ20_R1_[bp]\tQ20_R2_[bp]\tQ20_R1_[%]\tQ20_R2_[%]\tQ30_R1_[bp]\tQ30_R2_[bp]\tQ30_R1_[%]\tQ30_R2_[%]\tTotal_Sequenced_[bp]\tTotal_Sequenced_[reads]\n')
        Q20_Total = str(Q20_R1_bp + Q20_R2_bp)
        Q30_Total = str(Q30_R1_bp + Q30_R2_bp)
        Total_Sequenced_bp = str(raw_R1_bases + raw_R2_bases)
        Total_Sequenced_reads = str(raw_R1_reads + raw_R2_reads)
        Line = name + '\t' + str(raw_R1_reads) + '\t' + str(raw_R1_bases) + '\t' + str(raw_R2_reads) + '\t' + str(raw_R2_bases) + '\t' + Q20_Total + '\t' + Q30_Total + '\t' + str(Q20_R1_bp) + '\t' + str(Q20_R2_bp) + '\t' + str(Q20_R1_percent) + '\t' + str(Q20_R2_percent) + '\t' + str(Q30_R1_bp) + '\t' + str(Q30_R2_bp) + '\t' +str( Q30_R1_percent) + '\t' + str(Q30_R2_percent) + '\t' + Total_Sequenced_bp + '\t' + Total_Sequenced_reads
        f.write(Line)

def get_raw_stats(stats):
    with open(stats) as f:
        lines = f.readlines()
        raw_reads = int(lines[0].strip('\n').replace("total reads: ", ""))
        raw_bases = int(lines[1].strip('\n').replace("total bases: ", ""))
        Q20_bp = int(lines[2].strip('\n').replace("q20 bases: ", ""))
        Q20_percent = str(round(float(lines[4].strip('\n').replace("q20 percents: ", ""))/100, 4)) # need 0.XXXX notation for pipeline_stats_writer.sh
        Q30_bp = int(lines[3].strip('\n').replace("q30 bases: ", ""))
        Q30_percent = str(round(float(lines[5].strip('\n').replace("q30 percents: ", ""))/100, 4)) # need 0.XXXX notation for pipeline_stats_writer.sh
    return raw_reads, raw_bases, Q20_bp, Q20_percent, Q30_bp, Q30_percent

def all_raw_stats(r1_stats, r2_stats, name):
    """Makes a pre- and post-filtering output of trimmed info"""
    raw_output= name + "_raw_read_counts.txt"
    raw_R1_reads, raw_R1_bases, Q20_R1_bp, Q20_R1_percent, Q30_R1_bp, Q30_R1_percent = get_raw_stats(r1_stats)
    raw_R2_reads, raw_R2_bases, Q20_R2_bp, Q20_R2_percent, Q30_R2_bp, Q30_R2_percent = get_raw_stats(r2_stats)
    write_raw_stats(raw_R1_reads, raw_R1_bases, Q20_R1_bp, Q20_R1_percent, Q30_R1_bp, Q30_R1_percent, raw_R2_reads, raw_R2_bases, Q20_R2_bp, Q20_R2_percent, Q30_R2_bp, Q30_R2_percent, raw_output, name)

def main():
    args = parse_cmdline()
    all_raw_stats(args.r1_stats, args.r2_stats, args.name)

if __name__ == '__main__':
    main()