#!/usr/bin/env python3

import sys
import glob
import os
import os.path
from decimal import *
getcontext().prec = 4
import argparse

##Makes a summary Excel file when given a run folder from PhoeNiX
##Usage: >python Phoenix_Summary_Line_06-10-22.py -n Sequence_Name -t Trimmed_QC_Data_File -r Ratio_File -m MLST_File -q Quast_File -a AR_GAMMA_File -v Hypervirulence_GAMMA_File  -s synopsis_file -o Out_File
## Written by Rich Stanton (njr5@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary line')
    parser.add_argument('-n', '--name', required=True, help='sequence name')
    parser.add_argument('-t', '--trimmed', required=True, help='QC data file for trimmed reads')
    parser.add_argument('-r', '--ratio', required=True, help='assembly ratio file')
    parser.add_argument('-m', '--mlst', required=True, help='MLST file')
    parser.add_argument('-q', '--quast', required=True, help='QUAST file')
    parser.add_argument('-a', '--ar', required=True, help='AR GAMMA file')
    parser.add_argument('-v', '--vir', required=True, help='hypervirulence GAMMA file')
    parser.add_argument('-s', '--stats', dest="stats", required=True, help='Pipeline Stats file synopsis file')
    parser.add_argument('-o', '--out', required=True, help='output file name')
    return parser.parse_args()

def MLST_Info_Only(MLST_file):
    """Pulls MLST info from *_Scheme.mlst file"""
    f = open(MLST_file, 'r')
    String1 = f.readline()
    ST = String1.split()[2]
    Scheme = String1.split()[1]
    Out = ST + ' (' + Scheme + ')'
    return Out

def MLST_ST(MLST_file):
    """Pulls MLST info from *_Scheme.mlst file"""
    ST = 'Unknown'
    f = open(MLST_file, 'r')
    String1 = f.readline()
    ST = String1.split()[2]
    return ST

def MLST_Scheme(MLST_file):
    """Pulls MLST info from *_Scheme.mlst file"""
    Scheme = 'Unknown'
    f = open(MLST_file, 'r')
    String1 = f.readline()
    Scheme = String1.split()[1]
    return Scheme

def Contig_Count(input_quast):
    Contigs = '0'
    f = open(input_quast, 'r')
    String1 = f.readline()
    while String1 != '':
        if ('# contigs (>= 0 bp)' in String1):
            Contigs = String1.split()[-1]
            break
        String1 = f.readline()
    return Contigs

def Genome_Size(input_quast):
    Size = '0'
    f = open(input_quast, 'r')
    String1 = f.readline()
    while String1 != '':
        if ('Total length (>= 0 bp)' in String1):
            Size = String1.split()[-1]
            break
        String1 = f.readline()
    return Size

def N50_Length(input_quast):
    N50 = '0'
    f = open(input_quast, 'r')
    String1 = f.readline()
    while String1 != '':
        if ('N50' in String1):
            N50 = String1.split()[-1]
            break
        String1 = f.readline()
    return N50

def GC_Content(input_quast):
    NGC = '0'
    f = open(input_quast, 'r')
    String1 = f.readline()
    while String1 != '':
        if ('GC' in String1):
            GC = String1.split()[-1]
            break
        String1 = f.readline()
    return GC

def Assembly_Ratio(ratio_file):
    f = open(ratio_file, 'r')
    String1 = f.readline()
    Ratio = 'Unknown'
    SD = 'Unknown'
    while String1 != '':
        if ('Isolate_St.Devs' in String1):
            SD = String1.split()[1]
        elif ('Ratio' in String1):
            Ratio = String1.split()[1]
        String1 = f.readline()
    f.close()
    Out = Ratio + ' (' + SD + ')'
    return Out

def Assembly_Ratio_Length(ratio_file):
    f = open(ratio_file, 'r')
    String1 = f.readline()
    Length = 0
    while String1 != '':
        if ('Actual_length' in String1):
            Length = String1.split()[1]
        String1 = f.readline()
    f.close()
    Out = int(Length)
    return Out

def Assembly_Ratio_Species(ratio_file):
    f = open(ratio_file, 'r')
    String1 = f.readline()
    Species = 'Unknown'
    while String1 != '':
        if ('Tax:' in String1):
            Species = String1.split()[1:]
            Species = ' '.join(Species)
        String1 = f.readline()
    f.close()
    return Species

def Trimmed_BP(trimmed_counts_file):
    f = open(trimmed_counts_file, 'r')
    String1 = f.readline()
    String1 = f.readline()
    BP = String1.split()[-2]
    BP = int(BP)
    return BP

def Trim_Coverage(trimmed_counts_file, ratio_file):
    Length = Assembly_Ratio_Length(ratio_file)
    Trimmed = Trimmed_BP(trimmed_counts_file)
    Coverage = str(Decimal(Trimmed) / Decimal(Length))
    return Coverage

def Bla_Genes(input_gamma):
    f = open(input_gamma, 'r')
    Bla = []
    String1 = f.readline()
    for line in f:
        List1 = line.split()
        Cat = List1[0].split('__')[4]
        Gene = List1[0].split('__')[2]
        if ('LACTAM' in Cat.upper()):
            Bla.append(Gene)
    f.close()
    Bla.sort()
    return Bla

def Non_Bla_Genes(input_gamma):
    f = open(input_gamma, 'r')
    Non_Bla = []
    String1 = f.readline()
    for line in f:
        List1 = line.split()
        Cat = List1[0].split('__')[4]
        Gene = List1[0].split('__')[2]
        if ('LACTAM' in Cat.upper()) == False:
            Non_Bla.append(Gene)
    f.close()
    Non_Bla.sort()
    return Non_Bla

def HV_Genes(input_gamma):
    f = open(input_gamma, 'r')
    HV = []
    String1 = f.readline()
    for line in f:
        Gene = line.split()[0]
        HV.append(Gene)
    f.close()
    HV.sort()
    return HV

def QC_Pass(stats):
    f = open(stats, 'r')
    status = []
    reason = []
    String1 = f.readline()
    for line in f:
        print(line)
        if line.startswith("Auto Pass/FAIL"):
            line_status = line.split(":")[1]
            line_reason = line.split(":")[2]
            status.append(line_status.strip())
            reason.append(line_reason.strip())
    f.close()
    status= str(status[0])
    reason = str(reason[0])
    return status, reason

def Isolate_Line(ID, trimmed_counts, ratio_file, MLST_file, quast_file, gamma_ar, gamma_hv, stats):
    try:
        Coverage = Trim_Coverage(trimmed_counts, ratio_file)
    except:
        Coverage = 'Unknown'
    try:
        Genome_Length = Genome_Size(quast_file)
    except:
        Genome_Length = 'Unknown'
    try:
        Ratio = Assembly_Ratio(ratio_file)
    except:
        Ratio = 'Unknown'
    try:
        Contigs = Contig_Count(quast_file)
    except:
        Contigs = 'Unknown'
    try:
        GC = GC_Content(quast_file)
    except:
        GC = 'Unknown'
    try:
        Species = Assembly_Ratio_Species(ratio_file)
    except:
        Species = 'Unknown'
    try:
        ST = MLST_ST(MLST_file)
    except:
        ST = 'Unknown'
    try:
        Scheme = MLST_Scheme(MLST_file)
    except:
        Scheme = 'Unknown'
    try:
        Bla = Bla_Genes(gamma_ar)
        Bla = ','.join(Bla)
    except:
        Bla = 'Unknown'
    try:
        Non_Bla = Non_Bla_Genes(gamma_ar)
        Non_Bla = ','.join(Non_Bla)
    except:
        Non_Bla = 'Unknown'
    try:
        HV = HV_Genes(gamma_hv)
        HV = ','.join(HV)
    except:
        HV = 'Unknown'
    try:
        QC_Outcome, Reason = QC_Pass(stats)
    except:
        QC_Outcome = 'Unknown'
        Reason = ""
    Line = ID + '\t' + Coverage + '\t' + Genome_Length + '\t' + Ratio + '\t' + Contigs + '\t' + Species + '\t' + Scheme + '\t' + ST + '\t' + GC + '\t' + Bla + '\t' + Non_Bla + '\t' + HV + '\t' + QC_Outcome + '\t' + Reason
    return Line

def Isolate_Line_File(ID, trimmed_counts, ratio_file, MLST_file, quast_file, gamma_ar, gamma_hv, out_file, stats):
    Line = Isolate_Line(ID, trimmed_counts, ratio_file, MLST_file, quast_file, gamma_ar, gamma_hv, stats)
    Out = open(out_file, 'w')
    Out.write(Line)
    Out.close()

def main():
    args = parseArgs()
    Isolate_Line_File(args.name, args.trimmed, args.ratio, args.mlst, args.quast, args.ar, args.vir, args.out, args.stats)

if __name__ == '__main__':
    main()
