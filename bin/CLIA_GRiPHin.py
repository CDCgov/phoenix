#!/usr/bin/env python3

import sys
sys.dont_write_bytecode = True
import glob
import os
from decimal import *
import pandas as pd
import numpy as np
import argparse
import re
from re import search
import operator
import xlsxwriter as ws
from xlsxwriter.utility import xl_rowcol_to_cell
import csv
import string
from Bio import SeqIO
from itertools import chain
from species_specific_griphin import  transform_value, create_shiga_df, double_check_taxa_id, fill_taxa_id


##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python GRiPHin.py -s ./samplesheet.csv -a ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output --phoenix --scaffolds
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    # <3.0 uses gamma as AR output
    return "3.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-s', '--samplesheet', default=None, required=False, dest='samplesheet', help='PHoeNIx style samplesheet of sample,directory in csv format. Directory is expected to have PHoeNIx stype output.')
    parser.add_argument('-b', '--bldb', default=None, required=False, dest='bldb', help='If a directory is given rather than samplesheet GRiPHin will create one for all samples in the directory.')
    parser.add_argument('-d', '--directory', default=None, required=False, dest='directory', help='If a directory is given rather than samplesheet GRiPHin will create one for all samples in the directory.')
    parser.add_argument('-c', '--control_list', required=False, dest='control_list', help='CSV file with a list of sample_name,new_name. This option will output the new_name rather than the sample name to "blind" reports.')
    parser.add_argument('-a', '--ar_db', dest="ar_db", required=True, help='Pass the name of the amrfinder database used. Only used for documentation purposes')
    parser.add_argument('-o', '--output', default="", required=False, dest='output', help='Name of output file default is GRiPHin_Summary.xlsx.')
    parser.add_argument('--phx_version', default="Unknown", required=False, dest='phx_version', help='The version of phx used to produce GRiPHin_Summary row for the sample.')
    parser.add_argument('--coverage', default=30, required=False, dest='set_coverage', help='The coverage cut off default is 30x.')
    parser.add_argument('--shigapass', dest="shigapass", default=False, action='store_true', required=False, help='Use for when there are E. coli or Shigella isolates in samplesheet.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CYELLOW = '\033[93m'
CEND = '\033[0m'

def Get_Parent_Folder(directory):
    '''getting project and parent_folder info from the paths'''
    #Project - parent folder (first folder that is in the outdir)
    #relative/submission - rest of the path
    #first make sure we have an absolute path
    directory = os.path.abspath(directory)
    #handing if trailing backslash isn't in there.
    if directory[-1] != "/": 
        directory = directory + "/"
    # get project from directory path
    project = os.path.split(os.path.split(os.path.split(directory)[0])[0])[1]
    # get everything after CEMB
    parent_folder = os.path.split(os.path.split(os.path.split(os.path.split(directory)[0])[0])[0])[0]
    #parent_folder = os.path.split(cemb_path)[1].lstrip("/") # remove backslash on left side to make it clean  #this is only the last name of the folder not full path
    return project, parent_folder

def get_Q30(trim_stats, raw_stats):
    # trimmed data
    data_df = pd.read_csv(trim_stats, sep='\t', header=0)
    Trim_Q30_R1_percent = round(data_df["Q30_R1_[%]"]*100, 2)[0] #make percent and round to two decimals
    Trim_Q30_R2_percent = round(data_df["Q30_R2_[%]"]*100, 2)[0]
    Total_Trimmed_bp = data_df["Total_Sequenced_[bp]"][0]
    Total_Trimmed_reads = data_df["Total_Sequenced_[reads]"][0]
    Paired_Trimmed_reads = data_df["Paired_Sequenced_[reads]"][0]
    #do the same with raw
    data_raw_df = pd.read_csv(raw_stats, sep='\t', header=0)
    Q30_R1_percent = round(data_raw_df["Q30_R1_[%]"]*100, 2)[0] #make percent and round to two decimals
    Q30_R2_percent = round(data_raw_df["Q30_R2_[%]"]*100, 2)[0]
    Total_Raw_bp = data_raw_df["Total_Sequenced_[bp]"][0]
    Total_Raw_reads = data_raw_df["Total_Sequenced_[reads]"][0]
    return Q30_R1_percent, Q30_R2_percent, Total_Raw_bp, Total_Raw_reads, Total_Trimmed_bp, Paired_Trimmed_reads, Total_Trimmed_reads, Trim_Q30_R1_percent, Trim_Q30_R2_percent

def get_kraken_info(kraken_trim, kraken_wtasmbld, sample_name):
    try:
        with open(kraken_trim,"r") as f:
            for line in f:
                if line.startswith("U:"):
                    Trim_unclassified_percent = line.split(' ')[1].strip()
                elif line.startswith("G:"):
                    Trim_Genus_percent = line.split(' ')[1].strip()
                    Trim_Genus = line.split(' ')[2].strip()
                elif line.startswith("s:"):
                    Trim_Species_percent = line.split(' ')[1].strip()
                    Trim_Species = line.split(' ')[2].strip()
        Trim_kraken = Trim_Genus + " (" + Trim_Genus_percent + ") " + Trim_Species + " (" + Trim_Species_percent + ")"
    except FileNotFoundError:
        print("Warning: " + sample_name + ".trimd_summary.txt not found")
        Trim_kraken = 'Unknown'
        Trim_unclassified_percent = "Unknown"
        Trim_Genus_percent = "Unknown"
    try:
        with open(kraken_wtasmbld,"r") as f:
            for line in f:
                if line.startswith("U:"):
                    Asmbld_unclassified_percent = line.split(' ')[1].strip()
                elif line.startswith("G:"):
                    Asmbld_Genus_percent = line.split(' ')[1].strip()
                    Asmbld_Genus = line.split(' ')[2].strip()
                elif line.startswith("s:"):
                    Asmbld_Species_percent = line.split(' ')[1].strip()
                    Asmbld_Species = line.split(' ')[2].strip()
        Asmbld_kraken = Asmbld_Genus + " (" + Asmbld_Genus_percent + ") " + Asmbld_Species + " (" + Asmbld_Species_percent + ")"
    except FileNotFoundError:
        print("Warning: " + sample_name + ".wtasmbld_summary.txt not found")
        Asmbld_kraken = 'Unknown'
        Asmbld_unclassified_percent = "Unknown"
        Asmbld_Genus_percent = 0
    return Trim_kraken, Trim_Genus_percent, Asmbld_kraken, Asmbld_Genus_percent, Trim_unclassified_percent, Asmbld_unclassified_percent

def Calculate_Trim_Coverage(Total_Trimmed_bp, quast_report):
    """Taking the total sequenced bp from fastp and the assembly length from quast to calculate coverage."""
    with open(quast_report, 'r') as f:
        for line in f:
            if ('Total length' in line):
                Assembly_Length = int(line.split('\t')[1])
                break
    Coverage = float(round(Total_Trimmed_bp / Assembly_Length, 2))
    return Coverage, Assembly_Length

def Get_Assembly_Length(quast_report): #for -entry SCAFFOLDS or CDC_SCAFFOLDS
    """Taking the the assembly length from quast to get assembly length."""
    with open(quast_report, 'r') as f:
        for line in f:
            if ('Total length' in line):
                Assembly_Length = int(line.split('\t')[1])
                break
    return Assembly_Length

def get_scaffold_count(quast_report):
    scaffolds = '0'
    with open(quast_report, 'r') as f:
        for line in f:
            if ('# contigs (>= 0 bp)' in line):
                scaffolds = int(line.split()[-1])
                break
    return scaffolds

def Get_BUSCO_Gene_Count(busco_short_summary):
    """Parse BUSCO file to get the lineage and % of matching BUSCOs"""
    with open(busco_short_summary, 'r') as f:
        for line in f:
            if "The lineage dataset is:" in line:
                lineage = re.search(r': (.*) \(', line).group()
                lineage = re.sub(r":", "", lineage)
                lineage = re.sub(r"\(", "", lineage).strip()
            if "Complete BUSCOs" in line:
                found_buscos = int(re.search(r'\d+', line).group())
            if "Total BUSCO groups searched" in line:
                total_buscos = int(re.search(r'\d+', line).group())
    percent_busco = round(((found_buscos/total_buscos)*100),2)
    full_busco_line = str(percent_busco) + " (" + str(found_buscos) + "/" + str(total_buscos) +")"
    #busco_line = lineage + " (" + percent_busco + "%)" # old busco line
    busco_metrics = [lineage, percent_busco]
    return busco_metrics, full_busco_line

def get_gc_metrics(gc_file):
    '''Collects the gc and gc stdev.'''
    with open(gc_file, 'r') as f:
        for line in f:
            if "Species_GC_StDev:" in line:
                if "Not calculated on species with n<10 references" in line or "No Match Found" in line:
                    gc_stdev = "NA"
                    out_of_range_stdev = gc_stdev
                else:
                    gc_stdev = float((line.split("Species_GC_StDev: ",1)[1]).strip())
                    out_of_range_stdev = gc_stdev*2.58
            elif "Sample_GC_Percent:" in line:
                if "No Match Found" in line or "NA" in line:
                    sample_gc="NA"
                else:
                    sample_gc = float((line.split("Sample_GC_Percent: ",1)[1]).strip())
            elif "Species_GC_Mean:" in line:
                if "No Match Found" in line or "-" in line:
                    species_gc_mean="NA"
                else:
                    species_gc_mean = float((line.split("Species_GC_Mean: ",1)[1]).strip())
            else:
                pass
    return gc_stdev, sample_gc, out_of_range_stdev, species_gc_mean

def get_assembly_ratio(asmbld_ratio, tax_file):
    '''Collects the assembly ratio, stdev and method used to determine taxa.'''
    with open(asmbld_ratio, 'r') as f:
        for line in f:
            if "Tax: " in line:
                taxa = (line.split("Tax: ",1)[1]).strip()
            if "Isolate_St.Devs:" in line:
                stdev = (line.split("Isolate_St.Devs: ",1)[1]).strip()
                if stdev == 'NA' or stdev == 'No Match Found': #handling making stdev a float only if its a number
                    pass
                else:
                    stdev =  float(stdev)
            if "Ratio: " in line:
                ratio =  float((line.split("Ratio: ",1)[1]).strip())
    #assembly_ratio_line = ratio + " (stdev " + stdev + ") " + taxa #old way of reporting in one line
    with open(tax_file, 'r') as f:
        tax_method = f.readline().split("\t")[0]
    assembly_ratio_metrics = [ratio, stdev, tax_method]
    return assembly_ratio_metrics

def compile_alerts(coverage, assembly_stdev, gc_stdev):
    """
    No orphaned reads found after trimming
    <10 reference genomes for species identified so no stdev for assembly ratio or %GC content calculated
    >150x coverage or <40x coverage
    """
    alerts = []
    if coverage != "Unknown": # if its unknown it will fail already so skip
        if int(coverage) > 30 and int(coverage) < 40:
            alerts.append("Coverage between 30-40x("+ str(coverage) + "x)")
        elif int(coverage) > 100.00:
            alerts.append("Coverage >100x(" + str(coverage) + "x)")
    if str(assembly_stdev) == "NA":
        if str(gc_stdev) == "NA":
            alerts.append("Assembly ratio and GC% STDev are N/A <10 genomes as reference")
        else:
            alerts.append("Open Github issue assembly ratio STDev is N/A but not GC%. This shouldn't happen.")
    elif str(gc_stdev) == "NA":
        alerts.append("Open Github issue STDev is N/A for GC% and not assembly ratio. This shouldn't happen.")
    else:
        pass
    alerts = ', '.join(alerts)
    return alerts

def compile_warnings(Total_Trimmed_reads, Total_Raw_reads, Q30_R1_per, Q30_R2_per, Trim_Q30_R1_per, Trim_Q30_R2_per, scaffolds, gc_metrics, \
                     assembly_ratio_metrics, Trim_unclassified_percent, Wt_asmbld_unclassified_percent, kraken_trim_genus, kraken_wtasmbld_genus, Trim_Genus_percent, Asmbld_Genus_percent,\
                     fastani_warning, busco_id, FastANI_ID, FastANI_coverage, QC_reason):
    """
    <1,000,000 total reads for each raw and trimmed reads - Total_Sequenced_reads
    % raw and trimmed reads with Q30 average for R1 (<90%) and R2 (<70%) - Q30_R1_percent, Q30_R2_percent
    >200 scaffolds - scaffolds
    Checking that %GC content isn't >2.58 stdev away from the mean %GC content for the species determined - assembly_ratio_metrics
    Contamination check: >30% unclassified reads and confirm there is only 1 genera with >25% of assigned reads - Trim_kraken, wt_Asmbld_kraken
    """
    # check warnings
    warnings = []
    if Total_Trimmed_reads == "Unknown" or int(Total_Trimmed_reads) < int(1000000):
        warnings.append("<1,000,000 trimmed reads")
    if Total_Raw_reads == "Unknown" or int(Total_Raw_reads) < int(1000000):
        warnings.append("<1,000,000 raw reads")
    if Q30_R1_per == "Unknown" or float(Q30_R1_per) < float(90.00):
        warnings.append("Average Q30 of raw R1 reads <{:.2f}%".format(float(90.00)))
    if Q30_R2_per == "Unknown" or float(Q30_R2_per) < float(70.00):
        warnings.append("Average Q30 of raw R2 reads <{:.2f}%".format(int(70.00)))
    if Trim_Q30_R1_per == "Unknown" or float(Trim_Q30_R1_per) < float(90.00):
        try:
            warnings.append("Average Q30 of trimmed R1 reads <{:.2f}% ({:.2f}%)".format(float(90.00),float(Trim_Q30_R1_per)))
        except ValueError:
            warnings.append("Average Q30 of trimmed R1 reads <{:.2f}% ({})".format(float(90.00),Trim_Q30_R1_per))
    if Trim_Q30_R2_per == "Unknown" or float(Trim_Q30_R2_per) < float(70.00):
        try:
            warnings.append("Average Q30 of trimmed R2 reads <{:.2f}% ({:.2f}%)".format(float(70.00),float(Trim_Q30_R2_per)))
        except ValueError:
            warnings.append("Average Q30 of trimmed R2 reads <{:.2f}% ({})".format(float(70.00),Trim_Q30_R2_per))
    if Trim_unclassified_percent == "Unknown" or float(Trim_unclassified_percent) > float(30.00):
        warnings.append(">{:.2f}% unclassifed trimmed reads".format(int(30)))
    if len(kraken_trim_genus) >=2:
        warnings.append(">=2 genera had >{:.2f}% of reads assigned to them".format(int(25)))
    if Trim_Genus_percent == "Unknown" or float(Trim_Genus_percent) <float(70.00):
        try:
            warnings.append("<70% of reads assigned to top genera hit ({:.2f}%)".format(float(Trim_Genus_percent)))
        except ValueError:
            warnings.append("<70% of reads assigned to top genera hit ({})".format(Trim_Genus_percent))
    if gc_metrics[0] != "NA" and gc_metrics[0] != "Unknown":
        # sample_gc > (species_gc_mean + out_of_range_stdev)
        if float(gc_metrics[1]) > (float(gc_metrics[3])+float(gc_metrics[2])): #check that gc% is < 2.58 stdev away from mean gc of species
            warnings.append("GC% >2.58 stdev away from mean GC of {:.2f}%.".format(float(gc_metrics[3])))
    if scaffolds != "Unknown" and Wt_asmbld_unclassified_percent != "Unknown" and Asmbld_Genus_percent != "Unknown":
        if int(scaffolds) > int(200) and int(scaffolds) < int(500): # between 200-500 
            warnings.append("High scaffold count 200-500 ({}).".format(int(scaffolds)))
        if float(Wt_asmbld_unclassified_percent) > float(30.00):
            warnings.append(">{:.2f}% unclassifed weighted scaffolds.".format(int(30)))
        if float(Asmbld_Genus_percent) <float(70.00):
            warnings.append("<70% of weighted scaffolds assigned to top genera hit ({:.2f}%)".format(float(Asmbld_Genus_percent)))
    elif scaffolds == "Unknown" and Wt_asmbld_unclassified_percent == "Unknown" and Asmbld_Genus_percent == "Unknown":
        if "No assembly due to:" in QC_reason: # if there is already a QC reason for no assembly then don't add this warning
            pass
        else:
            warnings.append("No assembly file found possible SPAdes failure.")
    if len(kraken_wtasmbld_genus) >=2:
        warnings.append(">=2 genera had >{:.2f}% of wt scaffolds assigned to them.".format(int(25))) 
    if FastANI_ID != "Unknown":
        if float(FastANI_ID) < float(95.00):
            warnings.append("FastANI match is <95%.")
        if float(FastANI_coverage) < float(90.00):
            warnings.append("FastANI coverage is <90%.")
    if busco_id != "Unknown":
        if float(busco_id) < float(97.00):
            warnings.append("BUSCO match is <97%.")
    #add in fastani warning
    if fastani_warning != None:
        warnings.append(fastani_warning)
    # For spades failures, lack of reads after trimming or corruption we will simplify the warnings by supressing other warnings
    if "No assembly file found possible SPAdes failure." in warnings:
        warnings = "No assembly file found possible SPAdes failure."
    elif "is corrupt and is unable to be unzipped" in warnings:
        warnings = [item for item in warnings if "corrupt" in item]
    elif "The # of reads in raw R1/R2 files are NOT equal." in warnings:
        warnings = [item for item in warnings if "NOT equal" in item]
    # Reduce warnings if certain QC reasons are present
    if "The # of reads in raw R1/R2 files are NOT equal." in QC_reason:
        warnings = [item for item in warnings if "trimmed" not in item and "reads assigned" not in item]
        warnings.insert(0, "Skipped trimmed steps: unequal R1/R2 read counts.")
    elif "No reads remain after trimming" in QC_reason :
        warnings = [item for item in warnings if "trimmed" not in item and "reads assigned" not in item]
        warnings.insert(0, "Skipped trimmed steps: No reads remain after trimming.")
    elif "corrupt" in QC_reason :
        warnings = "Corrupted input FASTQ file(s): downstream steps skipped."
    if isinstance(warnings, list) and len(warnings) > 1:
        warnings = ', '.join(warnings).strip(", ")
    elif len(warnings) == 1:
        warnings = warnings[0]
    elif warnings == [""] or warnings == []:
        warnings = ""
    return warnings

def parse_kraken_report(kraken_trim_report, kraken_wtasmbld_report, sample_name):
    """Checking that only 1 genera with >25% of assigned reads."""
    kraken_trim_genus = []
    kraken_wtasmbld_genus = []
    try:
        #file is a VERY WERID so need some extra arguments
        with open(kraken_trim_report, mode='r', encoding='utf8', newline='\r') as f:
            for line in f: #weird file so its really just one big long line
                clean_line = line.replace('  ', '').strip('\n') #cleaning up
                split_line = clean_line.split('\n') #split so we can go line by line
                for thing in split_line:
                    if "\tG\t" in thing:
                        genus_percent = float(thing.split('\t')[0].replace(' ',''))
                        if genus_percent >= 25.00:
                            kraken_trim_genus.append(thing.replace(' ',''))
    except FileNotFoundError:
        print("Warning: " + sample_name + ".kraken2_trimd.summary.txt not found")
        kraken_trim_genus = 'Unknown'
    try:
        #file is a VERY WERID so need some extra arguments
        with open(kraken_wtasmbld_report, mode='r', encoding='utf8', newline='\r') as f:
            for line in f: #weird file so its really just one big long line
                clean_line = line.replace('  ', '').strip('\n') #cleaning up
                split_line = clean_line.split('\n') #split so we can go line by line
                if '\tunclassified'in split_line[0]:
                    unclass_amount = float(split_line[0].split('\t')[1].replace(' ',''))
                else:
                    unclass_amount = 0
                for thing in split_line:
                    if '\troot' in thing:
                        root_amount = float(thing.split('\t')[1].replace(' ',''))
                        total_amount = root_amount + unclass_amount
                    elif "\tG\t" in thing:
                        genus_amount = float(thing.split('\t')[1].replace(' ',''))
                        genus_percent = (genus_amount/total_amount)*100
                        if genus_percent >= 25.00:
                            kraken_wtasmbld_genus.append(thing.replace(' ',''))
    except FileNotFoundError:
        print("Warning: " + sample_name + ".kraken2_wtasmbld.summary.txt not found")
        kraken_wtasmbld_report = 'Unknown'
    return kraken_trim_genus, kraken_wtasmbld_genus

def Checking_auto_pass_fail(fairy_files, spades_fairy_file, coverage, length, assembly_stdev, asmbld_ratio, set_coverage, scaffolds, sample_name):
    """
    Checking auto pass fail conditions
    SPAdes failure would say: "run_failure,no_scaffolds,no_contigs"
    SPAdes  Success options: 
    1. "run_completed,scaffolds_created,contigs_created"
    2. "run_completed,no_scaffolds,contigs_created"
    """
    #assembly_stdev = assembly_ratio_line.split("(")[1].split(")")[0].split(" ")[1] # parse to get standard dev, old method
    QC_reason = []
    QC_result = [] # set as blank to begin with, need this to assign variable in QC_result == "FAIL": line
    QC_result.append("PASS") #set default as PASS
    #check output of fairy and determine if reads were corrupt or had unequal number of reads
    for fairy_file in fairy_files:
        with open(fairy_file, 'r') as f:
            for line in f:
                if ('FAILED CORRUPTION CHECK!' in line):
                    fastq_file_failure = str(line.split(' ')[10])
                    QC_result.append("FAIL")
                    QC_reason.append(str(fastq_file_failure) +" is corrupt and is unable to be unzipped")
                if ('FAILED: The number of reads in R1/R2 are NOT the same!' in line):
                    QC_result.append("FAIL")
                    QC_reason.append("The # of reads in raw R1/R2 files are NOT equal")
                if ('FAILED: There are 0 reads in' in line):
                    QC_result.append("FAIL")
                    QC_reason.append("No reads remain after trimming")
                if ('FAILED: No scaffolds in ' in line):
                    QC_result.append("FAIL")
                    QC_reason.append("No scaffolds were >500bp")
    # Check if spades_fairy_file exists and is not empty before trying to open it
    if spades_fairy_file and spades_fairy_file != "":
        try:
            with open(spades_fairy_file, 'r') as sp:
                for line in sp:
                        if ('no_scaffolds,no_contigs' in line):
                            QC_result.append("FAIL")
                            QC_reason.append("No scaffolds created by SPAdes")
                        elif ('no_scaffolds,contigs_created' in line):
                            QC_result.append("FAIL")
                            QC_reason.append("No scaffolds created by SPAdes only contigs")
        except FileNotFoundError:
            print(f"Warning: {spades_fairy_file} not found for sample {sample_name}")
    if coverage == "Unknown" or int(coverage) < int(set_coverage):
        QC_result.append("FAIL")
        if coverage == "Unknown": # if else really only needed so you don't end up with "unknownx"
            QC_reason.append("Coverage <"+ str(set_coverage) +"x (" + str(coverage) + ")")
        else:
            QC_reason.append("Coverage <"+ str(set_coverage) +"x (" + str(coverage) + "x)")
    if length == "Unknown" or int(length) <= 1000000:
        QC_result.append("FAIL")
        QC_reason.append("assembly <1,000,000bps (" + str(length) + ")")
    if str(assembly_stdev) != "NA": # have to have a second layer cuz you can't make NA a float, N/A means less than 10 genomes so no stdev calculated
        if str(asmbld_ratio) == "Unknown": # if there is no ratio file then fail the sample
            QC_result.append("FAIL")
            QC_reason.append("Assembly file not found")
        elif float(assembly_stdev) > 2.58:
            QC_result.append("FAIL")
            QC_reason.append("Assembly stdev >2.58 (" + str(assembly_stdev) + ")")
    if str(scaffolds) == "Unknown" or int(scaffolds) > int(500):
        QC_result.append("FAIL")
        QC_reason.append("High scaffold count >500 ({}).".format(str(scaffolds)))
    QC_reason = set(QC_reason)
    QC_reason = ', '.join(QC_reason)
    # Simplify error for when Assembly file not found
    check_QC_reason = QC_reason
    if "Assembly file not found" in check_QC_reason:
        QC_reason = "Assembly file not found"
        if "is corrupt and is unable to be unzipped" in check_QC_reason:
            new_QC_reason = [item for item in check_QC_reason.split(",") if "corrupt" in item ]
            if len(new_QC_reason) == 1:
                QC_reason = "No assembly due to: " + new_QC_reason[0] + "."
            elif len(new_QC_reason) == 2:
                QC_reason = "No assembly due to: " + sample_name + "_R1 and " + sample_name + "_R2 files being corrupted."
            else:
                print(f"This shouldn't happen, please open a github issue.")
        elif "The # of reads in raw R1/R2 files are NOT equal" in check_QC_reason:
            new_QC_reason = [item for item in check_QC_reason.split(",") if "NOT equal" in item]
            QC_reason = "No assembly due to: " + new_QC_reason[0] + "."
        elif "No scaffolds created by SPAdes" in check_QC_reason:
            new_QC_reason = [item for item in check_QC_reason.split(",") if "No scaffolds created by SPAdes" in item ]
            QC_reason = "No assembly due to: " + new_QC_reason[0] + "."
        elif "No reads" in check_QC_reason or "No scaffolds" in check_QC_reason and "No scaffolds created by SPAdes" not in check_QC_reason:
            new_QC_reason = [item for item in check_QC_reason.split(",") if "No reads" in item or "No scaffolds" in item]
            QC_reason = "No assembly due to: " + new_QC_reason[0] + "."
    #checking if it was a pass
    if any("FAIL" in sub for sub in QC_result):
        QC_result = "FAIL"
    else:
        QC_result = "PASS"
    return QC_result, QC_reason

def duplicate_column_clean(df):
    if len([x for x in list(df.columns) if list(df.columns).count(x) > 1]) > 0:
        #get column names that are duplicates
        dups = set([x for x in list(df.columns) if list(df.columns).count(x) > 1])
        # get dataframe for duplicate columns
        for dup in dups:
            new_col = df[dup].agg(';'.join, axis=1).astype(str).values[0]
            #drop old frame(s)
            df = df.drop(dup, axis=1)
            #add in new frame
            df[dup] = new_col
    return df

def parse_amrfinder_ar(amrfinder_file, sample_name, final_df, ar_db):
    """Parsing the AMRFinder file run on the antibiotic resistance database."""
    db = ar_db.replace(".tar.gz", "")
    amrfinder_df = pd.read_csv(amrfinder_file, sep='\t', header=0)
    # Drop rows where 'Scope' contains the string 'plus'
    amrfinder_df = amrfinder_df[~amrfinder_df['Scope'].str.contains('plus')]
    # Drop rows where 'Scope' contains the string 'EFFLUX'
    amrfinder_df = amrfinder_df[~amrfinder_df['Class'].str.contains('EFFLUX')]
    # Drop rows where 'Element_type' contains the string 'STRESS'
    amrfinder_df = amrfinder_df[~amrfinder_df['Element_type'].str.contains('STRESS')]
    amrfinder_df = amrfinder_df.sort_values(by='Contig_id') #sort by contig id
    #percent_BP_IDs = np.floor(amrfinder_df["%_Identity_to_reference_sequence"]).tolist() # round % to whole number
    percent_codon_IDs = np.floor(amrfinder_df["%_Identity_to_reference_sequence"]).tolist() # round % to whole number
    percent_lengths = np.floor(amrfinder_df["%_Coverage_of_reference_sequence"]).tolist() # round % to whole number
    #fix capitalization of resistance
    amrfinder_df["Class"] = amrfinder_df["Class"].apply(lambda x: x.capitalize())
    amrfinder_df["Subclass"] = amrfinder_df["Subclass"].apply(lambda x: x.capitalize())
    conferred_resistances = amrfinder_df["Class"] + "_" + amrfinder_df["Subclass"] #combine columns in file to get conferred resistance out of gene name
    contig_numbers = amrfinder_df["Contig_id"].str.replace(sample_name, "").str.split("_").str[1] #Parse "Contig" column in gamma file
    genes = amrfinder_df["Gene_symbol"] + "_" + amrfinder_df["Accession_of_closest_sequence"] #Parse "Gene_symbol" and "Accession_of_closest_sequence" column in gamma file to get gene name and accession
    method = amrfinder_df["Method"] #get method used to detect gene
    # loop through list of genes to combine with conferred resistance and make back into a pandas series
    column_name = ["{}_({})".format(gene, conferred_resistance) for gene, conferred_resistance in zip(genes, conferred_resistances)]
    match_type = ['NT' if x[-1] in ['X', 'N'] else 'AA' if x[-1] in ['P', 'M'] else '' for x in amrfinder_df['Method']]
    # loop through list of gamma info to combine ifdnto "code" for ID%/%cov:contig# and make back into a pandas series
    coverage = ["[{:.0f}{}/{:.0f}:#{}:{}]".format(percent_codon_ID, match_type, percent_length, contig_number, method) for percent_codon_ID, match_type, percent_length, contig_number, method in zip(percent_codon_IDs, match_type, percent_lengths, contig_numbers, method)]
    # Minimum % length required to be included in report, otherwise removed from list
    if bool([percent_length for percent_length in percent_lengths if int(percent_length) < 90]): 
        index_remove_postion = [ n for n,percent_length in enumerate(percent_lengths) if int(percent_length) < 90 ] # get index for value removed to remove from other lists (values less than 90)
        percent_lengths = [percent_length for percent_length in percent_lengths if int(percent_length) >= 90] # filter list to remove values below cutoff (keep those greater than or equal to 90)
        for index in sorted(index_remove_postion, reverse=True):
            del coverage[index]
            del column_name[index]
            del percent_codon_IDs[index]
    # Minimum % identity required to be included in report, otherwise removed from list
    if bool([percent_codon_ID for percent_codon_ID in percent_codon_IDs if int(percent_codon_ID) < 98]):
        index_remove_postion = [ n for n,percent_codon_ID in enumerate(percent_codon_IDs) if int(percent_codon_ID) < 98 ] # get index for value removed to remove from other lists (values less than 98)
        percent_codon_IDs = [percent_codon_ID for percent_codon_ID in percent_codon_IDs if int(percent_codon_ID) >= 98] # filter list to remove values below cutoff (keep those greater than or equal to 98)
        #loop through list of indexes to delete and remove them from the other lists so they all match
        for index in sorted(index_remove_postion, reverse=True):
            del coverage[index]
            del column_name[index]
            del percent_lengths[index]
    #building a new dataframe - create giant row
    df = pd.DataFrame(coverage).T
    if df.empty:
        df = pd.DataFrame({'WGS_ID':[sample_name], 'AR_Database':[db], 'No_AR_Genes_Found':['[-/-]'] })
        df.index = [sample_name]
    else:
        df.columns = column_name # add column names
        df["WGS_ID"] = sample_name
        df["AR_Database"] = db
        df["No_AR_Genes_Found"] = ""
        df.index = [sample_name]
    # Check for duplicate column names, multiple hits 
    #print(df["blaFOX-5_NG_049105.1(beta-lactam)"])
    df = duplicate_column_clean(df)
    final_df = pd.concat([final_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    #print(final_df["mph(D)_NC_017312(macrolide_lincosamide_streptogramin)"])
    return final_df

def parse_ani(fast_ani_file):
    """Parse ANI file to get format 99.98%ID-98.58%COV-Acinetobacter baumannii(Acinetobacter_baumannii_GCF_012935145.1_ASM1293514v1_genomic.fna.gz)."""
    with open(fast_ani_file) as f:
        first_line = f.readline().strip('\n')
        #try: #try to see if there is a second line.
        #    second_line = f.readlines()[0].strip('\n')
        #except IndexError:
        #    second_line = ""
    if "No MASH hit found" in first_line:
        FastANI_output_list = ['NA','NA','NA','NA']
        fastani_warning = "No MASH hit found."
    elif "No hits above an ANI value >=80%" in first_line:
        FastANI_output_list = ['NA','NA','NA','NA']
        fastani_warning = "No hits with >=80% ANI."
    else:
        fastani_warning = None
        ani_df = pd.read_csv(fast_ani_file, sep='\t', header=0) # should only be one line long.
        ID = ani_df["% ID"][0]
        coverage = ani_df["% Coverage"][0]
        organism = ani_df["Organism"][0]
        source_file = ani_df["Source File"][0]
        #Species_Support = str(ID) + "%ID-" + str(coverage) + "%COV-" + organism + "(" + source_file + ")" #old way of reporting
        FastANI_output_list = [source_file, ID, coverage, organism]
    return FastANI_output_list, fastani_warning

def Get_Metrics(set_coverage, ar_df, trim_stats, raw_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld_report, kraken_wtasmbld, quast_report, busco_short_summary, asmbld_ratio, gc_file, sample_name, fairy_files, spades_fairy_file, amrfinder_file, fast_ani_file, tax_file, ar_db):
    '''For each step to gather metrics try to find the file and if not then make all variables unknown'''
    try:
        Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Raw_reads, Total_Trimmed_bp, Paired_Trimmed_reads, Total_Trimmed_reads, Trim_Q30_R1_percent, Trim_Q30_R2_percent = get_Q30(trim_stats, raw_stats)
    except FileNotFoundError:
        print("Warning: " + sample_name + "_trimmed_read_counts.txt not found")
        Q30_R1_per = Q30_R2_per = Total_Raw_Seq_bp = Total_Raw_reads = Total_Trimmed_bp = Paired_Trimmed_reads = Total_Trimmed_reads = Trim_Q30_R1_percent = Trim_Q30_R2_percent = 'Unknown'
    # Try and except are in the get_kraken_info function to allow for cases where trimming was completed, but not assembly
    Trim_kraken, Trim_Genus_percent, Asmbld_kraken, Asmbld_Genus_percent, Trim_unclassified_percent, Wt_asmbld_unclassified_percent = get_kraken_info(kraken_trim, kraken_wtasmbld, sample_name)
    try:
        if Total_Trimmed_bp != "Unknown": # for -entry CDC_SCAFFOLDS where no reads are present to calculate. For all other scenerios run this portion
            Coverage, Assembly_Length = Calculate_Trim_Coverage(Total_Trimmed_bp, quast_report)
        elif Total_Trimmed_bp == "Unknown": # for -entry CDC_SCAFFOLDS and -entry SCAFFOLDS
            Assembly_Length = Get_Assembly_Length(quast_report)
            Coverage = 'Unknown'
        else:
            Coverage = Assembly_Length = 'Unknown'
    except FileNotFoundError:
        print("Warning: " + sample_name + "_summary.tsv not found")
        Coverage = Assembly_Length = 'Unknown'
    try:
        Scaffold_Count = get_scaffold_count(quast_report)
    except FileNotFoundError:
        print("Warning: " + sample_name + "_summary.tsv not found")
        Scaffold_Count = 'Unknown'
    try:
        busco_metrics, full_busco_line = Get_BUSCO_Gene_Count(busco_short_summary)
    except FileNotFoundError:
        lineage = percent_busco = full_busco_line = 'Unknown'
        busco_metrics = [lineage, percent_busco]
    try:
        assembly_ratio_metrics = get_assembly_ratio(asmbld_ratio, tax_file)
    except FileNotFoundError:
        print("Warning: " + sample_name + "_Assembly_ratio_*.txt or " + sample_name + ".tax not found.")
        ratio = stdev = tax_method = 'Unknown'
        assembly_ratio_metrics = [ratio, stdev, tax_method]
    try:
        gc_metrics = get_gc_metrics(gc_file)
    except FileNotFoundError:
        print("Warning: " + sample_name + "_GC_content_20230504.txt not found.")
        gc_stdev = sample_gc = out_of_range_stdev = species_gc_mean = 'Unknown'
        gc_metrics = [gc_stdev, sample_gc, out_of_range_stdev, species_gc_mean]
    try:
        QC_result, QC_reason = Checking_auto_pass_fail(fairy_files, spades_fairy_file, Coverage, Assembly_Length, assembly_ratio_metrics[1], assembly_ratio_metrics[0], set_coverage, Scaffold_Count, sample_name)
    
    except FileNotFoundError:
        print("Warning: Possibly coverage and assembly length was not calculated and/or "+ sample_name + "_Assembly_ratio_*.txt not found.")
        QC_result = QC_reason = 'Unknown'
    try:
        FastANI_output_list, fastani_warning = parse_ani(fast_ani_file)
    except FileNotFoundError: 
        print("Warning: " + sample_name + ".fastANI.txt not found")
        ani_source_file = fastani_ID = fastani_coverage = fastani_organism = 'Unknown'
        FastANI_output_list = [ani_source_file, fastani_ID, fastani_coverage, fastani_organism]
    try:
        ar_df = parse_amrfinder_ar(amrfinder_file, sample_name, ar_df, ar_db)
    except FileNotFoundError: 
        print("Warning: AMRFinder file " + sample_name + "_all_genes_blank.tsv not found.")
        df = pd.DataFrame({'WGS_ID':[sample_name], 'No_AR_Genes_Found':['File not found'], 'AR_Database':['AMRFinder file not found'] })
        df.index = [sample_name]
        ar_df = pd.concat([ar_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    try:
        alerts = compile_alerts(Coverage, assembly_ratio_metrics[1], gc_metrics[0])
    except:
        alerts = ""
    # try except in the function itself
    kraken_trim_genus, kraken_wtasmbld_genus = parse_kraken_report(kraken_trim_report, kraken_wtasmbld_report, sample_name)
    try:
        warnings = compile_warnings(Total_Trimmed_reads, Total_Raw_reads, Q30_R1_per, Q30_R2_per, Trim_Q30_R1_percent, Trim_Q30_R2_percent,\
                                    Scaffold_Count, gc_metrics, assembly_ratio_metrics, Trim_unclassified_percent, Wt_asmbld_unclassified_percent,\
                                    kraken_trim_genus, kraken_wtasmbld_genus, Trim_Genus_percent, Asmbld_Genus_percent, 
                                    fastani_warning, busco_metrics[1], FastANI_output_list[1], FastANI_output_list[2], QC_reason)
    except:
        warnings = ""
    return ar_df, Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Raw_reads, Paired_Trimmed_reads, Total_Trimmed_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, \
    Scaffold_Count, busco_metrics, assembly_ratio_metrics, gc_metrics, QC_result, QC_reason, full_busco_line

def Get_Files(directory, sample_name):
    '''Create file paths to collect files from sample folder.'''
    # if there is a trailing / remove it
    directory = directory.rstrip('/')
    # create file names
    trim_stats = directory + "/" + sample_name + "_trimmed_read_counts.txt"
    raw_stats = directory + "/" + sample_name + "_raw_read_counts.txt"
    kraken_trim = directory + "/" + sample_name + ".kraken2_trimd.top_kraken_hit.txt"
    kraken_trim_report = directory + "/" + sample_name + ".kraken2_trimd.summary.txt"
    kraken_wtasmbld = directory + "/" + sample_name + ".kraken2_wtasmbld.top_kraken_hit.txt"
    kraken_wtasmbld_report = directory + "/" + sample_name + ".kraken2_wtasmbld.summary.txt"
    quast_report = directory + "/" + sample_name + "_summary.tsv"
    # For glob patterns, try both directories
    fairy_files = glob.glob(directory + "/" + sample_name + "*_summary.txt")
    spades_fairy_file_1 = glob.glob(directory + "/" + sample_name + "*_spades_outcome.csv")
    spades_fairy_file = spades_fairy_file_1[0] if spades_fairy_file_1 else ""
    # For the remaining glob patterns, handle with try-except but attempt both directories
    try:
        amrfinder_file1 = glob.glob(directory + "/" + sample_name + "_all_genes_*.tsv")
        if amrfinder_file1:
            amrfinder_file = amrfinder_file1[0]
        else:
            amrfinder_file = directory + "/" + sample_name + "_all_genes_blank.tsv"
    except IndexError:
        amrfinder_file = directory + "/" + sample_name + "_all_genes_blank.tsv"
    try:
        busco_short_summary_1 = glob.glob(directory + "/short_summary.specific.*" + sample_name + ".filtered.scaffolds.fa.txt")
        if busco_short_summary_1:
            busco_short_summary = busco_short_summary_1[0]
        else:
            busco_short_summary = directory + "/short_summary.specific.blank" + sample_name + ".filtered.scaffolds.fa.txt"
    except IndexError:
        busco_short_summary = directory + "/short_summary.specific.blank" + sample_name + ".filtered.scaffolds.fa.txt"
    try:
        asmbld_ratio_1 = glob.glob(directory + "/" + sample_name + "_Assembly_ratio_*.txt")
        if asmbld_ratio_1:
            asmbld_ratio = asmbld_ratio_1[0]
        else:
            asmbld_ratio = directory + "/" + sample_name + "_Assembly_ratio_blank.txt"
    except IndexError:
        asmbld_ratio = directory + "/" + sample_name + "_Assembly_ratio_blank.txt"
    try:
        gc_1 = glob.glob(directory + "/" + sample_name + "_GC_content_*.txt")
        if gc_1:
            gc = gc_1[0]
        else:
            gc = directory + "/" + sample_name + "_GC_content_blank.txt"
    except IndexError:
        gc = directory + "/" + sample_name + "_GC_content_blank.txt"
    try:
        fast_ani_file_1 = glob.glob(directory + "/" + sample_name + "_REFSEQ_*.fastANI.txt")
        if fast_ani_file_1:
            fast_ani_file = fast_ani_file_1[0]
        else:
            fast_ani_file = directory + "/" + sample_name + ".fastANI.txt"
    except IndexError:
        fast_ani_file = directory + "/" + sample_name + ".fastANI.txt"
    # For regular paths, use the try_paths function
    tax_file = directory + "/" + sample_name + ".tax"
    return trim_stats, raw_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld_report, kraken_wtasmbld, quast_report, fairy_files, spades_fairy_file, busco_short_summary, asmbld_ratio, gc, amrfinder_file, fast_ani_file, tax_file

def Append_Lists(data_location, parent_folder, sample_name, Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Seq_reads, Paired_Trimmed_reads, Total_trim_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, \
            Scaffold_Count, busco_metrics, gc_metrics, assembly_ratio_metrics, QC_result, QC_reason, full_busco_line,
            data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L,Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, full_busco_line_L):
        data_location_L.append(data_location)
        parent_folder_L.append(parent_folder)
        Sample_Names.append(str(sample_name))
        Q30_R1_per_L.append(Q30_R1_per)
        Q30_R2_per_L.append(Q30_R2_per)
        #Total_Raw_Seq_bp_L.append(Total_Raw_Seq_bp)
        Paired_Trimmed_reads_L.append(Paired_Trimmed_reads)
        Total_Seq_reads_L.append(Total_Seq_reads)
        Total_trim_Seq_reads_L.append(Total_trim_Seq_reads)
        Trim_kraken_L.append(Trim_kraken)
        Asmbld_kraken_L.append(Asmbld_kraken)
        Coverage_L.append(Coverage)
        Assembly_Length_L.append(Assembly_Length)
        Species_Support_L.append(FastANI_output_list[0])
        fastani_organism_L.append(FastANI_output_list[3])
        fastani_ID_L.append(FastANI_output_list[1])
        fastani_coverage_L.append(FastANI_output_list[2])
        gc_L.append(gc_metrics[1])
        Scaffold_Count_L.append(Scaffold_Count)
        busco_lineage_L.append(busco_metrics[0])
        percent_busco_L.append(busco_metrics[1])
        assembly_ratio_L.append(assembly_ratio_metrics[0])
        assembly_stdev_L.append(assembly_ratio_metrics[1])
        tax_method_L.append(assembly_ratio_metrics[2])
        QC_result_L.append(QC_result)
        QC_reason_L.append(QC_reason)
        warnings_L.append(warnings)
        alerts_L.append(alerts)
        full_busco_line_L.append(full_busco_line)
        return data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
        Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, full_busco_line_L

def Create_df(phx_version, data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L,
Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, full_busco_line_L):
    phx_version_L = [str(phx_version)] * len(Sample_Names)
    #combine all metrics into a dataframe
    data = {'WGS_ID'             : Sample_Names,
    'Parent_Folder'              : parent_folder_L,
    'Data_Location'              : data_location_L,
    'PHX_Version'                : phx_version_L,
    'Minimum_QC_Check'           : QC_result_L,
    'Minimum_QC_Issues'          : QC_reason_L,
    'Warnings'                   : warnings_L,
    'Alerts'                     : alerts_L,
    'Raw_Q30_R1_[%]'             : Q30_R1_per_L,
    'Raw_Q30_R2_[%]'             : Q30_R2_per_L,
    #'Total_Raw_[bp]'             : Total_Raw_Seq_bp_L,
    'Total_Raw_[reads]'          : Total_Seq_reads_L,
    'Paired_Trimmed_[reads]'     : Paired_Trimmed_reads_L,
    'Total_Trimmed_[reads]'      : Total_trim_Seq_reads_L,
    'Estimated_Trimmed_Coverage' : Coverage_L,
    'GC[%]'                      : gc_L,
    'Scaffolds'                  : Scaffold_Count_L,
    'Assembly_Length'            : Assembly_Length_L,
    'Assembly_Ratio'             : assembly_ratio_L,
    'Assembly_StDev'             : assembly_stdev_L,
    'Final_Taxa_ID'              : "", # we will fill this later
    'Taxa_Source'                : tax_method_L,
    'BUSCO_Lineage'              : busco_lineage_L,
    'BUSCO_%Match'               : percent_busco_L,
    'Kraken_ID_Trimmed_Reads_%'  : Trim_kraken_L,
    'Kraken_ID_WtAssembly_%'     : Asmbld_kraken_L,
    'FastANI_Organism'           : fastani_organism_L, 
    'FastANI_%ID'                : fastani_ID_L, 
    'FastANI_%Coverage'          : fastani_coverage_L,
    'Species_Support_ANI'        : Species_Support_L}
    busco_data = {'WGS_ID'       : Sample_Names,
                'BUSCO'          : full_busco_line_L}
    df = pd.DataFrame(data)
    busco_df = pd.DataFrame(busco_data)
    return df, busco_df

def find_big_5(BLDB):
    df = pd.read_csv(BLDB)
    # Filter rows where 'Protein_name' contains any of the substrings
    big5_genes = ["KPC", "IMP", "NDM", "OXA", "VIM"]
    filtered_df = df[df["Protein name"].str.contains('|'.join(big5_genes), case=False, na=False)]
    # \xa0 is from hyperlinks as there is not normal spaces # Further filter for 'carbapenemase' or 'IR carbapenemase' in the assumed column (e.g., "Classification")
    final_df = filtered_df[filtered_df["Functional information"].isin(["carbapenemase", "IR carbapenemase", "carbapenemase\xa0view", "IR carbapenemase\xa0view"])]
    # Condition to check if "Protein_name" contains "OXA"
    oxa_condition = final_df["Protein name"].str.contains("OXA", case=False, na=False)
    # Condition to filter "Subfamily" only for rows where "Protein_name" contains "OXA"
    subfamily_condition = final_df["Subfamily"].isin(["OXA-48-like", "OXA-23-like", "OXA-24-like", "OXA-58-like", "OXA-143-like"])
    # Keep all rows where "Protein_name" does NOT contain "OXA"
    non_oxa_rows = final_df[~oxa_condition]
    # Keep only filtered rows where "Protein_name" contains "OXA" and "Subfamily" is in the list
    filtered_oxa_rows1 = final_df[oxa_condition & subfamily_condition]
    ###########print(filtered_oxa_rows1)
    filtered_oxa_rows = filtered_oxa_rows1[~(filtered_oxa_rows1["Natural (N) or Acquired (A)"].str.contains(r"N\s\(", na=False) & ~filtered_oxa_rows1["Subfamily"].str.contains("OXA-48-like", na=False))]
    # Combine both DataFrames
    filtered_final_df = pd.concat([non_oxa_rows, filtered_oxa_rows])
    # Select the relevant columns and drop complete duplicates
    unique_proteins = final_df[["Protein name", "Alternative protein names"]].drop_duplicates()
    # Flatten into a list and remove NaN values
    protein_list = unique_proteins.values.flatten()
    protein_list = [protein for protein in protein_list if pd.notna(protein)]  # Remove NaNs
    # Separate protein_list into two lists
    oxa_proteins = [protein for protein in protein_list if "OXA" in protein]
    non_oxa_proteins = [protein for protein in protein_list if "OXA" not in protein]
    return non_oxa_proteins, oxa_proteins

def big5_check(final_ar_df, BLDB):
    """"Function that will return list of columns to highlight if a sample has a hit for a big 5 gene."""
    columns_to_highlight = []
    final_ar_df = final_ar_df.drop(['AR_Database','WGS_ID'], axis=1)
    all_genes = final_ar_df.columns.tolist()
    big5_keep, big5_oxa_keep= find_big_5(BLDB)
    # loop through column names and check if they contain a gene we want highlighted. Then add to highlight list if they do. 
    for gene in all_genes: # loop through each gene in the dataframe of genes found in all isolates
        if gene == 'No_AR_Genes_Found':
            pass
        else:
            gene_name = gene.split('_(')[0] # remove drug name for matching genes
            drug = gene.split('_(')[1] # keep drug name to add back later
            if "-like" in gene_name:
                gene_name = gene_name.split('_bla')[0] # remove blaOXA-1-like name for matching genes -- just extra stuff that doesn't allow complete match
            # make sure we have a complete match for oxa 48/23/24/58/143 genes and oxa 48/23/24/58/143-like genes
            if "OXA" in gene_name: #check for complete blaOXA match
                [ columns_to_highlight.append(gene_name + "_(" + drug) for big5_oxa in big5_oxa_keep if gene_name == big5_oxa ]
            else: # for "blaIMP", "blaVIM", "blaNDM", and "blaKPC", this will take any thing with a matching substring to these
                [ columns_to_highlight.append(gene_name + "_(" + drug) for big5 in big5_keep if search(big5, gene_name) ]
    print(CYELLOW + "\nhighlighting colums:", columns_to_highlight, CEND)
    return columns_to_highlight

def column_letter(index):
    """Convert zero-based column index to Excel column letter."""
    letters = list(string.ascii_uppercase)
    if index < 26:
        return letters[index]
    else:
        return letters[index // 26 - 1] + letters[index % 26]  # Handle AA, AB, etc.

def write_to_excel(set_coverage, output, df, qc_max_col, ar_gene_count, columns_to_highlight, ar_df, ar_db):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    if output != "":
        writer = pd.ExcelWriter((output + '.xlsx'), engine='xlsxwriter')
    else:
        writer = pd.ExcelWriter(('GRiPHin_Summary.xlsx'), engine='xlsxwriter')
    # Convert the dataframe to an XlsxWriter Excel object.
    df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1)
    # Get the xlsxwriter workfbook worksheet objects for formating
    workbook = writer.book
    (max_row, max_col) = df.shape # Get the dimensions of the dataframe.
    worksheet = writer.sheets['Sheet1']
    # Setting columns to numbers so you can have commas that make it more human readable
    number_comma_format = workbook.add_format({'num_format': '#,##0'})
    ##worksheet.set_column('H:J', None, number_comma_format) # Total_seqs (raw and trimmed) Total_bp - use when python is >3.7.12
    ##worksheet.set_column('M:N', None, number_comma_format) # scaffolds and assembly_length - use when python is >3.7.12
    # set formating for python 3.7.12
    worksheet.conditional_format('O3:P' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_comma_format})
    worksheet.conditional_format('J3:L' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_comma_format})
    # Setting columns to float so its more human readable
    #number_dec_format = workbook.add_format({'num_format': '0.000'})
    #number_dec_format.set_align('left') - agh not working
    number_dec_2_format = workbook.add_format({'num_format': '0.00'})
    ##worksheet.set_column('K:L', None, number_dec_2_format) # Estimated_Trimmed_Coverage and GC% - use when python is >3.7.12
    ##worksheet.set_column('F:G', None, number_dec_2_format) # Q30%s - use when python is >3.7.12
    ##worksheet.set_column('O:P', None, number_dec_2_format) # Assembly Ratio and StDev - use when python is >3.7.12
    # set formating for python 3.7.12
    worksheet.conditional_format('M3:N' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
    worksheet.conditional_format('H3:I' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
    worksheet.conditional_format('Q3:R' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
    worksheet.conditional_format('Y3:Z' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
    #worksheet.set_column('W:X', None, number_dec_2_format) # FastANI ID and Coverage
    # getting values to set column widths automatically
    for idx, col in enumerate(df):  # loop through all columns
        series = df[col]
        #print(series)  # uncomment this to see what the issue is with the "mad" error
        if col == "Parent_Folder": #name is real long due to path so just making the width of the column header
            max_len = len("Parent_Folder")
        else:
            max_len = max((
            series.astype(str).map(len).max(),  # len of largest item
                len(str(series.name))  # len of column name/header
                )) + 1  # adding a little extra space
        worksheet.set_column(idx, idx, max_len)  # set column width
    # Setting colors for headers
    cell_format_light_blue = workbook.add_format({'bg_color': '#ADD8E6', 'font_color': '#000000', 'bold': True})
    cell_format_grey_blue = workbook.add_format({'bg_color': '#72A0C1', 'font_color': '#000000', 'bold': True})
    cell_format_green_blue = workbook.add_format({'bg_color': '#A3C1AD', 'font_color': '#000000', 'bold': True})
    cell_format_green = workbook.add_format({'bg_color': '#ACE1AF', 'font_color': '#000000', 'bold': True})
    cell_format_lightgrey = workbook.add_format({'bg_color': '#D5D8DC', 'font_color': '#000000', 'bold': True})
    cell_format_grey = workbook.add_format({'bg_color': '#AEB6BF', 'font_color': '#000000', 'bold': True})
    cell_format_darkgrey = workbook.add_format({'bg_color': '#808B96', 'font_color': '#000000', 'bold': True})
    # Headers
    #worksheet.set_column('', "PHoeNIx Summary", cell_format_light_blue)
    #worksheet.write('A1', "PHoeNIx Summary") #use for only 1 column in length
    #worksheet.set_column('A1:A1', None, cell_format_light_blue) #make summary column blue, #use for only 1 column in length
    worksheet.merge_range('A1:D1', "PHoeNIx Summary", cell_format_light_blue)
    worksheet.merge_range('E1:S1', "QC Metrics", cell_format_grey_blue)
    #taxa column 
    taxa_start_col  = column_letter(list(df.columns).index("Final_Taxa_ID"))  # Get index of start column
    taxa_end_col = column_letter(list(df.columns).index("Species_Support_ANI"))  # Get index of end column
    # Dynamically merge based on start and end column
    worksheet.merge_range(f"{taxa_start_col}1:{taxa_end_col}1", "Taxonomic Information", cell_format_green)
    worksheet.merge_range(0, qc_max_col, 0, (qc_max_col + ar_gene_count - 1), "Antibiotic Resistance Genes", cell_format_lightgrey)
    # making WGS IDs bold
    bold = workbook.add_format({'bold': True})
    worksheet.set_column('A3:A' + str(max_row + 2), None, bold)
    # Setting colors for footers and conditional formating
    yellow_format = workbook.add_format({'bg_color': '#FFEB9C', 'font_color': '#000000'}) # Light yellow fill with black text.
    red_format = workbook.add_format({'bg_color': '#F5B7B1', 'font_color': '#000000'}) # red fill with black text.
    orange_format = workbook.add_format({'bg_color': '#F5CBA7', 'font_color': '#000000', 'bold': True})
    orange_format.set_border(1) # add back border so it matches the rest of the column names
    orange_format_nb = workbook.add_format({'bg_color': '#F5CBA7', 'font_color': '#000000', 'bold': False})
    # Apply a conditional format for checking coverage is between set_coverage-100 in estimated coverage column. adding 2 to max row to account for headers
    worksheet.conditional_format('M3:M' + str(max_row + 2), {'type': 'cell', 'criteria': '<', 'value':  str(set_coverage), 'format': yellow_format})
    worksheet.conditional_format('M3:M' + str(max_row + 2), {'type': 'cell', 'criteria': '>', 'value':  100.00, 'format': yellow_format})
    # Apply a conditional format for auto pass/fail in Auto_PassFail coverage column.
    worksheet.conditional_format('D3:D' + str(max_row + 2), {'type': 'cell', 'criteria': 'equal to', 'value':  '"FAIL"', 'format': red_format})
    # conditional formating to highlight big 5 genes
    # Start iterating through the columns and the rows to apply the format
    column_count = 0
    for column in ar_df.columns:
        for gene in columns_to_highlight:
            if column == gene: # if the column is one of the big 5 genes to highlight
                col_adjustment = column_count + qc_max_col - 1 # adjust starting place to account for qc columns 
                cell = xl_rowcol_to_cell(1, col_adjustment)   # Gets the excel location like A1
                #cell_value = ar_df.iloc[2, column_count] # get the value in that cell
                worksheet.write(cell, column, orange_format)
        column_count = column_count + 1
    ##            for row in range(ar_df.shape[0]):
    ##                col_adjustment = column_count + qc_max_col - 1 # adjust starting place to account for qc columns 
    ##                row_adjustment = row + 2
    ##                cell = xl_rowcol_to_cell(row_adjustment, col_adjustment)   # Gets the excel location like A1
    ##                cell_value = ar_df.iloc[row, column_count] # get the value in that cell
    ##                if cell_value != "":
    ##                    worksheet.write(cell, cell_value, orange_format)
    ##    column_count = column_count + 1
    # Creating footers
    worksheet.write('A' + str(max_row + 4), 'Cells in YELLOW denote isolates outside of ' + str(set_coverage) + '-100X coverage', yellow_format)
    worksheet.write('A' + str(max_row + 5), 'Cells in ORANGE denote “Big 5” carbapenemase gene (i.e., blaKPC, blaNDM, blaOXA-48-like, blaVIM, and blaIMP) or an acquired blaOXA gene, please confirm what AR Lab Network HAI/AR WGS priority these meet.', orange_format_nb)
    worksheet.write('A' + str(max_row + 6), 'Cells in RED denote isolates that failed one or more auto failure triggers (cov < 30, assembly ratio stdev > 2.58, assembly length < 1Mbps)', red_format)
    # More footers - Disclaimer etc.
    # unbold
    no_bold = workbook.add_format({'bold': False})
    worksheet.write('A' + str(max_row + 7),"^Using Antibiotic Resistance Gene database " + ar_db + " AMRFinder Antimicrobial Resistance Reference Gene Database curated by NCBI using default thresholds see below. Output is shown as  [%Identity (NT) or (AA), %Coverage, Contig#, Method of Detection].", no_bold)
    worksheet.write('A' + str(max_row + 8),"AMRFinder Methods Explained:", no_bold)
    worksheet.write('A' + str(max_row + 9),"ALLELE - 100% sequence match over 100% of length to a protein named at the allele level in the AMRFinderPlus database.", no_bold)
    worksheet.write('A' + str(max_row + 10),"EXACT - 100% sequence match over 100% of length to a protein in the database that is not a named allele.", no_bold)
    worksheet.write('A' + str(max_row + 11),"BLAST - BLAST alignment is > 90% of length and > 90% identity to a protein in the AMRFinderPlus database.", no_bold)
    worksheet.write('A' + str(max_row + 12),"PARTIAL - BLAST alignment is > 50% of length, but < 90% of length and > 90% identity to the reference, and does not end at a contig boundary.", no_bold)
    worksheet.write('A' + str(max_row + 13),"PARTIAL_CONTIG_END - BLAST alignment is > 50% of length, but < 90% of length and > 90% identity to the reference, and the break occurrs at a contig boundary indicating that this gene is more likely to have been split by an assembly issue.", no_bold)
    worksheet.write('A' + str(max_row + 14),"HMM - HMM was hit above the cutoff, but there was not a BLAST hit that met standards for BLAST or PARTIAL. This does not have a suffix because only protein sequences are searched by HMM.", no_bold)
    worksheet.write('A' + str(max_row + 15),"INTERNAL_STOP - Translated blast reveals a stop codon that occurred before the end of the protein.", no_bold)
    worksheet.write('A' + str(max_row + 16),"POINT - Point mutation identified by blast.", no_bold)
    worksheet.write('A' + str(max_row + 17),"DISCLAIMER: These data are preliminary and subject to change. The identification methods used and the data summarized are for public health surveillance or investigational purposes only and must NOT be communicated to the patient, their care provider, or placed in the patient’s medical record. These results should NOT be used for diagnosis, treatment, or assessment of individual patient health or management.", bold)
    worksheet.write('A' + str(max_row + 18),"The suffix indicates whether the gene was identified in a protein (P) or nucleotide translated (X) input file. Nucleotide BLAST (BLASTN) is used to identify nucleotide point mutations and POINTN is used as the Method for those hits.", no_bold)

    #adding review and date info
    worksheet.write('A' + str(max_row + 20), "Reviewed by:", no_bold)
    worksheet.write('A' + str(max_row + 21), "Date:", no_bold)
    # add autofilter
    worksheet.autofilter(1, 0, max_row, max_col - 1)
    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

def blind_samples(final_df, control_file):
    """If you passed a file to -c this will swap out sample names to 'blind' the WGS_IDs in the final excel file."""
    with open(control_file, 'r') as controls:
        header = next(controls) # skip the first line of the samplesheet
        for line in controls:
            old_sample_name = line.split(",")[0]
            new_sample_name = line.split(",")[1].rstrip("\n")
            if new_sample_name != old_sample_name:
                final_df['WGS_ID'] = final_df['WGS_ID'].replace(old_sample_name, new_sample_name)
    # create new csv file
    return final_df

def create_samplesheet(directory):
    """Function will create a samplesheet from samples in a directory if -d argument passed."""
    directory = os.path.abspath(directory) # make sure we have an absolute path to start with
    with open("Directory_samplesheet.csv", "w") as samplesheet:
        samplesheet.write('sample,directory\n')
    dirs = sorted(os.listdir(directory))
    # If there are any new files added to the top directory they will need to be added here or you will get an error
    skip_list_a1 = glob.glob(directory + "/*_GRiPHin_Summary.*") # for if griphin is run on a folder that already has a report in it
    skip_list_a2 = glob.glob(directory + "/*_comparison") # for comparinator script
    skip_list_a3 = glob.glob(directory + "/*.html") # for clia files
    skip_list_a4 = glob.glob(directory + "/*.pdf") # for clia files
    skip_list_a = skip_list_a1 + skip_list_a2 + skip_list_a3 + skip_list_a4
    skip_list_a = [ gene.split('/')[-1] for gene in skip_list_a ]  # just get the excel name not the full path
    skip_list_b = ["BiosampleAttributes_Microbe.1.0.xlsx", "Sra_Microbe.1.0.xlsx", "Phoenix_Summary.tsv", "pipeline_info", "GRiPHin_Summary.xlsx", "multiqc", "samplesheet_converted.csv", "Directory_samplesheet.csv", "sra_samplesheet.csv"]
    skip_list = skip_list_a + skip_list_b
    #remove unwanted files
    dirs_cleaned = [item for item in dirs if item not in skip_list]
    try: #if there are numbers in the name then use that to sort
        dirs_sorted=sorted(dirs_cleaned, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    except: #if no numbers then use only alphabetically
        dirs_sorted=sorted(dirs_cleaned)
    for sample in dirs_sorted:
        with open("Directory_samplesheet.csv", "a") as samplesheet:
            if directory[-1] != "/": # if directory doesn't have trailing / add one
                directory = directory + "/"
            samplesheet.write(sample + "," + directory + sample + '\n')
    samplesheet = "Directory_samplesheet.csv"
    return samplesheet

def sort_samplesheet(samplesheet):
    df = pd.read_csv(samplesheet)
    #get list of sample ids to sort
    samples = df["sample"]
    try: #if there are numbers in the name then use that to sort
        samples_sorted=sorted(samples, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    except: #if no numbers then use only alphabetically
        samples_sorted=sorted(samples)
    df = df.set_index("sample")
    df = df.loc[samples_sorted]
    df.to_csv(samplesheet, sep=',', encoding='utf-8') #overwrite file

def convert_excel_to_tsv(output):
    '''Reads in the xlsx file that was just created, outputs as tsv version with first layer of headers removed'''
    if output != "":
        output_file = output
    else:
        output_file = 'GRiPHin_Summary'
    #Read excel file into a dataframe
    data_xlsx = pd.read_excel(output_file + '.xlsx', 'Sheet1', index_col=None, header=[1])
    #Replace all fields having line breaks with space
    #data_xlsx = data_xlsx.replace('\n', ' ',regex=True)
    #drop the footer information
    data_xlsx = data_xlsx.iloc[:-11] 
    #Write dataframe into csv
    data_xlsx.to_csv(output_file + '.tsv', sep='\t', encoding='utf-8',  index=False, lineterminator ='\n')

def clean_ar_df(df, ar_df, BLDB):
    ar_cols = list(ar_df)
    # move the column to head of list using index, pop and insert
    ar_cols.insert(0, ar_cols.pop(ar_cols.index('No_AR_Genes_Found')))
    ar_cols.insert(0, ar_cols.pop(ar_cols.index('AR_Database')))
    # use ix to reorder
    final_ar_df = ar_df.loc[:, ar_cols]
    ar_max_col = final_ar_df.shape[1] - 1 #remove one for the WGS_ID column
    # now we will check for the "big 5" genes for highlighting later.
    columns_to_highlight = big5_check(final_ar_df, BLDB)
    # combining all dataframes
    final_df = pd.merge(df, final_ar_df, how="left", on=["WGS_ID","WGS_ID"])
    #get database names and remove if file is not found in the database list
    ar_db = final_df['AR_Database'].unique().tolist()
    if 'AMRFinder file not found' in ar_db:
        ar_db.remove('AMRFinder file not found') #Don't want this reported as the ar_db
    ar_db = ",".join(ar_db)
    return final_df, ar_max_col, columns_to_highlight, final_ar_df, ar_db

def write_phoenix_summary(set_coverage, final_df, busco_df):
    #full_busco_line_L = full_busco_line_L.rename(columns=column_name_mapping)
    columns_to_remove = ['Parent_Folder', 'Data_Location', "No_AR_Genes_Found"]
    final_df = final_df.drop(columns=columns_to_remove)
    #add in busco line we want
    final_df = pd.merge(final_df, busco_df, on='WGS_ID')
    # Dictionary mapping old column names in GRiPHin to new column names in phoenix_summary.tsv
    column_name_mapping = {'WGS_ID':'ID','Minimum_QC_Check':'Auto_QC_Outcome', 'Estimated_Trimmed_Coverage':'Estimated_Coverage', 'Assembly_Length':'Genome_Length',
                            'GC[%]':'GC_%', 'BUSCO_Lineage':'BUSCO_DB', 'Assembly_StDev':'Assembly_Ratio_(STDev)', 'Scaffolds':'#_of_Scaffolds_>500bp',
                            'Kraken_ID_Trimmed_Reads_%':'Kraken2_Trimd', 'Kraken_ID_WtAssembly_%':'Kraken2_Weighted','Minimum_QC_Issues': 'Auto_QC_Failure_Reason','BUSCO_%Match':'BUSCO_Match'}
                            #'Primary_MLST_Scheme':'MLST_Scheme_1', 'Primary_MLST':'MLST_1','Secondary_MLST_Scheme':'MLST_2 }
    # Rename multiple columns
    final_df = final_df.rename(columns=column_name_mapping)
    column_list = final_df.columns.tolist()
    # get Resistance genes - ( in the column name and is not a digit
    cols_to_drop = ['ID', 'Auto_QC_Outcome', 'Auto_QC_Failure_Reason', 'Warnings', 'Alerts', 'Raw_Q30_R1_[%]', 'Raw_Q30_R2_[%]', 'Total_Raw_[reads]',
       'Paired_Trimmed_[reads]', 'Total_Trimmed_[reads]','Estimated_Coverage', 'GC_%', '#_of_Scaffolds_>500bp', 'Genome_Length',
       'Assembly_Ratio', 'Assembly_Ratio_(STDev)', 'Taxa_Source', 'BUSCO_DB', 'BUSCO_Match', 'BUSCO', 'Kraken2_Trimd', 'Kraken2_Weighted',
       'FastANI_Organism', 'FastANI_%ID', 'FastANI_%Coverage','Species_Support_ANI', 'AR_Database']
    # remove what we don't want
    #column_list = [item for item in column_list if item not in cols_to_drop]
    pattern = re.compile(r'\((.*?)\)$')  # Pattern to match the last string between "(" and ")" at the end of the string
    resistance_types = [pattern.search(col_name).group(1) for col_name in column_list if pattern.search(col_name)]
    #pattern = re.compile(r'\w{2}_\d+\.\d+')
    #resistance_types = [col_name.split(re.search(pattern, col_name).group(),1)[1].replace('_(', '').replace(')', '') for col_name in column_list ]
    #create ar db
    #ar_db_df = pd.DataFrame()
    #for resistance_type in resistance_types:
        #get all column names that match resistance type
    #    matching_col_resistance_type = [col for col in final_ar_df.columns if resistance_type in col]
        # Handle multiple matching columns and combining gene names with ,
    #    ar_db_df[resistance_type] = final_ar_df.apply(lambda row: ', '.join([col.split('_')[0] for col in matching_col_resistance_type if row[col] != ""]), axis=1)
    # Set make index into ID column
    #ar_db_df = ar_db_df.reset_index().rename(columns={'index':'ID'})
    # Apply the function to create the new Taxa_Confidence column
    final_df['Taxa_Confidence'] = final_df.apply(extract_confidence, axis=1)
    # Apply the function to create the new Species column
    final_df['Species'] = final_df.apply(extract_species, axis=1)
    # make column for warnings count
    final_df['Warning_Count'] = final_df['Warnings'].str.split(',').str.len()
    #make new dataframe
    final_df = final_df[['ID','PHX_Version','Auto_QC_Outcome','Warning_Count','Estimated_Coverage','Genome_Length','Assembly_Ratio_(STDev)','#_of_Scaffolds_>500bp','GC_%','BUSCO','BUSCO_DB',
                       'Final_Taxa_ID', 'Taxa_Source', 'FastANI_Organism','FastANI_%ID','FastANI_%Coverage','ShigaPass_Organism','Kraken2_Trimd','Kraken2_Weighted', 'Auto_QC_Failure_Reason']]#,'MLST_Scheme_1','MLST_1','MLST_Scheme_2','MLST_2']]
    #phx_df = pd.merge(phx_df, ar_db_df, on='ID')
    # add commas to make it more readable
    final_df["Genome_Length"] = final_df["Genome_Length"].apply(lambda x: "{:,.0f}".format(x) if pd.notna(x) and x != "Unknown" else x)
    #print out Phoenix_Summary.tsv
    final_df.to_csv('Phoenix_Summary.tsv', sep='\t', index=False)

def extract_species(row):
    """Define a function to extract the relevant information based on Taxa_Source."""
    if row['Taxa_Source'] == 'ANI_REFSEQ':
        return f"{row['FastANI_Organism']}"
    elif row['Taxa_Source'] == 'kraken2_wtasmbld':
        # Extract the first number in parentheses using regular expression
        scaffolds_match = re.search(r'\((\d+(\.\d+)?)', row['Kraken2_Weighted'])
        return f"{scaffolds_match.group(1)}"
    elif row['Taxa_Source'] == 'kraken2_trimmed':
        reads_match = re.search(r'\((\d+(\.\d+)?)', row['Kraken2_Trimd'])
        return f"{reads_match.group(1)}"
    elif row['Taxa_Source'] == 'Unknown':
        return "Unknown"
    else:
        return "This is a bug open a github issue"

def extract_confidence(row):
    """Define a function to extract the relevant information based on Taxa_Source."""
    if row['Taxa_Source'] == 'ANI_REFSEQ':
        return row['FastANI_%ID']
    elif row['Taxa_Source'] == 'kraken2_wtasmbld':
        return row['Kraken2_Weighted'].replace(r'\([^)]*\)', '')
    elif row['Taxa_Source'] == 'kraken2_trimmed':
        return row['Kraken2_Trimd'].replace(r'\([^)]*\)', '')
    elif row['Taxa_Source'] == 'Unknown':
        return "Unknown"
    else:
        return "This is a bug open a github issue"

def main():
    args = parseArgs()
    shiga_df = pd.DataFrame()
    # create empty lists to append to later
    Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, Scaffold_Count_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
    busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, data_location_L, parent_folder_L, full_busco_line_L = ([] for i in range(29))
    ar_df = pd.DataFrame() #create empty dataframe to fill later for AR genes
    # check if a directory or samplesheet was given
    if (args.samplesheet == None) and (args.directory == None): # if no directory give AND no sample sheet given exit
        sys.exit(CRED + "You MUST pass EITHER a samplesheet or a top directory of PHoeNIx output to create one.\n" + CEND)
    # If a directory is given then create a samplesheet from it if not use the samplesheet passed
    if args.directory !=None:
        samplesheet = create_samplesheet(args.directory)
    else:
        sort_samplesheet(args.samplesheet)
        samplesheet = args.samplesheet
    #input is a samplesheet that is "samplename,directory" where the directory is a phoenix like folder
    with open(samplesheet) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        header = next(csv_reader) # skip the first line of the samplesheet
        for row in csv_reader:
            sample_name = row[0]
            directory = "./GRiPHin/" + sample_name + "/" # all samples should be in this directory
            data_location, parent_folder = Get_Parent_Folder(directory)
            trim_stats, raw_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld_report, kraken_wtasmbld, quast_report, fairy_files, spades_fairy_file, busco_short_summary, asmbld_ratio, gc, amrfinder_file, fastani_file, tax_file = Get_Files(directory, sample_name)
            #Get the metrics for the sample
            ar_df, Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Seq_reads, Paired_Trimmed_reads, Total_trim_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, Scaffold_Count, busco_metrics, assembly_ratio_metrics, gc_metrics, QC_result, \
            QC_reason, full_busco_line = Get_Metrics(args.set_coverage, ar_df, trim_stats, raw_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld_report, kraken_wtasmbld, quast_report, busco_short_summary, asmbld_ratio, gc, sample_name, fairy_files, spades_fairy_file, amrfinder_file, fastani_file, tax_file, args.ar_db)
            #Collect this mess of variables into appeneded lists
            data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L , alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, full_busco_line_L = Append_Lists(data_location, parent_folder, sample_name, \
            Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Seq_reads, Paired_Trimmed_reads, Total_trim_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, Scaffold_Count, busco_metrics, \
            gc_metrics, assembly_ratio_metrics, QC_result, QC_reason, full_busco_line, \
            data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, full_busco_line_L)
            if args.shigapass == True:
                shiga_df = create_shiga_df(directory, sample_name, shiga_df, FastANI_output_list[3], "")
    # combine all lists into a dataframe
    df, busco_df = Create_df(args.phx_version, data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
    Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, full_busco_line_L)
    if args.shigapass == True:
        df = double_check_taxa_id(shiga_df, df)
    else:
        df['Final_Taxa_ID'] = df.apply(fill_taxa_id, axis=1)
    (qc_max_row, qc_max_col) = df.shape
    final_df, ar_max_col, columns_to_highlight, final_ar_df, ar_db_col = clean_ar_df(df, ar_df, args.bldb)
    # Checking if there was a control sheet submitted
    if args.control_list !=None:
        final_df = blind_samples(final_df, args.control_list)
    else:
        final_df = final_df
    write_to_excel(args.set_coverage, args.output, final_df, qc_max_col, ar_max_col, columns_to_highlight, final_ar_df, ar_db_col)
    convert_excel_to_tsv(args.output)
    write_phoenix_summary(args.set_coverage, final_df, busco_df)

if __name__ == '__main__':
    main()