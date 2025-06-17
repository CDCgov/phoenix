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
from pathlib import Path
from species_specific_griphin import clean_and_format_centar_dfs, create_centar_combined_df, transform_value, create_shiga_df, double_check_taxa_id, fill_taxa_id

# Set display options to show all rows and columns
pd.set_option('display.max_rows', None)  # Show all rows
pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.max_colwidth', None)  # Show all columns

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python GRiPHin.py -s ./samplesheet.csv -a ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output --phoenix --scaffolds
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "2.1.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-s', '--samplesheet', default=None, required=False, dest='samplesheet', help='PHoeNIx style samplesheet of sample,directory in csv format. Directory is expected to have PHoeNIx stype output.')
    parser.add_argument('-b', '--bldb', default=None, required=False, dest='bldb', help='If a directory is given rather than samplesheet GRiPHin will create one for all samples in the directory.')
    parser.add_argument('-d', '--directory', default=None, required=False, dest='directory', help='If a directory is given rather than samplesheet GRiPHin will create one for all samples in the directory.')
    parser.add_argument('-c', '--control_list', required=False, dest='control_list', help='CSV file with a list of sample_name,new_name. This option will output the new_name rather than the sample name to "blind" reports.')
    parser.add_argument('-a', '--ar_db', default=None, required=True, dest='ar_db', help='AR Gene Database file that is used to confirm srst2 gene names are the same as GAMMAs output.')
    parser.add_argument('-o', '--output', default="", required=False, dest='output', help='Name of output file default is GRiPHin_Summary.xlsx.')
    parser.add_argument('--phx_version', default="Unknown", required=False, dest='phx_version', help='The version of phx used to produce GRiPHin_Summary row for the sample.')
    parser.add_argument('--coverage', default=30, required=False, dest='set_coverage', help='The coverage cut off default is 30x.')
    parser.add_argument('--scaffolds', dest="scaffolds", default=False, action='store_true', help='Turn on with --scaffolds to keep samples from failing/warnings/alerts that are based on trimmed data. Default is off.')
    parser.add_argument('--updater', dest="updater", default=False, action='store_true', help='When passed files locations are checked in two locations on in -d and in the dir listed in samplesheet.valid.csv.')
    parser.add_argument('--phoenix', dest="phoenix", default=False, action='store_true', required=False, help='Use for -entry PHOENIX rather than CDC_PHOENIX, which is the default.')
    parser.add_argument('--shigapass', dest="shigapass", default=False, action='store_true', required=False, help='Use for when there are E. coli or Shigella isolates in samplesheet.')
    parser.add_argument('--centar', dest="centar", default=False, action='store_true', required=False, help='Use for when there are C. diff isolates in samplesheet.')
    parser.add_argument('--filter_samples', dest="filter_samples", default=False, action='store_true', required=False, help='Use for when there are C. diff isolates in samplesheet.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CYELLOW = '\033[93m'
CEND = '\033[0m'

def print_df(df_toprint, label, all):
    print(label+"  ------------------------------------------------------")
    print("Columns:", df_toprint.columns)
    if all == True:
        with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 1000, 'display.colheader_justify', 'center', 'display.precision', 2, 'display.max_colwidth', 100):  # more options can be specified also
            print(df_toprint)

def Get_Parent_Folder(directory):
    '''getting project and parent_folder info from the paths'''
    #Project - parent folder (first folder that is in the outdir)
    #relative/submission - rest of the path
    #first make sure we have an absolute path
    directory = os.path.abspath(directory)
    #handing if trailing backslash isn't in there.
    if directory[-1] == "\\":
        directory = directory[::-1]
    path=Path(directory)
    project = os.path.basename(os.path.dirname(directory))
    # get project from directory path
    #project = os.path.split(os.path.split(os.path.split(directory)[0])[0])[1]
    parent_folder = path.parent.parent.resolve() # added resolve to handle symlinks rather than using .absolute()
    # get everything after CEMB
    #parent_folder = os.path.split(os.path.split(os.path.split(os.path.split(directory)[0])[0])[0])[0]
    return project, parent_folder

def make_ar_dictionary(ar_db):
    seq_id_list = []
    ar_dic = {}
    for seq_record in SeqIO.parse(ar_db, "fasta"):
        seq_id_list.append(seq_record.id)
    gene_name_list = [seq_id.split("__")[2] for seq_id in seq_id_list]
    drug_list = [seq_id.split("__")[4]  for seq_id in seq_id_list]
    ar_dic = dict(zip(gene_name_list, drug_list))
    return ar_dic

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
        # Need to add check for cases when Krakens DB does not have Genus or Species level hits, leaving a blank spot that crashes downstream
        if len(Trim_Genus) == 0:
            Trim_Genus = "None"
        if len(Trim_Species) == 0:
            Trim_Species = "None"
        Trim_kraken = Trim_Genus + " (" + Trim_Genus_percent + ") " + Trim_Species + " (" + Trim_Species_percent + ")."
        #guess what the mlst scheme is to check later
        scheme_guess_kraken_trimd = Trim_Genus[0].lower() + Trim_Species[0:4]
    except FileNotFoundError:
        print("Warning: " + sample_name + ".trimd_summary.txt not found")
        Trim_kraken = 'Unknown'
        Trim_unclassified_percent = "Unknown"
        Trim_Genus_percent = "Unknown"
        scheme_guess_kraken_trimd = ""
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
        # Need to add check for cases when Krakens DB does not have Genus or Species level hits, leaving a blank spot that crashes downstream
        if len(Asmbld_Genus) == 0:
            Asmbld_Genus = "None"
        if len(Asmbld_Species) == 0:
            Asmbld_Species = "None"
        Asmbld_kraken = Asmbld_Genus + " (" + Asmbld_Genus_percent + ") " + Asmbld_Species + " (" + Asmbld_Species_percent + ")."
        #guess what the mlst scheme is to check later
        scheme_guess_kraken_wt = Asmbld_Genus[0].lower() + Asmbld_Species[0:4]
    except FileNotFoundError:
        print("Warning: " + sample_name + ".wtasmbld_summary.txt not found")
        Asmbld_kraken = 'Unknown'
        Asmbld_unclassified_percent = "Unknown"
        Asmbld_Genus_percent = "Unknown"
        scheme_guess_kraken_wt = ""
    return Trim_kraken, Trim_Genus_percent, Asmbld_kraken, Asmbld_Genus_percent, Trim_unclassified_percent, Asmbld_unclassified_percent, scheme_guess_kraken_wt, scheme_guess_kraken_trimd

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
    #busco_line = lineage + " (" + percent_busco + "%)" # old busco line
    busco_metrics =  [lineage, percent_busco]
    return busco_metrics

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
                if "No Match Found" in line:
                    sample_gc="NA"
                else:
                    sample_gc = float((line.split("Sample_GC_Percent: ",1)[1]).strip())
            elif "Species_GC_Mean:" in line:
                if "No Match Found" in line:
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
                taxa =  (line.split("Tax: ",1)[1]).strip()
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

def compile_alerts(scaffolds_entry, coverage, assembly_stdev, gc_stdev):
    """
    No orphaned reads found after trimming
    <10 reference genomes for species identified so no stdev for assembly ratio or %GC content calculated
    >150x coverage or <40x coverage
    """
    alerts = []
    if scaffolds_entry == False:
        if coverage != "Unknown": # if its unknown it will fail already so skip
            if int(coverage) > 30 and int(coverage) < 40:
                alerts.append("coverage between 30-40x("+ str(coverage) + "x)")
            elif int(coverage) > 100.00:
                alerts.append("coverage >100x(" + str(coverage) + "x)")
    else:
        pass
    if str(assembly_stdev) == "NA":
        if str(gc_stdev) == "NA":
            alerts.append("Assembly ratio and GC% STDev are N/A <10 genomes as reference.")
        else:
            alerts.append("Open Github issue assembly ratio STDev is N/A but not GC%. This shouldn't happen.")
    elif str(gc_stdev) == "NA":
        alerts.append("Open Github issue STDev is N/A for GC% and not assembly ratio. This shouldn't happen.")
    else:
        pass
    alerts = ', '.join(alerts)
    return alerts

def compile_warnings(scaffolds_entry, Total_Trimmed_reads, Total_Raw_reads, Q30_R1_per, Q30_R2_per, Trim_Q30_R1_per, Trim_Q30_R2_per, scaffolds, gc_metrics, \
                    assembly_ratio_metrics, Trim_unclassified_percent, Wt_asmbld_unclassified_percent, kraken_trim_genus, kraken_wtasmbld_genus, Trim_Genus_percent, Asmbld_Genus_percent,\
                    MLST_scheme_1, MLST_scheme_2, scheme_guess, genus, fastani_warning, busco_id, FastANI_ID, FastANI_coverage, srst2_warning, QC_reason):
    """
    <1,000,000 total reads for each raw and trimmed reads - Total_Sequenced_reads
    % raw and trimmed reads with Q30 average for R1 (<90%) and R2 (<70%) - Q30_R1_percent, Q30_R2_percent
    >200 scaffolds - scaffolds
    Checking that %GC content isn't >2.58 stdev away from the mean %GC content for the species determined - assembly_ratio_metrics
    Contamination check: >30% unclassified reads and confirm there is only 1 genera with >25% of assigned reads - Trim_kraken, wt_Asmbld_kraken
    """
    # check warnings
    warnings = []
    if scaffolds_entry == False:
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
                warnings.append("Average Q30 of trimmed R1 reads <{:.2f}% ({:.2f}%)".format(float(90.00),float(Trim_Q30_R2_per)))
            except ValueError:
                warnings.append("Average Q30 of trimmed R1 reads <{:.2f}% ({})".format(float(90.00),Trim_Q30_R2_per))
        if Trim_unclassified_percent == "Unknown" or float(Trim_unclassified_percent) > float(30.00):
            warnings.append(">{:.2f}% unclassifed trimmed reads.".format(int(30)))
        if len(kraken_trim_genus) >=2:
            warnings.append(">=2 genera had >{:.2f}% of reads assigned to them.".format(int(25)))
        if Trim_Genus_percent == "Unknown" or float(Trim_Genus_percent) <float(70.00):
            try:
                warnings.append("<70% of reads assigned to top genera hit ({:.2f}%)".format(float(Trim_Genus_percent)))
            except ValueError:
                warnings.append("<70% of reads assigned to top genera hit ({})".format(Trim_Genus_percent))
    else:
        pass
    if gc_metrics[0] != "NA" and gc_metrics[0] != "Unknown":
        # sample_gc > (species_gc_mean + out_of_range_stdev)
        if float(gc_metrics[1]) > (float(gc_metrics[3])+float(gc_metrics[2])): #check that gc% is < 2.58 stdev away from mean gc of species
            warnings.append("GC% >2.58 stdev away from mean GC of {:.2f}%".format(float(gc_metrics[3])))
    if scaffolds != "Unknown" and Wt_asmbld_unclassified_percent != "Unknown" and Asmbld_Genus_percent != "Unknown":
        if int(scaffolds) > int(200) and int(scaffolds) < int(500): # between 200-500 
            warnings.append("High scaffold count 200-500 ({})".format(int(scaffolds)))
        if float(Wt_asmbld_unclassified_percent) > float(30.00):
            warnings.append(">{:.2f}% unclassifed weighted scaffolds".format(int(30)))
        if float(Asmbld_Genus_percent) <float(70.00):
            warnings.append("<70% of weighted scaffolds assigned to top genera hit ({:.2f}%)".format(float(Asmbld_Genus_percent)))
    elif scaffolds == "Unknown" and Wt_asmbld_unclassified_percent == "Unknown" and Asmbld_Genus_percent == "Unknown":
        if "No assembly due to:" in QC_reason: # if there is already a QC reason for no assembly then don't add this warning
            pass
        else:
            warnings.append("No assembly file found possible SPAdes failure.")
    if len(kraken_wtasmbld_genus) >=2:
        warnings.append(">=2 genera had >{:.2f}% of wt scaffolds assigned to them".format(int(25))) 
    if MLST_scheme_1 != "-" and not MLST_scheme_1.startswith(scheme_guess):
        if genus == "Enterobacter" and MLST_scheme_1 == "ecloacae":
            pass
        elif genus == "Klebsiella" and MLST_scheme_1 == "koxytoca" or "kpneumoniae":
            pass
        elif genus == "Acinetobacter" and MLST_scheme_1.startswith("abaumannii"):
            pass
        elif genus == "Citrobacter" and MLST_scheme_1 == "cfreundii":
            pass
        elif genus == "Burkholderia" and MLST_scheme_1 == "bcepacia":
            pass
        elif (genus == "Shigella" or genus == "Escherichia") and MLST_scheme_1.startswith("ecoli"):
            pass
        else:
            warnings.append("Check 1st MLST scheme matches taxa IDed")
    if MLST_scheme_2 != "-" and not MLST_scheme_2.startswith(scheme_guess):
        if genus == "Enterobacter" and MLST_scheme_2 == "ecloacae":
            pass
        elif genus == "Klebsiella" and MLST_scheme_1 == "koxytoca" or "kpneumoniae":
            pass
        elif genus == "Acinetobacter" and MLST_scheme_2.startswith("abaumannii"):
            pass
        elif genus == "Citrobacter" and MLST_scheme_2 == "cfreundii":
            pass
        elif genus == "Burkholderia" and MLST_scheme_2 == "bcepacia":
            pass
        elif (genus == "Shigella" or genus == "Escherichia") and MLST_scheme_2.startswith("ecoli"):
            pass
        else:
            warnings.append("Check 2nd MLST scheme matches taxa IDed")
    if FastANI_ID != "Unknown":
        if float(FastANI_ID) < float(95.00):
            warnings.append("FastANI match is <95%")
        if float(FastANI_coverage) < float(90.00):
            warnings.append("FastANI coverage is <90%")
    if busco_id != "Unknown":
        if float(busco_id) < float(97.00):
            warnings.append("BUSCO match is <97%")
    #add in fastani warning
    if fastani_warning != None:
        warnings.append(fastani_warning)
    # add srst2 warning
    if srst2_warning != None:
        warnings.append(srst2_warning)
    warnings = ', '.join(warnings).strip(", ")
    # For spades failures, lack of reads after trimming or corruption we will simplify the warnings by supressing other warnings
    if "No assembly file found possible SPAdes failure." in warnings:
        warnings = "No assembly file found possible SPAdes failure."
    elif "is corrupt and is unable to be unzipped" in warnings:
        warnings = [item for item in warnings if "corrupt" in item]
    elif  "The # of reads in raw R1/R2 files are NOT equal." in warnings:
        warnings = [item for item in warnings if "NOT equal" in item]
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

def Checking_auto_pass_fail(fairy_files, scaffolds_entry, coverage, length, assembly_stdev, asmbld_ratio, set_coverage, scaffolds, sample_name):
    """Checking auto pass fail conditions"""
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
                    QC_reason.append(str(fastq_file_failure) + " is corrupt and is unable to be unzipped")
                if ('FAILED: The number of reads in R1/R2 are NOT the same!' in line):
                    QC_result.append("FAIL")
                    QC_reason.append("The # of reads in raw R1/R2 files are NOT equal")
                if ('FAILED: There are 0 reads in' in line):
                    QC_result.append("FAIL")
                    QC_reason.append("No reads remain after trimming")
                if ('FAILED: No scaffolds in ' in line):
                    QC_result.append("FAIL")
                    QC_reason.append("No scaffolds were >500bp")
        f.close()
    if scaffolds_entry == False: # if its not being used for scaffolds entry check estimated coverage otherwise don't
        if coverage == "Unknown" or int(coverage) < int(set_coverage):
            QC_result.append("FAIL")
            if coverage == "Unknown": # if else really only needed so you don't end up with "unknownx"
                QC_reason.append("coverage <"+ str(set_coverage) +"x (" + str(coverage) + ")")
            else:
                QC_reason.append("coverage <"+ str(set_coverage) +"x (" + str(coverage) + "x)")
    else:
        pass
    if length == "Unknown" or int(length) <= 1000000:
        QC_result.append("FAIL")
        QC_reason.append("assembly <1,000,000bps (" + str(length) + ")")
    if str(assembly_stdev) != "NA": # have to have a second layer cuz you can't make NA a float, N/A means less than 10 genomes so no stdev calculated
        if str(asmbld_ratio) == "Unknown": # if there is no ratio file then fail the sample
            QC_result.append("FAIL")
            QC_reason.append("Assembly file not found")
        elif float(assembly_stdev) > 2.58:
            QC_result.append("FAIL")
            QC_reason.append("assembly stdev >2.58 (" + str(assembly_stdev) + ")")
    if str(scaffolds) == "Unknown" or int(scaffolds) > int(500):
        QC_result.append("FAIL")
        QC_reason.append("High scaffold count >500 ({})".format(str(scaffolds)))
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
        elif "No reads" in check_QC_reason or "No scaffolds" in check_QC_reason:
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
            df.drop(dup, axis=1, inplace=True)
            #add in new frame
            df[dup] = new_col
    return df

def parse_gamma_ar(gamma_ar_file, sample_name, final_df):
    """Parsing the gamma file run on the antibiotic resistance database."""
    gamma_df = pd.read_csv(gamma_ar_file, sep='\t', header=0)
    DB = (gamma_ar_file.rsplit('/', 1)[-1]).replace(sample_name, "").rsplit('_')[1] + "_" + (gamma_ar_file.rsplit('/', 1)[-1]).replace(sample_name, "").rsplit('_')[2] + "([XNT/98AA/90]G:[98NT/90]S)"
    percent_BP_IDs = np.floor(gamma_df["BP_Percent"]*100).tolist() # round % to whole number
    percent_codon_IDs = np.floor(gamma_df["Codon_Percent"]*100).tolist() # round % to whole number
    percent_lengths = np.floor(gamma_df["Percent_Length"]*100).tolist() # round % to whole number
    conferred_resistances = gamma_df["Gene"].str.split("__").str[4] #parse "Gene" column in gamma file to get conferred resistance out of gene name
    contig_numbers = gamma_df["Contig"].str.replace(sample_name, "").str.split("_").str[1] #Parse "Contig" column in gamma file
    genes = gamma_df["Gene"].str.split("__").str[2] #Parse "Gene" column in gamma file to get gene name and accession
    # loop through list of genes to combine with conferred resistance and make back into a pandas series
    column_name = ["{}_({})".format(gene, conferred_resistance) for gene, conferred_resistance in zip(genes, conferred_resistances)]
    # loop through list of gamma info to combine into "code" for ID%/%cov:contig# and make back into a pandas series
    coverage = ["[{:.0f}NT/{:.0f}AA/{:.0f}:#{}]G".format(percent_BP_ID, percent_codon_ID, percent_length, contig_number) for percent_BP_ID, percent_codon_ID, percent_length, contig_number in zip(percent_BP_IDs, percent_codon_IDs, percent_lengths, contig_numbers)]
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
        df = pd.DataFrame({'WGS_ID':[sample_name], 'AR_Database':[DB], 'No_AR_Genes_Found':['[-/-]'] })
        df.index = [sample_name]
    else:
        df.columns = column_name # add column names
        df["WGS_ID"] = sample_name
        df["AR_Database"] = DB
        df["No_AR_Genes_Found"] = ""
        df.index = [sample_name]
    # Check for duplicate column names, multiple hits 
    df = duplicate_column_clean(df)
    final_df = pd.concat([final_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    return final_df

def parse_gamma_hv(gamma_hv_file, sample_name, final_df):
    """Parsing the gamma file run on the antibiotic resistance database."""
    gamma_df = pd.read_csv(gamma_hv_file, sep='\t', header=0)
    DB = (gamma_hv_file.rsplit('/', 1)[-1]).replace(sample_name, "").rsplit('_')[1] + "_" + (gamma_hv_file.rsplit('/', 1)[-1]).replace(sample_name, "").rsplit('_')[2].strip(".gamma")
    percent_BP_IDs = np.floor(gamma_df["BP_Percent"]*100).tolist() # round % to whole number
    percent_codon_IDs = np.floor(gamma_df["Codon_Percent"]*100).tolist() # round % to whole number
    percent_lengths = np.floor(gamma_df["Percent_Length"]*100).tolist() # round % to whole number
    conferred_resistances = gamma_df["Gene"].str.split("__").str[4] #parse "Gene" column in gamma file to get conferred resistance out of gene name
    contig_numbers = gamma_df["Contig"].str.replace(sample_name, "").str.split("_").str[1] #Parse "Contig" column in gamma file
    hv_column_name = gamma_df["Gene"] #Parse "Gene" column in gamma file to get gene name and accession
    # loop through list of gamma info to combine into "code" for ID%/%cov:contig# and make back into a pandas series
    coverage = ["[{:.0f}NT/{:.0f}AA/{:.0f}:#{}]G".format(percent_BP_ID, percent_codon_ID, percent_length, contig_number) for percent_BP_ID, percent_codon_ID, percent_length, contig_number in zip(percent_BP_IDs, percent_codon_IDs, percent_lengths, contig_numbers)]
    #building a new dataframe - create giant row
    df = pd.DataFrame(coverage).T
    if df.empty:
        df = pd.DataFrame({'WGS_ID':[sample_name], 'HV_Database':[DB], 'No_HVGs_Found':['[-/-]'] })
        df.index = [sample_name]
    else:
        df.columns = hv_column_name # add column names
        df["WGS_ID"] = sample_name
        df["HV_Database"] = DB
        df["No_HVGs_Found"] = ""
        df.index = [sample_name]
    # Check for duplicate column names
    df = duplicate_column_clean(df)
    final_df = pd.concat([final_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    return final_df

def parse_gamma_pf(gamma_pf_file, sample_name, pf_df):
    """Parsing the gamma file run on the plasmid marker database."""
    gamma_df = pd.read_csv(gamma_pf_file, sep='\t', header=0)
    DB = (gamma_pf_file.rsplit('/', 1)[-1]).replace(sample_name, "").rsplit('_')[1] + "_" + (gamma_pf_file.rsplit('/', 1)[-1]).replace(sample_name, "").rsplit('_')[2].strip(".gamma") + "([95NT/60]) "
    if DB == "":
        DB ="Unknown"
    percent_NT_IDs = np.floor(gamma_df["Match_Percent"]*100).tolist() # round % to whole number
    percent_lengths = np.floor(gamma_df["Length_Percent"]*100).tolist() # round % to whole number - this is the coverage
    contig_numbers = gamma_df["Contig"].str.replace(sample_name, "").str.split("_").str[1] #Parse "Contig" column in gamma file
    pf_column_name = gamma_df["Gene"] #Parse "Gene" column in gamma file to get gene name and accession
    # loop through list of gamma info to combine into "code" for ID%/%cov:contig# and make back into a pandas series
    pf_coverage = ["[{:.0f}NT/{:.0f}:#{}]G".format(percent_NT_ID, percent_length, contig_number) for percent_NT_ID, percent_length, contig_number in zip(percent_NT_IDs, percent_lengths, contig_numbers)]
    # Minimum % length required to be included in report, otherwise removed from list
    if bool([percent_length for percent_length in percent_lengths if int(percent_length) < 60]): 
        index_remove_postion = [ n for n,percent_length in enumerate(percent_lengths) if int(percent_length) < 60 ] # get index for value removed to remove from other lists (values less than 90)
        percent_lengths = [percent_length for percent_length in percent_lengths if int(percent_length) >= 60] # filter list to remove values below cutoff (keep those greater than or equal to 90)
        for index in sorted(index_remove_postion, reverse=True): #delete them in reverse order so that you don't throw off the subsequent indexes.
            del pf_coverage[index]
            del pf_column_name[index]
            del percent_NT_IDs[index]
    # Minimum % identity required to be included in report, otherwise removed from list
    if bool([percent_NT_ID for percent_NT_ID in percent_NT_IDs if int(percent_NT_ID) < 95]):
        index_remove_postion = [ n for n,percent_NT_ID in enumerate(percent_NT_IDs) if int(percent_NT_ID) < 95 ] # get index for value removed to remove from other lists (values less than 98)
        percent_NT_IDs = [percent_NT_ID for percent_NT_ID in percent_NT_IDs if int(percent_NT_ID) >= 95] # filter list to remove values below cutoff (keep those greater than or equal to 98)
        #reset indexes as if you deleted rows in the if loop above the index values will be off from index_remove_postion 
        pf_column_name.reset_index(drop=True, inplace=True)
        #loop through list of indexes to delete and remove them from the other lists so they all match
        for index in sorted(index_remove_postion, reverse=True): #delete them in reverse order so that you don't throw off the subsequent indexes.
            del pf_coverage[index]
            del pf_column_name[index]
            del percent_lengths[index]
    #building a new dataframe - create giant row
    df = pd.DataFrame(pf_coverage).T
    if df.empty:
        df = pd.DataFrame({'WGS_ID':[sample_name], 'Plasmid_Replicon_Database':[DB], 'No_Plasmid_Markers':['[-/-]'] })
        df.index = [sample_name]
    else:
        df.columns = pf_column_name # add column 'HV_Database':[DB], names
        df["WGS_ID"] = sample_name
        df["Plasmid_Replicon_Database"] = DB
        df['No_Plasmid_Markers'] = ""
        df.index = [sample_name]
    # Check for duplicate column names
    df = duplicate_column_clean(df)
    pf_df = pd.concat([pf_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    return pf_df

def parse_mlst(mlst_file, scheme_guess, sample_name):
    """Pulls MLST info from *_combined.tsv file."""
    Scheme_list = [[],[],[],[],[]] # create empty list to fill later this would be the MLST_Scheme_1	MLST_1	MLST_Scheme_2	MLST_2
    with open(mlst_file, 'r') as f:
        lines = f.readlines()
        lines.pop(0) # drop header line
        for rawline in lines:
            line=rawline.strip()
            split_line = line.split("\t")
            source = split_line[1]
            date = split_line[2]
            DB_ID = split_line[3] # scheme name (i.e Pasteur or Oxford etc)
            Scheme = str(split_line[4]) # scheme number
            if scheme_guess == "abaum" and "PARALOG" in Scheme: # supress Paralogs for Acinetobacter_baumannii
                print("Warning: surpressing " + Scheme + " in " + sample_name)
            else:
                #handle cases where the alleles are all -
                if (len(set(split_line[5:]))==1) and split_line[5:][0] == "-":
                    alleles = "-"
                else:
                    alleles = ".".join(split_line[5:]) # combine all alleles separated by -
                #exclusion list for later
                exclusion_list = [ "-", "Novel_allele", "Novel_profile", "Missing_allele", "Novel_allele-PARALOG", "Novel_profile-PARALOG", "Missing_allele-PARALOG" ]
                if DB_ID in Scheme_list[0]: # check if scheme name is already in the scheme list
                    for i in range(0,len(Scheme_list[0])): #loop through list of scheme names
                        if DB_ID == Scheme_list[0][i]: # looking for matching scheme name that was already in the list 
                            # If the scheme was already in the list then add the ST, alleles, source and data into that list within for that scheme
                            # Example: [['abaumannii(Pasteur)', 'abaumannii(Oxford)'], [['ST2'], ['ST195', 'ST1816-PARALOG']], [['cpn60(2)-fusA(2)-gltA(2)-pyrG(2)-recA(2)-rplB(2)-rpoB(2)'], ['gltA(1)-gyrB(3)-gdhB(3)-recA(2)-cpn60(2)-gpi(96)-rpoD(3)', 'gltA(1)-gyrB(3)-gdhB(189)-recA(2)-cpn60(2)-gpi(96)-rpoD(3)']], [['standard/srst2'], ['standard', 'standard/srst2']], [['2022-12-02'], ['2022-12-02', '2022-12-02']]]
                            if not any(x in Scheme for x in exclusion_list):
                                Scheme_list[1][i].append("ST"+str(Scheme)) # if there is a value add ST in front of number
                            else:
                                Scheme_list[1][i].append(Scheme) # just append what is there
                            Scheme_list[2][i].append(alleles)
                            Scheme_list[3][i].append(source)
                            Scheme_list[4][i].append(date)
                else: # if scheme name is not already in the scheme list add it
                    Scheme_list[0].append(DB_ID)
                    if not any(x in Scheme for x in exclusion_list):
                        Scheme_list[1].append(["ST"+str(Scheme)]) # if there is a value add ST in front of number
                    else:
                        Scheme_list[1].append([Scheme]) # just append what is there
                    Scheme_list[2].append([alleles])
                    Scheme_list[3].append([source])
                    Scheme_list[4].append([date])
    return Scheme_list

def parse_ani(fast_ani_file):
    """Parse ANI file to get format 99.98%ID-98.58%COV-Acinetobacter baumannii(Acinetobacter_baumannii_GCF_012935145.1_ASM1293514v1_genomic.fna.gz)."""
    with open(fast_ani_file) as f:
        first_line = f.readline().strip('\n')
        #try: #try to see if there is a second line.
        #    second_line = f.readlines()[0].strip('\n')
        #except IndexError:
        #    second_line = ""
    if "No MASH hit found" in first_line:
        FastANI_output_list = ['Unknown','Unknown','Unknown','Unknown']
        scheme_guess = "NA NA"
        fastani_warning = "No MASH hit found."
    elif "No hits above an ANI value >=80%" in first_line:
        FastANI_output_list = ['Unknown','Unknown','Unknown','Unknown']
        scheme_guess = "NA NA"
        fastani_warning = "No hits >=80% ANI see *.ani.txt for details."
    else:
        fastani_warning = ""
        ani_df = pd.read_csv(fast_ani_file, sep='\t', header=0) # should only be one line long.
        ID = ani_df["% ID"][0]
        coverage = ani_df["% Coverage"][0]
        organism = ani_df["Organism"][0].replace("-chromosome", "")
        organism = re.sub(r"-strain.*", "", organism) # remove everything after strain
        #if sp. is followed by anything other than a space add it.
        if "sp." in organism and re.search(r'sp\.\S', organism):
            organism = organism.replace("sp.","sp. ")
        source_file = ani_df["Source File"][0]
        #Species_Support = str(ID) + "%ID-" + str(coverage) + "%COV-" + organism + "(" + source_file + ")" #old way of reporting
        FastANI_output_list = [source_file, ID, coverage, organism]
        # get taxa to check mlst scheme
        scheme_guess = organism.split(' ')[0][0].lower() + organism.split(' ')[1][0:4]
    return FastANI_output_list, scheme_guess, fastani_warning

def parse_srst2_ar(srst2_file, ar_dic, final_srst2_df, sample_name):
#def parse_srst2_ar(srst2_file, final_srst2_df, sample_name):
    """Parsing the srst2 file run on the ar gene database."""
    srst2_df = pd.read_csv(srst2_file, sep='\t', header=0)
    percent_lengths = np.floor(srst2_df["coverage"]).tolist()
    genes = srst2_df["allele"].tolist()
    percent_BP_IDs = np.floor(100 - srst2_df["divergence"]).tolist()
    #ar_dic={}
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #    print(srst2_df)
    #for i in range(0, len(srst2_df.index)):
    #    print(i, srst2_df['annotation'].values[i])
    #    temp_allele=srst2_df['annotation'].values[i].split(':')[0].split(']')[1]
    #    temp_accession=srst2_df['annotation'].values[i].split(':')[2]
    #    temp_drug=srst2_df['annotation'].values[i].split(';')[-2]
    #    temp_id=temp_allele+'_'+temp_accession
    #    ar_dic[temp_id]=temp_drug
        
    # Since srst2 currently doesn't handle () in the gene names we will make a quick detour to fix this... now fixing annotations
    srst2_df.annotation = srst2_df.annotation.fillna(srst2_df.allele.map(ar_dic)) # this only fills in nas


    #### Add error to catch removed genes or possible database mismatching
    srst2_df['conferred_resistances'] = srst2_df['allele'].map(ar_dic)
    conferred_resistances = srst2_df['conferred_resistances'].tolist()
    ###!print(conferred_resistances)
    # loop through list of genes to combine with conferred resistance and make back into a pandas series
    column_name = ["{}_({})".format(gene, conferred_resistance) for gene, conferred_resistance in zip(genes, conferred_resistances)]
    # loop through list of srst2 info to combine into "code" for ID%/%cov:contig# and make back into a pandas series
    coverage = ["[{:.0f}NT/{:.0f}]S".format(percent_BP_ID, percent_length) for percent_BP_ID, percent_length in zip(percent_BP_IDs, percent_lengths)]
    # Minimum % length required to be included in report, otherwise removed from list
    if bool([percent_length for percent_length in percent_lengths if int(percent_length) < 90]): 
        index_remove_postion = [ n for n,percent_length in enumerate(percent_lengths) if int(percent_length) < 90 ] # get index for value removed to remove from other lists (values less than 90)
        percent_lengths = [percent_length for percent_length in percent_lengths if int(percent_length) >= 90] # filter list to remove values below cutoff (keep those greater than or equal to 90)
        for index in sorted(index_remove_postion, reverse=True):
            del coverage[index]
            del column_name[index]
            del percent_codon_IDs[index]
    # Minimum % identity required to be included in report, otherwise removed from list
    if bool([percent_BP_ID for percent_BP_ID in percent_BP_IDs if int(percent_BP_ID) < 98]):
        index_remove_postion = [ n for n,percent_BP_ID in enumerate(percent_BP_IDs) if int(percent_BP_ID) < 98 ] # get index for value removed to remove from other lists (values less than 98)
        percent_BP_IDs = [percent_BP_ID for percent_BP_ID in percent_BP_IDs if int(percent_BP_ID) >= 98] # filter list to remove values below cutoff (keep those greater than or equal to 98)
        #loop through list of indexes to delete and remove them from the other lists so they all match
        for index in sorted(index_remove_postion, reverse=True):
            del coverage[index]
            del column_name[index]
            del percent_lengths[index]
    #building a new dataframe - create giant row
    df = pd.DataFrame(coverage).T
    if df.empty: #check if its empty - which would be when nothing is found and/or no hits passed the filter
        df = pd.DataFrame({'WGS_ID':[sample_name]})
        df.index = [sample_name]
    else:
        df.columns = column_name # add column names
        df["WGS_ID"] = sample_name
        df.index = [sample_name]
    final_srst2_df = pd.concat([final_srst2_df, df], axis=0, sort=True, ignore_index=False).fillna("")

    return final_srst2_df

def Get_Metrics(phoenix_entry, scaffolds_entry, set_coverage, srst2_ar_df, pf_df, ar_df, hv_df, trim_stats, raw_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld_report, kraken_wtasmbld, quast_report, busco_short_summary, asmbld_ratio, gc_file, sample_name, mlst_file, fairy_file, gamma_ar_file, gamma_pf_file, gamma_hv_file, fast_ani_file, tax_file, srst2_file, ar_dic):
    '''For each step to gather metrics try to find the file and if not then make all variables unknown'''
    try:
        Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Raw_reads, Total_Trimmed_bp, Paired_Trimmed_reads, Total_Trimmed_reads, Trim_Q30_R1_percent, Trim_Q30_R2_percent = get_Q30(trim_stats, raw_stats)
    except FileNotFoundError:
        print("Warning: " + sample_name + "_trimmed_read_counts.txt not found")
        Q30_R1_per = Q30_R2_per = Total_Raw_Seq_bp = Total_Raw_reads = Total_Trimmed_bp = Paired_Trimmed_reads = Total_Trimmed_reads = Trim_Q30_R1_percent = Trim_Q30_R2_percent = 'Unknown'
    # Try and except are in the get_kraken_info function to allow for cases where trimming was completed, but not assembly
    Trim_kraken, Trim_Genus_percent, Asmbld_kraken, Asmbld_Genus_percent, Trim_unclassified_percent, Wt_asmbld_unclassified_percent, scheme_guess_kraken_wt, scheme_guess_kraken_trimd = get_kraken_info(kraken_trim, kraken_wtasmbld, sample_name)
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
        busco_metrics = Get_BUSCO_Gene_Count(busco_short_summary)
    except FileNotFoundError:
        if phoenix_entry == False: # suppress warning when busco is not run
            print("Warning: short_summary.specific." + sample_name + ".filtered.scaffolds.fa.txt not found.")
        lineage = percent_busco = 'Unknown'
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
        QC_result, QC_reason = Checking_auto_pass_fail(fairy_file, scaffolds_entry, Coverage, Assembly_Length, assembly_ratio_metrics[1], assembly_ratio_metrics[0], set_coverage, Scaffold_Count, sample_name)
    except FileNotFoundError: 
        print("Warning: Possibly coverage and assembly length was not calculated and/or "+ sample_name + "_Assembly_ratio_*.txt not found.")
        QC_result = QC_reason = 'Unknown'
    try:
        FastANI_output_list, scheme_guess_fastani, fastani_warning = parse_ani(fast_ani_file)
    except FileNotFoundError: 
        print("Warning: " + sample_name + ".fastANI.txt not found")
        ani_source_file = fastani_ID = fastani_coverage = fastani_organism = 'Unknown'
        FastANI_output_list = [ani_source_file, fastani_ID, fastani_coverage, fastani_organism]
        scheme_guess_fastani = fastani_warning = ""
    try:
        ar_df = parse_gamma_ar(gamma_ar_file, sample_name, ar_df)
    except FileNotFoundError: 
        print("Warning: Gamma file for ar database on " + sample_name + " not found")
        df = pd.DataFrame({'WGS_ID':[sample_name], 'No_AR_Genes_Found':['File not found'], 'AR_Database':['GAMMA file not found'] })
        df.index = [sample_name]
        ar_df = pd.concat([ar_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    try:
        pf_df = parse_gamma_pf(gamma_pf_file, sample_name, pf_df)
    except FileNotFoundError: 
        print("Warning: Gamma file for pf database on " + sample_name + " not found")
        df = pd.DataFrame({'WGS_ID':[sample_name], 'No_Plasmid_Markers':['File not found'], 'Plasmid_Replicon_Database':['GAMMA file not found'] })
        df.index = [sample_name]
        pf_df = pd.concat([pf_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    try:
        hv_df = parse_gamma_hv(gamma_hv_file, sample_name, hv_df)
    except FileNotFoundError: 
        print("Warning: Gamma file for hv database on " + sample_name + " not found")
        df = pd.DataFrame({'WGS_ID':[sample_name], 'No_HVGs_Found':['File not found'], 'HV_Database':['GAMMA file not found'] })
        df.index = [sample_name]
        hv_df = pd.concat([hv_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    try:
        #handling for if srst2 quietly failed
        srst2_failure_checks = ["failed gene detection","No AR genes found"]
        with open(srst2_file, 'r') as file:
            file_content = file.read()
        if any(x in file_content for x in srst2_failure_checks):
            if "No AR genes found" in file_content: # no warning needed, but do need to make dataframe empty to keep things moving
                srst2_warning = None
            else:
                srst2_warning = "failed srst2 gene detection."
            df = pd.DataFrame({'WGS_ID':[sample_name]})
            df.index = [sample_name]
            srst2_ar_df = pd.concat([srst2_ar_df, df], axis=0, sort=True, ignore_index=False).fillna("")
        else:
            srst2_ar_df = parse_srst2_ar(srst2_file, ar_dic, srst2_ar_df, sample_name)
            #srst2_ar_df = parse_srst2_ar(srst2_file, srst2_ar_df, sample_name)
            srst2_warning = None
    except (FileNotFoundError, pd.errors.EmptyDataError) : # second one for an empty dataframe - srst2 module creates a blank file
        if phoenix_entry == False: #supress warning when srst2 is not run
            print("Warning: " + sample_name + "__fullgenes__ResGANNCBI__*_srst2__results.txt not found")
        df = pd.DataFrame({'WGS_ID':[sample_name]})
        df.index = [sample_name]
        srst2_ar_df = pd.concat([srst2_ar_df, df], axis=0, sort=True, ignore_index=False).fillna("")
        srst2_warning = None
    try:
        alerts = compile_alerts(scaffolds_entry, Coverage, assembly_ratio_metrics[1], gc_metrics[0])
    except:
        alerts = ""
    # try except in the function itself
    kraken_trim_genus, kraken_wtasmbld_genus = parse_kraken_report(kraken_trim_report, kraken_wtasmbld_report, sample_name)
    #check which scheme guess to use - use fastANI > kraken wt > kraken_trimd
    if scheme_guess_fastani == "":
        if scheme_guess_kraken_wt == "":
            scheme_guess = scheme_guess_kraken_wt
            genus = kraken_wtasmbld_genus
        else: #use kraken trimmed
            scheme_guess = scheme_guess_kraken_trimd
            genus = kraken_trim_genus
    else:
        scheme_guess = scheme_guess_fastani
        genus = FastANI_output_list[3].split(" ")[0]
    try:
        Scheme_list = parse_mlst(mlst_file, scheme_guess, sample_name)
        if len(Scheme_list[0]) > 1: # If there is more than one scheme
            if Scheme_list[0][0] < Scheme_list[0][1]: # this if else is all just to make sure things are printing out in the same order.
                MLST_scheme_1 = Scheme_list[0][0] # get 1st scheme name from the list
                mlst_types_1=sorted(Scheme_list[1][0])[::-1]
                # Added to maintain ordering of later lists, since sorting possibly changed order of types
                order_1 = [Scheme_list[1][0].index(item) for item in mlst_types_1]
                sorted_alleles_1=[Scheme_list[2][0][i] for i in order_1]
                sorted_sources_1=[Scheme_list[3][0][i] for i in order_1]
                MLST_type_1 = ", ".join(mlst_types_1)
                MLST_alleles_1 = ",".join(sorted_alleles_1)
                MLST_source_1 = ",".join(sorted_sources_1)
                MLST_scheme_2 = Scheme_list[0][1] # get 2nd scheme name from the list
                mlst_types_2=sorted(Scheme_list[1][1])[::-1]
                # Added to maintain ordering of later lists, since sorting possibly changed order of types
                
                order_2 = [Scheme_list[1][1].index(item) for item in mlst_types_2]
                sorted_alleles_2=[Scheme_list[2][1][i] for i in order_2]
                sorted_sources_2=[Scheme_list[3][1][i] for i in order_2]
                MLST_type_2 = ", ".join(mlst_types_2)
                MLST_alleles_2 = ",".join(sorted_alleles_2)
                MLST_source_2 = ",".join(sorted_sources_2)
            else:
                MLST_scheme_1 = Scheme_list[0][1] # get 1st scheme name from the list, in this case its the 2nd element
                mlst_types_1=sorted(Scheme_list[1][1])[::-1]
                # Added to maintain ordering of later lists, since sorting possibly changed order of types
                order_1 = [Scheme_list[1][1].index(item) for item in mlst_types_1]
                sorted_alleles_1=[Scheme_list[2][1][i] for i in order_1]
                sorted_sources_1=[Scheme_list[3][1][i] for i in order_1]
                MLST_type_1 = ", ".join(mlst_types_1)
                MLST_alleles_1 = ",".join(sorted_alleles_1)
                MLST_source_1 = ",".join(sorted_sources_1)
                MLST_scheme_2 = Scheme_list[0][0] # get 2nd scheme name from the list, in this case its the first element
                mlst_types_2=sorted(Scheme_list[1][0])[::-1]
                # Added to maintain ordering of later lists, since sorting possibly changed order of types
                order_2 = [Scheme_list[1][0].index(item) for item in mlst_types_2]
                sorted_alleles_2=[Scheme_list[2][0][i] for i in order_2]
                sorted_sources_2=[Scheme_list[3][0][i] for i in order_2]
                MLST_type_2 = ", ".join(mlst_types_2)
                MLST_alleles_2 = ",".join(sorted_alleles_2)
                MLST_source_2 = ",".join(sorted_sources_2)
        else: # If there is only one scheme then the last scheme and type are just "-"
            MLST_scheme_1 = Scheme_list[0][0]
            MLST_type_1 = ", ".join(Scheme_list[1][0]) # join together the STs for this one scheme
            MLST_alleles_1 = ",".join(Scheme_list[2][0])
            MLST_source_1 = ",".join(Scheme_list[3][0])
            MLST_scheme_2 = ""
            MLST_type_2 = ""
            MLST_alleles_2 = ""
            MLST_source_2 = ""
    except FileNotFoundError: 
        print("Warning: " + sample_name + "_combined.tsv not found")
        MLST_scheme_1 = MLST_scheme_2 = MLST_type_1 = MLST_type_2 = MLST_alleles_1 = MLST_alleles_2 = MLST_source_1 = MLST_source_2 = 'Unknown'
    try:
        warnings = compile_warnings(scaffolds_entry, Total_Trimmed_reads, Total_Raw_reads, Q30_R1_per, Q30_R2_per, Trim_Q30_R1_percent, Trim_Q30_R2_percent,\
                                    Scaffold_Count, gc_metrics, assembly_ratio_metrics, Trim_unclassified_percent, Wt_asmbld_unclassified_percent,\
                                    kraken_trim_genus, kraken_wtasmbld_genus, Trim_Genus_percent, Asmbld_Genus_percent, MLST_scheme_1, MLST_scheme_2, scheme_guess,\
                                    genus, fastani_warning, busco_metrics[1], FastANI_output_list[1], FastANI_output_list[2], srst2_warning, QC_reason)
    except:
        warnings = ""
    return srst2_ar_df, pf_df, ar_df, hv_df, Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Raw_reads, Paired_Trimmed_reads, Total_Trimmed_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, \
    Scaffold_Count, busco_metrics, gc_metrics, assembly_ratio_metrics, QC_result, QC_reason, MLST_scheme_1, MLST_scheme_2, MLST_type_1, MLST_type_2, MLST_alleles_1, MLST_alleles_2, MLST_source_1, MLST_source_2
    

def try_paths(path1, path2):
    """Function to try the first path, then fall back to second path if file doesn't exist"""
    if os.path.exists(path1):
        return path1
    elif os.path.exists(path2):
        return path2
    else:
        return path1  # Return the first path even if it doesn't exist, for consistent error handling

def Get_Files(directory1, sample_name, directory2):
    print(sample_name)
    '''Create file paths to collect files from sample folder.'''
    # if there is a trailing / remove it
    directory1 = directory1.rstrip('/')
    directory2 = directory2.rstrip('/')
    # create file names
    trim_stats = try_paths( directory1 + "/qc_stats/" + sample_name + "_trimmed_read_counts.txt", directory2 + "/" + sample_name + "/qc_stats/" + sample_name + "_trimmed_read_counts.txt" )
    raw_stats = try_paths( directory1 + "/raw_stats/" + sample_name + "_raw_read_counts.txt", directory2 + "/" + sample_name + "/raw_stats/" + sample_name + "_raw_read_counts.txt" )
    kraken_trim = try_paths( directory1 + "/kraken2_trimd/" + sample_name + ".kraken2_trimd.top_kraken_hit.txt", directory2 + "/" + sample_name + "/kraken2_trimd/" + sample_name + ".kraken2_trimd.top_kraken_hit.txt" )
    kraken_trim_report = try_paths( directory1 + "/kraken2_trimd/" + sample_name + ".kraken2_trimd.summary.txt", directory2 + "/" + sample_name + "/kraken2_trimd/" + sample_name + ".kraken2_trimd.summary.txt" )
    kraken_wtasmbld = try_paths( directory1 + "/kraken2_asmbld_weighted/" + sample_name + ".kraken2_wtasmbld.top_kraken_hit.txt", directory2 + "/" + sample_name + "/kraken2_asmbld_weighted/" + sample_name + ".kraken2_wtasmbld.top_kraken_hit.txt" )
    kraken_wtasmbld_report = try_paths( directory1 + "/kraken2_asmbld_weighted/" + sample_name + ".kraken2_wtasmbld.summary.txt", directory2 + "/" + sample_name + "/kraken2_asmbld_weighted/" + sample_name + ".kraken2_wtasmbld.summary.txt" )
    quast_report = try_paths( directory1 + "/quast/" + sample_name + "_summary.tsv", directory2 + "/" + sample_name + "/quast/" + sample_name + "_summary.tsv" )
    mlst_file = try_paths( directory1 + "/mlst/" + sample_name + "_combined.tsv", directory2 + "/" + sample_name + "/mlst/" + sample_name + "_combined.tsv" )
    # For glob patterns, try both directories
    fairy_file_1 = glob.glob(directory1 + "/file_integrity/" + sample_name + "_*_summary.txt")
    fairy_file_2 = glob.glob(directory2 + "/" + sample_name + "/file_integrity/" + sample_name + "_*_summary.txt")
    fairy_file = fairy_file_1 if fairy_file_1 else fairy_file_2
    # For the remaining glob patterns, handle with try-except but attempt both directories
    try:
        busco_short_summary_1 = glob.glob(directory1 + "/BUSCO/short_summary.specific.*" + sample_name + ".filtered.scaffolds.fa.txt")
        if busco_short_summary_1:
            busco_short_summary = busco_short_summary_1[0]
        else:
            busco_short_summary_2 = glob.glob(directory2 + "/" + sample_name + "/BUSCO/short_summary.specific.*" + sample_name + ".filtered.scaffolds.fa.txt")
            if busco_short_summary_2:
                busco_short_summary = busco_short_summary_2[0]
            else:
                busco_short_summary = directory1 + "/BUSCO/short_summary.specific.blank" + sample_name + ".filtered.scaffolds.fa.txt"
    except IndexError:
        busco_short_summary = directory1 + "/BUSCO/short_summary.specific.blank" + sample_name + ".filtered.scaffolds.fa.txt"
    try:
        asmbld_ratio_1 = glob.glob(directory1 + "/" + sample_name + "_Assembly_ratio_*.txt")
        if asmbld_ratio_1:
            asmbld_ratio = asmbld_ratio_1[0]
        else:
            asmbld_ratio_2 = glob.glob(directory2 + "/" + sample_name + "/" + sample_name + "_Assembly_ratio_*.txt")
            if asmbld_ratio_2:
                asmbld_ratio = asmbld_ratio_2[0]
            else:
                asmbld_ratio = directory1 + "/" + sample_name + "_Assembly_ratio_blank.txt"
    except IndexError:
        asmbld_ratio = directory1 + "/" + sample_name + "_Assembly_ratio_blank.txt"
    try:
        gc_1 = glob.glob(directory1 + "/" + sample_name + "_GC_content_*.txt")
        if gc_1:
            gc = gc_1[0]
        else:
            gc_2 = glob.glob(directory2 + "/" + sample_name + "/" + sample_name + "_GC_content_*.txt")
            if gc_2:
                gc = gc_2[0]
            else:
                gc = directory1 + "/" + sample_name + "_GC_content_blank.txt"
    except IndexError:
        gc = directory1 + "/" + sample_name + "_GC_content_blank.txt"
    # Continue with similar pattern for remaining glob patterns
    try:
        gamma_ar_file_1 = glob.glob(directory1 + "/gamma_ar/" + sample_name + "_*.gamma")
        if gamma_ar_file_1:
            gamma_ar_file = gamma_ar_file_1[0]
        else:
            gamma_ar_file_2 = glob.glob(directory2 + "/" + sample_name + "/gamma_ar/" + sample_name + "_*.gamma")
            if gamma_ar_file_2:
                gamma_ar_file = gamma_ar_file_2[0]
            else:
                gamma_ar_file = directory1 + "/gamma_ar/" + sample_name + "_blank.gamma"
    except IndexError:
        gamma_ar_file = directory1 + "/gamma_ar/" + sample_name + "_blank.gamma"
    # Apply the same pattern for the remaining files
    try:
        gamma_pf_file_1 = glob.glob(directory1 + "/gamma_pf/" + sample_name + "_*.gamma")
        if gamma_pf_file_1:
            gamma_pf_file = gamma_pf_file_1[0]
        else:
            gamma_pf_file_2 = glob.glob(directory2 + "/" + sample_name + "/gamma_pf/" + sample_name + "_*.gamma")
            if gamma_pf_file_2:
                gamma_pf_file = gamma_pf_file_2[0]
            else:
                gamma_pf_file = directory1 + "/gamma_pf/" + sample_name + "_blank.gamma"
    except IndexError:
        gamma_pf_file = directory1 + "/gamma_pf/" + sample_name + "_blank.gamma"
    try: 
        gamma_hv_file_1 = glob.glob(directory1 + "/gamma_hv/" + sample_name + "_*.gamma")
        if gamma_hv_file_1:
            gamma_hv_file = gamma_hv_file_1[0]
        else:
            gamma_hv_file_2 = glob.glob(directory2 + "/" + sample_name + "/gamma_hv/" + sample_name + "_*.gamma")
            if gamma_hv_file_2:
                gamma_hv_file = gamma_hv_file_2[0]
            else:
                gamma_hv_file = directory1 + "/gamma_hv/" + sample_name + "_blank.gamma"
    except IndexError:
        gamma_hv_file = directory1 + "/gamma_hv/" + sample_name + "_blank.gamma"
    try:
        fast_ani_file_1 = glob.glob(directory1 + "/ANI/" + sample_name + "_REFSEQ_*.fastANI.txt")
        if fast_ani_file_1:
            fast_ani_file = fast_ani_file_1[0]
        else:
            fast_ani_file_2 = glob.glob(directory2 + "/" + sample_name + "/ANI/" + sample_name + "_REFSEQ_*.fastANI.txt")
            if fast_ani_file_2:
                fast_ani_file = fast_ani_file_2[0]
            else:
                fast_ani_file = directory1  + "/ANI/" + sample_name + ".fastANI.txt"
    except IndexError:
        fast_ani_file = directory1 + "/ANI/" + sample_name + ".fastANI.txt"
    # For regular paths, use the try_paths function
    tax_file = try_paths(  directory1 + "/" + sample_name + ".tax", directory2 + "/" + sample_name + "/" + sample_name + ".tax")
    try:
        srst2_file_1 = glob.glob(directory1 + "/srst2/" + sample_name + "__fullgenes__*_srst2__results.txt")
        if srst2_file_1:
            srst2_file = srst2_file_1[0]
        else:
            srst2_file_2 = glob.glob(directory2 + "/" + sample_name + "/srst2/" + sample_name + "__fullgenes__*_srst2__results.txt")
            if srst2_file_2:
                srst2_file = srst2_file_2[0]
            else:
                srst2_file = directory1 + "/srst2/" + sample_name + "__fullgenes__blank_srst2__results.txt"
    except IndexError:
        srst2_file = directory1 + "/srst2/" + sample_name + "__fullgenes__blank_srst2__results.txt"

    return trim_stats, raw_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld_report, kraken_wtasmbld, quast_report, mlst_file, fairy_file, busco_short_summary, asmbld_ratio, gc, gamma_ar_file, gamma_pf_file, gamma_hv_file, fast_ani_file, tax_file, srst2_file

def Append_Lists(data_location, parent_folder, sample_name, Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Seq_reads, Paired_Trimmed_reads, Total_trim_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, \
            Scaffold_Count, busco_metrics, gc_metrics, assembly_ratio_metrics, QC_result, QC_reason, MLST_scheme_1, MLST_scheme_2, MLST_type_1, MLST_type_2, MLST_alleles_1, MLST_alleles_2, MLST_source_1, MLST_source_2, \
            data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L,Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L):
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
        MLST_scheme_1_L.append(MLST_scheme_1)
        MLST_scheme_2_L.append(MLST_scheme_2)
        MLST_type_1_L.append(MLST_type_1),
        MLST_type_2_L.append(MLST_type_2)
        MLST_alleles_1_L.append(MLST_alleles_1)
        MLST_alleles_2_L.append(MLST_alleles_2)
        MLST_source_1_L.append(MLST_source_1)
        MLST_source_2_L.append(MLST_source_2)
        return data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
        Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L

def Create_df(phx_version, phoenix, data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L,
Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L):
    phx_version_L = [str(phx_version)] * len(Sample_Names)
    #combine all metrics into a dataframe
    if phoenix == True:
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
        'Kraken_ID_Trimmed_Reads_%'  : Trim_kraken_L,
        'Kraken_ID_WtAssembly_%'     : Asmbld_kraken_L,
        'FastANI_Organism'           : fastani_organism_L, 
        'FastANI_%ID'                : fastani_ID_L, 
        'FastANI_%Coverage'          : fastani_coverage_L,
        'Species_Support_ANI'        : Species_Support_L,
        'Primary_MLST_Scheme'        : MLST_scheme_1_L,
        'Primary_MLST_Source'        : MLST_source_1_L,
        'Primary_MLST'               : MLST_type_1_L,
        'Primary_MLST_Alleles'       : MLST_alleles_1_L,
        'Secondary_MLST_Scheme'      : MLST_scheme_2_L,
        'Secondary_MLST_Source'      : MLST_source_2_L,
        'Secondary_MLST'             : MLST_type_2_L,
        'Secondary_MLST_Alleles'     : MLST_alleles_2_L}
    else:
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
        'Species_Support_ANI'        : Species_Support_L,
        'Primary_MLST_Scheme'        : MLST_scheme_1_L,
        'Primary_MLST_Source'        : MLST_source_1_L,
        'Primary_MLST'               : MLST_type_1_L,
        'Primary_MLST_Alleles'       : MLST_alleles_1_L,
        'Secondary_MLST_Scheme'      : MLST_scheme_2_L,
        'Secondary_MLST_Source'      : MLST_source_2_L,
        'Secondary_MLST'             : MLST_type_2_L,
        'Secondary_MLST_Alleles'     : MLST_alleles_2_L}
    df = pd.DataFrame(data)
    return df

def srst2_dedup(srst2_ar_df, gamma_ar_df):
    ##### First, we will drop columns with "partial" in the name
    # Filter out columns that contain the substring
    columns_to_drop = [col for col in srst2_ar_df.columns if "partial" in col]
    # Drop the columns
    srst2_ar_df.drop(columns=columns_to_drop, inplace=True)
    ##### Second we will look for GAMMA + genes, but different alleles found my SRST2 and dedup them (AKA we will remove and not report them) #####
    # Iterate over each row in the first DataFrame -> this will give all the srst2 genes and alleles
    for idx, row in srst2_ar_df.iterrows():
        #create empty list for each row
        gene_list = []
        # For each column in the row find the srst2 hits
        for col_name in srst2_ar_df.columns:
            if pd.notna(row[col_name]) and row[col_name] != "":
                # split on _ to get just the gene/allele name for matching
                gene = col_name.split('_')[0].split("-")[0]
                gene_list.append(gene)
        # Extract unique genes in srst2 - positive SRST2 hits
        unique_gene_list = list(set(gene_list))
        # Check if the gene name (column name) exists in the GAMMA AR DataFrame -> if there is a match these would be GAMMA +, but a different allele and we want to remove these.
        for gene in unique_gene_list:
            if '(' in gene and ')' not in gene:
                ###!print("Fixing", gene, "to", gene+')')
                gene = gene + ')'
            if gamma_ar_df.columns.str.match(gene).any():
                # Check for partial string match in the column names of the gamma DataFrame
                matching_columns = gamma_ar_df.columns[gamma_ar_df.columns.str.contains(gene)].tolist()
                for column in matching_columns:
                    # check that for the column in question there is also a value found in the GAMMA column
                    if pd.notna(gamma_ar_df.at[str(idx), column]) and gamma_ar_df.at[str(idx), column] != "":
                        # Check if there is a value in the corresponding column in the second DataFrame
                        #srst2_ar_df.at[idx, srst2_ar_df.columns.str.contains(gene)] = ""
                        srst2_ar_df.loc[idx, srst2_ar_df.columns[srst2_ar_df.columns.str.contains(gene)]] = ""
                        print(f"sample {idx}: Value found in column '{gene}' of srst2_df and this matches the '{matching_columns}' of gamma_df (alleles with srst/gamma or only gamma positive) and was removed for deduplication purposes.")
    ##### Third, we will look at GAMMA- samples and check the # of gene alleles and filter to only have the top hits. ####
    # These are now the GAMMA neg hits and we will now check the number of alleles for each gene - first we do some dataframe rearranging to make it "easier"
    # Check if DataFrame is not empty
    if not srst2_ar_df.empty:
        gamma_neg_srst2 = srst2_ar_df.drop(srst2_ar_df.columns[srst2_ar_df.apply(lambda col: all(val == '' or pd.isna(val) for val in col))], axis=1)
        # check the number of alleles per gene so that we only have those that are singles and will be reported
        gamma_neg_genes = pd.DataFrame(gamma_neg_srst2.apply(lambda row: ['_'.join(column.split('_')[0:3]) for column in row.index if row[column]] or [""], axis=1))
        gamma_neg_genes.columns = ["neg_genes"] # a column name to make it easier
        # Iterate over each row
        for index, row in gamma_neg_genes.iterrows():
            gene_list = []
            multiple_occurrences = []
            count = 0
            for val in row["neg_genes"]:
                # Gonna need to fix this soon. There are genes that have dashes in the gene name that dont relate to allele
                if "(" in val:
                    if '-' in val:
                        if val.find('(') < val.find('-'):
                            gene_list.append(val.split(")")[0]+')') # In the cases where the gene name includes dashes inside of parentheses
                        else:
                            gene_list.append(val.split("-")[0]) # If the dash occurs before the parentheses...Not sure if this occurs, but attempting to catch it
                    else:
                        continue # This would break downstream anyway, but in the case that a gene does not have ANY dash in it
                else:
                    gene_list.append(val.split("-")[0]) # split to remove allele number and just have gene name
            for val in gene_list:
                count = count + 1 #only continue below if we have seen this gene before
                if gene_list.count(val) > 1 and val not in multiple_occurrences:
                    #get a dataframe of the gene in question
                    df = pd.DataFrame(gamma_neg_srst2.loc[row.name, gamma_neg_srst2.columns.str.contains(val)])
                    # Define the regex pattern to extract Percent_Match and Coverage and Extract Percent_Match and Coverage using str.extract
                    df[['Percent_Match', 'Coverage']] = df[row.name].str.extract(r'\[(\d+)NT/(\d+)\]S')
                    # Convert extracted values to integer type and keep NA values
                    df['Percent_Match'] = pd.to_numeric(df['Percent_Match'], errors='coerce')
                    df['Coverage'] = pd.to_numeric(df['Coverage'], errors='coerce')
                    df.dropna(subset=['Percent_Match'], inplace=True) #drop rows with no values
                    ### Filtering steps to get the top hit for srst2
                    # Step 1: Identify the Max Percent_Match
                    max_percent_match = df['Percent_Match'].max()

                    # Step 2: Filter rows to keep rows with the max Percent_Match
                    max_percent_match_rows = df[df['Percent_Match'] == max_percent_match]

                    # Step 3: Identify the max Coverage among rows with the max Percent_Match
                    max_coverage = max_percent_match_rows['Coverage'].max()

                    # Step 4: Drop rows with the lowest Percent_Match and, if needed, with the lowest Coverage
                    if len(max_percent_match_rows) > 1:  # Only consider Coverage if there are ties in Percent_Match
                        rows_to_drop = max_percent_match_rows[max_percent_match_rows['Coverage'] != max_coverage].index
                        # Drop the identified rows
                        df_cleaned = df.drop(rows_to_drop)
                    else:
                        max_percent_mismatch_rows = df[df['Percent_Match'] != max_percent_match]
                        df_cleaned = df.drop(max_percent_mismatch_rows.index)
                    # Now that we know what alleles we are keeping based on the top hits (these are the row names), we will get the inverse row names so we know what to drop
                    index_diff = df.index.difference(df_cleaned.index)
                    # For the sample in question remove the data from cell if the particular allele(s) we want to drop
                    srst2_ar_df.loc[index, index_diff.tolist()] = ""
                multiple_occurrences.append(val)
        srst2_ar_df.drop(srst2_ar_df.columns[srst2_ar_df.apply(lambda col: all(val == '' or pd.isna(val) for val in col))], axis=1, inplace=True)
    return srst2_ar_df

def order_ar_gene_columns(ar_combined_df, is_combine):
    ar_combined_ordered_df = pd.DataFrame() #create new dataframe to fill
    #fixing column orders
    if (is_combine):
        ar_combined_ordered_df = pd.concat([ar_combined_ordered_df, ar_combined_df[['AR_Database', 'UNI','No_AR_Genes_Found']]], axis=1, sort=False) # first adding back in ['AR_Database', 'WGS_ID']
    else:
        ar_combined_ordered_df = pd.concat([ar_combined_ordered_df, ar_combined_df[['AR_Database', 'WGS_ID','No_AR_Genes_Found']]], axis=1, sort=False) # first adding back in ['AR_Database', 'WGS_ID']
    ar_drugs_list = ar_combined_df.columns.str.extract('.*\\((.*)\\).*').values.tolist() # get all ar drug names form column names
    sorted_list = sorted(list(set([str(drug) for sublist in ar_drugs_list for drug in sublist]))) #get unique drug names (with set) and sort list
    sorted_drug_names = [x for x in sorted_list if x != 'nan'] #get unique drug names (with set) and drop nan that comes from WGS_ID column and sort
    # create list to add to
    all_column_names = []
    # loop over each gene with the same drug its name
    for drug in sorted_drug_names:
        drug = "(" + drug + ")"
        column_list = sorted([col for col in ar_combined_df.columns if drug in col]) # get column names filtered for each drug name
        ar_combined_ordered_df = pd.concat([ar_combined_ordered_df, ar_combined_df[column_list]], axis=1, sort=False) # setting column's order by combining dataframes
        all_column_names.append(column_list)
    # unnest list
    all_column_names = list(chain(*all_column_names))
    #add back AR_DB and WGS_ID to front of list
    all_column_names.insert(0, "No_AR_Genes_Found")
    if (is_combine):
        all_column_names.insert(0, "UNI")
    else:
        all_column_names.insert(0, "WGS_ID")
    all_column_names.insert(0, "AR_Database")
    # reorder columns - should be alphabetical by drug name and within drug name genes are alphabetically listed
    ar_combined_ordered_df = ar_combined_ordered_df.reindex(all_column_names, axis=1)
    return ar_combined_ordered_df

def add_srst2(ar_df, srst2_ar_df, is_combine):
    if srst2_ar_df.empty:
        return ar_df
    ar_combined_df = pd.DataFrame() #create new dataframe to fill
    common_cols = ar_df.columns.intersection(srst2_ar_df.columns) #get column names that are in both dataframes --> These are GAMMA +
    # Combine values in cells for columns that are in both dataframes as these would be the same gene alleles for GAMMA and SRST2
    for col in common_cols:
        if col != "WGS_ID" and col != 'UNI':
            ar_combined_df[col] = (srst2_ar_df[col].map(str) + ":" + ar_df[col]).replace(':', "")
            ar_combined_df[col] = ar_combined_df[col].map(lambda x: str(x).lstrip(':').rstrip(':')) # clean up : for cases where there isn't a gamma and srst2 for all rows
            ar_combined_df = ar_combined_df.copy() #defragment to correct "PerformanceWarning: DataFrame is highly fragmented."
        else:
            ar_combined_df[col] = srst2_ar_df[col]
    # check if you missed any rows, if there is a sample in ar_db, that is not in the srst2 then you will have it have NA in rows when joined
    # drop columns from srst2 dataframe that are in common in the ar_db as these are already in ar_combined_df
    srst2_ar_df.drop(common_cols, axis = 1, inplace=True) # This will leave GAMMA- samples
    ar_df.drop(common_cols, axis = 1, inplace=True)
    # Add cols that are unique to srst2
    ############### DEDUPING FOR SRST2 ###############
    srst2_ar_df = srst2_dedup(srst2_ar_df, ar_df.join(ar_combined_df))
    ar_combined_df = ar_combined_df.join(srst2_ar_df)
    # Add cols that are unique to gamma ar_df
    ar_combined_df = ar_combined_df.join(ar_df)
    #fixing column orders
    ar_combined_ordered_df = order_ar_gene_columns(ar_combined_df, is_combine)
    return ar_combined_ordered_df

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

def big5_check(final_ar_df, is_combine, BLDB):
    """"Function that will return list of columns to highlight if a sample has a hit for a big 5 gene."""
    columns_to_highlight = []
    if (is_combine):
        final_ar_df = final_ar_df.drop(['AR_Database','UNI'], axis=1)
    else:
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

def Combine_dfs(df, ar_df, pf_df, hv_df, srst2_ar_df, phoenix, is_combine, BLDB):
    hv_cols = list(hv_df)
    pf_cols = list(pf_df)
    ar_cols = list(ar_df)
    # move the column to head of list using index, pop and insert
    pf_cols.insert(0, pf_cols.pop(pf_cols.index('No_Plasmid_Markers')))
    pf_cols.insert(0, pf_cols.pop(pf_cols.index('Plasmid_Replicon_Database')))
    hv_cols.insert(0, hv_cols.pop(hv_cols.index('No_HVGs_Found')))
    hv_cols.insert(0, hv_cols.pop(hv_cols.index('HV_Database')))
    ar_cols.insert(0, ar_cols.pop(ar_cols.index('No_AR_Genes_Found')))
    ar_cols.insert(0, ar_cols.pop(ar_cols.index('AR_Database')))
    # use ix to reorder
    pf_df = pf_df.loc[:, pf_cols]
    hv_df = hv_df.loc[:, hv_cols]
    ar_df = ar_df.loc[:, ar_cols]
    # if we run -entry PHOENIX then skip
    if phoenix == True:
        final_ar_df = ar_df
    else:
        # combining srst2 and gamma ar dataframes
        final_ar_df = add_srst2(ar_df, srst2_ar_df, is_combine)
    ar_max_col = final_ar_df.shape[1] - 1 #remove one for the WGS_ID column
    # now we will check for the "big 5" genes for highlighting later.
    columns_to_highlight = big5_check(final_ar_df, is_combine, BLDB)
    # combining all dataframes
    if (is_combine):
        final_df = pd.merge(df, final_ar_df, how="left", on=["UNI","UNI"])
        final_df = pd.merge(final_df, hv_df, how="left", on=["UNI","UNI"])
        final_df = pd.merge(final_df, pf_df, how="left", on=["UNI","UNI"])
    else:
        final_df = pd.merge(df, final_ar_df, how="left", on=["WGS_ID","WGS_ID"])
        final_df = pd.merge(final_df, hv_df, how="left", on=["WGS_ID","WGS_ID"])
        final_df = pd.merge(final_df, pf_df, how="left", on=["WGS_ID","WGS_ID"])
    #get database names and remove if file is not found in the database list
    ar_db = final_df['AR_Database'].unique().tolist()
    if 'GAMMA file not found' in ar_db:
        ar_db.remove('GAMMA file not found') #Don't want this reported as the ar_db
    ar_db = ",".join(ar_db)
    hv_db = final_df['HV_Database'].unique().tolist()
    if 'GAMMA file not found' in hv_db:
        hv_db.remove("GAMMA file not found")
    hv_db = ",".join(hv_db)
    pf_db = final_df['Plasmid_Replicon_Database'].unique().tolist()
    if 'GAMMA file not found' in pf_db:
        pf_db.remove("GAMMA file not found")
    pf_db = ",".join(pf_db)
    return final_df, ar_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db

def column_letter(index):
    """Convert zero-based column index to Excel column letter."""
    letters = list(string.ascii_uppercase)
    if index < 26:
        return letters[index]
    else:
        return letters[index // 26 - 1] + letters[index % 26]  # Handle AA, AB, etc.

def write_to_excel(set_coverage, output, df, qc_max_col, ar_gene_count, pf_gene_count, hv_gene_count, columns_to_highlight, ar_df, pf_db, ar_db, hv_db, phoenix, shigapass, centar, centar_df_lens):
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
    if phoenix == True:
        worksheet.conditional_format('W3:X' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
        #worksheet.set_column('U:V', None, number_dec_2_format) # FastANI ID and Coverage
    else:
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
    ## colors for centar
    cell_format_p1 = workbook.add_format({'bg_color': '#DB7093', 'font_color': '#000000', 'bold': True})
    cell_format_p2 = workbook.add_format({'bg_color': '#FF69B4', 'font_color': '#000000', 'bold': True})
    cell_format_p3 = workbook.add_format({'bg_color': '#FFB6C1', 'font_color': '#000000', 'bold': True})
    cell_format_p4 = workbook.add_format({'bg_color': '#FFC0CB', 'font_color': '#000000', 'bold': True})
    # Headers
    #worksheet.set_column('A1:A1', None, cell_format_light_blue) #make summary column blue, #use for only 1 column in length
    if "PHX_Version" in df.columns:
        worksheet.merge_range('A1:D1', "PHoeNIx Summary", cell_format_light_blue)
        worksheet.merge_range('E1:R1', "QC Metrics", cell_format_grey_blue)
    else: # allow for backward compatibility with versions <2.2.0
        worksheet.merge_range('A1:C1', "PHoeNIx Summary", cell_format_light_blue)
        worksheet.merge_range('D1:S1', "QC Metrics", cell_format_grey_blue)
    #taxa columns 
    # Find start and end column letters
    # to allow backwards compatability with v2.1.1 we need a little try and catch...  
    if "Final_Taxa_ID" in df.columns:
        taxa_start_col  = column_letter(list(df.columns).index("Final_Taxa_ID"))  # Get index of start column
    elif "Taxa_Source" in df.columns:
        taxa_start_col  = column_letter(list(df.columns).index("Taxa_Source"))  # Get index of start column
    else:
        raise ValueError("Final_Taxa_ID and Taxa_Source in the created dataframe. Something went wrong, please open a github issue to report the problem.")
    taxa_end_col = column_letter(list(df.columns).index("Species_Support_ANI"))  # Get index of end column
    # Dynamically merge based on start and end column
    worksheet.merge_range(f"{taxa_start_col}1:{taxa_end_col}1", "Taxonomic Information", cell_format_green)
    #MLST columns 
    # Define start and end column based on centar condition
    mlst_start_col = column_letter(list(df.columns).index("Primary_MLST_Scheme"))  # Get index of start column
    if centar:
        mlst_end_col = column_letter(list(df.columns).index("MLST Clade"))  # Use "MLST Clade" if centar is True
    else:
        mlst_end_col = column_letter(list(df.columns).index("Secondary_MLST_Alleles"))  # Otherwise, use "Secondary_MLST_Alleles"
    # Dynamically merge based on start and end column
    #worksheet.merge_range(f"{mlst_start_col}1:{mlst_end_col}1", "MLST Information", cell_format_green)
    if centar == True:
        # qc_max_col centar columns to make merging easier so we need to substract the total number of centar columns from the qc_max_col to get the right starting point
        # as part of combine_GRiPHins.py organism specifc columns are in qc_max_col
        qc_minus_centar = qc_max_col - sum(centar_df_lens)
        if centar_df_lens[0] <= 1:
            # Just write to the single cell
            worksheet.write(0, qc_minus_centar, "Toxin A/B Variants", cell_format_p4)
        else:
            # Safe to merge multiple cells
            worksheet.merge_range(0, (qc_minus_centar), 0, (qc_minus_centar + centar_df_lens[0] - 1), "Toxin A/B Variants", cell_format_p4)
        worksheet.merge_range(0, (qc_minus_centar + centar_df_lens[0]), 0, (qc_minus_centar + centar_df_lens[0] + centar_df_lens[1] - 1), "Other Toxins", cell_format_p3) #-1 is to account for MLST clade being in the MLST columns, but in the centar dataframe
        worksheet.merge_range(0, (qc_minus_centar + centar_df_lens[0] + centar_df_lens[1]), 0, (qc_minus_centar + centar_df_lens[0] + centar_df_lens[1] + centar_df_lens[2] - 1), "C. difficile Specific AR Mutations", cell_format_p2)
        worksheet.merge_range(0, (qc_minus_centar + centar_df_lens[0] + centar_df_lens[1] + centar_df_lens[2]), 0, (qc_minus_centar + sum(centar_df_lens) - 1), "ML Predicted Ribotype", cell_format_p1)
        worksheet.merge_range(0, (qc_max_col), 0, (qc_max_col + ar_gene_count)-1, "Antibiotic Resistance Genes", cell_format_lightgrey) 
        worksheet.merge_range(0, (qc_max_col + ar_gene_count), 0 ,(qc_max_col + ar_gene_count + hv_gene_count - 1), "Hypervirulence Genes^^", cell_format_grey)
        worksheet.merge_range(0, (qc_max_col + ar_gene_count + hv_gene_count), 0, (qc_max_col + ar_gene_count + pf_gene_count + hv_gene_count - 1), "Plasmid Incompatibility Replicons^^^", cell_format_darkgrey)
        # needed this for anoter set of samples... not sure what the differences are  -- cdc_phx with --centar
        #worksheet.merge_range(0, (qc_max_col), 0, (qc_max_col + centar_df_lens[0] - 2), "Toxin A/B Variants", cell_format_p4)
        #worksheet.merge_range(0, (qc_max_col + centar_df_lens[0] - 1), 0, (qc_max_col + centar_df_lens[0] + centar_df_lens[1] - 2), "Other Toxins", cell_format_p3)
        #worksheet.merge_range(0, (qc_max_col + centar_df_lens[0] + centar_df_lens[1] - 1), 0, (qc_max_col + centar_df_lens[0] + centar_df_lens[1] + centar_df_lens[2] - 2), "C. difficile Specific AR Mutations", cell_format_p2)
        #worksheet.merge_range(0, (qc_max_col + centar_df_lens[0] + centar_df_lens[1] + centar_df_lens[2] - 1), 0, (qc_max_col + sum(centar_df_lens) - 2), "ML Predicted Ribotype", cell_format_p1)
        #worksheet.merge_range(0, (qc_max_col + sum(centar_df_lens) - 1 ), 0, (qc_max_col + sum(centar_df_lens) + ar_gene_count - 2), "Antibiotic Resistance Genes", cell_format_lightgrey) #-1 is to account for MLST clade being in the MLST columns, but in the centar dataframe
        #worksheet.merge_range(0, (qc_max_col + sum(centar_df_lens) - 1 + ar_gene_count), 0 ,(qc_max_col + sum(centar_df_lens) + ar_gene_count + hv_gene_count - 2), "Hypervirulence Genes^^", cell_format_grey)
        #worksheet.merge_range(0, (qc_max_col + sum(centar_df_lens) - 1 + ar_gene_count + hv_gene_count), 0, (qc_max_col + sum(centar_df_lens) + ar_gene_count + pf_gene_count + hv_gene_count - 2), "Plasmid Incompatibility Replicons^^^", cell_format_darkgrey)
    else:
        worksheet.merge_range(0, qc_max_col, 0, (qc_max_col + ar_gene_count - 1), "Antibiotic Resistance Genes", cell_format_lightgrey)
        worksheet.merge_range(0, (qc_max_col + ar_gene_count), 0 ,(qc_max_col + ar_gene_count + hv_gene_count - 1), "Hypervirulence Genes^^", cell_format_grey)
        worksheet.merge_range(0, (qc_max_col + ar_gene_count + hv_gene_count), 0, (qc_max_col + ar_gene_count + pf_gene_count + hv_gene_count - 1), "Plasmid Incompatibility Replicons^^^", cell_format_darkgrey)
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
    worksheet.write('A' + str(max_row + 5), 'Cells in ORANGE denote “Big 5” carbapenemase gene (i.e., blaKPC, blaNDM, blaOXA48-like, blaVIM, and blaIMP) or an acquired blaOXA gene, please confirm what AR Lab Network HAI/AR WGS priority these meet.', orange_format_nb)
    worksheet.write('A' + str(max_row + 6), 'Cells in RED denote isolates that failed one or more auto failure triggers (cov < 30, assembly ratio stdev > 2.58, assembly length < 1Mbps, >500 scaffolds)', red_format)
    # More footers - Disclaimer etc.
    # unbold
    no_bold = workbook.add_format({'bold': False})
    worksheet.write_url('A' + str(max_row + 7), 'https://github.com/CDCgov/phoenix/wiki/Pipeline-Overview#mlst-allele-symbols', string="Click for a full explaination of symbols used in MLST allele markers. The source of the MLST determination can be 'assembly' (MLST), 'reads' (SRST2) or both 'assembly/reads'.")
    worksheet.write('A' + str(max_row + 8),"^Using Antibiotic Resistance Gene database " + ar_db + " (ResFinder, ARG-ANNOT, NCBI Bacterial Antimicrobial Resistance Reference Gene Database) using output thresholds ([98AA/90]G:[98NT/90]S); gene matches from S:(SRST2) with [%Nuc_Identity, %Coverage], or from G:(GAMMA) with [%Nuc_Identity, %AA_Identity,  %Coverage]; GAMMA gene matches indicate associated contig.", no_bold)
    worksheet.write('A' + str(max_row + 9),"^^Using CDC-compiled iroB, iucA, peg-344, rmpA, and rmpA2 hypervirulence gene database ( " + hv_db + " ); gene matches noted with [%Nuc_Identity, %AA_Identity,  %Coverage].", no_bold)
    worksheet.write('A' + str(max_row + 10),"^^^Using the plasmid incompatibility replicons plasmidFinder database ( " + pf_db + " ) using output thresholds [95NT/60]; replicon matches noted with [%Nuc_Identity, %Coverage].", no_bold)
    worksheet.write('A' + str(max_row + 11),"DISCLAIMER: These data are preliminary and subject to change. The identification methods used and the data summarized are for public health surveillance or investigational purposes only and must NOT be communicated to the patient, their care provider, or placed in the patient’s medical record. These results should NOT be used for diagnosis, treatment, or assessment of individual patient health or management.", bold)
    #adding review and date info
    worksheet.write('A' + str(max_row + 13), "Reviewed by:", no_bold)
    worksheet.write('D' + str(max_row + 13), "Date:")
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

def create_samplesheet(input_directory, scaffolds_entry):
    """Function will create a samplesheet from samples in a directory if -d argument passed."""
    directory = os.path.abspath(input_directory) # make sure we have an absolute path to start with
    with open("Directory_samplesheet.csv", "w") as samplesheet:
        samplesheet.write('sample,directory\n')
    dirs = sorted(os.listdir(input_directory))
    # Filter directories based on the presence of *_1.trim.fastq.gz files
    valid_directories = [ directory for directory in dirs if glob.glob(os.path.join(input_directory, directory, "*_summaryline.tsv")) ]
    files_glob = "*_summaryline.tsv"
    # Identify and warn about excluded directories
    excluded_dirs = [excluded_dir for excluded_dir in dirs if excluded_dir not in valid_directories]
    print(f"\n\033[93m Warning: The following directories '{excluded_dirs}' were excluded from analysis because no '{files_glob}' files weren't found in these locations.\033[0m\n")
    try: #if there are numbers in the name then use that to sort
        dirs_sorted=sorted(valid_directories, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    except: #if no numbers then use only alphabetically
        dirs_sorted=sorted(valid_directories)
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
    data_xlsx = data_xlsx.iloc[:-10] 
    #Write dataframe into csv
    data_xlsx.to_csv(output_file + '.tsv', sep='\t', encoding='utf-8',  index=False, lineterminator ='\n')

def get_second_dir(samplesheet, sample_name):
    """Function to get the second directory for updater."""
    with open(samplesheet, 'r') as f:
        # Skip header
        next(f)
        for line in f:
            # Check if line starts with sample name
            if line.startswith(sample_name + ','):
                # Split on comma and take the second item (directory)
                directory = line.split(',')[1].strip()
                # Pattern to capture the directory name immediately before /{sample_name}
                pattern = rf'/([^/]+)/{re.escape(sample_name)}'
                # Split on sample name and take the first part to get just the directory path
                return "./" + re.search(pattern, directory).group(1)
    # If sample not found, return empty string
    print(f"\033[93mWarning: Sample '{sample_name}' not found in samplesheet.\033[0m")
    return ""

def main():
    args = parseArgs()
    # create empty lists to append to later
    Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, Scaffold_Count_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
    busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L, data_location_L, parent_folder_L= ([] for i in range(36))
    ar_df = pd.DataFrame() #create empty dataframe to fill later for AR genes
    pf_df = pd.DataFrame() #create another empty dataframe to fill later for Plasmid markers
    hv_df = pd.DataFrame() #create another empty dataframe to fill later for hypervirulence genes
    srst2_ar_df = pd.DataFrame()
    shiga_df = pd.DataFrame()
    centar_dfs = []
    # Since srst2 currently doesn't handle () in the gene names we will make a quick detour to fix this... first making a dictionary
    ar_dic = make_ar_dictionary(args.ar_db)
    # check if a directory or samplesheet was given
    if (args.samplesheet == None) and (args.directory == None): # if no directory give AND no sample sheet given exit
        sys.exit(CRED + "You MUST pass EITHER a samplesheet or a top directory of PHoeNIx output to create one.\n" + CEND)
    # If a directory is given then create a samplesheet from it if not use the samplesheet passed
    if args.directory != None:
        samplesheet = create_samplesheet(args.directory, args.scaffolds)
    else:
        sort_samplesheet(args.samplesheet)
        samplesheet = args.samplesheet
    if args.updater == True:
        input_samplesheet_df = pd.read_csv(args.samplesheet)
        samples_to_run = input_samplesheet_df["sample"].tolist()
    if args.centar == True and args.samplesheet != None and args.filter_samples == True: 
        # When using species specific pipelines and --samplesheet is  given this means we need to make sure only samples in samplesheet are run
        input_samplesheet_df = pd.read_csv(args.samplesheet)
        output_dir_string = str(args.output).replace("_GRiPHin_Summary","").replace("_GRiPHin","")
        input_samplesheet_df = input_samplesheet_df[input_samplesheet_df["directory"].str.contains(fr"/{str(output_dir_string)}", na=False, regex=True)]
        samples_to_run = input_samplesheet_df["sample"].tolist()
    #input is a samplesheet that is "samplename,directory" where the directory is a phoenix like folder
    with open(samplesheet) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        header = next(csv_reader) # skip the first line of the samplesheet
        csv_rows = list(csv_reader)  # Convert the iterator to a list to reuse it
        if args.centar == True or args.updater == True and args.samplesheet != None and args.filter_samples == True:
            filtered_out_samples = [row[0] for row in csv_rows if any(sample not in row[0] for sample in samples_to_run)]
            csv_rows = [row for row in csv_rows if row[0] in samples_to_run]
            print("\n\033[93m Warning: The following sample(s) are not in samplesheet and were filtered out of reporting in griphin: {}\033[0m\n".format(list(set(filtered_out_samples) - set(samples_to_run))))
        for row in csv_rows:
            sample_name = row[0]
            if args.updater == True:
                # If updater is true then we need to get the second directory for updater
                directory2 = get_second_dir(args.samplesheet, sample_name)
            else:
                directory2 = ""
            directory = row[1]
            # check if species specific information is present
            data_location, parent_folder = Get_Parent_Folder(directory)
            trim_stats, raw_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld_report, kraken_wtasmbld, quast_report, mlst_file, fairy_file, busco_short_summary, asmbld_ratio, gc, gamma_ar_file, gamma_pf_file, gamma_hv_file, fast_ani_file, tax_file, srst2_file = Get_Files(directory, sample_name, directory2)
            #Get the metrics for the sample
            srst2_ar_df, pf_df, ar_df, hv_df, Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Seq_reads, Paired_Trimmed_reads, Total_trim_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, Scaffold_Count, busco_metrics, gc_metrics, assembly_ratio_metrics, QC_result, \
            QC_reason, MLST_scheme_1, MLST_scheme_2, MLST_type_1, MLST_type_2, MLST_alleles_1, MLST_alleles_2, MLST_source_1, MLST_source_2 = Get_Metrics(args.phoenix, args.scaffolds, args.set_coverage, srst2_ar_df, pf_df, ar_df, hv_df, trim_stats, raw_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld_report, kraken_wtasmbld, quast_report, busco_short_summary, asmbld_ratio, gc, sample_name, mlst_file, fairy_file, gamma_ar_file, gamma_pf_file, gamma_hv_file, fast_ani_file, tax_file, srst2_file, ar_dic)
            #Collect this mess of variables into appeneded lists
            data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L , alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L , MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L = Append_Lists(data_location, parent_folder, sample_name, \
            Q30_R1_per, Q30_R2_per, Total_Raw_Seq_bp, Total_Seq_reads, Paired_Trimmed_reads, Total_trim_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, Scaffold_Count, busco_metrics, \
            gc_metrics, assembly_ratio_metrics, QC_result, QC_reason, MLST_scheme_1, MLST_scheme_2, MLST_type_1, MLST_type_2, MLST_alleles_1, MLST_alleles_2, MLST_source_1, MLST_source_2, \
            data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L)
            if args.shigapass == True:
                shiga_df = create_shiga_df(directory, sample_name, shiga_df, FastANI_output_list[3])
            if args.centar == True:
                centar_df = create_centar_combined_df(directory, sample_name)
                centar_dfs.append(centar_df)
    # combine all lists into a dataframe
    df = Create_df(args.phx_version, args.phoenix, data_location_L, parent_folder_L, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Raw_Seq_bp_L, Total_Seq_reads_L, Paired_Trimmed_reads_L, Total_trim_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
    Scaffold_Count_L, busco_lineage_L, percent_busco_L, gc_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L , MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L)
    if args.shigapass == True:
        df = double_check_taxa_id(shiga_df, df)
    else:
        df['Final_Taxa_ID'] = df.apply(fill_taxa_id, axis=1)
    if args.centar == True:
        try:
            full_centar_df = pd.concat(centar_dfs, ignore_index=True) # combine rows of c diff samples into one c diff df
        except ValueError as e:
            if "No objects to concatenate" in str(e):
                print(CYELLOW + "\nThere was a ValueError: 'No objects to concatenate'. Check that --output is the same as the phx dir its used to get a path with --centar and --samplesheet. Search for 'output_dir_string' to find origin of the error.\n" + CEND)
            else:
                raise  # re-raise if it's not the one you expected
        ordered_centar_df, A_B_Tox_len, other_Tox_len, mutant_len, RB_type_len = clean_and_format_centar_dfs(full_centar_df)
        centar_df_lens = [ A_B_Tox_len, other_Tox_len, mutant_len, RB_type_len ]
        # combing centar with phx qc information
        df = pd.concat([df, ordered_centar_df], axis=1)
    else:
        ordered_centar_df = pd.DataFrame()
        A_B_Tox_len = other_Tox_len = mutant_len = RB_type_len = 0
        centar_df_lens = [0,0,0,0] # we sum later so needs to be a list of numbers
    (qc_max_row, qc_max_col) = df.shape
    pf_max_col = pf_df.shape[1] - 1 #remove one for the WGS_ID column
    hv_max_col = hv_df.shape[1] - 1 #remove one for the WGS_ID column
    final_df, ar_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db = Combine_dfs(df, ar_df, pf_df, hv_df, srst2_ar_df, args.phoenix, False, args.bldb)
    # Checking if there was a control sheet submitted
    if args.control_list !=None:
        final_df = blind_samples(final_df, args.control_list)
    else:
        final_df = final_df
    write_to_excel(args.set_coverage, args.output, final_df, qc_max_col, ar_max_col, pf_max_col, hv_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db, args.phoenix, args.shigapass, args.centar, centar_df_lens)
    convert_excel_to_tsv(args.output)

if __name__ == '__main__':
    main()
    
