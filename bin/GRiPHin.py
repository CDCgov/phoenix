#!/usr/bin/env python3

import sys
import glob
import os
from decimal import *
import pandas as pd
import numpy as np
import argparse
import json
import re
from re import search
import operator
import xlsxwriter as ws
from xlsxwriter.utility import xl_rowcol_to_cell
import csv
from Bio import SeqIO

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python GRiPHin.py -s ./samplesheet.csv -a ../PHX/phoenix/assets/databases/ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output
## Written by Jill Hagey (qpk9@cdc.gov)

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-s', '--samplesheet', default=None, required=False, dest='samplesheet', help='PHoeNIx style samplesheet of sample,directory in csv format. Directory is expected to have PHoeNIx stype output.')
    parser.add_argument('-d', '--directory', default=None, required=False, dest='directory', help='If a directory is given rather than samplesheet GRiPHin will create one for all samples in the directory.')
    parser.add_argument('-c', '--control_list', required=False, dest='control_list', help='CSV file with a list of sample_name,new_name. This option will output the new_name rather than the sample name to "blind" reports.')
    parser.add_argument('-a', '--ar_db', default=None, required=True, dest='ar_db', help='AR Gene Database file that is used to confirm srst2 gene names are the same as GAMMAs output.')
    parser.add_argument('-o', '--output', default="", required=False, dest='output', help='Name of output file default is GRiPHin_Report.xlsx.')
    parser.add_argument('-p', '--platform', default=None, required=False, dest='platform', help='String for the sequencing platform used.') #kinda dumping this and just calling it parent_folder
    parser.add_argument('--coverage', default=30, required=False, dest='set_coverage', help='The coverage cut off default is 30s.')
    parser.add_argument('--phoenix', dest="phoenix", default=False, action='store_true', required=False, help='Use for -entry PHOENIX rather than CDC_PHOENIX which is the default.')
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'

def Get_Parent_Folder(directory, platform):
    '''getting project and platform info from the paths'''
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
    cemb_path = os.path.split(os.path.split(os.path.split(os.path.split(directory)[0])[0])[0])[0]
    if platform == None:
        platform = os.path.split(cemb_path)[1].lstrip("/") # remove backslash on left side to make it clean
    return project, platform

def make_ar_dictionary(ar_db):
    seq_id_list = []
    ar_dic = {}
    for seq_record in SeqIO.parse(ar_db, "fasta"):
        seq_id_list.append(seq_record.id)
    gene_name_list = [seq_id.split("__")[2] for seq_id in seq_id_list]
    drug_list = [seq_id.split("__")[4]  for seq_id in seq_id_list]
    ar_dic = dict(zip(gene_name_list, drug_list))
    return ar_dic

def get_Q30(trim_stats):
    data_df = pd.read_csv(trim_stats, sep='\t', header=0)
    Q30_R1_percent = round(data_df["Q30_R1_[%]"]*100, 2)[0] #make percent and round to two decimals
    Q30_R2_percent = round(data_df["Q30_R2_[%]"]*100, 2)[0]
    Total_Sequenced_bp = data_df["Total_Sequenced_[bp]"][0]
    Total_Sequenced_reads = data_df["Total_Sequenced_[reads]"][0]
    return Q30_R1_percent, Q30_R2_percent, Total_Sequenced_bp, Total_Sequenced_reads

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
    return Trim_kraken, Asmbld_kraken, Trim_unclassified_percent, Asmbld_unclassified_percent

def Calculate_Trim_Coverage(Total_Sequenced_bp, quast_report):
    """Taking the total sequenced bp from fastp and the assembly length from quast to calculate coverage."""
    with open(quast_report, 'r') as f:
        for line in f:
            if ('Total length' in line):
                Assembly_Length = int(line.split('\t')[1])
                break
    Coverage = float(round(Total_Sequenced_bp / Assembly_Length, 2))
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

def get_assembly_ratio(asmbld_ratio, tax_file):
    '''Collects the assembly ratio, stdev and method used to determine taxa.'''
    with open(asmbld_ratio, 'r') as f:
        for line in f:
            if "Tax: " in line:
                taxa =  (line.split("Tax: ",1)[1]).strip()
            if "Isolate_St.Devs:" in line:
                stdev = (line.split("Isolate_St.Devs: ",1)[1]).strip()
                if stdev == 'N/A': #handling making stdev a float only if its a number
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

def compile_alerts(coverage, assembly_stdev):
    """
    No orphaned reads found after trimming
    <10 reference genomes for species identified so no stdev for assembly ratio or %GC content calculated
    >150x coverage or <40x coverage
    """
    alerts = []
    if str(assembly_stdev) == "N/A":
        alerts.append("STDev was N/A, <10 genomes as reference")
    if coverage != "Unknown": # if its unknown it will fail already so skip
        if int(coverage) > 30 and int(coverage) < 40:
            alerts.append("coverage between 30-40x("+ str(coverage) + "x)")
        elif int(coverage) > 100.00:
            alerts.append("coverage above 100x(" + str(coverage) + "x)")
    alerts = ', '.join(alerts)
    return alerts

def compile_warnings(Total_Seq_reads, Q30_R1_per, Q30_R2_per, scaffolds, assembly_ratio_metrics, Trim_unclassified_percent, Asmbld_unclassified_percent, kraken_trim_genus):
    """
    <1,000,000 total reads for each raw and trimmed reads - Total_Sequenced_reads
    % reads with Q30 average for R1 (<90%) and R2 (<70%) - Q30_R1_percent, Q30_R2_percent
    >200 scaffolds - scaffolds
    Checking that %GC content is within 2.58 stdev away from the mean %GC content for the species determined - assembly_ratio_metrics
    Contamination check: >30% unclassified reads and confirm there is only 1 genera with >25% of assigned reads - Trim_kraken, Asmbld_kraken
    """
    warnings = []
    if Total_Seq_reads == "Unknown" or int(Total_Seq_reads) < int(1000000):
        warnings.append("<1,000,000 reads")
    if Q30_R1_per == "Unknown" or float(Q30_R1_per) < float(90.00):
        warnings.append("<{:.2f}% reads".format(float(90.00)))
    if Q30_R2_per == "Unknown" or float(Q30_R1_per) < float(70.00):
        warnings.append("<{:.2f}% reads".format(int(70.00)))
    if scaffolds == "Unknown" or int(scaffolds) > int(200):
        warnings.append(">200 scaffolds".format(int(70)))
    if Trim_unclassified_percent == "Unknown" or float(Trim_unclassified_percent) > float(30.00):
        warnings.append(">{:.2f}% unclassifed trimmed reads".format(int(30)))
    if Asmbld_unclassified_percent == "Unknown" or float(Asmbld_unclassified_percent) > float(30.00):
        warnings.append(">{:.2f}% unclassifed scaffolds".format(int(30)))
    if len(kraken_trim_genus) >=2:
        warnings.append(">=2 genus had >{:.2f}% of reads assigned to them.".format(int(25)))
    warnings = ', '.join(warnings)
    return warnings

def parse_kraken_report(kraken_trim_report):
    """Checking that only 1 genera with >25% of assigned reads."""
    kraken_trim_genus = []
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
        print("Warning: " + sample_name + ".kraken2_trimd.report.txt not found")
        kraken_trim_genus = 'Unknown'
    return kraken_trim_genus

def Checking_auto_pass_fail(coverage, length, assembly_stdev, asmbld_ratio, set_coverage):
    """Checking auto pass fail conditions"""
    #assembly_stdev = assembly_ratio_line.split("(")[1].split(")")[0].split(" ")[1] # parse to get standard dev, old method
    if coverage == "Unknown" or int(coverage) < int(set_coverage):
        QC_result = "FAIL"
        QC_reason = "coverage below "+ str(set_coverage) +"x(" + str(coverage) + "x)"
    elif int(length) <= 1000000 or length == "Unknown":
        QC_result = "FAIL"
        QC_reason = "smaller than 1,000,000bps(" + str(length) + ")"
    elif str(assembly_stdev) != "N/A": # have to have a second layer cuz you can't make NA a float, N/A means less than 10 genomes so no stdev calculated
        if str(asmbld_ratio) == "Unknown": # if there is no ratio file then fail the sample
            QC_result = "FAIL"
            QC_reason="Assembly file not Found"
        elif float(assembly_stdev) > 2.58:
            QC_result = "FAIL"
            QC_reason="assmebly stdev above 2.58(" + str(assembly_stdev) + ")"
        else:
            QC_result = "PASS"
            QC_reason = ""
    else:
        QC_result = "PASS"
        QC_reason = ""
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
    #print(df["blaFOX-5_NG_049105.1(beta-lactam)"])
    df = duplicate_column_clean(df)
    final_df = pd.concat([final_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    #print(final_df["mph(D)_NC_017312(macrolide_lincosamide_streptogramin)"])
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

def parse_mlst(mlst_file):
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
            #print(split_line[5:])
            alleles = "-".join(split_line[5:]) # combine all alleles separated by -
            if DB_ID in Scheme_list[0]: # check if scheme name is already in the scheme list
                for i in range(0,len(Scheme_list[0])): #loop through list of scheme names
                    if DB_ID == Scheme_list[0][i]: # looking for matching scheme name that was already in the list 
                        # If the scheme was already in the list then add the ST, alleles, source and data into that list within for that scheme
                        # Example: [['abaumannii(Pasteur)', 'abaumannii(Oxford)'], [['ST2'], ['ST195', 'ST1816-PARALOG']], [['cpn60(2)-fusA(2)-gltA(2)-pyrG(2)-recA(2)-rplB(2)-rpoB(2)'], ['gltA(1)-gyrB(3)-gdhB(3)-recA(2)-cpn60(2)-gpi(96)-rpoD(3)', 'gltA(1)-gyrB(3)-gdhB(189)-recA(2)-cpn60(2)-gpi(96)-rpoD(3)']], [['standard/srst2'], ['standard', 'standard/srst2']], [['2022-12-02'], ['2022-12-02', '2022-12-02']]]
                        if Scheme != "-" and Scheme != "Novel_allele" and Scheme != "Novel_profile":
                            Scheme_list[1][i].append("ST"+str(Scheme)) # if there is a value add ST in front of number
                        else:
                            Scheme_list[1][i].append(Scheme) # just append what is there
                        Scheme_list[2][i].append(alleles)
                        Scheme_list[3][i].append(source)
                        Scheme_list[4][i].append(date)
            else: # if scheme name is not already in the scheme list add it
                Scheme_list[0].append(DB_ID) 
                if Scheme != "-" and Scheme != "Novel_allele" and Scheme != "Novel_profile":
                    Scheme_list[1].append(["ST"+Scheme]) # if there is a value add ST in front of number
                else:
                    Scheme_list[1].append([Scheme]) # just append what is there
                Scheme_list[2].append([alleles])
                Scheme_list[3].append([source])
                Scheme_list[4].append([date])
    #print(Scheme_list)
    return Scheme_list

def parse_ani(fast_ani_file):
    """Parse ANI file to get format 99.98%ID-98.58%COV-Acinetobacter baumannii(Acinetobacter_baumannii_GCF_012935145.1_ASM1293514v1_genomic.fna.gz)."""
    ani_df = pd.read_csv(fast_ani_file, sep='\t', header=0) # should only be one line long.
    ID = ani_df["% ID"][0]
    coverage = ani_df["% Coverage"][0]
    organism = ani_df["Organism"][0]
    source_file = ani_df["Source File"][0]
    #Species_Support = str(ID) + "%ID-" + str(coverage) + "%COV-" + organism + "(" + source_file + ")" #old way of reporting
    FastANI_output_list = [source_file, ID, coverage, organism]
    return FastANI_output_list

def parse_srst2_ar(srst2_file, ar_dic, final_srst2_df, sample_name):
    """Parsing the srst2 file run on the ar gene database."""
    srst2_df = pd.read_csv(srst2_file, sep='\t', header=0)
    percent_lengths = np.floor(srst2_df["coverage"]).tolist()
    genes = srst2_df["allele"].tolist()
    percent_BP_IDs = np.floor(100 - srst2_df["divergence"]).tolist()
    # Since srst2 currently doesn't handle () in the gene names we will make a quick detour to fix this... now fixing annotations
    #srst2_df.annotation = srst2_df.annotation.fillna(srst2_df.allele.map(ar_dic)) # this only fills in nas
    srst2_df['conferred_resistances'] = srst2_df['allele'].map(ar_dic)
    conferred_resistances = srst2_df['conferred_resistances'].tolist()
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

def Get_Metrics(set_coverage, srst2_ar_df, pf_df, ar_df, hv_df, trim_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld, quast_report, busco_short_summary, asmbld_ratio, sample_name, mlst_file, gamma_ar_file, gamma_pf_file, gamma_hv_file, fast_ani_file, tax_file, srst2_file, ar_dic):
    '''For each step to gather metrics try to find the file and if not then make all variables unknown'''
    try:
        Q30_R1_per, Q30_R2_per, Total_Seq_bp, Total_Seq_reads = get_Q30(trim_stats)
    except FileNotFoundError:
        print("Warning: " + sample_name + "_trimmed_read_counts.txt not found")
        Q30_R1_per = Q30_R2_per = Total_Seq_bp = Total_Seq_reads = 'Unknown'
    # Try and except are in the get_kraken_info function to allow for cases where trimming was completed, but not assembly
    Trim_kraken, Asmbld_kraken, Trim_unclassified_percent, Asmbld_unclassified_percent = get_kraken_info(kraken_trim, kraken_wtasmbld, sample_name)
    try:
        if Total_Seq_bp != "Unknown": # for -entry CDC_SCAFFOLDS where no reads are present to calculate. For all other scenerios run this portion
            Coverage, Assembly_Length = Calculate_Trim_Coverage(Total_Seq_bp, quast_report)
        elif Total_Seq_bp == "Unknown": # for -entry CDC_SCAFFOLDS and -entry SCAFFOLDS
            Assembly_Length = Get_Assembly_Length(quast_report)
            Coverage = 'Unknown'
        else:
            Coverage = Assembly_Length = 'Unknown'
    except FileNotFoundError:
        print("Warning: " + sample_name + "_report.tsv not found")
        Coverage = Assembly_Length = 'Unknown'
    try:
        Scaffold_Count = get_scaffold_count(quast_report)
    except FileNotFoundError:
        print("Warning: " + sample_name + "_report.tsv not found")
        Scaffold_Count = 'Unknown'
    try:
        busco_metrics = Get_BUSCO_Gene_Count(busco_short_summary)
    except FileNotFoundError:
        print("Warning: short_summary.specific." + sample_name + ".filtered.scaffolds.fa.txt not found.")
        lineage = percent_busco = 'Unknown'
        busco_metrics = [lineage, percent_busco]
    try:
        assembly_ratio_metrics = get_assembly_ratio(asmbld_ratio, tax_file)
    except FileNotFoundError:
        print("Warning: " + sample_name + "_Assembly_ratio_*.txt not found.")
        ratio = stdev = tax_method = 'Unknown'
        assembly_ratio_metrics = [ratio, stdev, tax_method]
    try:
        QC_result, QC_reason = Checking_auto_pass_fail(Coverage, Assembly_Length, assembly_ratio_metrics[1], assembly_ratio_metrics[0], set_coverage)
    except FileNotFoundError: 
        print("Warning: Possibly coverage and assembly length was not calculated and/or"+ sample_name + "_Assembly_ratio_*.txt not found.")
        QC_result = QC_reason = 'Unknown'
    try:
        FastANI_output_list = parse_ani(fast_ani_file)
    except FileNotFoundError: 
        print("Warning: " + sample_name + ".fastANI.txt not found")
        ani_source_file = fastani_ID = fastani_coverage = fastani_organism = 'Unknown'
        FastANI_output_list = [ani_source_file, fastani_ID, fastani_coverage, fastani_organism]
    try:
        Scheme_list = parse_mlst(mlst_file)
        if len(Scheme_list[0]) > 1: # If there is more than one scheme
            if Scheme_list[0][0] < Scheme_list[0][1]: # this if else is all just to make sure things are printing out in the same order.
                MLST_scheme_1 = Scheme_list[0][0] # get 1st scheme name from the list
                mlst_types_1=sorted(Scheme_list[1][0])[::-1]
                MLST_type_1 = ", ".join(mlst_types_1)
                MLST_alleles_1 = ",".join(Scheme_list[2][0])
                MLST_source_1 = Scheme_list[3][0]
                MLST_scheme_2 = Scheme_list[0][1] # get 2nd scheme name from the list
                mlst_types_2=sorted(Scheme_list[1][1])[::-1]
                MLST_type_2 = ", ".join(mlst_types_2)
                MLST_alleles_2 = ",".join(Scheme_list[2][1])
                MLST_source_2 = Scheme_list[3][1]
            else:
                MLST_scheme_1 = Scheme_list[0][1] # get 1st scheme name from the list, in this case its the 2nd element
                mlst_types_1=sorted(Scheme_list[1][1])[::-1]
                MLST_type_1 = ", ".join(mlst_types_1)
                MLST_alleles_1 = ",".join(Scheme_list[2][1])
                MLST_source_1 = Scheme_list[3][1]
                MLST_scheme_2 = Scheme_list[0][0] # get 2nd scheme name from the list, in this case its the first element
                mlst_types_2=sorted(Scheme_list[1][0])[::-1]
                MLST_type_2 = ", ".join(mlst_types_2)
                MLST_alleles_2 = ",".join(Scheme_list[2][0])
                MLST_source_2 = Scheme_list[3][0]
        else: # If there is only one scheme then the last scheme and type are just "-"
            MLST_scheme_1 = Scheme_list[0][0]
            MLST_type_1 = ", ".join(Scheme_list[1][0]) # join together the STs for this one scheme
            MLST_alleles_1 = ",".join(Scheme_list[2][0])
            MLST_source_1 = Scheme_list[3][0]
            MLST_scheme_2 = "-"
            MLST_type_2 = "-"
            MLST_alleles_2 = "-"
            MLST_source_2 = "-"
    except FileNotFoundError: 
        print("Warning: " + sample_name + "_combined.tsv not found")
        MLST_scheme_1 = MLST_scheme_2 = MLST_type_1 = MLST_type_2 = MLST_alleles_1 = MLST_alleles_2 = MLST_source_1 = MLST_source_2 = 'Unknown'
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
        srst2_ar_df = parse_srst2_ar(srst2_file, ar_dic, srst2_ar_df, sample_name)
    except (FileNotFoundError, pd.errors.EmptyDataError) : # second one for an empty dataframe - srst2 module creates a blank file 
        print("Warning: " + sample_name + "__fullgenes__ResGANNCBI__*_srst2__results.txt not found")
        df = pd.DataFrame({'WGS_ID':[sample_name]})
        df.index = [sample_name]
        srst2_ar_df = pd.concat([srst2_ar_df, df], axis=0, sort=True, ignore_index=False).fillna("")
    #try:
    alerts = compile_alerts(Coverage, assembly_ratio_metrics[1])
    # try except in the function itself
    kraken_trim_genus = parse_kraken_report(kraken_trim_report)
    try:
        warnings = compile_warnings(Total_Seq_reads, Q30_R1_per, Q30_R2_per, Scaffold_Count, assembly_ratio_metrics, Trim_unclassified_percent, Asmbld_unclassified_percent, kraken_trim_genus)
    except:
        warnings = ""
    return srst2_ar_df, pf_df, ar_df, hv_df, Q30_R1_per, Q30_R2_per, Total_Seq_bp, Total_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, \
    Scaffold_Count, busco_metrics, assembly_ratio_metrics, QC_result, QC_reason, MLST_scheme_1, MLST_scheme_2, MLST_type_1, MLST_type_2, MLST_alleles_1, MLST_alleles_2, MLST_source_1, MLST_source_2

def Get_Files(directory, sample_name):
    '''Create file paths to collect files from sample folder.'''
    # if there is a trailing / remove it
    directory = directory.rstrip('/')
    # create file names
    trim_stats = directory + "/qc_stats/" + sample_name + "_trimmed_read_counts.txt"
    kraken_trim = directory + "/kraken2_trimd/" + sample_name + ".trimd_summary.txt"
    kraken_trim_report = directory + "/kraken2_trimd/" + sample_name + ".kraken2_trimd.report.txt"
    kraken_wtasmbld = directory + "/kraken2_asmbld_weighted/" + sample_name + ".wtasmbld_summary.txt"
    quast_report = directory + "/quast/" + sample_name + "_report.tsv"
    mlst_file = directory + "/mlst/" + sample_name + "_combined.tsv"
    # This creates blank files for if no file exists. Varibles will be made into "Unknown" in the Get_Metrics function. Need to only do this for files determined by glob
    # You only need this for glob because glob will throw an index error if not.
    try:
        busco_short_summary =  glob.glob(directory + "/BUSCO/short_summary.specific.*" + sample_name + ".filtered.scaffolds.fa.txt")[0]
    except IndexError:
        busco_short_summary =  directory + "/BUSCO/short_summary.specific.blank" + sample_name + ".filtered.scaffolds.fa.txt"
    try:
        asmbld_ratio = glob.glob(directory + "/" + sample_name + "_Assembly_ratio_*.txt")[0]
    except IndexError:
        asmbld_ratio = directory + "/" + sample_name + "_Assembly_ratio_blank.txt"
    try:
        gamma_ar_file = glob.glob(directory + "/gamma_ar/" + sample_name + "_*.gamma")[0]
    except IndexError:
        gamma_ar_file = directory + "/gamma_ar/" + sample_name + "_blank.gamma"
    try:
        gamma_pf_file = glob.glob(directory + "/gamma_pf/" + sample_name + "_*.gamma")[0]
    except IndexError:
        gamma_pf_file = directory + "/gamma_pf/" + sample_name + "_blank.gamma"
    try: 
        gamma_hv_file = glob.glob(directory + "/gamma_hv/" + sample_name + "_*.gamma")[0]
    except IndexError:
        gamma_hv_file = directory + "/gamma_hv/" + sample_name + "_blank.gamma"
    fast_ani_file = directory + "/ANI/" + sample_name + ".fastANI.txt"
    tax_file = directory + "/" + sample_name + ".tax" # this file will tell you if kraken2 wtassembly, kraken2 trimmed (reads) or fastani determined the taxa
    try:
        srst2_file = glob.glob(directory + "/srst2/" + sample_name + "__fullgenes__*_srst2__results.txt")[0]
    except IndexError:
        srst2_file = directory + "/srst2/" + sample_name + "__fullgenes__blank_srst2__results.txt"
    return trim_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld, quast_report, mlst_file, busco_short_summary, asmbld_ratio, gamma_ar_file, gamma_pf_file, gamma_hv_file, fast_ani_file, tax_file, srst2_file

def Append_Lists(data_location, platform, sample_name, Q30_R1_per, Q30_R2_per, Total_Seq_bp, Total_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, \
            Scaffold_Count, busco_metrics, assembly_ratio_metrics, QC_result, QC_reason, MLST_scheme_1, MLST_scheme_2, MLST_type_1, MLST_type_2, MLST_alleles_1, MLST_alleles_2, MLST_source_1, MLST_source_2, \
            Data_Locations, Platforms, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Seq_bp_L, Total_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L):
        Data_Locations.append(data_location)
        Platforms.append(platform)
        Sample_Names.append(str(sample_name))
        Q30_R1_per_L.append(Q30_R1_per)
        Q30_R2_per_L.append(Q30_R2_per)
        Total_Seq_bp_L.append(Total_Seq_bp)
        Total_Seq_reads_L.append(Total_Seq_reads)
        Trim_kraken_L.append(Trim_kraken)
        Asmbld_kraken_L.append(Asmbld_kraken)
        Coverage_L.append(Coverage)
        Assembly_Length_L.append(Assembly_Length)
        Species_Support_L.append(FastANI_output_list[0])
        fastani_organism_L.append(FastANI_output_list[3])
        fastani_ID_L.append(FastANI_output_list[1])
        fastani_coverage_L.append(FastANI_output_list[2])
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
        return Data_Locations, Platforms, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Seq_bp_L, Total_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
        Scaffold_Count_L, busco_lineage_L, percent_busco_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L

def Create_df(phoenix, Data_Locations, Platforms, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Seq_bp_L, Total_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L,
Scaffold_Count_L, busco_lineage_L, percent_busco_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L):
    #combine all metrics into a dataframe
    if phoenix == True:
        data = {'WGS_ID'             : Sample_Names,
        'Parent_Folder'              : Platforms,
        'Data_Location'              : Data_Locations,
        'Minimum_QC_Check'           : QC_result_L,
        'QC_Issues'                  : QC_reason_L,
        'Warnings'                   : warnings_L,
        'Alerts'                     : alerts_L,
        'Q30_R1_[%]'                 : Q30_R1_per_L,
        'Q30_R2_[%]'                 : Q30_R2_per_L,
        'Total_Sequenced_[bp]'       : Total_Seq_bp_L,
        'Total_Sequenced_[reads]'    : Total_Seq_reads_L,
        'Estimated_Trimmed_Coverage' : Coverage_L,
        'Scaffolds'                  : Scaffold_Count_L,
        'Assembly_Length'            : Assembly_Length_L,
        'Assembly_Ratio'             : assembly_ratio_L,
        'Assembly_StDev'             : assembly_stdev_L,
        'Taxa_Source'                : tax_method_L,
        'Kraken_ID_Raw_Reads_%'      : Trim_kraken_L,
        'Kraken_ID_WtAssembly_%'     : Asmbld_kraken_L,
        'FastANI_Organism'           : fastani_organism_L, 
        'FastANI_%ID'                : fastani_ID_L, 
        'FastANI_%Coverage'          : fastani_coverage_L,
        'Species_Support_ANI'        : Species_Support_L,
        'Primary_MLST_Scheme_Name'   : MLST_scheme_1_L,
        'Primary_MLST'               : MLST_type_1_L,
        'Primary_MLST_alleles'       : MLST_alleles_1_L,
        'Secondary_MLST_Scheme_Name' : MLST_scheme_2_L,
        'Secondary_MLST'             : MLST_type_2_L,
        'Secondary_MLST_alleles'     : MLST_alleles_2_L}
    else:
        data = {'WGS_ID'                 : Sample_Names,
        'Parent_Folder'              : Platforms,
        'Data_Location'              : Data_Locations,
        'Minimum_QC_Check'           : QC_result_L,
        'QC_Issues'                  : QC_reason_L,
        'Warnings'                   : warnings_L,
        'Alerts'                     : alerts_L,
        'Q30_R1_[%]'                 : Q30_R1_per_L,
        'Q30_R2_[%]'                 : Q30_R2_per_L,
        'Total_Sequenced_[bp]'       : Total_Seq_bp_L,
        'Total_Sequenced_[reads]'    : Total_Seq_reads_L,
        'Estimated_Trimmed_Coverage' : Coverage_L,
        'Scaffolds'                  : Scaffold_Count_L,
        'Assembly_Length'            : Assembly_Length_L,
        'Assembly_Ratio'             : assembly_ratio_L,
        'Assembly_StDev'             : assembly_stdev_L,
        'Taxa_Source'                : tax_method_L,
        'BUSCO_Lineage'              : busco_lineage_L,
        'BUSCO_%Match'               : percent_busco_L,
        'Kraken_ID_Raw_Reads_%'      : Trim_kraken_L,
        'Kraken_ID_WtAssembly_%'     : Asmbld_kraken_L,
        'FastANI_Organism'           : fastani_organism_L, 
        'FastANI_%ID'                : fastani_ID_L, 
        'FastANI_%Coverage'          : fastani_coverage_L,
        'Species_Support_ANI'        : Species_Support_L,
        'Primary_MLST_Scheme_Name'   : MLST_scheme_1_L,
        'Primary_MLST'               : MLST_type_1_L,
        'Primary_MLST_alleles'       : MLST_alleles_1_L,
        'Secondary_MLST_Scheme_Name' : MLST_scheme_2_L,
        'Secondary_MLST'             : MLST_type_2_L,
        'Secondary_MLST_alleles'     : MLST_alleles_2_L}
    df = pd.DataFrame(data)
    return df

def add_srst2(ar_df, srst2_ar_df):
    ar_combined_df = pd.DataFrame() #create new dataframe to fill
    ar_combined_ordered_df = pd.DataFrame() #create new dataframe to fill
    common_cols = ar_df.columns.intersection(srst2_ar_df.columns) #get column names that are in both dataframes
    # Combine values in cells for columns that are in both dataframes
    for col in common_cols:
        if col != "WGS_ID":
            ar_combined_df[col] = (srst2_ar_df[col].map(str) + ":" + ar_df[col]).replace(':', "")
            ar_combined_df[col] = ar_combined_df[col].map(lambda x: str(x).lstrip(':').rstrip(':')) # clean up : for cases where there isn't a gamma and srst2 for all rows
            ar_combined_df = ar_combined_df.copy() #defragment to correct "PerformanceWarning: DataFrame is highly fragmented."
        else:
            ar_combined_df[col] = srst2_ar_df[col]
    # check if you missed any rows, if there is a sample in ar_db, that is not in the srst2 then you will have it have NA in rows when joined
    # drop columns from srst2 dataframe that are in common in the ar_db as these are already in ar_combined_df
    srst2_ar_df.drop(common_cols, axis = 1, inplace=True)
    ar_df.drop(common_cols, axis = 1, inplace=True)
    # Add cols that are unique to srst2
    ar_combined_df = ar_combined_df.join(srst2_ar_df)
    # Add cols that are unique to gamma ar_df
    ar_combined_df = ar_combined_df.join(ar_df)
    #fixing column orders
    ar_combined_ordered_df = pd.concat([ar_combined_ordered_df, ar_combined_df[['AR_Database', 'WGS_ID']]], axis=1, sort=False) # first adding back in ['AR_Database', 'WGS_ID']
    ar_drugs_list = ar_combined_df.columns.str.extract('.*\((.*)\).*').values.tolist() # get all ar drug names form column names
    sorted_list = sorted(list(set([str(drug) for sublist in ar_drugs_list for drug in sublist]))) #get unique drug names (with set) and sort list
    sorted_drug_names = [x for x in sorted_list if x != 'nan'] #get unique drug names (with set) and drop nan that comes from WGS_ID column and sort
    #sorted_drug_names = sorted(list(set([str(drug) for sublist in ar_drugs_list for drug in sublist]))[1:])
    # loop over each gene with the same drug its name
    for drug in sorted_drug_names:
        drug = "(" + drug + ")"
        column_list = sorted([col for col in ar_combined_df.columns if drug in col]) # get column names filtered for each drug name
        ar_combined_ordered_df = pd.concat([ar_combined_ordered_df, ar_combined_df[column_list]], axis=1, sort=False) # setting column's order by combining dataframes
    return ar_combined_ordered_df

def big5_check(final_ar_df):
    """"Function that will return list of columns to highlight if a sample has a hit for a big 5 gene."""
    columns_to_highlight = []
    final_ar_df = final_ar_df.drop(['AR_Database','WGS_ID'], axis=1)
    all_genes = final_ar_df.columns.tolist()
    big5_keep = [ "blaIMP", "blaVIM", "blaNDM", "blaKPC"] # list of genes to highlight
    blaOXA_48_like = [ "blaOXA-48", "blaOXA-54", "blaOXA-162", "blaOXA-181", "blaOXA-199", "blaOXA-204", "blaOXA-232", "blaOXA-244", "blaOXA-245", "blaOXA-247", "blaOXA-252", "blaOXA-370", "blaOXA-416", "blaOXA-436", \
    "blaOXA-438", "blaOXA-439", "blaOXA-484", "blaOXA-505", "blaOXA-514", "blaOXA-515", "blaOXA-517", "blaOXA-519", "blaOXA-535", "blaOXA-538", "blaOXA-546", "blaOXA-547", "blaOXA-566", "blaOXA-567", "blaOXA-731", \
    "blaOXA-788", "blaOXA-793", "blaOXA-833", "blaOXA-894", "blaOXA-918", "blaOXA-920", "blaOXA-922", "blaOXA-923", "blaOXA-924", "blaOXA-929", "blaOXA-933", "blaOXA-934", "blaOXA-1038", "blaOXA-1039", "blaOXA-1055", "blaOXA-1119", "blaOXA-1146" ]
    # combine lists of all genes we want to highlight
    #all_big5_keep = big5_keep + blaOXA_48_like
    # remove list of genes that look like big 5 but don't have activity
    big5_drop = [ "blaKPC-62", "blaKPC-63", "blaKPC-64", "blaKPC-65", "blaKPC-66", "blaKPC-72", "blaKPC-73", "blaOXA-163", "blaOXA-405"]
    # loop through column names and check if they contain a gene we want highlighted. Then add to highlight list if they do. 
    for gene in all_genes: # loop through each gene in the dataframe of genes found in all isolates
        if gene == 'No_AR_Genes_Found':
            pass
        else:
            gene_name = gene.split('_(')[0] # remove drug name for matching genes
            drug = gene.split('_(')[1] # keep drug name to add back later
            # make sure we have a complete match for blaOXA-48 and blaOXA-48-like genes
            if gene_name.startswith("blaOXA"): #check for complete blaOXA match
                [ columns_to_highlight.append(gene_name + "(" + drug) for big5_keep_gene in blaOXA_48_like if gene_name == big5_keep_gene ]
            else: # for "blaIMP", "blaVIM", "blaNDM", and "blaKPC", this will take any thing with a matching substring to these
                for big5 in big5_keep:
                    if search(big5, gene_name): #search for big5 gene substring in the gene name
                        columns_to_highlight.append(gene_name + "(" + drug)
    #loop through list of genes to drop and removed if they are in the highlight list
    for bad_gene in big5_drop:
        #search for big5 gene substring in the gene name and remove if it is
        [columns_to_highlight.remove(gene) for gene in columns_to_highlight if bad_gene in gene]
    return columns_to_highlight

def Combine_dfs(df, ar_df, pf_df, hv_df, srst2_ar_df, phoenix):
    hv_cols = list(hv_df)
    pf_cols = list(pf_df)
    ar_cols = list(ar_df)
    # move the column to head of list using index, pop and insert
    pf_cols.insert(0, pf_cols.pop(pf_cols.index('No_Plasmid_Markers')))
    pf_cols.insert(0, pf_cols.pop(pf_cols.index('Plasmid_Replicon_Database')))
    hv_cols.insert(0, hv_cols.pop(hv_cols.index('No_HVGs_Found')))
    hv_cols.insert(0, hv_cols.pop(hv_cols.index('HV_Database')))
    ar_cols.insert(0, ar_cols.pop(ar_cols.index('AR_Database')))
    # use ix to reorder
    pf_df = pf_df.loc[:, pf_cols]
    hv_df = hv_df.loc[:, hv_cols]
    ar_df = ar_df.loc[:, ar_cols]
    # if we run -entry PHOENIX then skip
    if phoenix == True: 
        # combining srst2 and gamma ar dataframes
        final_ar_df = add_srst2(ar_df, srst2_ar_df)
    else:
        final_ar_df = ar_df
    ar_max_col = final_ar_df.shape[1] - 1 #remove one for the WGS_ID column
    # now we will check for the "big 5" genes for highlighting later.
    columns_to_highlight = big5_check(final_ar_df)
    # combining all dataframes
    final_df = pd.merge(df, final_ar_df, how="left", on=["WGS_ID","WGS_ID"])
    final_df = pd.merge(final_df, hv_df, how="left", on=["WGS_ID","WGS_ID"])
    final_df = pd.merge(final_df, pf_df, how="left", on=["WGS_ID","WGS_ID"])
    #get database names
    ar_db = final_df['AR_Database'].tolist()[0]
    hv_db = final_df['HV_Database'].tolist()[0]
    pf_db = final_df['Plasmid_Replicon_Database'].tolist()[0]
    return final_df, ar_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db,

def write_to_excel(set_coverage, output, df, qc_max_col, ar_gene_count, pf_gene_count, hv_gene_count, columns_to_highlight, ar_df, pf_db, ar_db, hv_db, phoenix):
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    if output != "":
        writer = pd.ExcelWriter((output + '_GRiPHin_Report.xlsx'), engine='xlsxwriter')
    else:
        writer = pd.ExcelWriter(('GRiPHin_Report.xlsx'), engine='xlsxwriter')
    # Convert the dataframe to an XlsxWriter Excel object.
    df.to_excel(writer, sheet_name='Sheet1', index=False, startrow=1)
    # Get the xlsxwriter workfbook worksheet objects for formating
    workbook = writer.book
    (max_row, max_col) = df.shape # Get the dimensions of the dataframe.
    worksheet = writer.sheets['Sheet1']
    # Setting columns to numbers so you can have commas that make it more human readable
    number_comma_format = workbook.add_format({'num_format': '#,##0'})
    number_comma_format.set_align('left')
    worksheet.set_column('J:K', None, number_comma_format) # Total_seqs Total_bp
    worksheet.set_column('M:N', None, number_comma_format) # scaffolds and assembly_length
    # Setting columns to float so its more human readable
    #number_dec_format = workbook.add_format({'num_format': '0.000'})
    #number_dec_format.set_align('left')
    number_dec_2_format = workbook.add_format({'num_format': '0.00'})
    number_dec_2_format.set_align('left')
    worksheet.set_column('L:L', None, number_dec_2_format) # Estimated_Trimmed_Coverage
    worksheet.set_column('H:I', None, number_dec_2_format) # Q30%s
    worksheet.set_column('O:P', None, number_dec_2_format) # Assembly Ratio and StDev
    if phoenix == True:
        worksheet.set_column('U:V', None, number_dec_2_format) # FastANI ID and Coverage
    else:
        worksheet.set_column('W:X', None, number_dec_2_format) # FastANI ID and Coverage
    # getting values to set column widths automatically
    for idx, col in enumerate(df):  # loop through all columns
        series = df[col]
        #print(series)  # uncomment this to see what the issue is with the "mad" error
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
    worksheet.merge_range('A1:C1', "PHoeNIx Summary", cell_format_light_blue)
    worksheet.merge_range('D1:P1', "QC Metrics", cell_format_grey_blue)
    worksheet.merge_range('Q1:Y1', "Taxonomic Information", cell_format_green)
    worksheet.merge_range('Z1:AE1', "MLST Schemes", cell_format_green_blue)
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
    # Apply a conditional format for checking coverage is between 40-100 in estimated coverage column. adding 2 to max row to account for headers
    worksheet.conditional_format('L3:L' + str(max_row + 2), {'type': 'cell', 'criteria': '<', 'value':  40.00, 'format': yellow_format})
    worksheet.conditional_format('L3:L' + str(max_row + 2), {'type': 'cell', 'criteria': '>', 'value':  100.00, 'format': yellow_format})
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
    worksheet.write('A' + str(max_row + 5), 'Cells in ORANGE denote “Big 5” carbapenemase gene (i.e., blaKPC, blaNDM, blaoxa-48-like, blaVIM, and blaIMP) or an acquired blaOXA gene, please confirm what AR Lab Network HAI/AR WGS priority these meet.', orange_format_nb)
    worksheet.write('A' + str(max_row + 6), 'Cells in RED denote isolates that failed one or more auto failure triggers (cov < 30, stdev > 2.58, assembly length < 1Mbps)', red_format)
    # More footers - Disclaimer etc.
    # unbold
    no_bold = workbook.add_format({'bold': False})
    worksheet.write('A' + str(max_row + 7),"^Using Antibiotic Resistance Gene database " + ar_db + " (ResFinder, ARG-ANNOT, NCBI Bacterial Antimicrobial Resistance Reference Gene Database) using output thresholds ([98AA/90]G:[98NT/90]S); gene matches from S:(SRST2) with [%Nuc_Identity, %Coverage], or from G:(GAMMA) with [%Nuc_Identity, %AA_Identity,  %Coverage]; GAMMA gene matches indicate associated contig.", no_bold)
    worksheet.write('A' + str(max_row + 8),"^^Using CDC-compiled iroB, iucA, peg-344, rmpA, and rmpA2 hypervirulence gene database ( " + hv_db + " ); gene matches noted with [%Nuc_Identity, %AA_Identity,  %Coverage].", no_bold)
    worksheet.write('A' + str(max_row + 9),"^^^Using the plasmid incompatibility replicons plasmidFinder database ( " + pf_db + " ) using output thresholds [95NT/60]; replicon matches noted with [%Nuc_Identity, %Coverage].", no_bold)
    worksheet.write('A' + str(max_row + 10),"DISCLAIMER: These data are preliminary and subject to change. The identification methods used and the data summarized are for public health surveillance or investigational purposes only and must NOT be communicated to the patient, their care provider, or placed in the patient’s medical record. These results should NOT be used for diagnosis, treatment, or assessment of individual patient health or management.", bold)
    #adding review and date info
    worksheet.write('A' + str(max_row + 12), "Reviewed by:", no_bold)
    worksheet.write('D' + str(max_row + 12), "Date:")
    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

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
    with open("GRiPHin_samplesheet.csv", "w") as samplesheet:
        samplesheet.write('sample,directory\n')
    dirs = sorted(os.listdir(directory))
    # If there are any new files added to the top directory they will need to be added here or you will get an error
    skip_list_a = glob.glob(directory + "/*_GRiPHin_Report.xlsx") # for if griphin is run on a folder that already has a report in it
    skip_list_a = [ gene.split('/')[-1] for gene in skip_list_a ]  # just get the excel name not the full path
    skip_list_b = ["Phoenix_Output_Report.tsv", "pipeline_info", "GRiPHin_Report.xlsx", "multiqc", "samplesheet_converted.csv", "GRiPHin_samplesheet.csv", "sra_samplesheet.csv"]
    skip_list = skip_list_a + skip_list_b
    for sample in dirs:
        if sample not in skip_list:
            with open("GRiPHin_samplesheet.csv", "a") as samplesheet:
                if directory[-1] != "/": # if directory doesn't have trailing / add one
                    directory = directory + "/"
                samplesheet.write(sample + "," + directory + sample + '\n')
    samplesheet = "GRiPHin_samplesheet.csv"
    return samplesheet

def main():
    args = parseArgs()
    # create empty lists to append to later
    Data_Locations, Platforms, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Seq_bp_L, Total_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, Scaffold_Count_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
    busco_lineage_L, percent_busco_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L = ([] for i in range(33))
    ar_df = pd.DataFrame() #create empty dataframe to fill later for AR genes
    pf_df = pd.DataFrame() #create another empty dataframe to fill later for Plasmid markers
    hv_df = pd.DataFrame() #create another empty dataframe to fill later for hypervirulence genes
    srst2_ar_df = pd.DataFrame()
    # Since srst2 currently doesn't handle () in the gene names we will make a quick detour to fix this... first making a dictionary
    ar_dic = make_ar_dictionary(args.ar_db)
    # check if a directory or samplesheet was given
    if (args.samplesheet == None) and (args.directory == None): # if no directory give AND no sample sheet given exit
        sys.exit(CRED + "You MUST pass EITHER a samplesheet or a top directory of PHoeNIx output to create one.\n" + CEND)
    # If a directory is given then create a samplesheet from it if not use the samplesheet passed
    if args.directory !=None:
        samplesheet = create_samplesheet(args.directory)
    else:
        samplesheet = args.samplesheet
    #input is a samplesheet that is "samplename,directory" where the directory is a phoenix like folder
    with open(samplesheet) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        header = next(csv_reader) # skip the first line of the samplesheet
        csv_reader = sorted(csv_reader, key=operator.itemgetter(1), reverse=False) # sort first so we have an order to the output. 
        for row in csv_reader:
            sample_name = row[0]
            directory = row[1]
            data_location, platform = Get_Parent_Folder(directory, args.platform)
            trim_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld, quast_report, mlst_file, busco_short_summary, asmbld_ratio, gamma_ar_file, gamma_pf_file, gamma_hv_file, fast_ani_file, tax_file, srst2_file = Get_Files(directory, sample_name)
            #Get the metrics for the sample
            srst2_ar_df, pf_df, ar_df, hv_df, Q30_R1_per, Q30_R2_per, Total_Seq_bp, Total_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, Scaffold_Count, busco_metrics, assembly_ratio_metrics, QC_result, \
            QC_reason, MLST_scheme_1, MLST_scheme_2, MLST_type_1, MLST_type_2, MLST_alleles_1, MLST_alleles_2, MLST_source_1, MLST_source_2 = Get_Metrics(args.set_coverage, srst2_ar_df, pf_df, ar_df, hv_df, trim_stats, kraken_trim, kraken_trim_report, kraken_wtasmbld, quast_report, busco_short_summary, asmbld_ratio, sample_name, mlst_file, gamma_ar_file, gamma_pf_file, gamma_hv_file, fast_ani_file, tax_file, srst2_file, ar_dic)
            #Collect this mess of variables into appeneded lists
            Data_Locations, Platforms, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Seq_bp_L, Total_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L , alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L , MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L = Append_Lists(data_location, \
            platform, sample_name, Q30_R1_per, Q30_R2_per, Total_Seq_bp, Total_Seq_reads, Trim_kraken, Asmbld_kraken, Coverage, Assembly_Length, FastANI_output_list, warnings, alerts, Scaffold_Count, busco_metrics, \
            assembly_ratio_metrics, QC_result, QC_reason, MLST_scheme_1, MLST_scheme_2, MLST_type_1, MLST_type_2, MLST_alleles_1, MLST_alleles_2, MLST_source_1, MLST_source_2, \
            Data_Locations, Platforms, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Seq_bp_L, Total_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
            Scaffold_Count_L, busco_lineage_L, percent_busco_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L, MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L)
    # combine all lists into a dataframe
    df = Create_df(args.phoenix, Data_Locations, Platforms, Sample_Names, Q30_R1_per_L, Q30_R2_per_L, Total_Seq_bp_L, Total_Seq_reads_L, Trim_kraken_L, Asmbld_kraken_L, Coverage_L, Assembly_Length_L, Species_Support_L, fastani_organism_L, fastani_ID_L, fastani_coverage_L, warnings_L, alerts_L, \
    Scaffold_Count_L, busco_lineage_L, percent_busco_L, assembly_ratio_L, assembly_stdev_L, tax_method_L, QC_result_L, QC_reason_L, MLST_scheme_1_L, MLST_scheme_2_L, MLST_type_1_L, MLST_type_2_L, MLST_alleles_1_L , MLST_alleles_2_L, MLST_source_1_L, MLST_source_2_L)
    (qc_max_row, qc_max_col) = df.shape
    pf_max_col = pf_df.shape[1] - 1 #remove one for the WGS_ID column
    hv_max_col = hv_df.shape[1] - 1 #remove one for the WGS_ID column
    final_df, ar_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db = Combine_dfs(df, ar_df, pf_df, hv_df, srst2_ar_df, args.phoenix)
    # Checking if there was a control sheet submitted
    if args.control_list !=None:
        final_df = blind_samples(final_df, args.control_list)
    else:
        final_df = final_df
    write_to_excel(args.set_coverage, args.output, final_df, qc_max_col, ar_max_col, pf_max_col, hv_max_col, columns_to_highlight, final_ar_df, pf_db, ar_db, hv_db, args.phoenix)

if __name__ == '__main__':
    main()
