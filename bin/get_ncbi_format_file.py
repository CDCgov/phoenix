#!/usr/bin/env python3

__author__ = 'Yueli Zheng'

#disable cache usage in the Python so __pycache__ isn't formed. If you don't do this using 'nextflow run cdcgov/phoenix...' a second time will causes and error
import sys
sys.dont_write_bytecode = True
import argparse
import os
import time
import pandas as pd
import retrieve_taxo_mlst
import tools
import glob
from load_files import FileLoader

# Function to get the script version
def get_version():
    return "2.0.0"

def parseArgs(args=None):
    """"""
    parser = argparse.ArgumentParser(description='Script to generate template NCBI excel sheets for samples in PHoeNIx run')
    parser.add_argument('-i','--input', type=str, default=None, required=False, help='A GRiPHin_samplesheet. If not passed then you need to pass a directory.')
    parser.add_argument('--biosample-type', type=str, required=True)
    parser.add_argument('-o','--output', type=str, required=True)
    parser.add_argument('-m', '--microbe', dest='microbe_example', type=str, required=True)
    parser.add_argument('-s', '--sra', dest='sra_metadata', type=str, required=True)
    parser.add_argument('-g', '--griphin_summary', dest='griphin_summary', type=str, required=True)
    parser.add_argument('-b', '--osii_bioprojects', dest='osii_bioprojects', type=str, required=True)
    parser.add_argument('-d', '--directory', default=None, dest='directory', type=str, required=False, help='Pass a PHoeNIx output directory rather than a GRiPHin_samplesheet.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

def get_isolate_dirs(directory, griphin_summary):
    # If there are any new files added to the top directory they will need to be added here or you will get an error
    directories = glob.glob(directory + "/*") # for if griphin is run on a folder that already has a report in it
    ext = [".tsv", ".csv", ".xlsx" ]
    filtered_dirs = [val for val in directories if not val.endswith((".tsv", ".csv", ".xlsx"))]
    skip_list_dirs = ["pipeline_info", "multiqc"]
    # check that directories have the files needed.
    df = pd.read_csv(griphin_summary, header=0, sep='\t')
    df_fails = df[df['Minimum_QC_Check'].str.contains('FAIL')]
    failed_id_list = df_fails["WGS_ID"].tolist()
    skip_list = skip_list_dirs + failed_id_list
    # remove bad directories from list
    dirs_cleaned = [val for val in filtered_dirs if not any(str(bad_val) in str(val) for bad_val in skip_list)] #remove unwanted files
    try: #if there are numbers in the name then use that to sort
        dirs_sorted=sorted(dirs_cleaned, key=lambda x: int("".join([i for i in x if i.isdigit()])))
    except: #if no numbers then use only alphabetically
        dirs_sorted=sorted(dirs_cleaned)
    return dirs_sorted

# get biosample format file
def load_bio_projects(sample_type, isolate_list_path, microbe_example):
    isolate_names = tools.extract_string(isolate_list_path, "/")
    metafile = {}
    if "Microbe" in sample_type or "microbe" in sample_type:
        #biosample_example_path = os.path.join(os.path.dirname(__file__), microbe_example)
        df = pd.read_excel(microbe_example, header=0, sheet_name='Microbe.1.0', comment='#')
        #index = df.index[df.apply(lambda x: all(keyword in str(x.iloc[:]) for keyword in ["sample_name", "sample_title"]), axis=1)]
        #columns = list(df.iloc[index.values].values[0])
        # Remove newline characters from column names
        columns = df.columns
        for isolate in isolate_names:
            sample_id = isolate
            sample = {}
            for column_name in columns:
                sample[column_name] = ""
            sample_meta = metainfo(sample)
            metafile[sample_id] = sample_meta
    return metafile, columns

def ncbi_excel_loader(biosample_example_path, isolate_list_path):
    isolate_names = tools.extract_string(isolate_list_path, "/")
    metafile = {}
    df = pd.read_excel(biosample_example_path, header=None, sheet_name='Microbe.1.0')
    index = df.index[df.apply(lambda x: all(keyword in str(x.iloc[:]) for keyword in ["sample_name", "sample_title"]), axis=1)]
    columns = list(df.iloc[index.values].values[0])
    for isolate in isolate_names:
        sample_id = isolate
        sample = {}
        for column_name in columns:
            sample[column_name] = ""
        sample_meta = metainfo(sample)
        metafile[sample_id] = sample_meta
    return metafile

class metainfo:
    def __init__(self, samplecontent):
        self.sampleContent = samplecontent

def load_sra(ncbi_project_type, isolate_list_path, sra_metadata):
    isolate_names = tools.extract_string(isolate_list_path, "/")
    srameta = {}
    if "Microbe" in ncbi_project_type or "microbe" in ncbi_project_type:
        #sra_format_path = os.path.join(os.path.dirname(__file__), sra_metadata)
        columns = pd.read_excel(sra_metadata, header=None, sheet_name='SRA_data').iloc[0].tolist()
        for isolate in isolate_names:
            sample_id = isolate
            sample = {}
            for column_name in columns:
                sample[column_name] = ""
            sample_meta = metainfo(sample)
            srameta[sample_id] = sample_meta
    return srameta, columns

def check_project(one_isoalte_taxo, bioprojects_taxo):
    flag = False
    # check P, then G, then G + S
    if (one_isoalte_taxo['P'] + '(P)') in bioprojects_taxo.keys():
        flag = True
        return one_isoalte_taxo['P'] + '(P)'
    if (one_isoalte_taxo['G'] + '(G)') in bioprojects_taxo.keys():
        flag = True
        return one_isoalte_taxo['G'] + '(G)'
    if (one_isoalte_taxo['G'] + '(G)' + " " + one_isoalte_taxo['s'] + '(s)') in bioprojects_taxo.keys():
        flag = True
        return one_isoalte_taxo['G'] + '(G)' + " " + one_isoalte_taxo['s'] + '(s)'
    if (one_isoalte_taxo['G'] + '(G)') == "Staphylococcus(G)" and (one_isoalte_taxo['s'] + '(s)') != "aureus(s)":
        flag = True
        return one_isoalte_taxo['G'] + '(G)' + " " + "non-aureus species(s)"
    else:
        return ""

def fill_meta_values(meta_info_biosample, isolate_taxs, mlst_info_isolate, bioprojects_taxo):
    project = ""
    for isolate_name in list(isolate_taxs.keys()):
        meta_info_biosample[isolate_name].sampleContent['*sample_name'] = isolate_name
        if check_project(isolate_taxs[isolate_name], bioprojects_taxo) != "":
            project = check_project(isolate_taxs[isolate_name], bioprojects_taxo)
            meta_info_biosample[isolate_name].sampleContent['bioproject_accession'] = bioprojects_taxo[project]
        else:
            meta_info_biosample[isolate_name].sampleContent['bioproject_accession'] = ""
        meta_info_biosample[isolate_name].sampleContent['*organism'] = isolate_taxs[isolate_name]['G'] + " " + isolate_taxs[isolate_name]['s']
        meta_info_biosample[isolate_name].sampleContent['strain'] = ""
        meta_info_biosample[isolate_name].sampleContent['isolate'] = isolate_name
        meta_info_biosample[isolate_name].sampleContent['host'] = ""
        meta_info_biosample[isolate_name].sampleContent['isolation_source'] = ""
        meta_info_biosample[isolate_name].sampleContent['*collection_date'] = ""
        meta_info_biosample[isolate_name].sampleContent['*geo_loc_name'] = "USA"
        meta_info_biosample[isolate_name].sampleContent['*sample_type'] = "whole organism"
        meta_info_biosample[isolate_name].sampleContent['MLST'] = mlst_info_isolate[isolate_name]
    return meta_info_biosample

def fill_sra(sra, seq_machine_isolate):
    for key, value in sra.items():
        sra[key].sampleContent['sample_name'] = key
        sra[key].sampleContent['library_ID'] = key
        sra[key].sampleContent['title'] = "Illumina whole-genome sequencing of " + key
        sra[key].sampleContent['library_strategy'] = "WGS"
        sra[key].sampleContent['library_source'] = "GENOMIC"
        sra[key].sampleContent['library_selection'] = "RANDOM"
        sra[key].sampleContent['library_layout'] = "paired"
        sra[key].sampleContent['platform'] = "ILLUMINA"
        sra[key].sampleContent['instrument_model'] = ""
        sra[key].sampleContent['design_description'] = ""
        sra[key].sampleContent['filetype'] = "fastq"
        sra[key].sampleContent['filename'] = key + "_R1_001.fastq.gz"
        sra[key].sampleContent['filename2'] = key + "_R2_001.fastq.gz"
        sra[key].sampleContent['filename3'] = ""
        sra[key].sampleContent['filename4'] = ""
        sra[key].sampleContent['assembly'] = ""
        sra[key].sampleContent['fasta_file'] = ""
    return sra

def check_input(input):
    flag = ""
    df = pd.read_csv(input)
    first_row = df.iloc[0]
    if "sample" in df.columns and "directory" in df.columns:
        flag = "csv"
        return flag
    if "/" in first_row.values:
        flag = "txt"
        return flag
    else:
        return flag

def base_output(output_path, sra, sra_columns, bio_attribute, biosample_columns):
    path = output_path
    if not os.path.exists(path):
        # Create the directory
        os.makedirs(path)
    path_bio = path + "/BiosampleAttributes_Microbe.1.0.xlsx"
    path_sra = path + "/Sra_Microbe.1.0.xlsx"
    df_biosample = pd.DataFrame(tools.purify_dict(bio_attribute)).T.reset_index(drop=True)
    #make sure the column order is correct
    df_biosample[biosample_columns]
    df_biosample.to_excel(path_bio, index=False)
    df_sra = pd.DataFrame(tools.purify_dict(sra)).T.reset_index(drop=True)
    #make sure the column order is correct
    df_sra[sra_columns]
    df_sra.to_excel(path_sra, index=False)
    #add disclaimer
    add_disclaimer(df_biosample, path_bio, "Sheet1")
    add_disclaimer(df_sra, path_sra, "SRA_data")

def add_disclaimer(df, input_excel, input_sheet_name):
    # Load the Excel file using pandas ExcelWriter
    with pd.ExcelWriter(input_excel, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name=input_sheet_name, index=False)
        workbook  = writer.book
        worksheet = writer.sheets[input_sheet_name]
        # getting values to set column widths automatically
        for idx, col in enumerate(df):  # loop through all columns
            series = df[col]
            max_len = max((
            series.astype(str).map(len).max(),  # len of largest item
                len(str(series.name))  # len of column name/header
                )) + 1  # adding a little extra space
            worksheet.set_column(idx, idx, max_len)  # set column width
        # Add text to the cell
        # Change the font color to red
        red_format = workbook.add_format({'color': 'red', 'bold': True, 'text_wrap': True})
        # Change the font color to orange
        orange_format = workbook.add_format({'color': 'orange', 'bold': True, 'text_wrap': True})
        biosample_delete_warning = """Do the following before upload:
1. Delete this row and the rows below!
2. At minimum fill out the following columns: 
    - Host: e.g., Homo sapiens, animal, environmental, other
    - Collection Date: Specimen collection year only
    - Isolate Source: Specimen source. 
        * Homo sapiens or animal: blood, urine, etc
        * environmental: device or surface """
        sra_delete_warning = """Do the following before upload:
1. Delete this row and the rows below!
2. Fill out 'design_description' column with a short description of our library prep info and any other pertinent information. Ex: Sequenced using Nextera XT library prep kit, 2 x 250.
3. Fill out 'instrument_model' column with your illumina model type and number (if it has one). Ex: Illumina HiSeq 1500."""
        disclaimer_text = """As a reminder, please do not submit raw sequencing data to the CDC HAI-Seq BioProject (531911) that is auto populated in this sheet unless you are a state public health laboratory, a CDC partner or have been directed to do so by DHQP. The BioProject accession IDs in this file are specifically designated for domestic HAI bacterial pathogen sequencing data, \
including from the Antimicrobial Resistance Laboratory Network (AR Lab Network), state public health labs, surveillance programs, and outbreaks. For inquiries about the appropriate BioProject location for your data, please contact HAISeq@cdc.gov."""
        # Determine the number of rows already filled
        num_rows = df.shape[0] + 2
        # Add text to the cell
        if input_sheet_name == "SRA_data":
            worksheet.merge_range('A' + str(num_rows+1) + ':J' + str(num_rows+4), sra_delete_warning, orange_format)
            worksheet.merge_range('A' + str(num_rows+5) + ':J' + str(num_rows+8), disclaimer_text, red_format)
        else:
            worksheet.merge_range('A' + str(num_rows+1) + ':J' + str(num_rows+8), biosample_delete_warning, orange_format)
            worksheet.merge_range('A' + str(num_rows+9) + ':J' + str(num_rows+12), disclaimer_text, red_format)


def base_function(isolate_full_path, sample_type, output, microbe_example, sra_metadata, osii_bioprojects, directory, griphin_summary):
    bioprojects_taxo = FileLoader().load_bioproject(osii_bioprojects)
    if directory !=None:
        determined_isolate_full_path = get_isolate_dirs(directory, griphin_summary)
    else:
        determined_isolate_full_path = isolate_full_path
    biosample, biosample_columns = load_bio_projects(sample_type, determined_isolate_full_path, microbe_example)
    tax_info_isolate = retrieve_taxo_mlst.retrieve_taxo(determined_isolate_full_path)
    mlst_info_isolate = retrieve_taxo_mlst.retrieve_mlst_nonovel(determined_isolate_full_path)
    seq_machine_isolate = retrieve_taxo_mlst.retrieve_instrument_model(determined_isolate_full_path)
    fill_meta_values(biosample, tax_info_isolate, mlst_info_isolate, bioprojects_taxo)
    sra, sra_columns = load_sra(sample_type, determined_isolate_full_path, sra_metadata)
    fill_sra(sra, seq_machine_isolate)
    base_output(output, sra, sra_columns, biosample, biosample_columns)

def manage_functions(input_file, sample_type, output, microbe_example, sra_metadata, osii_bioprojects, directory, griphin_summary):
    # case 1: providing the full path of isolates in csv file
    if input_file != None: # If a griphin file is not passed then assume the files are in the directory where script is being run.
        if directory != None:
            exit("Pick Either -i or -d to pass. You can't have both.")
        else:
            if check_input(input_file) == "csv":
                isolate_full_path_griphin = FileLoader().get_full_path_griphin(input_file)
                base_function(isolate_full_path_griphin, sample_type, output, microbe_example, sra_metadata, osii_bioprojects, directory, griphin_summary)

            # case 2: providing the path of isolates in txt file
            if check_input(input_file) == "txt":
                tax_files_path = FileLoader().import_isolate_list(input_file)
                isolate_full_path = FileLoader().get_full_path(tax_files_path)
                base_function(isolate_full_path, sample_type, output, microbe_example, sra_metadata, osii_bioprojects, directory, griphin_summary)
                print("base path and isolate name file")
            if check_input(input_file) != "csv" and check_input(input_file) != "txt":
                print("The input format is not correct")
    else:
        isolate_full_path = None
        base_function(isolate_full_path, sample_type, output, microbe_example, sra_metadata, osii_bioprojects, directory, griphin_summary)

def main():
    start_time = time.time()
    args = parseArgs()
    manage_functions(args.input, args.biosample_type, args.output, args.microbe_example, args.sra_metadata, args.osii_bioprojects, args.directory, args.griphin_summary)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("time used: ", elapsed_time)

if __name__ == '__main__':
    main()