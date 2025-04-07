#!/usr/bin/env python3

#disable cache usage in the Python so __pycache__ isn't formed. If you don't do this using 'nextflow run cdcgov/phoenix...' a second time will causes and error
import sys
sys.dont_write_bytecode = True
import glob
import os
import pandas as pd
import numpy as np
import argparse
import xlsxwriter as ws
from xlsxwriter.utility import xl_rowcol_to_cell
from openpyxl import load_workbook
from itertools import chain
from GRiPHin import order_ar_gene_columns, Combine_dfs, write_to_excel, convert_excel_to_tsv, print_df

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
##Usage: >python GRiPHin.py -g1 ./samplesheet.csv -a ResGANNCBI_20220915_srst2.fasta -c control_file.csv -o output --phoenix --scaffolds
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='''Script to create new griphin excel sheet by combining two griphin summaries. The -g2 is considered the "new" file and thus when 
    samples with data in the -g1 file will have values overwritten to be the values in the -g2 file.''')
    parser.add_argument('-g1', '--old_griphin', default=None, required=False, dest='griphin_old', help='The first griphin excel file to combine.')
    parser.add_argument('-g2', '--new_griphin', default=None, required=False, dest='griphin_new', help='The second griphin excel file to combine.')
    parser.add_argument('-o', '--output', required=False, default=None, dest='output', help='Name of output file default is GRiPHin_Summary.xlsx.')
    parser.add_argument('-b', '--bldb', required=True, default=None, dest='bldb', help='Name of output file default is GRiPHin_Summary.xlsx.')
    parser.add_argument('-s', '--samplesheet', required=False, default=None, dest='samplesheet', help='samplesheet with sample,directory columns. Used to doublecheck sample names.')
    parser.add_argument('--griphin_list', required=False, default=None, type=str, dest='griphin_list', help='pass instead of -g1/-g2 when you want to combine more than 2 griphins. If you just pass --griphin_list the script assumes you have multiple griphin_summary.xlsx files in the current dir. You can also pass a csv that just has the full paths to the griphin files you want to combine.')
    parser.add_argument('--coverage', default=30, required=False, dest='set_coverage', help='The coverage cut off default is 30x.')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CYELLOW = '\033[93m'
CEND = '\033[0m'

def sort_columns_to_primary_ungrouped(df_to_sort):
    default_df = [ "UNI", "WGS_ID", "Parent_Folder", "Data_Location", "Minimum_QC_Check", "Minimum_QC_Issues", "Warnings", "Alerts", "Raw_Q30_R1_[%]", "Raw_Q30_R2_[%]", "Total_Raw_[reads]", "Paired_Trimmed_[reads]", "Total_Trimmed_[reads]", "Estimated_Trimmed_Coverage", "GC[%]", "Scaffolds", "Assembly_Length", "Assembly_Ratio", "Assembly_StDev", "Final_Taxa_ID", "Taxa_Source", "BUSCO_Lineage", "BUSCO_%Match", "Kraken_ID_Raw_Reads_%", "Kraken_ID_WtAssembly_%", "ShigaPass_Organism", "FastANI_Organism", "FastANI_%ID", "FastANI_%Coverage", "Species_Support_ANI", "Primary_MLST_Scheme", "Primary_MLST_Source", "Primary_MLST", "Primary_MLST_Alleles", "Secondary_MLST_Scheme", "Secondary_MLST_Source", "Secondary_MLST", "Secondary_MLST_Alleles", "MLST Clade", "Toxinotype", "Toxin-A_sub-type", "tcdA", "Toxin-B_sub-type", "tcdB", "tcdC_Variant", "tcdC other mutations", "tcdC", "tcdR", "tcdE", "PaLoc_NonTox_Variant", "PaLoc_NonTox other mutations", "PaLoc_NonTox", "cdtA", "cdtB", "cdtR_Variant", "cdtR other mutations", "cdtR", "cdtAB1", "cdtAB2", "gyrA known mutations", "gyrA other mutations", "gyrA", "gyrB known mutations", "gyrB other mutations", "gyrB", "dacS known mutations", "dacS other mutations", "dacS", "feoB known mutations", "feoB other mutations", "feoB", "fur known mutations", "fur other mutations", "fur", "gdpP known mutations", "gdpP other mutations", "gdpP", "glyC known mutations", "glyC other mutations", "glyC", "hemN known mutations", "hemN other mutations", "hemN", "hsmA known mutations", "hsmA other mutations", "hsmA","lscR known mutations", "lscR other mutations", "lscR", "marR known mutations", "marR other mutations", "marR", "murG known mutations", "murG other mutations", "murG", "nifJ known mutations", "nifJ other mutations", "nifJ", "PNimB known mutations", "PNimB other mutations", "PNimB", "rpoB known mutations", "rpoB other mutations", "rpoB", "rpoC known mutations", "rpoC other mutations", "rpoC", "sdaB known mutations", "sdaB other mutations", "sdaB", "thiH known mutations", "thiH other mutations", "thiH", "vanR known mutations", "vanR other mutations", "vanR", "vanS known mutations", "vanS other mutations", "vanS", "CEMB RT Crosswalk", "Inferred RT", "Probability", "ML Note", "Plasmid Info", "AR_Database", "No_AR_Genes_Found", "HV_Database", "No_HVGs_Found", "Plasmid_Replicon_Database", "No_Plasmid_Markers" ]
    sorted_list = [col for col in default_df if col in df_to_sort.columns ]
    sorted_df = df_to_sort[ sorted_list ]
    return sorted_df

def read_excels(file_path1, file_path2, samplesheet):
    #get number of footer lines for each file to skip when reading in as a pd
    footer_lines1 = detect_footer_lines(file_path1)
    footer_lines2 = detect_footer_lines(file_path2)
    # Read the Excel file, skipping the first row and using the second row as the header
    try: #check that this is an excel file
        df_1 = pd.read_excel(file_path1,
            skiprows=1,  # Skip the first header row
            header=0,    # Use the second row as the header
            skipfooter=footer_lines1,engine='openpyxl')
        df_1.insert(0, 'UNI', df_1.apply(lambda x:'%s/%s/%s' % (x['Parent_Folder'],x['Data_Location'],x['WGS_ID']),axis=1))
    except Exception as e:
        raise ValueError(f"The input file is not a valid Excel file: {file_path1}")
    try: #check that this is an excel file
        df_2 = pd.read_excel(file_path2,
            skiprows=1,  # Skip the first header row
            header=0,    # Use the second row as the header
            skipfooter=footer_lines2,engine='openpyxl')
        df_2.insert(0, 'UNI', df_2.apply(lambda x:'%s/%s/%s' % (x['Parent_Folder'],x['Data_Location'],x['WGS_ID']),axis=1))
    except Exception as e:
        raise ValueError(f"The input file is not a valid Excel file: {file_path2}")
    # because the dir is pulled into the module when making a griphin file some samples that aren't in the samplesheet can end up in the griphin file so we will remove those
    #print_df(df_1, "DF_1-precheck")
    #print_df(df_2, "DF_2-precheck")
    if samplesheet != None:
        samplesheet = pd.read_csv(samplesheet, header=0)
        if 'directory' in samplesheet.columns:
            samples = samplesheet['directory'].to_list()
        elif 'sample' in samplesheet.columns:
            samples =  samplesheet['sample'].tolist()
        df_1 = df_1[df_1['UNI'].isin(samples)]
        df_2 = df_2[df_2['UNI'].isin(samples)]
    # check that we have the files from the same entry
    phoenix, shiga, centar_1, centar_2, all_centar = check_column_presence(df_1, file_path1, df_2, file_path2)
    if centar_1 == True:
        #ordered_centar_df_1, centar_df_lens_1,centar_df_column_names_1 = split_centar_df(centar_1, file_path1)
        df1_qc, df1_gene_hits_with_centar = split_dataframe(df_1,"Toxinotype")
        df1_centar, df1_gene_hits = split_dataframe(df1_gene_hits_with_centar,"AR_Database")
        df1_ar, df1_pf_hv = split_dataframe(df1_gene_hits,"HV_Database")
        df1_hv, df1_pf = split_dataframe(df1_pf_hv,"Plasmid_Replicon_Database")
        pd.set_option('display.max_colwidth', None) 
        pd.set_option('display.max_rows', None)
        print_df(df1_qc, "DF1-QC", False)
        print_df(df1_centar, "DF1-CENTAR", False)
        print_df(df1_ar, "DF1-AR", False)
        print_df(df1_hv, "DF1-HV", False)
        print_df(df1_pf, "DF1-PF", False)
        print(df1_qc['UNI'])
        print(df1_centar['UNI'])
        print(df1_ar['UNI'])
        print(df1_hv['UNI'])
        print(df1_pf['UNI'])
    else:
        df1_qc, df1_gene_hits = split_dataframe(df_1,"AR_Database")
        df1_ar, df1_pf_hv = split_dataframe(df1_gene_hits,"HV_Database")
        df1_hv, df1_pf = split_dataframe(df1_pf_hv,"Plasmid_Replicon_Database")
    if centar_2 == True:
        #ordered_centar_df_2, centar_df_lens_2, centar_df_column_names_2 = split_centar_df(centar_2, file_path2)
        df2_qc, df2_gene_hits_with_centar = split_dataframe(df_2,"Toxinotype")
        df2_centar, df2_gene_hits = split_dataframe(df2_gene_hits_with_centar,"AR_Database")
        df2_ar, df2_pf_hv = split_dataframe(df2_gene_hits,"HV_Database")
        df2_hv, df2_pf = split_dataframe(df2_pf_hv,"Plasmid_Replicon_Database")
        centar_headlines=[['UNI', 'WGS_ID', 'Toxinotype', 'Toxin-A_sub-type', 'tcdA', 'Toxin-B_sub-type', 'tcdB'], ['tcdC_Variant', 'tcdC other mutations', 'tcdC', 'tcdR', 'tcdE', 'cdtA', 'cdtB', 'cdtR_Variant', 'cdtR other mutations', 'cdtR', 'cdtAB1', 'cdtAB2', 'PaLoc_NonTox_Variant', 'PaLoc_NonTox other mutations', 'PaLoc_NonTox'], ['gyrA known mutations', 'gyrA other mutations', 'gyrA', 'gyrB known mutations', 'gyrB other mutations', 'gyrB', 'dacS known mutations', 'dacS other mutations', 'dacS', 'feoB known mutations', 'feoB other mutations', 'feoB', 'fur known mutations', 'fur other mutations', 'fur', 'gdpP known mutations', 'gdpP other mutations', 'gdpP', 'glyC known mutations', 'glyC other mutations', 'glyC', 'hemN known mutations', 'hemN other mutations', 'hemN', 'hsmA known mutations', 'hsmA other mutations', 'hsmA', 'lscR known mutations', 'lscR other mutations', 'lscR', 'marR known mutations', 'marR other mutations', 'marR', 'murG known mutations', 'murG other mutations', 'murG', 'nifJ known mutations', 'nifJ other mutations', 'nifJ', 'PNimB known mutations', 'PNimB other mutations', 'PNimB', 'PNimB |','rpoB known mutations', 'rpoB other mutations', 'rpoB', 'rpoC known mutations', 'rpoC other mutations', 'rpoC', 'sdaB known mutations', 'sdaB other mutations', 'sdaB', 'thiH known mutations', 'thiH other mutations', 'thiH', 'vanR known mutations', 'vanR other mutations', 'vanR', 'vanS known mutations', 'vanS other mutations', 'vanS'], ['CEMB RT Crosswalk', 'Inferred RT', 'Probability', 'ML Note', 'Plasmid Info']]
        if centar_1:
            found_headlines = []
            found_headlines2 = [[],[],[],[]]
            centar_df_lens=[0,0,0,0]
            for section in range(0,len(centar_headlines)):
                for liner in centar_headlines[section]:
                    if liner in df1_centar.columns.tolist() or liner in df2_centar.columns.tolist():
                        found_headlines.append(liner)
                        if liner != 'WGS_ID' and liner != 'UNI':
                            found_headlines2[section].append(liner)
                            centar_df_lens[section]+=1
            combined_centar_df = pd.concat([df1_centar, df2_centar], axis = 0)
            combined_ordered_centar_df = combined_centar_df[found_headlines]
            centar_df_column_names = found_headlines2
        else:
            found_headlines = []
            found_headlines2 = [[],[],[],[]]
            centar_df_lens=[0,0,0,0]
            for section in range(0,len(centar_headlines)):
                for liner in centar_headlines[section]:
                    if liner in df2_centar.columns.tolist():
                        found_headlines.append(liner)
                        if liner != 'WGS_ID' and liner != 'UNI':
                            found_headlines2[section].append(liner)
                            centar_df_lens[section]+=1
            combined_ordered_centar_df = df2_centar
    else:
        df2_qc, df2_gene_hits = split_dataframe(df_2,"AR_Database")
        df2_ar, df2_pf_hv = split_dataframe(df2_gene_hits,"HV_Database")
        df2_hv, df2_pf = split_dataframe(df2_pf_hv,"Plasmid_Replicon_Database")
        if centar_1:
            for section in range(0,len(centar_headlines)):
                for liner in centar_headlines[section]:
                    if liner in df1_centar.columns.tolist():
                        found_headlines.append(liner)
                        found_headlines2[section].append(liner)
                        centar_df_lens[section]+=1
            combined_ordered_centar_df = df1_centar
        else:
            combined_ordered_centar_df=pd.DataFrame()
            centar_df_lens=[0,0,0,0]
            centar_df_column_names=[]
    #combine qc columns
    combined_df_qc = combine_qc_dataframes(df1_qc, df2_qc)
    #combine ar columns
    #print_df(df1_ar, "A1")
    #print_df(df1_hv, "A2")
    #print_df(df1_pf, "A3")
    #print_df(df2_ar, "B1")
    #print_df(df2_hv, "B2")
    #print_df(df2_pf, "B3")
    combined_df_ar, samples_to_add = combine_gene_dataframes(df1_ar, df2_ar)
    #print_df(combined_df_ar, "C1")
    #make sure the order of the ar genes is correct
    order_combined_ar_df = order_ar_gene_columns(combined_df_ar, True)
    print(CYELLOW + "\nAdding sample(s) to the GRiPHin summary:", samples_to_add.tolist(),"\n" + CEND)
    # combine pf and hv dataframes
    combined_df_pf, samples_to_add = combine_gene_dataframes(df1_pf, df2_pf)
    combined_df_hv, samples_to_add = combine_gene_dataframes(df1_hv, df2_hv)
    return combined_df_qc, order_combined_ar_df, combined_df_pf, combined_df_hv, phoenix, shiga, all_centar, combined_ordered_centar_df, centar_df_lens, centar_df_column_names
 
def combine_centar(ordered_centar_df_1, centar_df_lens_1, centar_df_column_names_1, ordered_centar_df_2, centar_df_lens_2, centar_df_column_names_2):
    # Combine the ordered_centar_df DataFrames
    reordered_centar_df = pd.concat([ordered_centar_df_1, ordered_centar_df_2], ignore_index=False, sort=False)
    # Combine the centar_df_lens lists
    centar_df_lens = [max(lens) for lens in zip(centar_df_lens_1, centar_df_lens_2)]
    #separate dataframes
    RB_type = [ "CEMB RT Crosswalk", "Inferred RT", "Probability", "ML Note", "Plasmid Info" ]
    RB_type_col = [col for col in reordered_centar_df.columns if any(substring in col for substring in RB_type) ]
    A_B_Tox = [ "Toxinotype", "Toxin-A_sub-type", "tcdA", "Toxin-B_sub-type", "tcdB"]
    A_B_Tox_col = [col for col in reordered_centar_df.columns if any(substring in col for substring in A_B_Tox) ]
    other_Tox = [ "tcdC", "tcdR", "tcdE", "cdtA", "cdtB", "cdtR", "cdtAB1", "cdtAB2", "non-tox", "PaLoc" ]
    other_Tox_col = [col for col in reordered_centar_df.columns if any(substring in col for substring in other_Tox) ]
    mutants = [ 'gyr','dac','feo','fur','gdp','gly','hem','hsm','Isc','mur', 'mur','nifJ','PNim','rpo','sda','thi','Van','mutations' ]
    mutations_col = [col for col in reordered_centar_df.columns if any(substring in col for substring in mutants) ]
    # List of mutation names to remove
    mutations_to_remove = ['tcdC other mutations', 'cdtR other mutations', 'PaLoc_NonTox other mutations']
    # Remove each mutation name if it exists in mutations_col
    mutations_col = [mutation for mutation in mutations_col if mutation not in mutations_to_remove]
    existing_columns_in_order = ['WGS_ID'] + A_B_Tox_col + other_Tox_col + mutations_col + RB_type_col
    ordered_centar_df = reordered_centar_df[existing_columns_in_order]
    centar_column_packages = []
    centar_column_packages.append(A_B_Tox_col)
    centar_column_packages.append(other_Tox_col)
    centar_column_packages.append(mutations_col)
    centar_column_packages.append(RB_type_col)
    # Combine the centar_df_column_names lists
    return ordered_centar_df, centar_df_lens, centar_column_packages

def split_centar_df(centar, excel):
    footer_lines1 = detect_footer_lines(excel)
    df = pd.read_excel(excel,
        header=[0, 1],    # Use the second row as the header
        skipfooter=footer_lines1,engine='openpyxl')
    df.insert(0, 'UNI', df.apply(lambda x:'%s/%s/%s' % (x[('PHoeNIx Summary','Parent_Folder')],x[('PHoeNIx Summary', 'Data_Location')],x[('PHoeNIx Summary', 'WGS_ID')]),axis=1))
    if centar == True:
        df_1 = pd.read_excel(excel,
            skiprows=1,  # Skip the first header row
            header=0,    # Use the second row as the header
            skipfooter=footer_lines1,engine='openpyxl')
        df_1.insert(0, 'UNI', df_1.apply(lambda x:'%s/%s/%s' % (x['Parent_Folder'],x['Data_Location'],x['WGS_ID']),axis=1))
        # Get the column names
        columns = df_1.columns.tolist()
        # Find the indices of the specified columns
        ## Might need to fiddle with MLST Clade here eventually...might not if original griphins get made right
        start_index = columns.index('Toxinotype')
        end_index = columns.index('AR_Database')
        # Get the columns between the specified columns
        columns_between = [columns[0]] + columns[start_index:end_index]
        # Subset the DataFrame
        ordered_centar_df = df_1[columns_between]
        ordered_centar_df = sort_columns_to_primary_ungrouped(ordered_centar_df)
        #get length of sub dataframes
        centar_df_lens = [ len(df['Toxin A/B Variants'].columns), len(df['Other Toxins'].columns), len(df['C. difficile Specific AR Mutations'].columns), len(df['ML Predicted Ribotype'].columns) ] # have to have it in a list for sum() later
        centar_df_column_names = [ df['Toxin A/B Variants'].columns, df['Other Toxins'].columns, df['C. difficile Specific AR Mutations'].columns, df['ML Predicted Ribotype'].columns ] # have to have it in a list for sum() later
    else:
        centar_df_lens = [0,0,0,0]
        centar_df_column_names = [[],[],[],[]]
        ordered_centar_df = pd.DataFrame()
    #print_df(ordered_centar_df, "SCDF", False)
    return ordered_centar_df, centar_df_lens, centar_df_column_names

def combine_gene_dataframes(old_df, new_df):
    # Ensure the first column is 'WGS_ID'
    #if new_df.columns[0] != 'WGS_ID' or old_df.columns[0] != 'WGS_ID':
    if new_df.columns[0] != 'UNI' and old_df.columns[0] != 'UNI':
        raise ValueError("The first column in both dataframes must be either 'WGS_ID' or 'UNI'")
    df1_gene = old_df.set_index('UNI')
    df2_gene = new_df.set_index('UNI')
    # Identify samples to add and print them
    samples_to_add = df2_gene.index.difference(df1_gene.index)
    # Combine dataframes, prioritizing new_df values and aligning columns
    combined_df = df1_gene.combine_first(df2_gene)
    ###!print_df(combined_df, "E3")
    # Add missing columns from new_df to old_df
    columns_to_add = [col for col in df2_gene.columns if col not in df1_gene.columns and col not in ["AR_Database", "HV_Database","Plasmid_Replicon_Database", "No_AR_Genes_Found", "No_HVGs_Found", "No_Plasmid_Markers", "WGS_ID", "UNI"]]
    for col in columns_to_add:
        combined_df[col] = combined_df[col].fillna(df2_gene[col])
    # Revert index to column
    combined_df.reset_index(inplace=True)
    return combined_df, samples_to_add

def combine_qc_dataframes(df1_qc, df2_qc):
    # Convert the 'WGS_ID' columns to the index to facilitate updating
    df1 = df1_qc.set_index('UNI')
    df2 = df2_qc.set_index('UNI')
    # Update df1_qc with matching rows from df2_qc
    df1.update(df2)
    # Append non-matching rows from df2_qc to df1_qc using pd.concat
    combined_df = pd.concat([df1, df2[~df2.index.isin(df1.index)]])
    # Reset the index to restore the 'UNI' column
    combined_df.reset_index(inplace=True)
    return combined_df

def add_blank_centar_columns(df, reference_df):
    # Get the list of column names from the reference dataframe
    columns = reference_df.columns.tolist()
    # Find the indices of the "Toxinotype" and "AR_Database" columns
    try:
        start_index = columns.index("Toxinotype")
        end_index = columns.index("AR_Database")
    except ValueError:
        raise ValueError("One of the required columns ('Toxinotype' or 'AR_Database') is missing in the reference dataframe.")

    # Get all the columns between "Toxinotype" (inclusive) and "AR_Database" (exclusive)
    centar_columns = columns[start_index:end_index]
    # Add blank columns for each centar-related column
    for col in centar_columns:
        if col not in df.columns:
            df[col] = pd.NA  # Fill with NaNs or blank
    return df

def check_column_presence(df1, df1_path, df2, df2_path):
    df1_has_column = "BUSCO_Lineage" in df1.columns
    df2_has_column = "BUSCO_Lineage" in df2.columns
    if df1_has_column and df2_has_column: # set that its a CDC PHOENIX run
        phoenix = False
    elif not df1_has_column and not df2_has_column: # set that its a PHOENIX run
        phoenix = True
    else:
        if df1_has_column:
            raise ValueError(f"{CRED}The old griphin file was produced from -entry CDC_PHOENIX and the new griphin summary wasn't. These files aren't compatible.{CEND}")
        else:
            raise ValueError(f"{CRED}The new griphin file was produced from -entry CDC_PHOENIX and the old griphin summary wasn't. These files aren't compatible.{CEND}")
    #check for centar
    if "Toxin-A_sub-type" in df1.columns or "Toxin-A_sub-type" in df2.columns:
        if "Toxin-A_sub-type" not in df1.columns:
            print("Adding centar blank columns to " + df1_path)
            df1 = add_blank_centar_columns(df1, df2)
            centar_1 = False
        else:
            centar_1 = True
        if "Toxin-A_sub-type" not in df2.columns:
            print("Adding centar blank columns to " + df2_path)
            df2 = add_blank_centar_columns(df2, df1)
            centar_2 = False
        else:
            centar_2 = True
    else:
        centar_1 = centar_2 = False
    if centar_1 == True or centar_2 == True:
        all_centar = True
    else:
        all_centar = False
    #check for shigapass
    if "ShigaPass_Organism" in df1.columns or "ShigaPass_Organism" in df2.columns:
        shiga = True
    else:
        shiga = False
    return phoenix, shiga, centar_1, centar_2, all_centar

def split_dataframe(df, split_column):
    # Ensure the first column is 'UNI'
    if df.columns[0] != 'UNI':
         raise ValueError("The first column must be 'UNI'")
    # Find the index of the split_column
    split_index = df.columns.get_loc(split_column)
    # Create two DataFrames based on the split index, including 'UNI' at the beginning FOR BOTH
    df_before = df.iloc[:, [0] + list(range(1, split_index))]  # All columns before the split column + 'WGS_ID'
    #df_before = df.iloc[:, list(range(1, split_index))]  # All columns before the split column + 'WGS_ID'
    df_after = df.iloc[:, [0] + list(range(split_index, df.shape[1]))]  # The split column and all columns after it + 'WGS_ID'
    return df_before, df_after

##### May need to revisit this, VERY different from merge version
def add_and_concatenate(df1, df2):
    # Step 1: Preserve column order from df1 and append missing columns from df2
    all_columns = list(df1.columns) + [col for col in df2.columns if col not in df1.columns]
    # Step 1: Find the union of columns in both DataFrames
    #all_columns = set(df1.columns).union(set(df2.columns))
    # Step 2: Add missing columns to each DataFrame at once by reindexing
    df1 = df1.reindex(columns=all_columns)
    df2 = df2.reindex(columns=all_columns)
    # Step 3: Concatenate the two DataFrames row-wise
    result = pd.concat([df1, df2], ignore_index=True, sort=False)
    # Step 4: Reset the index and return the final DataFrame
    result.reset_index(drop=True, inplace=True)
    return result

def detect_footer_lines(file_path, sheet_name=0):
    # Load the workbook
    wb = load_workbook(file_path, read_only=True, data_only=True)
    ws = wb[sheet_name if isinstance(sheet_name, str) else wb.sheetnames[sheet_name]]
    footer_lines = 3 #set to 3 because review by row has data in another cell and there are two all blank spacing rows
    is_footer_section = False
    # Count footer lines
    for row in reversed(tuple(ws.rows)):
        first_cell = row[0].value
        other_cells = [cell.value for cell in row[1:]]
        # Check if the first cell has text and all other cells are empty
        if isinstance(first_cell, str) and all(cell is None for cell in other_cells):
            footer_lines += 1
            is_footer_section = True
        elif is_footer_section:
            break
    return footer_lines

def read_excel(file_path, old_phoenix, reference_qc_df, reference_centar_df, samplesheet):
    #get number of footer lines for each file to skip when reading in as a pd
    pd.set_option('display.max_colwidth', None)
    footer_lines = detect_footer_lines(file_path)
    # Read the Excel file, skipping the first row and using the second row as the header
    try: #check that this is an excel file
        df = pd.read_excel(file_path,
            skiprows=1,  # Skip the first header row
            header=0,    # Use the second row as the header
            skipfooter=footer_lines,engine='openpyxl')
        df.insert(0, 'UNI', df.apply(lambda x:'%s/%s/%s' % (x['Parent_Folder'],x['Data_Location'],x['WGS_ID']),axis=1))
    except Exception as e:
        raise ValueError(f"The input file is not a valid Excel file: {file_path}")
    # check that we have the files from the same entry
    if "BUSCO_Lineage" in df.columns: # set that its a CDC PHOENIX run
        phoenix = False
    else:
        phoenix = False
    if old_phoenix != phoenix: 
        raise ValueError(f"{CRED}The one griphin file in your set was produced from -entry CDC_PHOENIX and some were not wasn't. These files aren't compatible to combine. {CEND}")
    #check for centar
    if "Toxin-A_sub-type" not in df.columns:
        print("Adding centar blank columns to " + file_path)
        for col in reference_centar_df.columns:
            if col not in df.columns:
                df[col] = pd.NA  # Fill with NaNs or blank
        centar = False
    else:
        centar = True
    #check for shigapass
    if "ShigaPass_Organism" in df.columns:
        shiga = True
    else:
        shiga = False
    #parse old griphin df
    if samplesheet != None:
        samplesheet = pd.read_csv(samplesheet, header=0)
        if 'directory' in samplesheet.columns:
            samples = samplesheet['directory'].to_list()
        elif 'sample' in samplesheet.columns:
            samples =  samplesheet['sample'].tolist()
        #filters the DataFrame df_1 to include only the rows where the values in the 'WGS_ID' column are present in the samples list.
        df = df[df['UNI'].isin(samples)]
    # check that we have the files from the same entry
    ordered_centar_df, centar_df_lens, centar_df_column_names = split_centar_df(centar, file_path)
    if centar == True:
        df_qc, df_gene_hits_with_centar = split_dataframe(df,"Toxinotype")
        df_centar, df_gene_hits = split_dataframe(df_gene_hits_with_centar,"AR_Database")
        df_ar, df_pf_hv = split_dataframe(df_gene_hits,"HV_Database")
        df_hv, df_pf = split_dataframe(df_pf_hv,"Plasmid_Replicon_Database")
        # set standard order they need to be in
        centar_headlines=[['UNI', 'WGS_ID', 'Toxinotype', 'Toxin-A_sub-type', 'tcdA', 'Toxin-B_sub-type', 'tcdB'], ['tcdC_Variant', 'tcdC other mutations', 'tcdC', 'tcdR', 'tcdE', 'cdtA', 'cdtB', 'cdtR_Variant', 'cdtR other mutations', 'cdtR', 'cdtAB1', 'cdtAB2', 'PaLoc_NonTox_Variant', 'PaLoc_NonTox other mutations', 'PaLoc_NonTox'], ['gyrA known mutations', 'gyrA other mutations', 'gyrA', 'gyrB known mutations', 'gyrB other mutations', 'gyrB', 'dacS known mutations', 'dacS other mutations', 'dacS', 'feoB known mutations', 'feoB other mutations', 'feoB', 'fur known mutations', 'fur other mutations', 'fur', 'gdpP known mutations', 'gdpP other mutations', 'gdpP', 'glyC known mutations', 'glyC other mutations', 'glyC', 'hemN known mutations', 'hemN other mutations', 'hemN', 'hsmA known mutations', 'hsmA other mutations', 'hsmA', 'lscR known mutations', 'lscR other mutations', 'lscR', 'marR known mutations', 'marR other mutations', 'marR', 'murG known mutations', 'murG other mutations', 'murG', 'nifJ known mutations', 'nifJ other mutations', 'nifJ', 'PNimB known mutations', 'PNimB other mutations', 'PNimB', 'PNimB |','rpoB known mutations', 'rpoB other mutations', 'rpoB', 'rpoC known mutations', 'rpoC other mutations', 'rpoC', 'sdaB known mutations', 'sdaB other mutations', 'sdaB', 'thiH known mutations', 'thiH other mutations', 'thiH', 'vanR known mutations', 'vanR other mutations', 'vanR', 'vanS known mutations', 'vanS other mutations', 'vanS'], ['CEMB RT Crosswalk', 'Inferred RT', 'Probability', 'ML Note', 'Plasmid Info']]
        #double length of dataframes for centar_df_lens to ensure merging columns for excelsheet works correctly - sometimes UNI is there sometimes WGS_ID
        found_headlines = []
        found_headlines2 = [[],[],[],[]]
        centar_df_lens=[0,0,0,0]
        for section in range(0,len(centar_headlines)):
            for header in centar_headlines[section]:
                if header in df_centar.columns.tolist():
                    found_headlines.append(header)
                    if header != 'WGS_ID' and header != 'UNI':
                        found_headlines2[section].append(header)
                        centar_df_lens[section]+=1
            ordered_centar_df = df_centar[found_headlines]
            centar_df_column_names = found_headlines2
    else:
        df_qc, df_gene_hits = split_dataframe(df,"AR_Database")
        df_ar, df_pf_hv = split_dataframe(df_gene_hits,"HV_Database")
        df_hv, df_pf = split_dataframe(df_pf_hv,"Plasmid_Replicon_Database")
    # Ensure the first column is 'UNI'
    if (df_qc.columns[0] != 'UNI' and reference_qc_df.columns[0] != 'UNI'):
       raise ValueError("The first column in both dataframes must be 'UNI'")
    # Set WGS_ID as index
    # Set WGS_ID as index
    df1_qc = reference_qc_df.set_index('UNI')
    ref_qc = df_qc.set_index('UNI')
    # Identify samples to add and print them
    samples_to_add = ref_qc.index.difference(df1_qc.index)
    #make sure the order of the ar genes is correct
    order_ar_df = order_ar_gene_columns(df_ar, True)
    return df_qc, order_ar_df, df_pf, df_hv, phoenix, shiga, centar, ordered_centar_df, centar_df_lens, centar_df_column_names

def update_centar_columns(centar_df_column_names_final, centar_df_column_names):
    # Iterate over each corresponding pair of lists in centar_df_lens_final and centar_df_lens
    for i, (final_list, new_list) in enumerate(zip(centar_df_column_names_final, centar_df_column_names)):
        # Convert both lists to sets to get unique column names
        unique_columns = set(final_list).union(set(new_list))
        # Update centar_df_lens_final with the unique columns
        centar_df_column_names_final[i] = list(unique_columns)
    # After processing all lists, calculate the total number of unique names in each list
    total_unique_columns = [len(columns) for columns in centar_df_column_names_final]
    return total_unique_columns

def main():
    # looping through excel files
    args = parseArgs()
    #figure out the name of the output file 
    if args.output != None:
        output_file = args.output
    else:
        # Derive output file name from input file name
        output_file = args.griphin_new.replace("_GRiPHin_Summary.xlsx", "")
    # checking what the input type is
    if args.griphin_old != None and args.griphin_new != None:
        combined_df_qc_final, combined_df_ar_final, combined_df_pf_final, combined_df_hv_final, phoenix_final, shiga_final, centar_final, ordered_centar_df_final, centar_df_lens_final, centar_df_column_names_final = read_excels(args.griphin_new, args.griphin_old, args.samplesheet)
    else:
        if args.griphin_list == None:
            griphin_files = glob.glob("*_GRiPHin_Summary.xlsx")
            if len(griphin_files) < 2:
                raise ValueError(f"{CRED}Need at least two GRiPHin files for combination when using --griphin_list.{CEND}")
        else:
            # A file was provided, read it and get the paths from the file
            griphin_files = []
            with open(args.griphin_list, 'r') as file:
                for line in file:
                    line = line.strip()
                    if line:  # Skip empty lines
                        if os.path.isfile(line):  # Ensure the path exists
                            griphin_files.append(line)
                        else:
                            print(f"Warning: File '{line}' does not exist.")
        # add centar_df to the almost_final_ar_df
        # Initialize the combination with the first file
        base_file = griphin_files.pop(0)
        #combine first two files
        combined_df_qc_final, combined_df_ar_final, combined_df_pf_final, combined_df_hv_final, phoenix_final, shiga_final, centar_final, ordered_centar_df_final, centar_df_lens_final, centar_df_column_names_final = read_excels(base_file, griphin_files[0], args.samplesheet)
        combined_dataframes_final = [ combined_df_qc_final, combined_df_ar_final, combined_df_pf_final, combined_df_hv_final, ordered_centar_df_final ]
        # Iterate over remaining files, progressively combining them with the base
        for count, next_file in enumerate(griphin_files[1:], start=1):
            # checking
            combined_df_qc_next, combined_df_ar_next, combined_df_pf_next, combined_df_hv_next, phoenix, shiga, centar, ordered_centar_df_next, centar_df_lens_next, centar_df_column_names_next = read_excel(next_file, phoenix_final, combined_df_qc_final, ordered_centar_df_final, args.samplesheet)
            combined_dataframes_next = [ combined_df_qc_next, combined_df_ar_next, combined_df_pf_next, combined_df_hv_next, ordered_centar_df_next ]
            # Update flags
            phoenix_final = phoenix_final or phoenix
            shiga_final = shiga_final or shiga
            centar_final = centar_final or centar
            combined_dataframes_final = [ add_and_concatenate(df1, df2) for df1, df2 in zip(combined_dataframes_final, combined_dataframes_next) ]
            if centar_final == True:
                # Update centar_df_lens1 with the maximum of each corresponding number
                centar_df_lens_final = update_centar_columns(centar_df_column_names_final, centar_df_column_names_next)
            # Unpack the tuple into individual DataFrames
            combined_df_qc_final, combined_df_ar_final, combined_df_pf_final, combined_df_hv_final, ordered_centar_df_final = combined_dataframes_final
    # add centar_df to the almost_final_ar_df
    if centar_final == True:
        # Reset the indices of both DataFrames during concat so they are aligned and we don't get NaNs in the final row of the dataframe
        ordered_centar_df_final = sort_columns_to_primary_ungrouped(ordered_centar_df_final)
        #combined_df_qc = pd.concat([combined_df_qc, ordered_centar_df], axis=1)
        combined_df_qc_final = pd.merge(combined_df_qc_final, ordered_centar_df_final, how="left", on = ['UNI', 'UNI'])
    # call function from griphin script to combine all dfs
    final_df, ar_max_col, columns_to_highlight, final_ar_df, final_pf_db, final_ar_db, final_hv_db = Combine_dfs(combined_df_qc_final, combined_df_ar_final, combined_df_pf_final, combined_df_hv_final, pd.DataFrame(), phoenix_final, True, args.bldb)
    #get other information for excel writing
    combined_df_qc_final = combined_df_qc_final.drop('UNI', axis = 1)
    (qc_max_row, qc_max_col) = combined_df_qc_final.shape

    pf_max_col = combined_df_pf_final.shape[1] - 1 #remove one for the UNI column
    hv_max_col = combined_df_hv_final.shape[1] - 1 #remove one for the UNI column
    #write excel sheet
    final_df = final_df.drop('UNI', axis=1)
    write_to_excel(args.set_coverage, output_file, final_df, qc_max_col, ar_max_col, pf_max_col, hv_max_col, columns_to_highlight, final_ar_df, final_pf_db, final_ar_db, final_hv_db, phoenix_final, shiga_final, centar_final, centar_df_lens_final)
    #write tsv from excel
    convert_excel_to_tsv(output_file)

if __name__ == '__main__':
    main()