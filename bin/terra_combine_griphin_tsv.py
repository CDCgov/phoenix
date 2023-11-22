#!/usr/bin/env python3

# importing the required modules
import glob
import pandas as pd
import argparse
import re
from re import search
from itertools import chain

##Makes a summary tsv file when given a series of griphin tsv files
##Usage: >python Terra_combine_griphin_tsv.py -o Output_Report.tsv
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "v1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a combined GRiPHin summary excel sheet')
    parser.add_argument('-o', '--out', dest='output_file', required=False, default=None, help='output file name')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    parser.add_argument('files', nargs=argparse.REMAINDER)
    return parser.parse_args()
 
def combine_tsvs(file_list):
    # create a new dataframe to store the merged excel file.
    excl_merged = pd.DataFrame()
    count = 1
    # pd.read_excel(file_path) reads the excel data into pandas dataframe.
    for file in file_list:
        #read in file and skip header
        df = pd.read_csv(file, header=0, sep='\t')
        #drop rows that have all NA
        df = df.dropna(axis=0,how='all')
        if count == 1: # for the first file
            #add dataframe together
            excl_merged = pd.concat([excl_merged, df])
        else: #next files
            sorted_cols = separate_column_type(excl_merged, df)
            excl_merged = pd.concat([excl_merged, df], ignore_index=True)
            # reorder df cols
            excl_merged = excl_merged.reindex(sorted_cols, axis=1)
        count = count + 1
    return excl_merged

def separate_column_type(excl_merged, df):
    #separate out dataframes for each type
    # get qc/taxa number of columns
    qc_df = df.loc[:,'WGS_ID':'Secondary_MLST_Alleles']
    qc_col_list = list(qc_df.columns)
    # get ar number of columns - df
    ar_df = df.loc[:,'AR_Database':'HV_Database']
    ar_df = ar_df.drop(['HV_Database'], axis=1)
    # get ar number of columns - excl_merged
    ar_df_merged = excl_merged.loc[:,'AR_Database':'HV_Database']
    ar_df_merged = ar_df_merged.drop(['HV_Database'], axis=1)
    # get all unique columns
    col_list = set(list(ar_df.columns) + list(ar_df_merged.columns))
    ar_col_list = fix_ar_col_order(col_list)
    #make sure correct columns are at the beginning
    ar_col_list.insert(0, "AR_Database")
    ar_col_list.insert(1, "No_AR_Genes_Found")
    # get hv number of cols - df
    hv_df = df.loc[:,'HV_Database':'Plasmid_Replicon_Database']
    hv_df = hv_df.drop(['Plasmid_Replicon_Database'], axis=1)
    # get hv number of columns - excl_merged
    hv_df_merged = excl_merged.loc[:,'HV_Database':'Plasmid_Replicon_Database']
    hv_df_merged = hv_df_merged.drop(['Plasmid_Replicon_Database'], axis=1)
    # get and order HV columns
    hv_col_list = sorted(set(list(hv_df.columns) + list(hv_df_merged.columns)))
    #make sure the correct columns are at the beginning
    hv_col_list.insert(0, hv_col_list.pop(hv_col_list.index("HV_Database")))
    hv_col_list.insert(1, hv_col_list.pop(hv_col_list.index("No_HVGs_Found")))
    # get pf number of cols - df
    pf_df = df.loc[:,'Plasmid_Replicon_Database':]
    # get hv number of columns - excl_merged
    pf_df_merged = excl_merged.loc[:,'Plasmid_Replicon_Database':]
    # get and order pf columns
    pf_col_list = sorted(set(list(pf_df.columns) + list(pf_df_merged.columns)))
    pf_col_list.insert(0, pf_col_list.pop(pf_col_list.index("Plasmid_Replicon_Database")))
    pf_col_list.insert(1, pf_col_list.pop(pf_col_list.index("No_Plasmid_Markers")))
    # add all columns together
    sorted_col_list = qc_col_list + ar_col_list + hv_col_list + pf_col_list
    return sorted_col_list

def fix_ar_col_order(col_list):
    ar_drugs_list = [re.findall('.*\((.*)\).*', col) for col in col_list]
    ar_drugs_list = sorted(set(list(chain.from_iterable(ar_drugs_list))))
    final_ar_list = []
    # loop over each gene with the same drug its name
    for drug in ar_drugs_list:
        drug = "(" + drug + ")"
        ar_column_list = sorted([col for col in col_list if drug in col]) # get column names filtered for each drug name
        final_ar_list.append(ar_column_list)
    # un-nest list
    final_ar_list = list(chain.from_iterable(final_ar_list))
    return final_ar_list

def get_variables(file_list):
    #create list to check that file types are all the same (CDC vs NOT)
    check_list = []
    coverage_list = []
    for file in file_list:
        df = pd.read_csv(file, header=0, sep='\t')
        df = df.dropna(axis=0,how='all')
        if "BUSCO_Lineage" in df.columns:
            phoenix = False
            check_list.append(phoenix)
        else:
            phoenix = True
            check_list.append(phoenix)
        #get coverage information
        try: #if there are no failures due to coverage then the coverage is unknown
            qc_col = df['Minimum_QC_Issues']
            qc_col = qc_col.dropna(axis=0,how='all')
            cell_val = qc_col[qc_col.str.contains("coverage <")]
            coverage = re.findall(r'coverage <\d+x', cell_val[0])
            coverage_list.append(coverage[0])
            #check that coverage is the same in all files
            if len(set(coverage_list)) != 1:
                print("Error: There are different coverage cut offs in files!")
                exit()
        except AttributeError:
            print("Warning: the coverage used to run {} could not be determined!".format(file))
    #check the values are all the same
    if len(set(check_list)) != 1:
        print("Error: Files are a mix of CDC and Not CDC versions of PHX!")
        exit()
    return phoenix

def write_combined_tsv(df, output):
    if output != None:
        output_file = output + '_GRiPHin_Summary.tsv'
    else:
        output_file = 'GRiPHin_Summary.tsv'
    #Write dataframe into csv
    df.to_csv(output_file, sep='\t', index=False, line_terminator='\n')

def main():
    args = parseArgs()
    # get files in the path
    file_list = glob.glob("*GRiPHin_*_Summary.tsv")
    # check that the file_list isn't empty
    if len(file_list) == 0:
        print("Error: No GRiPHin_Summary.tsv files were found using *GRiPHin_Summary.tsv!")
        exit()
    phoenix = get_variables(file_list)
    df = combine_tsvs(file_list)
    write_combined_tsv(df, args.output_file)

if __name__ == '__main__':
    main()