#!/usr/bin/env python3

# importing the required modules
import glob
import pandas as pd
import argparse
import xlsxwriter as ws
from xlsxwriter.utility import xl_rowcol_to_cell
import re
from re import search
from itertools import chain

##Makes a summary Excel file when given a series of griphin xlsx files
##Usage: >python terra_combine_griphin.py -o Output_Report.xlsx
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a combined GRiPHin summary excel sheet')
    parser.add_argument('-o', '--out', dest='output_file', required=False, default=None, help='output file name')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    parser.add_argument('files', nargs=argparse.REMAINDER)
    return parser.parse_args()
 
def combine_excels(file_list):
    # create a new dataframe to store the merged excel file.
    excl_merged = pd.DataFrame()
    count = 1
    # pd.read_excel(file_path) reads the excel data into pandas dataframe.
    for file in file_list:
        #read in file and skip header
        df = pd.read_excel(file, skiprows=0, header=1)
        #drop rows that have all NA
        df = df.dropna(axis=0,how='all')
        if count == 1: # for the first file
            #remove footer info
            df.drop(df.tail(8).index, inplace = True)
            #add dataframe together
            excl_merged = pd.concat([excl_merged, df])
        else: #next files
            sorted_cols = separate_column_type(excl_merged, df)
            #remove footer info
            df.drop(df.tail(8).index, inplace = True)
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

def write_excel(output_file, df, set_coverage, phoenix, qc_max_col, ar_gene_count, pf_gene_count, hv_gene_count, columns_to_highlight, ar_df, pf_db, ar_db, hv_db):
    # exports the dataframe into excel file with specified name.
    # Create a Pandas Excel writer using XlsxWriter as the engine.
    if output_file != None and output_file != "GRiPHin_Summary.xlsx":#check that its not "GRiPHin_Summary.xlsx"
        writer = pd.ExcelWriter((output_file + '_GRiPHin_Summary.xlsx'), engine='xlsxwriter')
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
    # set formating for python 3.7.12
    worksheet.conditional_format('O3:P' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_comma_format})
    worksheet.conditional_format('J3:L' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_comma_format})
    # Setting columns to float so its more human readable
    number_dec_2_format = workbook.add_format({'num_format': '0.00'})
    # set formating for python 3.7.12
    worksheet.conditional_format('M3:N' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
    worksheet.conditional_format('H3:I' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
    worksheet.conditional_format('Q3:R' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
    if phoenix == True:
        worksheet.conditional_format('W3:X' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
    else:
        worksheet.conditional_format('Y3:Z' + str(max_row + 2), {'type': 'cell', 'criteria': 'not equal to', 'value': '"Unknown"', 'format': number_dec_2_format})
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
    worksheet.merge_range('A1:C1', "PHoeNIx Summary", cell_format_light_blue)
    worksheet.merge_range('D1:R1', "QC Metrics", cell_format_grey_blue)
    if phoenix == True: #for non-CDC entry points
        worksheet.merge_range('S1:Y1', "Taxonomic Information", cell_format_green)
    else:
        worksheet.merge_range('S1:AA1', "Taxonomic Information", cell_format_green)
    if phoenix == True: #for non-CDC entry points
        worksheet.merge_range('Z1:AG1', "MLST Schemes", cell_format_green_blue)
    else:
        worksheet.merge_range('AB1:AI1', "MLST Schemes", cell_format_green_blue)
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
                worksheet.write(cell, column, orange_format)
        column_count = column_count + 1
    # Creating footers
    worksheet.write('A' + str(max_row + 4), 'Cells in YELLOW denote isolates outside of ' + str(set_coverage) + '-100X coverage', yellow_format)
    worksheet.write('A' + str(max_row + 5), 'Cells in ORANGE denote “Big 5” carbapenemase gene (i.e., blaKPC, blaNDM, blaOXA-48-like, blaVIM, and blaIMP) or an acquired blaOXA gene, please confirm what AR Lab Network HAI/AR WGS priority these meet.', orange_format_nb)
    worksheet.write('A' + str(max_row + 6), 'Cells in RED denote isolates that failed one or more auto failure triggers (cov < 30, assembly ratio stdev > 2.58, assembly length < 1Mbps)', red_format)
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
    # add autofilter
    worksheet.autofilter(1, 0, max_row, max_col - 1)
    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

def big5_check(final_ar_df):
    """"Function that will return list of columns to highlight if a sample has a hit for a big 5 gene."""
    columns_to_highlight = []
    final_ar_df = final_ar_df.drop(['AR_Database'], axis=1)
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
                [ columns_to_highlight.append(gene_name + "_(" + drug) for big5_keep_gene in blaOXA_48_like if gene_name == big5_keep_gene ]
            else: # for "blaIMP", "blaVIM", "blaNDM", and "blaKPC", this will take any thing with a matching substring to these
                for big5 in big5_keep:
                    if search(big5, gene_name): #search for big5 gene substring in the gene name
                        columns_to_highlight.append(gene_name + "_(" + drug)
    #loop through list of genes to drop and removed if they are in the highlight list
    for bad_gene in big5_drop:
        #search for big5 gene substring in the gene name and remove if it is
        [columns_to_highlight.remove(gene) for gene in columns_to_highlight if bad_gene in gene]
    return columns_to_highlight

def get_variables(file_list):
    #create list to check that file types are all the same (CDC vs NOT)
    check_list = []
    coverage_list = []
    for file in file_list:
        df = pd.read_excel(file, header=1)
        df = df.dropna(axis=0,how='all')
        if "BUSCO_Lineage" in df.columns:
            phoenix = False
            check_list.append(phoenix)
        else:
            phoenix = True
            check_list.append(phoenix)
        #get coverage information
        cell_val = df[df['WGS_ID'].str.contains("Cells in YELLOW")].iloc[0]['WGS_ID']
        coverage = re.findall(r'\d+-100X', cell_val)[0]
        coverage = re.sub('-100X', "" , str(coverage))
        coverage_list.append(coverage)
        #check that coverage is the same in all files
        if len(set(coverage_list)) != 1:
            print("Error: There are different coverage cut offs in files!")
            exit()
    #check the values are all the same
    if len(set(check_list)) != 1:
        print("Error: Files are a mix of CDC and Not CDC versions of PHX!")
        exit()
    return coverage, phoenix

def get_column_counts(df):
    # get qc number of columns
    qc_df = df.loc[:,'WGS_ID':'Secondary_MLST_Alleles']
    qc_max_col = int(len(qc_df.columns))
    # get ar number of columns
    ar_df = df.loc[:,'AR_Database':'HV_Database']
    ar_df = ar_df.drop(['HV_Database'], axis=1)
    ar_max_col = int(len(ar_df.columns))
    ar_db = list(ar_df['AR_Database'].unique())
    try:
        ar_db.remove("GAMMA file not found")
    except ValueError:
        pass
    ar_db = ','.join(ar_db)
    # get hv number of columns
    hv_df = df.loc[:,'HV_Database':'Plasmid_Replicon_Database']
    hv_df = hv_df.drop(['Plasmid_Replicon_Database'], axis=1)
    hv_max_col = int(len(hv_df.columns))
    hv_db = list(hv_df['HV_Database'].unique())
    try:
        hv_db.remove("GAMMA file not found")
    except ValueError:
        pass
    hv_db = ','.join(hv_db)
    # get hv number of columns
    pf_df = df.loc[:,'Plasmid_Replicon_Database':]
    pf_max_col = int(len(pf_df.columns))
    pf_db = list(pf_df['Plasmid_Replicon_Database'].unique())
    try:
        pf_db.remove("GAMMA file not found")
    except ValueError:
        pass
    pf_db = ','.join(pf_db)
    return qc_max_col, ar_max_col, pf_max_col, hv_max_col, ar_df, pf_db, ar_db, hv_db

def main():
    args = parseArgs()
    # get files in the path
    file_list = glob.glob("*_Summary.xlsx")
    # check that the file_list isn't empty
    if len(file_list) == 0:
        print("Error: No GRiPHin_Summary.xlsx files were found using *_Summary.xlsx!")
        exit()
    set_coverage, phoenix = get_variables(file_list)
    df = combine_excels(file_list)
    qc_max_col, ar_gene_count, pf_gene_count, hv_gene_count, ar_df, pf_db, ar_db, hv_db = get_column_counts(df)
    columns_to_highlight = big5_check(ar_df)
    write_excel(args.output_file, df, set_coverage, phoenix, qc_max_col, ar_gene_count, pf_gene_count, hv_gene_count, columns_to_highlight, ar_df, pf_db, ar_db, hv_db)

if __name__ == '__main__':
    main()