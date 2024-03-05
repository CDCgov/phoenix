#!/usr/bin/env python3
"""
CLIA WGS Run Summary generating script for CLIA Phoenix
12-14-2023
@author: Frank Bao
email: dyp9@cdc.gov
"""

# Function to get the script version
def get_version():
    return "1.0.0"

import argparse
import pandas as pd
from weasyprint import HTML, CSS
from datetime import date
import logging
import os, glob

# Minimum cutoff for Error
COVERAGE_CUTOFF = 50
ASSEMBLY_CUTOFF = 2.58
ASSEMBLY_LENGTH = 1000000
SCAFFOLD_NUM = 500

#Cutoff for warnings
GC = 2.58
RAW_READ  = 1000000
TRIMMED_READ = 1000000
R1_Q30_RAW = 90
R2_Q30_RAW = 70
R1_Q30_TRIM = 90
R2_Q30_TRIM = 70
SCAFFOLD = 200
BUSCO = 97
FASTANI_MATCH = 95
FASTANI_COVERAGE = 90
KRAKEN_ASSEMBLY_GENUS = 70

note_text1 = "Failed QC - Numbers below thresholds will be highlighted in red: <30x coverage; Assembly stdev >2.58; Min assembly length <1,000,000 bp; >500 scaffolds"
note_text2 = "QC Warning  - Numbers below thresholds will be highlighted in yellow: <1,000,000 total reads for raw and trimmed; % trimmed reads with Q30 average for R1 (<90%) and R2 (<70%)"
note_text4 = "Taxa ID Failed - Only Taxa IDed with FastANI are acceptable, all other methods will be highlighted in red"
note_text3 = "Taxa ID QC Warning  - Numbers below thresholds will be highlighted in yellow: BUSCO ID <97%; FastANI match <95%; FastANI coverage <90%"
note_text5 = "Assembly ratio and GC% STDev are NaN when there are <10 genomes as reference."
note_text6 = "AR genes were determined based on: ."

def check_time(start_time):
    from datetime import datetime
    # Parse the timestamp
    timestamp_obj = datetime.fromisoformat(start_time.replace("Z", "+00:00"))

    # Format the timestamp in the desired format
    formatted_start_time = timestamp_obj.strftime("%m/%d/%Y-%H:%M")
    end_time = datetime.now()
    formatted_end_time = end_time.strftime("%m/%d/%Y-%H:%M")
    return formatted_start_time, formatted_end_time

def color_string_red(val):
    if val == "Unknown":
        return f'<font style="background-color: lightcoral">{val}</font>'
    elif val != "ANI_REFSEQ":
        return f'<font style="background-color: lightcoral">{val}</font>'
    elif val == "ANI_REFSEQ":
        return val

def color_red(val, cutoff, flag):
    if val == "Unknown":
        return f'<font style="background-color: lightcoral">{val}</font>'
    elif flag == ">":
        if float(str(val).replace(",","")) > float(cutoff):
            return f'<font style="background-color: lightcoral">{val}</font>'
        else:
            return val
    else:
        if float(str(val).replace(",","")) < float(cutoff):
            return f'<font style="background-color: lightcoral">{val}</font>'
        else:
            return val
        
def color_unknown_red(val):
    if val == "Unknown":
        return f'<font style="background-color: lightcoral">{val}</font>'
    elif val == "FAIL":
        return f'<font style="background-color: lightcoral">{val}</font>'
    else:
        return val

def color_yellow(val, cutoff, flag):
    if val == "Unknown":
        return f'<font style="background-color: lightcoral">{val}</font>'
    elif flag == ">":
        if float(str(val).replace(",","")) > float(cutoff):
            return f'<font style="background-color: lightgoldenrodyellow">{val}</font>'
        else:
            return val
    else:
        if float(str(val).replace(",","")) < float(cutoff):
            return f'<font style="background-color: lightgoldenrodyellow">{val}</font>'
        else:
            return val

# Define a function to join non-NaN values with comma
def combine_strings(row):
    values = [str(value) for value in row if pd.notnull(value)]
    return ', '.join(values)

def clia_report(args):
    project_dir = args.project_dir
    ar_database = args.ar_database.replace("amrfinderdb_", "").replace(".tar.gz", "")
    run_start_time, run_end_time = check_time(args.start_time)

    summary = os.path.join(project_dir,'Phoenix_Summary.tsv')
    if os.path.exists(summary):
        summary_df = pd.read_csv(summary,sep='\t')
        #summary_df.drop(columns=["Auto_QC_Outcome","Warning_Count","Estimated_Coverage","Genome_Length","Assembly_Ratio_(STDev)","#_of_Scaffolds_>500bp","GC_%","BUSCO"], inplace=True)
        summary_df.drop(summary_df.columns[[i for i in range(1,10)] + [i for i in range(11,16)]], axis=1, inplace=True)
        # Combine strings in columns containing "Beta-lactam"
        beta_lactam_columns = [col for col in summary_df.columns if 'Beta-lactam' in col]
        summary_df['Beta-lactam Resistance Genes'] = summary_df[beta_lactam_columns].apply(lambda row: combine_strings(row), axis=1)
        # Combine strings in non "Beta-lactam" columns
        non_beta_lactam_columns = [col for col in list(summary_df.columns) if col not in (beta_lactam_columns + ['ID', 'Species','Beta-lactam Resistance Genes'])]
        summary_df['Other Resistance Genes'] = summary_df[non_beta_lactam_columns].apply(lambda row: combine_strings(row), axis=1)
        #drop individual combined columns
        summary_df.drop(columns=beta_lactam_columns, inplace=True)
        summary_df.drop(columns=non_beta_lactam_columns, inplace=True)
        summary_df = summary_df.rename(columns={'Species':'FastANI Species'})
    else:
        logging.info(f"No summary file: {summary}")
    
    griphin_summary = f"{project_dir}/*GRiPHin_Summary.tsv"
    matching_files = glob.glob(griphin_summary)
    if len(matching_files) == 1:
        columns_to_read = ['WGS_ID','Minimum_QC_Check','Raw_Q30_R1_[%]','Raw_Q30_R2_[%]','Total_Raw_[reads]','Total_Trimmed_[reads]','Estimated_Trimmed_Coverage', 'GC[%]','Scaffolds','Assembly_Length','Assembly_Ratio','Assembly_StDev']
        griphin_df = pd.read_csv(matching_files[0],sep='\t', usecols=columns_to_read)
        griphin_df = griphin_df.rename(columns={'WGS_ID':'ID','Minimum_QC_Check':'Minimum QC','Raw_Q30_R1_[%]':'R1 Q30 (%)', 'Raw_Q30_R2_[%]':'R2 Q30 (%)','Total_Raw_[reads]':'Total Raw', 'GC[%]': 'GC (%)',
                                                'Total_Trimmed_[reads]':'Total Trimmed','Estimated_Trimmed_Coverage':'Estimated Coverage','Assembly_Length':'Assembly Length','Assembly_Ratio':'Assembly Ratio','Assembly_StDev':'Assembly StDev'})
        #griphin_df['Assembly Length'] = griphin_df['Assembly Length'].apply(lambda x: f'{int(x.replace(",", "")):,}' if x != "Unknown" and x == x else x).astype(str)
        griphin_df[['Total Raw','Total Trimmed','Scaffolds']] = griphin_df[['Total Raw','Total Trimmed','Scaffolds']].applymap(lambda x: f'{int(x.replace(",", "")):,}' if x != "Unknown" and x == x else x).astype(str)
        #griphin_df['Total Trimmed'] = griphin_df['Total Trimmed'].apply(lambda x: f'{int(x.replace(",", "")):,}' if x != "Unknown" and x == x else x).astype(str)
        # Round to 2 decimals
        griphin_df[['Assembly Ratio', 'Assembly StDev']] = griphin_df[['Assembly Ratio', 'Assembly StDev']].applymap(lambda x: round(float(x), 2) if x != "Unknown" and x == x else x)
        #unknowns
        griphin_df[['GC (%)','Assembly Ratio']] = griphin_df[['GC (%)','Assembly Ratio']].applymap(lambda x: color_unknown_red(x))
        #Failures
        griphin_df['Minimum QC'] = griphin_df['Minimum QC'].apply(lambda x: color_unknown_red(x))
        # Failures
        griphin_df['Estimated Coverage'] = griphin_df['Estimated Coverage'].apply(lambda x: color_red(x, COVERAGE_CUTOFF, "<"))
        griphin_df['Assembly StDev'] = griphin_df['Assembly StDev'].apply(lambda x: color_red(x, ASSEMBLY_CUTOFF, ">"))
        griphin_df['Assembly Length'] = griphin_df['Assembly Length'].apply(lambda x: color_red(x, ASSEMBLY_LENGTH, "<"))
        griphin_df['Scaffolds'] = griphin_df['Scaffolds'].apply(lambda x: color_red(x, SCAFFOLD_NUM, ">"))
        # Warnings
        griphin_df['R1 Q30 (%)'] = griphin_df['R1 Q30 (%)'].apply(lambda x: color_yellow(x, R1_Q30_RAW, "<"))
        griphin_df['R2 Q30 (%)'] = griphin_df['R2 Q30 (%)'].apply(lambda x: color_yellow(x, R2_Q30_RAW, "<"))
        griphin_df['Total Raw'] = griphin_df['Total Raw'].apply(lambda x: color_yellow(x, RAW_READ, "<"))
        griphin_df['Total Trimmed'] = griphin_df['Total Trimmed'].apply(lambda x: color_yellow(x, TRIMMED_READ, "<"))
        #griphin_df['Scaffolds'] = griphin_df['Scaffolds'].apply(lambda x: color_red(x, SCAFFOLD_NUM, ">"))

        taxa_columns_to_read = ['WGS_ID','Taxa_Source','BUSCO_Lineage','BUSCO_%Match','Kraken_ID_Raw_Reads_%','Kraken_ID_WtAssembly_%','FastANI_%ID','FastANI_%Coverage','Taxa_Source']
        taxa_df = pd.read_csv(matching_files[0],sep='\t', usecols=taxa_columns_to_read)
        taxa_df = taxa_df.rename(columns={'WGS_ID':'ID','BUSCO_Lineage':'BUSCO Lineage', 'Taxa_Source':'Taxonomic ID Source', 'FastANI_%ID':'FastANI Match (%)','FastANI_%Coverage':'Bases Aligned to FastANI Taxon (%)',
                                            'BUSCO_%Match':'BUSCO Match (%)','Kraken_ID_Raw_Reads_%':'Kraken Raw Reads (%)','Kraken_ID_WtAssembly_%':'Kraken Assembly (%)'})
        #unknowns
        taxa_df[['BUSCO Lineage','Kraken Raw Reads (%)','Kraken Assembly (%)']] = taxa_df[['BUSCO Lineage','Kraken Raw Reads (%)','Kraken Assembly (%)']].applymap(lambda x: color_unknown_red(x))
        # Failures
        taxa_df['Taxonomic ID Source'] = taxa_df['Taxonomic ID Source'].apply(lambda x: color_string_red(x))
        # Warnings
        taxa_df['FastANI Match (%)'] = taxa_df['FastANI Match (%)'].apply(lambda x: color_yellow(x, FASTANI_MATCH, "<"))
        taxa_df['Bases Aligned to FastANI Taxon (%)'] = taxa_df['Bases Aligned to FastANI Taxon (%)'].apply(lambda x: color_yellow(x, FASTANI_COVERAGE, "<"))
        #taxa_df['Genome Size'] = taxa_df['Genome Size'].apply(lambda x: color_red(x, ASSEMBLY_LENGTH, "<"))
        taxa_df['BUSCO Match (%)'] = taxa_df['BUSCO Match (%)'].apply(lambda x: color_yellow(x, BUSCO, "<"))
    else:
        logging.info(f"No GRiPHin file: {matching_files[0]}")
    
    page_title = "WGS Run Summary"
    title = "WGS Run Summary"
    stitle1 = "Run Information"
    stitle2 = "Data Quality"
    stitle4 = "Taxonomic Identification Quality"
    stitle3 = "Predicted Organism and Resistance Genes"
    today = date.today().strftime("%m/%d/%Y")
    today_file_name = date.today().strftime("%Y-%m-%d")
    footer = "Research use only."

    # HTML #
    html = f'''
    <!DOCTYPE html>
    <html lang="en-us">
        <head>
            <meta charset="UTF-8">
            <title>{page_title}</title>
            <style>
            body {{
            font-family: Arial, Helvetica, sans-serif;
            }}
            h1 {{
            font-family: Arial, Helvetica, sans-serif;
            text-align: center;
            color: #1e5c85;
            }}
            h2 {{
            font-family: Arial, Helvetica, sans-serif;
            margin-top: 1em;
            margin-bottom: 0.4em;
            text-align: left;
            color: #D6672E;
            }}
            table {{
            border-collapse: collapse;
            border: none;
            font-size: 0.9em;
            }}
            th, td {{
            text-align: left;
            padding-left: 5px;
            padding-right: 5px;
            padding-top: 1px;
            padding-bottom: 1px;
            border: none;
            }}

            table.dataframe tr:nth-child(even) {{background-color: #f2f2f2;}}

            table.dataframe th {{
            background-color: #1e5c85;
            color: white;
            }}
            </style>
        </head>
        <body>
        <header>
        <h1>{title}</h1>
        </header>
        <hr>
        <article>
        <h2>{stitle1}</h2>
        <table>
        <tr>
        <td>Lab:</td>
        <td>Clinical and Environmental Microbiology Branch</td>
        </tr>
        <tr>
        <td>Run Operator:</td>
        <td>Frank Bao</td>
        </tr>
        <tr>
        <td>Run Start date:</td>
        <td>{run_start_time}</td>
        </tr>
        <tr>
        <td>Run End date:</td>
        <td>{run_end_time}</td>
        </tr>
        <tr>
        <td>Report date:</td>
        <td>{today}</td>
        </tr>
        </table>
        <h2>{stitle2}</h2>
        {griphin_df.to_html(index=False, escape=False, justify="left", classes = 'table table-bordered', border=1)}
        <p style='font-size: 0.8em'>1. {note_text1}</br>
        2. {note_text2}</br>
        3. {note_text5}</p>
        <h2>{stitle4}</h2>
        {taxa_df.to_html(index=False, escape=False, justify="left", classes = 'table table-bordered', border=1)}
        <p style='font-size: 0.8em'>1. {note_text4}</br>
        2. {note_text3}</p>
        <h2>{stitle3}</h2>
        {summary_df.to_html(index=False, justify="left")}
        <p style='font-size: 0.8em'>1. {note_text6}</br>
        <p></p>
        </article>
        <hr>
        <footer>
        <p><i>This report was created using the PHoeNIx <a href="https://github.com/CDCgov/phoenix">{args.phx_version}</a> bioinformatics pipeline CLIA entry point.</br>
        AR genes were identified with AMRFinderPlus <a href="https://github.com/ncbi/amr">v{args.amrfinder_version}</a> using database version {ar_database}</br>
        {footer}</i></p>
        </footer>
    </body>
    </html>
    '''

    # font is good for html, but too big for pdf so adjusting for pdf
    html_css = '''
    body {
        font-size:1em; 
    }
    h1 {
        font-size: 1.5em;
    }
    h2 {
        font-size: 1.2em;
    }
    table {
        font-size: 0.9em;
    }
    footer {
        font-size: 0.9em;
    }
    '''
    # Define the PDF styling
    pdf_css = CSS(string=
    '''
    @page {
        size: letter;
        margin: 0in 0.44in 0.2in 0.44in;
    }
    body {
        font-family: Arial;
        font-size: 10px;
    }
    table {
        font-size: 8px;
    }
    h1 {
        font-family: Arial;
        font-size: 16px;
    }
    h2 {
        font-family: Arial;
        font-size: 14px;
    }
    footer {
        font-size: 9px;
    }
    ''')

    html_file = f"WGS_Run_Summary_report_{today_file_name}.html"
    pdf_file = f"WGS_Run_Summary_report_{today_file_name}.pdf"
    with open(html_file, "w") as f:
        f.write(html_css + html)
    # Convert HTML to PDF
    HTML(string=html).write_pdf(pdf_file, stylesheets=[pdf_css])
    logging.info("Report generation is completed.")

def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawTextHelpFormatter,
        description = "Summary report generating script for CLIA Phoenix"
    )
    parser.add_argument('-p', '--project_dir', type=str, help='Phoenix output folder', required=True)
    parser.add_argument('-t', '--start_time', dest="start_time", type=str, help='Time the pipeline started.', required=False)
    parser.add_argument('--phx_version',dest="phx_version", type=str, help='Phoenix version', required=True)
    parser.add_argument('--amrfinder_version',dest="amrfinder_version", type=str, help='AMRFinderPlus version', required=True)
    parser.add_argument('--ar_database',dest="ar_database", type=str, help='AR Databased used with AMRFinderPlus.', required=True)
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    opts = parser.parse_args()
    return opts

def main():
    logging.basicConfig(filename=f"report_generate_{date.today()}.log",
                        level=logging.INFO,
                    format='[%(asctime)s %(levelname)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w')
    opts = parse_args()
    clia_report(opts)


if __name__ == "__main__":
    main()
