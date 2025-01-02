#!/usr/bin/env python3
import argparse

##Makes a summary Excel file when given a series of output summary line files from PhoeNiX
## Written by Jill Hagey (qpk9@cdc.gov)

# Function to get the script version
def get_version():
    return "1.0.0"

def parseArgs(args=None):
    parser = argparse.ArgumentParser(description='Script to generate a PhoeNix summary excel sheet.')
    parser.add_argument('-s', '--summary_line', default=None, required=True, dest='summary_line_file', help='PHoeNIx style samplesheet of sample,directory in csv format. Directory is expected to have PHoeNIx stype output.')
    parser.add_argument('-m', '--mlst', default=None, required=True, dest='mlst_file', help='If a directory is given rather than samplesheet GRiPHin will create one for all samples in the directory.')
    parser.add_argument('-g', '--gamma', required=True, dest='gamma_report', help='CSV file with a list of sample_name,new_name. This option will output the new_name rather than the sample name to "blind" reports.')
    parser.add_argument('-o', '--output', default="", required=False, dest='output', help='Name of output file default is GRiPHin_Summary.xlsx.')
    parser.add_argument('-a', '--amrfinder', dest="amrfinder_report", required=True, help='Pass the name of the amrfinder database used. Only used for documentation purposes')
    parser.add_argument('--version', action='version', version=get_version())# Add an argument to display the version
    return parser.parse_args()

#set colors for warnings so they are seen
CRED = '\033[91m'+'\nWarning: '
CEND = '\033[0m'


def MLST_Scheme(MLST_file):
    """Pulls MLST info from *_Scheme.mlst file"""
    Scheme_list = [[],[],[],[],[]]
    with open(MLST_file, 'r') as f:
        lines = f.readlines()
        lines.pop(0)
        #print(len(lines), lines)
        for rawline in lines:
            line=rawline.strip()
            split_line = line.split("\t")
            #print("\n".join(split_line))
            source = split_line[1]
            date = split_line[2]
            DB_ID = split_line[3]
            Scheme = str(split_line[4])
            alleles = "-".join(split_line[5:])
            if DB_ID in Scheme_list[0]:
                print("In Scheme_list[0]")
                print(Scheme_list[0])
                for i in range(0,len(Scheme_list[0])):
                    if DB_ID == Scheme_list[0][i]:
                        print("Adding to", Scheme_list[0][i], i)
                        if Scheme != "-" and "Novel" not in Scheme: #if Scheme != "-" and Scheme != "Novel_allele" and Scheme != "Novel_profile":
                            Scheme_list[1][i].append("ST"+str(Scheme))
                        else:
                            Scheme_list[1][i].append(Scheme)
                        Scheme_list[2][i].append(alleles)
                        Scheme_list[3][i].append(source)
                        Scheme_list[4][i].append(date)
                        print(Scheme_list)
            else:
                print("NOT in Scheme_list[0]")
                print(Scheme_list[0], Scheme_list[1], Scheme_list[2], Scheme_list[3], Scheme_list[4])
                Scheme_list[0].append(DB_ID)
                if Scheme != "-" and Scheme != "Novel_allele" and Scheme != "Novel_profile":
                    Scheme_list[1].append(["ST"+Scheme])
                else:
                    Scheme_list[1].append([Scheme])
                Scheme_list[2].append([alleles])
                Scheme_list[3].append([source])
                Scheme_list[4].append([date])
                print(Scheme_list[0], Scheme_list[1], Scheme_list[2], Scheme_list[3], Scheme_list[4])
            for i in Scheme_list:
                for j in i:
                    print(j)
            #print("\n".join(Scheme_list))
    print(Scheme_list)
    return Scheme_list

def convert_MLST_scheme(MLST_file):
    try:
        Scheme = MLST_Scheme(MLST_file)
        if len(Scheme[0]) > 1:
            if Scheme[0][0].lower() < Scheme[0][1].lower():
                MLST_scheme_1 = Scheme[0][0]
                mlst_types_1=sorted(Scheme[1][0])[::-1]
                MLST_type_1 = ",".join(mlst_types_1)
                if Scheme[3][0] == "srst2":
                    MLST_type_1 += '^'
                MLST_scheme_2 = Scheme[0][1]
                mlst_types_2=sorted(Scheme[1][1])[::-1]
                MLST_type_2 = ",".join(mlst_types_2)
                if Scheme[3][1] == "srst2":
                    MLST_type_2 += '^'
            else:
                MLST_scheme_1 = Scheme[0][1]
                mlst_types_1=sorted(Scheme[1][1])[::-1]
                MLST_type_1 = ",".join(mlst_types_1)
                if Scheme[3][1] == "srst2":
                    MLST_type_1 += '^'
                MLST_scheme_2 = Scheme[0][0]
                mlst_types_2=sorted(Scheme[1][0])[::-1]
                MLST_type_2 = ",".join(mlst_types_2)
                if Scheme[3][0] == "srst2":
                    MLST_type_2 += '^'
        else:
            MLST_scheme_1 = Scheme[0][0]
            MLST_type_1 = ",".join(Scheme[1][0])
            MLST_scheme_2 = "-"
            MLST_type_2 = "-"
    except:
        MLST_scheme_1 = 'Unknown'
        MLST_scheme_2 = 'Unknown'
        MLST_type_1 = 'Unknown'
        MLST_type_2 = 'Unknown'


def replace_value_in_column(input_tsv, output_tsv, column_name, old_value, new_value):
    # Open the input TSV file for reading
    with open(input_tsv, 'r', newline='') as infile:
        reader = csv.DictReader(infile, delimiter='\t')
        
        # Get the fieldnames from the input file
        fieldnames = reader.fieldnames

        # Open the output TSV file for writing
        with open(output_tsv, 'w', newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            
            # Process each row
            for row in reader:
                if row[column_name] == old_value:
                    row[column_name] = new_value
                writer.writerow(row)

def main():
    args = parseArgs()
    MLST_Scheme(args.mlst_file)
    replace_value_in_column(args.summary_line_file, args.output)


if __name__ == '__main__':
    main()