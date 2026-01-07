__author__ = 'Yueli Zheng'

import re
import pandas as pd
import tools


# isolate full path
# retrieve the tax files for each isolate


def retrieve_taxo(fullpath):
    isolate_taxs = {}
    for path in fullpath:
        file_dict = {}
        try: # handling for cases where no .tax file is created
            tax_file = path + "/" + path.split("/")[len(path.split("/")) - 1] + ".tax"
            file_content = pd.read_csv(tax_file, header=None).drop(0)
            for tax_level in list(file_content.values):
                file_dict[tax_level[0].split(":")[0]] = tax_level[0].split("\t")[1].strip()
            isolate_taxs[path.split("/")[len(path.split("/")) - 1]] = file_dict
        except FileNotFoundError:
            file_dict = {'D': '', 'P': '', 'C': '', 'O': '', 'F': '', 'G': '', 's': ''}
            isolate_taxs[path.split("/")[len(path.split("/")) - 1]] = file_dict
    return isolate_taxs


def on_bad_lines(bad_line):
    """This helps with mlst import when there are variable numbers of columns. use for newer versions of python."""
    # Remove the extra tab or field from the bad line
    bad_line = bad_line[:-1]
    # Return the bad line, so that it is not skipped
    return bad_line

def retrieve_mlst(fullpath):
    isolate_mlst = {}
    for path in fullpath: # if a griphin samplesheet was passed
        try: # handling for cases where no _combined.tsv file is created
            mlst_file = path + "/mlst/" + path.split("/")[len(path.split("/")) - 1] + "_combined.tsv"
            file_content = pd.read_csv(mlst_file, delimiter='\t', header=0, usecols = ["WGS_ID","Source", "Pulled_on", "Database", "ST"])
            if file_content.empty:
                print(path.split("/")[len(path.split("/")) - 1] + ": No mlst info found")
            else:
                # check "PARALOG"
                result = file_content.astype(str).apply(lambda x: x.str.contains("PARALOG", case=False)).values.flatten()
                if result.any():
                    file_content = file_content.drop(file_content.index[file_content["ST"].str.contains("PARALOG")]).reindex()
                # check "Oxford" "Pasteur"
                indices = file_content.index[file_content["Database"].str.contains("Oxford") |
                                            file_content["Database"].str.contains("Pasteur") | file_content["Database"].str.contains("Achtman")]
                tmpdict = tools.rearrange_oxford_pasteur(file_content, indices)
                isolate_mlst[file_content["WGS_ID"][0].strip()] = tmpdict[file_content["WGS_ID"][0].strip()]
        except FileNotFoundError:
            isolate_mlst[file_content["WGS_ID"][0].strip()] = ""
    return isolate_mlst


def retrieve_mlst_nonovel(fullpath):
    isolate_mlst = {}
    for path in fullpath:
        try:
            #print("start new isolate: " + path.split("/")[len(path.split("/")) - 1] + "_combined.tsv")
            isolate_name = path.split("/")[len(path.split("/")) - 1].replace('.filtered.scaffolds.fa', '')
            mlst_file = path + "/mlst/" + path.split("/")[len(path.split("/")) - 1] + "_combined.tsv"
            file_content = pd.read_csv(mlst_file, delimiter='\t', header=0, usecols = ["WGS_ID", "Source", "Pulled_on", "Database", "ST"])
            #clean file names
            file_content["WGS_ID"] = file_content["WGS_ID"].str.strip().str.replace('.filtered.scaffolds.fa', '')
            if file_content.empty:
                isolate_mlst[isolate_name] = ""
            else:
                if file_content["Database"][0] == '-':
                    file_content = file_content.drop(file_content[(file_content["Database"] == '-')].index)
                # check "Novel", "-", "Missing"
                result = file_content.astype(str).apply(lambda x: x.str.contains("Novel|Missing|-", case=False)).values.flatten()
                if result.any():
                    file_content["ST"] = file_content["ST"].astype(str)
                    file_content = file_content.drop(file_content.index[file_content["ST"].str.contains("Novel|Missing|-", case=False)]).reset_index(drop=True)
                if file_content.shape[0] == 0:
                    isolate_mlst[isolate_name] = ""
                else:
                    # check "Oxford" "Pasteur"
                    indices = file_content.index[file_content["Database"].str.contains("Oxford") |
                                                file_content["Database"].str.contains("Pasteur") | file_content["Database"].str.contains("Achtman")]
                    tmp_dict = tools.rearrange_oxford_pasteur(file_content, indices)
                    isolate_mlst[str(file_content["WGS_ID"][0]).strip()] = tmp_dict[str(file_content["WGS_ID"][0]).strip()]
        except FileNotFoundError:
            isolate_mlst[isolate_name] = ""
    return isolate_mlst

def retrieve_instrument_model(fullpath):
    seqmachine_model_list = {}
    seqmachine_model_isolate = {}
    isolate_no_model = []
    for path in fullpath:
        isolate_name = path.split("/")[len(path.split("/")) - 1]
        fastq_path = path + "/fastp_trimd/" + isolate_name + "_1.trim.fastq.gz"
        seq_machine_model = pd.read_csv(fastq_path, compression='gzip', nrows=1).columns.values[0].split(":")[0].replace("@", "").strip()
        seqmachine_model_list[isolate_name] = seq_machine_model
        instrument = search_instrument(seq_machine_model)[0]
        if instrument == "":
            isolate_no_model.append(isolate_name)
        seqmachine_model_isolate[isolate_name] = instrument
    print('The following isolates do not have sequence machine model information provided:')
    for item in isolate_no_model:
        print(item)
    return seqmachine_model_isolate


def search_instrument(fastq_model):
    for key in InstrumentIDs:
        if re.search(key, fastq_model):
            return InstrumentIDs[key]
    return [""]


InstrumentIDs = {"HWI-M[0-9]{4}$": ["MiSeq"],
                 "HWUSI": ["Genome Analyzer IIx"],
                 "M[0-9]{5}$": ["MiSeq"],
                 "HWI-C[0-9]{5}$": ["HiSeq"],
                 "C[0-9]{5}$": ["HiSeq"],
                 "HWI-D[0-9]{5}$": ["HiSeq"],
                 "D[0-9]{5}$": ["HiSeq"],
                 "J[0-9]{5}$": ["HiSeq"],
                 "K[0-9]{5}$": ["HiSeq"],
                 "E[0-9]{5}$": ["HiSeq X"],
                 "NB[0-9]{6}$": ["NextSeq"],
                 "NS[0-9]{6}$": ["NextSeq"],
                 "MN[0-9]{5}$": ["MiniSeq"],
                 "FS[0-9]{5}$": ["ISeq"]}