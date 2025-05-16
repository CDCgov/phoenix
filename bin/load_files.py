__author__ = "Yueli Zheng"

import os
import sys
import yaml
import pandas as pd
import tools


class FileLoader:
    def __init__(self):
        """"""

    # This list of isolates are required to submit to NCBI after the review or the publication
    def import_isolate_list(self, isolatelist):
        isolate_file_content = pd.read_csv(isolatelist, sep="\t", header=None)
        pos_miseq = tools.find_index(
            list(isolate_file_content[0]), "MiSeqAnalysisFiles"
        )[0]
        pos_nov = tools.find_index(
            list(isolate_file_content[0]), "NovaSeqAnalysisFiles"
        )[0]
        pos_haiseq = tools.find_index(list(isolate_file_content[0]), "NCBI_HAISeq")[0]
        miseq_list = list(isolate_file_content[0])[pos_miseq:pos_nov]
        nov_list = list(isolate_file_content[0])[pos_nov:pos_haiseq]
        hai_list = list(isolate_file_content[0])[
            pos_haiseq : list(isolate_file_content[0]).__len__()
        ]
        pos_list_miseq = self.find_pos_sep_project(miseq_list)
        pos_list_haiseq = self.find_pos_sep_project(hai_list)
        pos_list_novseq = self.find_pos_sep_project(nov_list)
        path_base_runname = []
        path_base_runname.append(self.prep_path(miseq_list, pos_list_miseq))
        path_base_runname.append(self.prep_path(nov_list, pos_list_novseq))
        path_base_runname.append(self.prep_path(hai_list, pos_list_haiseq))
        return path_base_runname

    def find_pos_sep_project(self, isolatelist):
        pos = []
        for iso in isolatelist:
            if "project" in iso:
                pos.append(isolatelist.index(iso))
            if "runname" in iso:
                pos.append(isolatelist.index(iso))
        return pos

    # prepare the isolates' paths to get the taxo info of the isolate list.
    def prep_path(self, isolate_list, pos_list):
        path_base_run_name = {}
        path_base = isolate_list[0]
        run_name_list = []
        for i in range(0, len(pos_list), 1):
            sample_list = []
            sample_dict = {}
            if i + 1 < len(pos_list):
                runname = isolate_list[pos_list[i]].split(":")[1].strip()
                for j in range(pos_list[i], pos_list[i + 1], 1):
                    if (
                        "project" not in isolate_list[j]
                        and "runname" not in isolate_list[j]
                    ):
                        sample_list.append(isolate_list[j])
                sample_dict[runname] = sample_list
            if i + 1 >= len(pos_list):
                runname = isolate_list[pos_list[i]].split(":")[1].strip()
                for j in range(pos_list[i], len(isolate_list), 1):
                    if (
                        "project" not in isolate_list[j]
                        and "runname" not in isolate_list[j]
                    ):
                        sample_list.append(isolate_list[j])
                sample_dict[runname] = sample_list
            run_name_list.append(sample_dict)
        path_base_run_name[path_base] = run_name_list
        return path_base_run_name

    def get_full_path(self, path_base_run_name):
        path = []
        for base_path in path_base_run_name:
            for runname_list in list(base_path.values()):
                for run_name in runname_list:
                    runname = list(run_name.keys())[0]
                    for asample in list(run_name.values())[0]:
                        apath = (
                            list(base_path.keys())[0] + "/" + runname + "/" + asample
                        )
                        path.append(apath)
        return path

    def get_full_path_griphin(self, griphin):
        path = pd.read_csv(griphin, sep="\t", header=None).drop(0)
        path_list = []
        for value in list(path.values):
            path_list.append(value[0].split(",")[1])
        return path_list

    def get_isolate_name(self, griphin):
        path = pd.read_csv(griphin, header=None).drop(0)
        isolate_name = []
        for value in list(path.values):
            isolate_name.append(value[0].split(",")[0])
        return isolate_name

    # load the bioprojects
    def load_bioproject(self, bioproject_file):
        # osii_bioproject_file = os.path.join(os.path.dirname(__file__), bioproject_file)
        # bioprojects = pd.read_csv(bioproject_file, sep='\t', header=None)
        with open(bioproject_file, "r") as file:
            bioprojects_taxo = yaml.safe_load(file)
        bioprojects_taxo_clean = {
            key: value
            for key, value in bioprojects_taxo.items()
            if value is not None and value != "blank"
        }
        return bioprojects_taxo_clean


def main():
    metafile = sys.argv[1]
    inputfile2 = sys.argv[2]
    isolate_list = sys.argv[3]
    griffin_csv = sys.argv[4]
    ob = FileLoader()
