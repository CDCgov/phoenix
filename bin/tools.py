__author__ = 'Yueli Zheng'


def find_index(alist, search_string):
    indices = [i for i, item in enumerate(alist) if search_string in item]
    if indices.__len__() == 0:
        print("No " + search_string + " found.")
    return indices


def extract_string(my_list, special_string):
    extracted_elements = [item.split(special_string)[-1] for item in my_list]
    return extracted_elements


def rearrange_oxford_pasteur(df, indices):
    isolate_mlst = {}
    if len(indices) == 0:
        tmp = "MLST" + str(df["ST"][0]).strip() + "_" + df["Database"][0].strip() + ", "
        isolate_mlst[str(df["WGS_ID"][0]).strip()] = ",".join(tmp.rsplit(",", 1)[:-1])
    if len(indices) == 1:
        tmp = "MLST" + str(df["ST"][indices[0]]).strip() + "_" + df["Database"][indices[0]].strip().split("(")[1].strip(")")
        isolate_mlst[str(df["WGS_ID"][0]).strip()] = ", ".join(tmp.split(", ", 1)[::-1])
    if len(indices) == 2:
        tmp = ""
        count = 1
        for index in indices:
            if count < 2:
                ending = ", "
            else:
                ending = ""
            tmp = tmp + "MLST" + str(df["ST"][index]).strip() + "_" + df["Database"][index].strip().split("(")[1].strip(")") + ending
            count = count + 1
        isolate_mlst[str(df["WGS_ID"][0]).strip()] = ", ".join(tmp.split(", ")[::-1])
    if len(indices) == 3:
        tmp = ""
        count = 1
        for index in indices:
            if count < 3:
                ending = ", "
            else:
                ending = ""
            tmp = tmp + "MLST" + str(df["ST"][index]).strip() + "_" + df["Database"][index].strip().split("(")[1].strip(")") + ending
            count = count + 1
        isolate_mlst[str(df["WGS_ID"][0]).strip()] = ",".join(tmp.split(", ")[::-1])
    if len(indices) >= 3:
        isolate_mlst[str(df["WGS_ID"][0]).strip()] = "ERROR: Unexpected MLST results. Review ~/mlst/*_combined.tsv file!"
    return isolate_mlst

def purify_dict(dict_container):
    pure_dict = {}
    for key, value in dict_container.items():
        pure_dict[key] = value.sampleContent
    return pure_dict