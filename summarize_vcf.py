import subprocess as sp
import re
import json
import argparse

def get_files(input_json, process):
    '''
    (str) -> list[str]

    Extracts the file paths from input_json

    Parameters
    ----------
    - input_json (str): path to the file that contains the vcf files

    '''
    data = {}

    with open(input_json) as file:
        data = json.load(file)
    file.close()

    #get mutect2 files
    vcf_files = {}

    for case_id, case in data["cases"].items():
        if process in case.keys():
            vcf_files[case_id] = case[process]["file"]

    return vcf_files


def get_num_calls(file):
    '''
    (str) -> int

    Returns the number of calls in the vcf file

    Parameters
    ----------
    - file (str): path to the vcf file 

    '''
    return int(
        sp.check_output(
            "zcat " + file + "| grep -v \"#\" | wc -l",
            shell=True
            )
            .decode('ascii')
            .strip()
        )


def get_num_pass(file):
    '''
    (str) -> int

    Returns the number of PASS calls in the vcf file

    Parameters
    ----------
    - file (str): path to the vcf file 

    '''
    return int(
        sp.check_output(
            "zcat " + file + " | awk \'!/^#/ {count[$7]++} END {print count[\"PASS\"]}\'",
            shell=True
            )
            .decode('ascii')
            .strip()
        )


def get_key(regex, full_name):
    '''
    (str, str) -> str

    Returns the portion of full_name that matches the given regex.
    Will return full_name if no matches are found.

    Parameters
    ----------
    - regex (str): regex specifying the abbreviated name
    - full_name (str): the full name of the file

    '''
    match = re.search(regex, full_name)    
    return match.group() if match is not None else full_name


def get_sn(file):
    '''
    (str) -> dict

    Selects statistics from the SN section of the bcftools summary generated from file

    Parameters
    ----------
    - file (str): path to the vcf file to be analyzed

    '''
    sn_stats = {}
    target_stats = ["records", "SNPs", "MNPs", "indels"]

    regex = "^SN.*" + target_stats[0]
    for stat in target_stats[1:]:
        regex = regex + "|" + "^SN.*"+ stat

    data = sp.check_output(
        "bcftools stats -f \"PASS\" " + file + "|grep -E '" + regex + "' |awk '{print $(NF)}'",
        shell=True
        ).decode('ascii').strip()
    data = data.split("\n")
    
    for count, stat in enumerate(target_stats):
        sn_stats[stat] = int(data[count])

    return sn_stats


def get_json_data(input_json, file_re):
    '''
    (str, str) -> str

    Returns a json object with keys as the file names abbreviated to the regex
    specified by file_re, containing summary data.

    Parameters
    ----------
    - input_json (str): file containing the paths to the vcf files to be analyzed
    - file_re (str): regex specifying the abbreviated name

    '''
    processes = ["mutations", "wg_structual_variants"]
    data = {}

    with open(input_json) as file:
        data = json.load(file)
    file.close()

    cases = { case_id : {} for case_id in data["cases"].keys()}
    for case in cases.keys():
        for process in processes:
            cases[case][process] = {}

    counter = 0

    for process in processes:
        vcf_files = get_files(input_json, process)
        for case_id, vcf_file in vcf_files.items():
            print(counter)
            cases[case_id][process]["num_calls"] = get_num_calls(vcf_file)
            cases[case_id][process]["num_pass"] = get_num_pass(vcf_file)
            cases[case_id][process]["bcf_PASS_summary"] = get_sn(vcf_file)
            counter = counter + 1

    with open('output.json', 'w', encoding='utf-8') as file:
        json.dump(cases, file, ensure_ascii=False, indent=4)
    

if __name__ == '__main__':
    
    # create top-level parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str, required=True)

    args = parser.parse_args()

    mutect_vcf_file_re = "PANX[^\/]*(?=\.mutect2)"
    get_json_data(args.infile, mutect_vcf_file_re)
