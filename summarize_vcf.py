import subprocess as sp
import re
import json
import argparse
import pandas

def get_files(data, processes):
    '''
    (dict, list[str]) -> list[str]

    Extracts file paths for each process in processes from data

    Parameters
    ----------
    - data (dict): dictionary containing parsed data
    - processes (list[str]): list of processes

    '''
    files = { case: {} for case in data["cases"].keys() }

    for case_id, case in data["cases"].items():
        for process in processes:
            try:
                files[case_id][process] = case[process]["file"]
            except:
                print(f"Analysis report generation: {case_id} is missing {process} file link")

    return files


def get_cli_output(command):
    '''
    (str) -> int
    Returns the number of calls specified by command

    Parameters
    ----------
    - command (str): command for command line

    '''
    return int(
        sp.check_output(
            command,
            shell=True
            )
            .decode('ascii')
            .strip()
        )


def vcf_get_sn(file):
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

    data =  sp.check_output(
            "bcftools stats -f \"PASS\" " + file + "|grep -E '" + regex + "' |awk '{print $(NF)}'",
            shell=True
            ).decode('ascii').strip()
    data = data.split("\n")
    
    for count, stat in enumerate(target_stats):
        sn_stats[stat] = int(data[count])

    return sn_stats

def get_mutations_data(file):
    '''
    (str) -> dict

    Selects statistics for the mutations process from file

    Parameters
    ----------
    - file (str): path to file to be parsed

    '''
    data = {}
    data["num_calls"] = get_cli_output(f"zcat {file} | grep -v \"#\" | wc -l")
    data["num_pass"] = get_cli_output("zcat " + file + " | awk \'!/^#/ {count[$7]++} END {print count[\"PASS\"]}\'")
    data["bcf_PASS_summary"] = vcf_get_sn(file)
    return data

def get_wgsv_data(file):
    '''
    (str) -> dict

    Selects statistics for the wg_structural_variants process from file

    Parameters
    ----------
    - file (str): path to file to be parsed

    '''
    #essentially a repeat of get_mutations_data, but repeated since parsing may diverge
    data = {}
    data["num_calls"] = get_cli_output(f"zcat {file} | grep -v \"#\" | wc -l")
    data["num_pass"] = get_cli_output("zcat " + file + " | awk \'!/^#/ {count[$7]++} END {print count[\"PASS\"]}\'")
    data["bcf_PASS_summary"] = vcf_get_sn(file)
    return data

def get_msv_data(file):
    '''
    (str) -> dict

    Selects statistics for the mavis_structural_variants process from file

    Parameters
    ----------
    - file (str): path to file to be parsed

    '''
    data = {}
    data["num_fusions"] = get_cli_output(f"cat {file}| grep -v \"#\" | wc -l")
    return data

def get_fusions_data(file):
    '''
    (str) -> dict

    Selects statistics for the fusions process from file

    Parameters
    ----------
    - file (str): path to file to be parsed

    '''
    data = {}
    data["num_fusions"] = get_cli_output(f"cat {file}| grep -v \"#\" | wc -l")
    return data

def get_cns_data(file):
    '''
    (str) -> dict

    Selects statistics for the copynumber_segmentation process from file

    Parameters
    ----------
    - file (str): path to file to be parsed

    '''
    data = {}
    data["num_cnv"] = get_cli_output(f"tail -n +2 {file} | wc -l")
    return data

def format_dataframe(data_list):
    df = pandas.json_normalize(data_list)
    df.set_index("id", inplace=True)
    df.rename(columns=lambda x:x.replace(".", "_"), inplace=True)

    return df

def parse_record(input_json):
    '''
    (dict) -> df

    Returns a pandas dataframe object containing summary data.

    Parameters
    ----------
    - input_json (str): file containing the paths to the vcf files to be analyzed

    '''
    processes = [
        "mutations",
        "wg_structual_variants",
        "mavis_structual_variants",
        "fusions",
        "copynumber_segmentation"
    ]

    process_parsers = {
        "mutations": get_mutations_data,
        "wg_structual_variants": get_wgsv_data,
        "mavis_structual_variants": get_msv_data,
        "fusions": get_fusions_data,
        "copynumber_segmentation": get_cns_data,
    }
    
    data = {}

    with open(input_json) as file:
        data = json.load(file)
    file.close()

    caches = { process: [] for process in processes }

    files = get_files(data, processes)

    for case in data["cases"].keys():
        for process in files[case].keys():
            try:
                record = process_parsers[process](files[case][process])
                record["id"] = case
                caches[process].append(record)
            except:
                print(f"Error: invalid file given for {case} in {process}")

    with open('output.json', 'w', encoding='utf-8') as file:
        json.dump(caches, file, ensure_ascii=False, indent=4)

    dataframes = { process: format_dataframe(caches[process]) for process in processes }

    # for v in dataframes.values():
    #     print(v.to_string())

    return dataframes
    

if __name__ == '__main__':
    
    # create top-level parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str, required=True)

    args = parser.parse_args()

    parse_record(args.infile)
