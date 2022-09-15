import subprocess as sp
import re
import json
import argparse


def get_file_names(input_file):
    '''
    (str) -> list[str]

    Returns the file names specified in input_file as a list

    Parameters
    ----------
    - input_file (str): path to the file that contains the vcf file names 

    '''
    with open(input_file) as file:
        file_names = file.readlines()
        file_names = [line.rstrip() for line in file_names]
    file.close()
    return file_names


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



def get_json_data(input_file, file_re):
    '''
    (str, str) -> str

    Returns a json object with keys as the file names abbreviated to the regex
    specified by file_re, containing the number of calls and the number of PASS calls.

    Parameters
    ----------
    - input_file (str): file containing the paths to the vcf files to be analyzed
    - file_re (str): regex specifying the abbreviated name

    '''

    vcf_file_names = get_file_names(input_file)

    data_dict = {}
    for file in vcf_file_names:
        if get_key(file_re, file) in data_dict:
            raise Exception("Duplicate file name, exiting...")
        else:
            data_dict[get_key(file_re, file)] = { "data" : {} }
            data_dict[get_key(file_re, file)]["data"]["num_calls"] = get_num_calls(file)
            data_dict[get_key(file_re, file)]["data"]["num_pass"] = get_num_pass(file)

    with open('output.json', 'w', encoding='utf-8') as file:
        json.dump(data_dict, file, ensure_ascii=False, indent=4)
    

if __name__ == '__main__':
    
    # create top-level parser
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=str, required=True)

    args = parser.parse_args()

    vcf_file_re = "PANX[^\/]*(?=\.mutect2)"
    get_json_data(args.infile, vcf_file_re)
