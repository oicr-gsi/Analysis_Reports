import subprocess as sp
import re
import json


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


def get_short_name(regex, full_name):
    '''
    (str, str) -> str

    Returns the portion of full_name that matches the given regex.
    Will return the full name if no matches are found.

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

    Returns a json object with the file name abbreviated to the regex
    specified by file_re and number of calls and the number of PASS calls.

    Parameters
    ----------
    - input_file (str): file containing the paths to the vcf files to be analyzed
    - file_re (str): regex specifying the abbreviated name

    '''

    vcf_file_names = get_file_names(input_file)

    data_dict = { i : {"file_name" : {}, "data" : {} } for i in range(len(vcf_file_names))}

    for count, file in enumerate(vcf_file_names):
        data_dict[count]["file_name"] = get_short_name(file_re, file)
        data_dict[count]["data"]["num_calls"] = get_num_calls(file)
        data_dict[count]["data"]["num_pass"] = get_num_pass(file)
    
    json_object = json.dumps(data_dict, indent=4)
    return json_object


#driver code
input_file = "in_files.txt"
vcf_file_re = "PANX[^/]*vep\.vcf\.gz$"
print(get_json_data(input_file, vcf_file_re))
