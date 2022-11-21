import json
from jinja2 import Environment, FileSystemLoader
import pandas as pd
import matplotlib.pyplot as plt
import sqlite3

import argparse
import os
import numpy as np
import time
import math
import requests
import gzip
import sys
import pathlib
from weasyprint import HTML
from weasyprint import CSS
from typing import Dict, Type, List, Callable, Union, Tuple, Set, Any

db = f"/scratch2/groups/gsi/staging/qcetl_v1/"


class Process:
    data = {}
    sample_ids = {}
    context = {}
    processes = {
        "mutect2": "mutect2_matched_by_tumor_group",
        "delly": "delly_matched_by_tumor_group",
        "sequenza": "sequenza_by_tumor_group",
        "rsem": "rsem",
        "starfusion": "starfusion",
        "metadata": "mutect2_matched_by_tumor_group"
    }

    columns: Dict[str, int]
    headers: Dict[str, str]
    cache_name: str
    table_name: str
    name: str
    blurb: str
    section: str

    def __init__(self):
        pass

    def get_data(self):
        con = sqlite3.connect(db + f"{self.cache_name}/latest")
        cur = con.cursor()

        data =  {id: {} for id in self.sample_ids}

        for id in self.sample_ids:
            wfr = self.data[id][self.processes[self.section]]["wfrun"]

            for row in cur.execute( # forloop automatically ignores empty SELECTs
                f"""
                select * from {self.table_name} where "Workflow Run SWID" like '%{wfr}%';
                """):
                data[id]["sample_id"] = id
                for column in self.columns.keys():
                    data[id][column] = (
                        row[self.columns[column]]
                        if isinstance(row[self.columns[column]], str)
                        else round(row[self.columns[column]], 2)
                    )
        
        cur.close()
        con.close()
        return data

    def load_context(self):
        Process.context[self.section + "_data"] = self.get_data()
        Process.context[self.section + "_headers"] = self.headers
        Process.context[self.section + "_blurb"] = self.blurb


class StarFusion(Process):
    def __init__(self):
        Process.__init__(self)
        self.cache_name = "analysis_starfusion"
        self.table_name = "analysis_starfusion_analysis_starfusion_1"
        self.section= "starfusion"
        self.columns = {
            "num_records": 8,
        }

        self.headers = {
            "sample_id": "Sample ID",
            "num_records": "Number of Records",
        }
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.Curabitur sit amet commodo lorem. Vivamus mollis mauris urna, eget volutpat ex pharetra vel. Suspendisse ultrices commodo feugiat. Morbi aliquam semper ligula eget interdum. Cras vestibulum tempus rutrum.
            """

class RSEM(Process):
    def __init__(self):
        Process.__init__(self)
        self.cache_name = "analysis_rsem"
        self.table_name = "analysis_rsem_analysis_rsem_1"
        self.section = "rsem"
        self.columns = {
            "total": 8,
            "pct_non_zero": 9,
            "Q0": 10,
            "Q0.05": 11,
            "Q0.1": 12,
            "Q0.25": 13,
            "Q0.5": 14,
            "Q0.75": 15,
            "Q0.9": 16,
            "Q0.95": 17,
            "Q1": 18,
        }

        self.headers = {
            "sample_id": "Sample ID",
            "total": "Total",
            "pct_non_zero": "Percent of Non-Zero Records",
            "Q0": "Q0",
            "Q0.05": "Q0.05",
            "Q0.1": "Q0.1",
            "Q0.25": "Q0.25",
            "Q0.5": "Q0.5",
            "Q0.75": "Q0.75",
            "Q0.9": "Q0.9",
            "Q0.95": "Q0.95",
            "Q1": "Q1",
        }

        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.Curabitur sit amet commodo lorem. Vivamus mollis mauris urna, eget volutpat ex pharetra vel. Suspendisse ultrices commodo feugiat. Morbi aliquam semper ligula eget interdum. Cras vestibulum tempus rutrum.
            """

class Mutect2(Process):
    def __init__(self):
        Process.__init__(self)
        self.cache_name = "analysis_mutect2"
        self.table_name = "analysis_mutect2_analysis_mutect2_1"
        self.section = "mutect2"
        self.columns = {
            "num_calls": 8,
            "num_PASS": 9,
            "num_SNPs": 10,
            "num_multi_SNPs": 11,
            "num_indels": 12,
            "titv_ratio": 13,
            "pct_ti": 15,
            "pct_tv": 16,
            "num_MNPs": 14,
        }

        self.headers = {
            "sample_id": "Sample ID",
            "num_calls": "Number of Calls",
            "num_PASS": "Number of PASS calls",
            "num_SNPs": "Number of SNPs",
            "num_multi_SNPs": "Number of multi-SNPs",
            "num_indels": "Number of indels",
            "titv_ratio": "TI/TV ratio",
            "pct_ti": "Percent TI",
            "pct_tv": "Percent of TV",
            "num_MNPs": "Number of MNPs",

        }
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.Curabitur sit amet commodo lorem. Vivamus mollis mauris urna, eget volutpat ex pharetra vel. Suspendisse ultrices commodo feugiat. Morbi aliquam semper ligula eget interdum. Cras vestibulum tempus rutrum.
            """

class Metadata(Process):
    def __init__(self):
        Process.__init__(self)
        self.cache_name = "analysis_mutect2"
        self.table_name = "analysis_mutect2_analysis_mutect2_1"
        self.section = "metadata"
        self.columns = {
            "group_id": 1,
            "tissue_type": 5,
            "tissue_origin": 6,
            "library_design": 2,
        }
        self.headers = {
            "sample_id": "Sample ID",
            "group_id": "Group ID",
            "tissue_type": "Tissue Type",
            "tissue_origin": "Tissue Origin",
            "library_design": "Library Design"
        }
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem. Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat. Morbi aliquam semper ligula eget interdum.
            Cras vestibulum tempus rutrum.
            """




def makepdf(html, outputfile):
    """
    (str) -> None
    
    Generates a PDF file from a string of HTML
   
    Parameters
    ----------
    - html (str) String of formated HTML
    - outputfile (str): Name of the output PDF file
    """
    
    htmldoc = HTML(string=html, base_url=__file__)
    htmldoc.write_pdf(outputfile, stylesheets=[CSS('./static/css/style.css')], presentational_hints=True)

def generate_report():
    with open("IRIS.json") as f:
        Process.data = json.load(f)
    #     with open("input.json", "w") as fo:
    #         json.dump(data, fo, indent=4)

    Process.sample_ids = list(Process.data.keys())

    classes = {
        "mutect2": Mutect2(),
        "rsem": RSEM(),
        "starfusion": StarFusion(),
        "metadata": Metadata(),
    }

    for p in classes.keys():
        classes[p].load_context()

    
    with open('meta_context.json', 'w', encoding='utf-8') as file:
        json.dump(Process.context, file, ensure_ascii=False, indent=4)

    environment = Environment(loader=FileSystemLoader("templates/"))
    results_template = environment.get_template("metadata.html")

    contents = results_template.render(Process.context)

    with open("meta_report.html", "w", encoding="utf-8") as results:
        results.write(contents)
        print("wrote to sample_report.html")
    
    makepdf(contents, "my_report.pdf")


if __name__ == "__main__":
    generate_report()
