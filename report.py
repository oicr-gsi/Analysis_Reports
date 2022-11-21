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

db = f"/scratch2/groups/gsi/staging/qcetl_v1/"


class Proccess:
    columns: {}
    headers: {}
    get_data: {}
    




processes = {
    "mutect2_matched_by_tumour_group" : "mutect2",
    "delly_matched_by_tumour_group" : "delly",
    "sequenza_by_tumour_group": "sequenza",
    "rsem": "rsem",
    "starfusion": "starfusion"
}

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

def get_starfusion_data(sample_ids, data, columns):
    con = sqlite3.connect(db + "analysis_starfusion/latest")
    cur = con.cursor()

    starfusion_data =  {id: {} for id in sample_ids}

    for id in sample_ids:
        wfr = data[id]["starfusion"]["wfrun"]

        for row in cur.execute( # forloop automatically ignores empty SELECTs
            f"""
            select * from analysis_starfusion_analysis_starfusion_1 where "Workflow Run SWID" like '%{wfr}%';
            """):
            starfusion_data[id]["sample_id"] = id
            for column in columns.keys():
                starfusion_data[id][column] = row[columns[column]]
    
    cur.close()
    con.close()
    return starfusion_data


def generate_report():
    data = {}
    context = {}
    with open("IRIS.json") as f:
        data = json.load(f)
    #     with open("input.json", "w") as fo:
    #         json.dump(data, fo, indent=4)

    sample_ids = list(data.keys())
    metadata = {id: {} for id in sample_ids}

    context["metadata_headers"] = {
        "sample_id": "Sample ID",
        "group_id": "Group ID",
        "tissue_type": "Tissue Type",
        "tissue_origin": "Tissue Origin",
        "library_design": "Library Design"
    }

    context["rsem_headers"] = {
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

    metadata_columns = {
        "group_id": 1,
        "tissue_type": 5,
        "tissue_origin": 6,
        "library_design": 2,
    }

    rsem_columns = {
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

    con = sqlite3.connect(db + "analysis_mutect2/latest")
    cur = con.cursor()

    for id in sample_ids:
        wfr = data[id]["mutect2_matched_by_tumor_group"]["wfrun"]

        for row in cur.execute( # forloop automatically ignores empty SELECTs
            f"""
            select * from analysis_mutect2_analysis_mutect2_1 where "Workflow Run SWID" like '%{wfr}%';
            """):
            metadata[id]["sample_id"] = id
            for column in metadata_columns.keys():
                metadata[id][column] = row[metadata_columns[column]]
    context["metadata"] = metadata

    cur.close()
    con.close()
    
    rsem_data = {id: {} for id in sample_ids}


    con = sqlite3.connect(db + "analysis_rsem/latest")
    cur = con.cursor()
    for id in sample_ids:
        wfr = data[id]["rsem"]["wfrun"]

        for row in cur.execute( # forloop automatically ignores empty SELECTs
            f"""
            select * from analysis_rsem_analysis_rsem_1 where "Workflow Run SWID" like '%{wfr}%';
            """):
            rsem_data[id]["sample_id"] = id
            for column in rsem_columns.keys():
                rsem_data[id][column] = round(row[rsem_columns[column]], 2)

    context["rsem_data"] = rsem_data


    cur.close()
    con.close()

    starfusion_columns = {
        "num_records": 8,
    }

    context["starfusion_data"] = get_starfusion_data(sample_ids, data, starfusion_columns)
    context["starfusion_columns"] = starfusion_columns
    context["starfusion_headers"] = {
        "sample_id": "Sample ID",
        "num_records": "Number of Records",
    }

    with open('meta_context.json', 'w', encoding='utf-8') as file:
        json.dump(context, file, ensure_ascii=False, indent=4)

    environment = Environment(loader=FileSystemLoader("templates/"))
    results_template = environment.get_template("metadata.html")

    contents = results_template.render(context)

    with open("meta_report.html", "w", encoding="utf-8") as results:
        results.write(contents)
        print("wrote to sample_report.html")
    
    makepdf(contents, "my_report.pdf")


if __name__ == "__main__":
    generate_report()