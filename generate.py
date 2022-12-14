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
from section import (
    RSEMSection,
    SequenzaSection,
    DellySection,
    StarFusionSection,
    Mutect2Section,
    RawSeqDataSection,
    CasesSection,
    CallReadyAlignmentsSection,
    HeaderSection,
)

from tables import Table


db = f"/.mounts/labs/gsiprojects/gsi/gsiusers/jqian/gsi-qc-etl/gsiqcetl/"


class Report:
    def __init__(self, project, release):
        self.data = {}
        self.sample_ids = []
        self.context = {"sections":{}, "header": {}}
        self.header = HeaderSection(project, release)
        self.sections = [
            CasesSection(),
            RawSeqDataSection(),
            CallReadyAlignmentsSection(),
            Mutect2Section(),
            SequenzaSection(),
            DellySection(),
            RSEMSection(),
            StarFusionSection(),
        ]

    def load_context(self):
        self.context["header"] = self.header.load_context()
        for section in self.sections:
            self.context["sections"][section.name] = section.load_context()

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

def generate_report(input, output, use_stage):
    infile = input if input else "IRIS.json"
    outfile = output if output else "Analysis_Report.pdf"
    table = Table(infile, use_stage) #initializing table data
    print(table.data.keys())
    report = Report(table.project, table.release)

    report.load_context()

    with open('meta_context.json', 'w', encoding='utf-8') as file:
        json.dump(report.context, file, ensure_ascii=False, indent=4)
    
    # with open("IRIS-3.json") as f:
    #     with open("input.json", "w") as input:
    #         json.dump(json.load(f), input, indent=4)

    environment = Environment(loader=FileSystemLoader("templates/"))
    results_template = environment.get_template("metadata.html")

    contents = results_template.render(report.context)

    with open("meta_report.html", "w", encoding="utf-8") as results:
        results.write(contents)
        print("wrote to sample_report.html")
    
    makepdf(contents, outfile)
    print(f"created report {outfile}")


if __name__ == "__main__":
    print("starting")

    parser = argparse.ArgumentParser(
        description="Generates a MOH Data Release Report"
    )

    parser.add_argument(
        '-i',
        '--infile',
        type=str,
        required=False,
        help="Name of the input file. Default looks for IRIS.json"
    )
    parser.add_argument(
        '-o',
        '--outfile',
        type=str,
        required=False,
        help="Name of output file. Default names pdf Analysis_Report.pdf"
    )
    parser.add_argument(
        '--stage',
        '--staging',
        action="store_true",
        help="Use qcetl data from stage",
    )
    args = parser.parse_args()

    print(args)

    generate_report(input=args.infile, output=args.outfile, use_stage=args.stage)
