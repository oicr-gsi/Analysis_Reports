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
    MetadataSection,
    CallReadyAlignmentsSection,
    HeaderSection,
)

from tables import Table


db = f"/.mounts/labs/gsiprojects/gsi/gsiusers/jqian/gsi-qc-etl/gsiqcetl/"


class Report:
    data = {}
    sample_ids = []
    context = {"sections":{}, "header": {}}
    header = HeaderSection()
    sections = [
        MetadataSection(),
        RawSeqDataSection(),
        CallReadyAlignmentsSection(),
        Mutect2Section(),
        SequenzaSection(),
        DellySection(),
        RSEMSection(),
        StarFusionSection(),
    ]

    def load_context(self):
        Report.context["header"] = self.header.load_context()
        for section in self.sections:
            Report.context["sections"][section.name] = section.load_context()

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


def thousands_separator(value):
    return f"{value:,.2f}"

def generate_report():
    report = Report()
    table = Table("IRIS-3.json") #initializing table data
    report.load_context()

    with open('meta_context.json', 'w', encoding='utf-8') as file:
        json.dump(report.context, file, ensure_ascii=False, indent=4)
    
    # with open("IRIS-3.json") as f:
    #     with open("input.json", "w") as input:
    #         json.dump(json.load(f), input, indent=4)

    environment = Environment(loader=FileSystemLoader("templates/"))
    # environment.filters["thousands_separator"] = thousands_separator
    results_template = environment.get_template("metadata.html")

    contents = results_template.render(report.context)

    # with open("meta_report.html", "w", encoding="utf-8") as results:
    #     results.write(contents)
    #     print("wrote to sample_report.html")
    
    makepdf(contents, "my_report.pdf")
    print("created report my_report.pdf")


if __name__ == "__main__":
    print("starting")

    generate_report()
