import os
import sys
import json
import argparse
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from weasyprint import CSS
from section import ( 
    HeaderSection, 
    CasesSection,
    DellySection, 
    Mutect2Section, 
    RSEMSection,
    StarFusionSection,
)

# Report class outlines the structure and order or a report
class Report:
    def __init__(self, workflow_ids, base_db_path):
        self.workflow_ids = workflow_ids
        self.base_db_path = base_db_path
        self.header = HeaderSection()  
        self.sections = [
            CasesSection(),
            DellySection(),
            Mutect2Section(),
            RSEMSection(),
            StarFusionSection(),
        ]

    def load_context(self):
        report_context = {
            "header": self.header.load_context(),  
            "sections": {}
        }
        for section in self.sections:
            section_context = section.load_context(self.workflow_ids, self.base_db_path)
            report_context["sections"][section.name] = section_context

        return report_context

def extract_workflow_ids(data):
    '''
    Extracts workflow ids from a JSON file and returns it as a list. 
    Note: The function first check to see if the data is in the form of a
          dictionary or a list. 

    Parameters
    ----------
    - The dataset in the form of a dictionary or a list. 
    '''
    workflow_ids = []

    if isinstance(data, dict):
        for key, value in data.items():
            if key == 'workflow_id':
                workflow_ids.append(value)
            elif value is not None:
                workflow_ids.extend(extract_workflow_ids(value))
    elif isinstance(data, list):
        for item in data:
            if item is not None:
                workflow_ids.extend(extract_workflow_ids(item))

    return workflow_ids


def makepdf(html, outputfile):
    '''
    (str) -> None
    
    Generates a PDF file from a string of HTML
   
    Parameters
    ----------
    - html (str) String of formated HTML
    - outputfile (str): Name of the output PDF file
    '''
    css_file = os.path.join(os.path.dirname(__file__), './static/css/style.css')

    htmldoc = HTML(string=html, base_url=__file__)
    htmldoc.write_pdf(outputfile, stylesheets=[CSS(css_file)], presentational_hints=True)


def generate_report(input, output, use_stage):
    '''
    (str, str, bool) -> None
    
    Generates a report using data from input file to output file. use_stage indicates if
    data should be pulled from production or stage
      
    Parameters
    ----------
    - input (str): name of input file
    - output (str): name of the output PDF file
    - use_stage: set to True if using data from staging
    '''
    infile = input if input else "ar_input.json"
    outfile = output if output else "Analysis_Report.pdf"

    # Read the input JSON file
    with open(infile, 'r') as file:
        data = json.load(file)

    workflow_ids = extract_workflow_ids(data)
    base_db_path = "/scratch2/groups/gsi/staging/qcetl_v1/" if use_stage else "/scratch2/groups/gsi/production/qcetl_v1/"


    report = Report(workflow_ids, base_db_path)
    report_context = report.load_context()

    # Generate HTML content using Jinja2 templates
    template_dir = os.path.join(os.path.dirname(__file__), './templates')
    environment = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
    results_template = environment.get_template("base.html")

    contents = results_template.render(report_context)

    makepdf(contents, outfile)
    print(f"Created report {outfile}")


if __name__ == "__main__":
    # Create parser for command line args
    parser = argparse.ArgumentParser(
        description="Generates an Analysis Data Release Report"
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

    print(f"Reading input from {args.infile}")
    
    generate_report(input=args.infile, output=args.outfile, use_stage=args.stage)
