import os
import sys
import json
import argparse
import logging
import re
import gzip
import shutil
import pandas as pd
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from weasyprint import CSS
from section import ( 
    HeaderSection, 
    CasesSection,
    RawSeqDataSection,
    CallReadyAlignmentsSection,
    Mutect2Section, 
    DellySection, 
    PurpleSection,
    MrdSection,
    RSEMSection,
    StarFusionSection,
)

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Report class outlines the structure and order or a report
class Report:
    def __init__(self, cases_data, project, workflow_ids, assay):
        logger.debug("Initializing Report object.")
        self.cases_data = cases_data
        self.workflow_ids = workflow_ids
        self.assay = assay
        self.header = HeaderSection(project, assay) 
        self.sections = self._get_sections() 
    
    def _get_sections(self):
        '''
        Returns a list of sections to be included in the report based on the assay type

        Parameters: None
        Returns:
        - A list of section objects to be included in the report
        '''
        if self.assay == "WGTS":
            return [
                CasesSection(self.assay),
                #RawSeqDataSection(self.assay),
                CallReadyAlignmentsSection(self.assay),
                Mutect2Section(self.assay),
                DellySection(self.assay),
                PurpleSection(self.assay),
                RSEMSection(self.assay),
                StarFusionSection(self.assay),
            ]
        elif self.assay == "WGS":
            return [
                CasesSection(self.assay),
                #RawSeqDataSection(self.assay),
                CallReadyAlignmentsSection(self.assay),
                Mutect2Section(self.assay),
                DellySection(self.assay),
                PurpleSection(self.assay),
            ]
        elif self.assay == "pWGS+WGS":
            return [
                CasesSection(self.assay),
                #RawSeqDataSection(self.assay),
                CallReadyAlignmentsSection(self.assay),
                Mutect2Section(self.assay),
                DellySection(self.assay),
                PurpleSection(self.assay),
                MrdSection(self.assay),
            ]
        elif self.assay == "pWGS+WGTS":
            return [
                CasesSection(self.assay),
                #RawSeqDataSection(self.assay),
                CallReadyAlignmentsSection(self.assay),
                Mutect2Section(self.assay),
                DellySection(self.assay),
                PurpleSection(self.assay),
                MrdSection(self.assay),
                RSEMSection(self.assay),
                StarFusionSection(self.assay),
            ]
        else:
            logger.error(f"Unsupported assay type: {self.assay}")
            raise ValueError(f"Assay type {self.assay} not recognized. Supported assays are WGTS, WGS and pWGS.")   

    def load_context(self):
        report_context = {
            "header": self.header.load_context(),  
            "sections": {}
        }
        for section in self.sections:
            section_context = section.load_context(self.cases_data, self.workflow_ids, self.assay)
            if section_context:
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

def get_fp_records(provenance, workflow_ids):
    '''
    Extracts records from the file_provenance_path that match the workflow_ids.
    Parameters
    ----------
    - file_provenance_path (str): Path to the file provenance file.
    - workflow_ids (list): List of workflow IDs to filter by.
    Returns
    -------
    - list: List of records matching the workflow IDs.
    '''
    records = []

    if is_gzipped(provenance):
        infile = gzip.open(provenance, 'rt', errors='ignore')
    else:
        infile = open(provenance, 'r', errors='ignore')

    header = infile.readline().strip().split('\t')
    if 'Workflow Run SWID' not in header:
        infile.close()
        logger.error("'Workflow Run SWID' column not found in header.")
        raise ValueError("'Workflow Run SWID' column not found in header.")

    swid_idx = header.index('Workflow Run SWID')
    workflow_ids_set = set(workflow_ids)
    for line in infile:
        row = line.strip().split('\t')
        if len(row) > swid_idx and row[swid_idx] in workflow_ids_set:
            records.append(row)

    infile.close()
    return header, records

def parse_fp_records(header, records):
    '''
    Parses the file provenance records and returns a DataFrame and the project name.
    Parameters
    ----------
    - header (list): List of column names.
    - records (list): List of fp rows to parse.
    Returns
    -------
    - DataFrame: Parsed DataFrame with relevant columns.
    '''
    col_indices = {
        'Workflow Run ID': header.index('Workflow Run SWID'),
        'Root Sample Name': header.index('Root Sample Name'),
        'Sample Attributes': header.index('Sample Attributes'),
        'LIMS ID': header.index('LIMS ID'),
        'Study Title': header.index('Study Title'),
        'Sample Name': header.index('Sample Name'),
    }

    case = []
    lims_dict = {}
    study_titles = set()

    for row in records:
        workflow_swid = row[col_indices['Workflow Run ID']]
        donor = row[col_indices['Root Sample Name']]
        sample_attributes = row[col_indices['Sample Attributes']]
        lims_id = row[col_indices['LIMS ID']]
        study_title = row[col_indices['Study Title']]
        sample_name = row[col_indices['Sample Name']]
        study_titles.add(study_title)

        if workflow_swid not in lims_dict:
            lims_dict[workflow_swid] = set()
        lims_dict[workflow_swid].add(lims_id)

        metadata_dict = {}
        for item in sample_attributes.split(';'):
            if '=' in item:
                key, value = item.split('=', 1)
                metadata_dict[key] = value

        case.append({
            'Donor': donor,
            'Group ID': metadata_dict.get('geo_group_id'),
            'Library Type': metadata_dict.get('geo_library_source_template_type'),
            'Tissue Type': metadata_dict.get('geo_tissue_type'),
            'Tissue Origin': metadata_dict.get('geo_tissue_origin'),
            'Tissue Preparation': metadata_dict.get('geo_tissue_preparation'),
            'External ID': metadata_dict.get('geo_external_name'),
            'Workflow Run ID': workflow_swid,
            'Sample Name': sample_name,
        })

    cases = pd.DataFrame(case).drop_duplicates()
    cases['LIMS ID'] = cases['Workflow Run ID'].map(
        lambda swid: ','.join(lims_dict.get(swid, [])))
    cases['SampleID'] = cases.apply(
        lambda row: f"{row['Donor']}_{row['Tissue Origin']}_{row['Tissue Type']}_{row['Library Type']}_{row['Group ID']}",
        axis=1)

    project = study_titles.pop() if len(study_titles) == 1 else ""

    return cases, project


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

    try:
        htmldoc = HTML(string=html, base_url=__file__)
        htmldoc.write_pdf(outputfile, stylesheets=[CSS(css_file)], presentational_hints=True)
    except Exception as e:
        logger.error(f"Failed to generate PDF: {e}")
        raise


def generate_report(input, output, temp_dir, assay):
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

    provenance = "/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz"
    workflow_ids = extract_workflow_ids(data)
    header, records = get_fp_records(provenance, workflow_ids)
    cases_data, project = parse_fp_records(header, records)

    report = Report(cases_data, project, workflow_ids, assay)
    report_context = report.load_context()

    # Generate HTML content using Jinja2 templates
    template_dir = os.path.join(os.path.dirname(__file__), './templates')
    environment = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
    results_template = environment.get_template("base.html")

    contents = results_template.render(report_context)

    makepdf(contents, outfile)

    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
        logger.info(f"Temporary directory {temp_dir} cleaned up.")

    logger.info(f"Report generation complete: {outfile}")

def is_gzipped(file):
    '''
    (str) -> bool

    Return True if file is gzipped

    Parameters
    ----------
    - file (str): Path to file
    '''
    # open file in rb mode
    infile = open(file, 'rb')
    header = infile.readline()
    infile.close()
    if header.startswith(b'\x1f\x8b\x08'):
        return True
    else:
        return False

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
    '-a',
    '--assay',
    type=str,
    required=True,
    help="Specify assay type (e.g., WGTS, WGS, pWGS)"
    )
    args = parser.parse_args()

    logger.info(f"Reading input from {args.infile}")
    temp_dir = 'temp'
    generate_report(input=args.infile, output=args.outfile, temp_dir=temp_dir, assay=args.assay)
