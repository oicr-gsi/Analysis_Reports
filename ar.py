import os
import argparse
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML
from weasyprint import CSS
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

# Report class outlines the structure and order or a report
class Report:
    def __init__(self, project, release):
        self.context = {"sections":{}, "header": {}} #context to be passed to jinja2 templating
        self.header = HeaderSection(project, release) 
        self.sections = [ #sections will appear in the following order in the report
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
        """
        None -> None
    
        Get the data and load it into a context dict for jinja2 to generate html
        """
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

    css_file = os.path.join(os.path.dirname(__file__), './static/css/style.css')
    
    htmldoc = HTML(string=html, base_url=__file__)
    htmldoc.write_pdf(outputfile, stylesheets=[CSS(css_file)], presentational_hints=True)


def generate_report(input, output, use_stage):
    """
    (str, str, bool) -> None
    
    Generates a report using data from input file to output file. use_stage indicates if
    data should be pulled from production or stage
      
    Parameters
    ----------
    - input (str): name of input file
    - output (str): name of the output PDF file
    - use_stage: set to True if using data from staging
    """
    infile = input if input else "ar_input.json"
    outfile = output if output else "Analysis_Report.pdf"
    table = Table(infile, use_stage) #initializing table data
    report = Report(table.project, table.release) #initialize report structure

    report.load_context()

    # used to debug issues with context
    # with open('ar_context.json', 'w', encoding='utf-8') as file:
    #     json.dump(report.context, file, ensure_ascii=False, indent=4)

    template_dir = os.path.join(os.path.dirname(__file__), './templates')
    environment = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
    results_template = environment.get_template("base.html")

    contents = results_template.render(report.context)
    
    makepdf(contents, outfile)
    print(f"Created report {outfile}")


if __name__ == "__main__":
    #create parser for command line args
    parser = argparse.ArgumentParser(
        description="Generates a Analysis Data Release Report"
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
