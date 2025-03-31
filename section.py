from tables import (
    CasesTable,
    RSEMTable,
    DellyTable,
    Mutect2Table,
    StarFusionTable,
)
from typing import List, Any
from datetime import date

# Section class defines a section of the report
class Section:
    title: str  # title of section
    blurb: str  # blurb of section
    tables: List[Any]  # tables of section
    name: str  # name of section, used as key in context for jinja2 templating

    def load_context(self, workflow_ids, base_db_path):
        '''
        Returns a dict containing the context of the section, including its tables.

        Parameters:
        - workflow_ids: List of workflow IDs extracted from the input data
        - base_db_path: Path to the base directory for databases

        Returns:
        - A dictionary with the section's context, including title, blurb, and tables
        '''
        context = {
            "title": self.title,
            "blurb": self.blurb,
            "tables": {},
        }
        # For each table in the section, load its context
        for tcount, table in enumerate(self.tables):
            table_context = table.load_context(workflow_ids, base_db_path)  # Pass workflow_ids and base_db_path to table
            context["tables"][tcount] = table_context  # Store the table context in the section context
        return context

#CasesSection class defines the section for cases
class CasesSection(Section):
    def __init__(self):
        self.title = "Donors"
        self.blurb = '''
        The following donors are included in this release.
        '''
        self.name = "cases"
        self.tables = [
            CasesTable(),
        ]

    def load_context(self, workflow_ids, base_db_path):
        context = super().load_context(workflow_ids, base_db_path)
        return context

# DellySection class defines the section for delly workflow
class DellySection(Section):
    def __init__(self):
        self.title = "Genomic Structural Variants"
        self.blurb = '''
        Summary metrics for structural variants generated from the WG Tumour/Normal pairs.
        Calls are generated with delly.
        '''
        self.name = "delly"
        self.tables = [
            DellyTable(),
        ]
    
    def load_context(self, workflow_ids, base_db_path):
        context = super().load_context(workflow_ids, base_db_path)
        return context

# HeaderSection class defines the header of the report
class HeaderSection(Section):
    def __init__(self):
        self.title = "IRIS" 
        self.name = "header"
        self.blurb = '''
        The data release report summarizes a variety of metrics generated from our 
        quality control and analysis workflows. All IRIS cases are processed through 
        our WGTS (Whole Genome, Transcriptome) sequencing and analysis pipelines, and
        include a single tumour sample with a matched normal.
        '''

    def load_context(self):
        '''
        Returns the context of the header section

        Parameters: None
        Returns:
        - A dictionary containing the title, current date, and the blurb for the header section
        '''
        context = {
            "title": self.title,
            "date": date.today().strftime("%Y-%m-%d"), 
            "blurb": self.blurb,
        }
        return context

# Mutect2Section class defines the section for the mutect2 workflow
class Mutect2Section(Section):
    def __init__(self):
        self.title = "Mutations"
        self.blurb = '''
        Summary metrics for somatic mutations (snvs + indels) generated from the WG Tumour/Normal pairs.
        Calls are generated with mutect2, and annotated with variant effect predictor.
        '''
        self.name = "mutect2"
        self.tables = [
            Mutect2Table(),
        ]
    
    def load_context(self, workflow_ids, base_db_path):
        context = super().load_context(workflow_ids, base_db_path)
        return context

# RSEMSection class defines the section for RSEM workflow
class RSEMSection(Section):
    def __init__(self):
        self.title = "Gene Expression"
        self.blurb = '''
        Summary metrics for normalized expression (TPM, transcripts per million) for genes in gencode release 31 
        (https://www.gencodegenes.org/human/release_31.html). \n
        Expression values are generated with RSEM.
        '''
        self.name = "rsem"
        self.tables = [
            RSEMTable(),
        ]
    
    def load_context(self, workflow_ids, base_db_path):
        context = super().load_context(workflow_ids, base_db_path)
        return context

#StarFusionSection class defines the section for StarFusion
class StarFusionSection(Section):
    def __init__(self):
        self.title = "Gene Fusions"
        self.blurb = '''
        Summary metrics for identified gene fusions.
        Fusions are detected with STAR-fusion.
        '''
        self.name = "starfusion"
        self.tables = [
            StarFusionTable(),
        ]

    def load_context(self, workflow_ids, base_db_path):
        context = super().load_context(workflow_ids, base_db_path)
        return context