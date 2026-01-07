from tables import (
    CasesTable,
    WGLaneLevelTable,
    WTLaneLevelTable,
    WGCallReadyTable,
    WTCallReadyTable,
    Mutect2Table,
    DellyTable,
    PurpleTable,
    MrdTable,
    RSEMTable,
    StarFusionTable,
)
from typing import List, Any
from datetime import date
import json

def load_section_config(filepath='./templates/blurb.json'):
    with open(filepath, "r") as f:
        return json.load(f)

# Section class defines a section of the report
class Section:
    title: str  # title of section
    blurb: str  # blurb of section
    tables: List[Any]  # tables of section
    name: str  # name of section, used as key in context for jinja2 templating
    assay: str  # assay type, e.g., 'WGS', 'RNA-Seq', etc.

    def _get_blurb(self):
        config = load_section_config()
        return config.get(self.assay, {}).get(self.name, {}).get("blurb", "")

    def load_context(self, cases_data, workflow_ids, assay=None):
        '''
        Returns a dict containing the context of the section, including its tables.

        Parameters:
        - cases_data: DataFrame containing the all the information related to the Donors like Sample ID, Tissue Type, External ID, etc.

        Returns:
        - A dictionary with the section's context, including title, blurb, and tables
        '''
        if assay:
            self.assay = assay
            self.blurb = self._get_blurb()

        context = {
            "title": self.title,
            "blurb": self.blurb,
            "tables": {},
        }
        has_data = False 

        # For each table in the section, load its context
        for tcount, table in enumerate(self.tables):
            table_context = table.load_context(cases_data, workflow_ids)  
            
            if table_context["data"]:
                context["tables"][tcount] = table_context 
                has_data = True

                if 'plots' in table_context:
                    context['tables'][tcount]['plots'] = table_context['plots']

        if not has_data:
            return None

        return context

# HeaderSection class defines the header of the report
class HeaderSection(Section):
    def __init__(self, project, assay):
        self.project = project
        self.assay = assay
        self.title = project 
        self.name = "header"
        self.blurb = self._get_blurb()

    def load_context(self):
        '''
        Returns the context of the header section

        Parameters: None
        Returns:
        - A dictionary containing the title, current date, and the blurb for the header section
        '''
        context = {
            "project": self.project,
            "title": self.title,
            "date": date.today().strftime("%Y-%m-%d"), 
            "blurb": self.blurb,
        }
        return context

#CasesSection class defines the section for cases
class CasesSection(Section):
    def __init__(self, assay):
        self.title = "Cases"
        self.assay = assay
        self.name = "cases"
        self.blurb = self._get_blurb()
        self.tables = [
            CasesTable(),
        ]

# RawSeqDataSection class defines the section for lane level data
class RawSeqDataSection(Section):
    def __init__(self, assay):
        self.title = "Raw Sequence Data"
        self.assay = assay
        self.name = "raw_seq_data"
        self.blurb = self._get_blurb()
        self.tables = [
            WGLaneLevelTable(),
            WTLaneLevelTable(),
        ]

#CallReadyAlignmentsSection class defines the section for call ready alignments
class CallReadyAlignmentsSection(Section):
    def __init__(self, assay):
        self.title = "Aligned Sequence Data"
        self.assay = assay
        self.name = "call_ready"
        self.blurb = self._get_blurb()
        self.tables = [
            WGCallReadyTable(),
            WTCallReadyTable(),
        ]
    
# Mutect2Section class defines the section for the mutect2 workflow
class Mutect2Section(Section):
    def __init__(self, assay):
        self.title = "Mutations"
        self.assay = assay
        self.name = "mutect2"
        self.blurb = self._get_blurb()
        self.tables = [
            Mutect2Table(),
        ]

# DellySection class defines the section for delly workflow
class DellySection(Section):
    def __init__(self, assay):
        self.title = "Genomic Structural Variants"
        self.assay = assay
        self.name = "delly"
        self.blurb = self._get_blurb()
        self.tables = [
            DellyTable(),
        ]

# PurpleSection class defines the section for purple workflow
class PurpleSection(Section):
    def __init__(self, assay):
        self.title = "Purity and Ploidy Estimation"
        self.assay = assay
        self.name = "purple"
        self.blurb = self._get_blurb()
        self.tables = [
            PurpleTable(),
        ]

# MrdSection class defines the section for mrdetect workflow
class MrdSection(Section):
    def __init__(self, assay):
        self.title = "Minimal Residual Disease Detection"
        self.assay = assay
        self.name = "mrdetect"
        self.blurb = self._get_blurb()
        self.tables = [
            MrdTable(),
        ]

# RSEMSection class defines the section for RSEM workflow
class RSEMSection(Section):
    def __init__(self, assay):
        self.title = "Gene Expression"
        self.assay = assay
        self.name = "rsem"
        self.blurb = self._get_blurb()
        self.tables = [
            RSEMTable(),
        ]

#StarFusionSection class defines the section for StarFusion
class StarFusionSection(Section):
    def __init__(self, assay):
        self.title = "Gene Fusions"
        self.assay = assay
        self.name = "starfusion"
        self.blurb = self._get_blurb()
        self.tables = [
            StarFusionTable(),
        ]