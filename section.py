from tables import (
    CasesTable,
    WGLaneLevelTable,
    WTLaneLevelTable,
    WGCallReadyTable,
    WTCallReadyTable,
    Mutect2Table,
    DellyTable,
    PurpleTable,
    RSEMTable,
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

    def load_context(self, cases_data, workflow_ids):
        '''
        Returns a dict containing the context of the section, including its tables.

        Parameters:
        - cases_data: DataFrame containing the all the information related to the Donors like Sample ID, Tissue Type, External ID, etc.

        Returns:
        - A dictionary with the section's context, including title, blurb, and tables
        '''
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
    def __init__(self, project):
        self.project = project
        self.title = project 
        self.name = "header"
        self.blurb = '''
        <p>This data release report summarizes metrics generated from OICR's 
        quality control and analysis workflows. All cases in this report are processed through 
        the Whole Genome and Transcriptome (WGTS) sequencing and analysis pipeline. </p>
        
        <p><strong>Publications resulting from this data are requested to include the following acknowledgement statement:</strong></p>
        
        <p style="font-style: italic; color: grey;">
        This study was conducted with the support of the Ontario Institute for Cancer Research's Genomics Program 
        (genomics.oicr.on.ca) through funding provided by the Government of Ontario.
        </p>
        '''

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
    def __init__(self):
        self.title = "Cases"
        self.blurb = '''
        The following cases are included in this release. Each case 
        includes two samples, a tumor with a matched normal/reference.
        Whole genome (WG) libraries are generated from the tumour/normal pair. 
        Whole transcriptome (WT) libraries are generated from only the tumour sample. 
        '''
        self.name = "cases"
        self.tables = [
            CasesTable(),
        ]

    def load_context(self, cases_data, workflow_ids):
        context = super().load_context(cases_data, workflow_ids)
        return context

# RawSeqDataSection class defines the section for lane level data
class RawSeqDataSection(Section):
    def __init__(self):
        self.title = "Raw Sequence Data"
        self.blurb = """
        All libraries are sequenced on the illumina Novaseq X Plus platform to generate demultiplexed FASTQ files.
        """
        self.name = "raw_seq_data"
        self.tables = [
            WGLaneLevelTable(),
            WTLaneLevelTable(),
        ]
    def load_context(self, cases_data, workflow_ids):
        context = super().load_context(cases_data, workflow_ids)
        return context

#CallReadyAlignmentsSection class defines the section for call ready alignments
class CallReadyAlignmentsSection(Section):
    def __init__(self):
        self.title = "Call Ready Alignments"
        self.blurb = """
        Raw sequence data (fastq) is trimmed to remove adapter sequence and aligned to the hg38 genomic reference.
        Each sample may have multiple bam files depending on how many lanes of sequence data has been generated. 
        The lane level alignments are merged and processed to a call ready state.
        """
        self.name = "call_ready"
        self.tables = [
            WGCallReadyTable(),
            WTCallReadyTable(),
        ]
    def load_context(self, cases_data, workflow_ids):
        context = super().load_context(cases_data, workflow_ids)
        return context

# Mutect2Section class defines the section for the mutect2 workflow
class Mutect2Section(Section):
    def __init__(self):
        self.title = "Mutations"
        self.blurb = '''
        Call ready alignments from a tumor/normal pair are used to generate somatic variants (snvs + indels). 
        Variants are generated with mutect2, and annotated with variant effect predictor.
        '''
        self.name = "mutect2"
        self.tables = [
            Mutect2Table(),
        ]
    
    def load_context(self, cases_data, workflow_ids):
        context = super().load_context(cases_data, workflow_ids)
        return context

# DellySection class defines the section for delly workflow
class DellySection(Section):
    def __init__(self):
        self.title = "Genomic Structural Variants"
        self.blurb = '''
        Call ready alignments from a tumour/normal pair are analyzed with delly 
        to generate somatic structural variants(deletions, duplications, inversion, insertions, translocations).
        '''
        self.name = "delly"
        self.tables = [
            DellyTable(),
        ]
    
    def load_context(self, cases_data, workflow_ids):
        context = super().load_context(cases_data, workflow_ids)
        return context

# PurpleSection class defines the section for purple workflow
class PurpleSection(Section):
    def __init__(self):
        self.title = "Purity and Ploidy Estimation"
        self.blurb = '''
        Call ready alignments, somatic and structural variants are analyzed with purple to estimate the purity and copy number
        profile of the tumour sample.
        '''
        self.name = "purple"
        self.tables = [
            PurpleTable(),
        ]
    
    def load_context(self, cases_data, workflow_ids):
        context = super().load_context(cases_data, workflow_ids)
        return context

# RSEMSection class defines the section for RSEM workflow
class RSEMSection(Section):
    def __init__(self):
        self.title = "Gene Expression"
        self.blurb = '''
        Aligned whole transcriptome data is analyzed with RSEM to generate expression calls.
        '''
        self.name = "rsem"
        self.tables = [
            RSEMTable(),
        ]
    
    def load_context(self, cases_data, workflow_ids):
        context = super().load_context(cases_data, workflow_ids)
        return context

#StarFusionSection class defines the section for StarFusion
class StarFusionSection(Section):
    def __init__(self):
        self.title = "Gene Fusions"
        self.blurb = '''
        Aligned whole transcriptome data is analyzed with STAR-fusion and arriba to generate gene fusion calls.
        '''
        self.name = "starfusion"
        self.tables = [
            StarFusionTable(),
        ]

    def load_context(self, cases_data, workflow_ids):
        context = super().load_context(cases_data, workflow_ids)
        return context