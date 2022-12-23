from tables import (
    RSEMTable, 
    SequenzaTable,
    DellyTable,
    StarFusionTable,
    Mutect2Table,
    CasesTable,
    WTLaneLevelTable,
    WGCallReadyTable,
    WTCallReadyTable,
    WGLaneLevelTable,
)
from typing import List, Any
from datetime import date

# Section class defines a section of the report
class Section:
    title: str  # title of section
    blurb: str  # blurb of section
    tables: List[Any] # tables of section
    name: str   # name of section, used as key in context for jinja2 templating

    def load_context(self):
        """
        None -> dict

        Returns a dict of the context of the section and its tables

        """
        context = {
            "title": self.title,
            "blurb": self.blurb,
            "tables": {},
            "plots": {},
        }
        for tcount, table in enumerate(self.tables):
            context["tables"][tcount] = table.load_context()
            context["tables"][tcount]["plots"] = {}
            for pcount, (column, plot) in enumerate(table.plots.items()):
                context["tables"][tcount]["plots"][pcount] = plot.load_context(f"{table.process}_{column}")
        return context

#CallReadyAlignmentsSection class defines the section for call ready alignments
class CallReadyAlignmentsSection(Section):
    def __init__(self):
        self.title = "Call Ready Alignments"
        self.blurb = """
        All data from each sample is merged and processed to a call ready state.
        """
        self.name = "call_ready"
        self.tables = [
            WGCallReadyTable(),
            WTCallReadyTable(),
        ]

#CasesSection class defines the section for cases
class CasesSection(Section):
    def __init__(self):
        self.title = "Cases"
        self.blurb = """
        The following cases are included in this release.
        """
        self.name = "cases"
        self.tables = [
            CasesTable(),
        ]

# DellySection class defines the section for delly workflow
class DellySection(Section):
    def __init__(self):
        self.title = "Genomic Structural Variants"
        self.blurb = """
        Summary metrics for structural variants  generated from the WG Tumour/Normal pairs.
        Calls are generated with delly.
        """
        self.name = "delly"
        self.tables = [
            DellyTable(),
        ]

# HeaderSection class defines the header of the report
class HeaderSection(Section):
    def __init__(self, project, release):
        self.project = project
        self.release = release
        self.title = "Marathon of Hope"
        self.name = "header"
        self.blurb = """
        The data release report summarizes a variety of metrics generated from our 
        quality control and analysis workflows. All MOH cases are processed through 
        our WGTS (Whole Genome, Transcriptome) sequencing and analysis pipelines, and
        include a single tumour sample with a matched normal.
        """
        
    def load_context(self):
        """
        None -> dict

        Returns a dict of the context of the section

        """
        context = {
            "project": self.project,
            "date": date.today().strftime("%Y-%m-%d"),
            "title": self.title,
            "blurb": self.blurb,
            "release": self.release,
        }
        return context
    
#Mutect2Section class defines the section for the mutect2 workflow
class Mutect2Section(Section):
    def __init__(self):
        self.title = "Mutations"
        self.blurb = """
        Summary metrics for somatic mutations (snvs + indels) generated from the WG Tumour/Normal pairs.
        Calls are generated with mutect2, and annotated with variant effect predictor.
        """
        self.name = "mutect2"
        self.tables = [
            Mutect2Table(),
        ]

#RawSeqDataSection class defines the section for raw sequence data
class RawSeqDataSection(Section):
    def __init__(self):
        self.title = "Raw Sequence Data"
        self.blurb = """
        Samples were sequenced on one or more sequencing runs.
        """
        self.name = "raw_seq_data"
        self.tables = [
            WGLaneLevelTable(),
            WTLaneLevelTable(),
        ]

#RSEMSection class defines the section for RSEM workflow
class RSEMSection(Section):
    def __init__(self):
        self.title = "Gene Expression"
        self.blurb = """
        Summary metrics for normalized expression (TPM, transcripts per million) for genes in gencode release 31 
        (https://www.gencodegenes.org/human/release_31.html). \n
        Expression values are generated with RSEM.
        """
        self.name = "rsem"
        self.tables = [
            RSEMTable(),
        ]

#SequenzaSection class defines the section for Sequenza
class SequenzaSection(Section):
    def __init__(self):
        self.title = "Copy Number Alterations"
        self.blurb = """
        Summary metrics for somatic copy number alterations generated from the WG Tumour/Normal pairs.
        \nInitial calls are generated with varscan, then processed with Sequenza.
        Multiple copy number profiles are generated and released over a range of tunable gamma settings.
        Metrics displayed are for gamma = 500.
        """
        self.name = "sequenza"
        self.tables = [
            SequenzaTable(),
        ]

#StarFusionSection class defines the section for StarFusion
class StarFusionSection(Section):
    def __init__(self):
        self.title = "Gene Fusions"
        self.blurb = """
        Summary metrics for identified gene fusions.
        Fusions are detected with STAR-fusion.
        """
        self.name = "starfusion"
        self.tables = [
            StarFusionTable(),
        ]
