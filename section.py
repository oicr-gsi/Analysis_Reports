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
from typing import Dict, Type, List, Callable, Union, Tuple, Set, Any
from datetime import date


class Section:
    title: str
    blurb: str
    tables: List[Any]
    name: str

    def __init__(self):
        pass

    def load_context(self):
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

class CallReadyAlignmentsSection(Section):
    def __init__(self):
        self.title = "Call Ready Alignments"
        self.blurb = """
        All data from each sample is merged and processed to a call ready state.
        """
        self.name = "call_ready"
        self.tables = [
            WGCallReadyTable(
                "Whole Genome, Normal Samples",
                "WG",
                "Normal",
            ),
            WGCallReadyTable(
                "Whole Genome, Tumour Samples",
                "WG",
                "Tumour",
            ),
            WTCallReadyTable(),
        ]


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

class HeaderSection(Section):
    def __init__(self):
        self.project = "IRIS"
        self.title = "Marathon of Hope"
        self.name = "header"
        
    def load_context(self):
        context = {
            "project": self.project,
            "date": date.today().strftime("%Y-%m-%d"),
            "title": self.title
        }
        return context
        

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
    
class RSEMSection(Section):
    def __init__(self):
        self.title = "Gene Expression"
        self.blurb = """
        Summary metrics for normalized expression (TPM, gencode release 31)
        \nExpression values are generated with RSEM
        """
        self.name = "rsem"
        self.tables = [
            RSEMTable(),
        ]
    
class SequenzaSection(Section):
    def __init__(self):
        self.title = "Copy Number Alterations"
        self.blurb = """
        Summary metrics for somatic copy number alterations generated from the WG Tumour/Normal pairs.
        Initial calls are generated with varscan, then processed with Sequenza
        """
        self.name = "sequenza"
        self.tables = [
            SequenzaTable(),
        ]

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

class RawSeqDataSection(Section):
    def __init__(self):
        self.title = "Raw Sequence Data"
        self.blurb = """
        Samples were sequenced on one or more sequencing runs.
        """
        self.name = "raw_seq_data"
        self.tables = [
            WGLaneLevelTable(
                "Whole Genome, Normal Sample",
                "WG",
                "Normal",
            ),
            WGLaneLevelTable(
                "Whole Genome, Tumour Sample",
                "WG",
                "Tumour",
            ),
            WTLaneLevelTable(),
        ]
