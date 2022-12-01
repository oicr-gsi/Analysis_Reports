from tables import (
    RSEMTable, 
    SequenzaAltSolnTable,
    SequenzaFGATable,
    DellyTable,
    StarFusionTable,
    Mutect2Table,
    Mutect2TITVStatsTable,
    BAMQC4Table,
    MetadataTable,
    RNASeqQCTable,
    BAMQC4MergedTable,
    RNASeqQCMergedTable,
)
from typing import Dict, Type, List, Callable, Union, Tuple, Set, Any


class Section:
    title: str
    blurb: str
    tables: List[Any]
    name: str

    def __init__(self):
        pass

    def load_context(self):
        context = {}
        context["title"] = self.title
        context["blurb"] = self.blurb
        context["tables"] = {}
        for count, table in enumerate(self.tables):
            context["tables"][count] = table.load_context()
        return context

class RSEMSection(Section):
    def __init__(self):
        self.title = "RSEM"
        self.blurb = """
        It was the best of times, it was the worst of times, it was the age of
        wisdom, it was the age of foolishness, it was the epoch of belief, it was
        the epoch of incredulity...
        """
        self.name = "rsem"
        self.tables = [
            RSEMTable(),
        ]
    
class SequenzaSection(Section):
    def __init__(self):
        self.title = "Sequenza"
        self.blurb = """
        It was the best of times, it was the worst of times, it was the age of
        wisdom, it was the age of foolishness, it was the epoch of belief, it was
        the epoch of incredulity...
        """
        self.name = "sequenza"
        self.tables = [
            SequenzaAltSolnTable(),
            SequenzaFGATable(),
        ]

class DellySection(Section):
    def __init__(self):
        self.title = "Delly"
        self.blurb = """
        It was the best of times, it was the worst of times, it was the age of
        wisdom, it was the age of foolishness, it was the epoch of belief, it was
        the epoch of incredulity...
        """
        self.name = "delly"
        self.tables = [
            DellyTable(),
        ]

class StarFusionSection(Section):
    def __init__(self):
        self.title = "StarFusion"
        self.blurb = """
        It was the best of times, it was the worst of times, it was the age of
        wisdom, it was the age of foolishness, it was the epoch of belief, it was
        the epoch of incredulity...
        """
        self.name = "starfusion"
        self.tables = [
            StarFusionTable(),
        ]

class Mutect2Section(Section):
    def __init__(self):
        self.title = "Mutect2"
        self.blurb = """
        It was the best of times, it was the worst of times, it was the age of
        wisdom, it was the age of foolishness, it was the epoch of belief, it was
        the epoch of incredulity...
        """
        self.name = "mutect2"
        self.tables = [
            Mutect2Table(),
            Mutect2TITVStatsTable(),
        ]

class MetadataSection(Section):
    def __init__(self):
        self.title = "Metadata"
        self.blurb = """
        It was the best of times, it was the worst of times, it was the age of
        wisdom, it was the age of foolishness, it was the epoch of belief, it was
        the epoch of incredulity...
        """
        self.name = "metadata"
        self.tables = [
            MetadataTable(),
        ]

class RawSeqDataSection(Section):
    def __init__(self):
        self.title = "Raw Sequence Data"
        self.blurb = """
        It was the best of times, it was the worst of times, it was the age of
        wisdom, it was the age of foolishness, it was the epoch of belief, it was
        the epoch of incredulity...
        """
        self.name = "raw_data"
        self.tables = [
            # BAMQC4Table(),
            # RNASeqQCTable(),
            # BAMQC4MergedTable(),
            RNASeqQCMergedTable(),
        ]
