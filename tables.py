from typing import Dict, Type, List, Callable, Union, Tuple, Set, Any
import sqlite3
import json
from table_columns import (
    RSEMTableColumns,
    SequenzaTableColumns,
    DellyTableColumns,
    Mutect2TableColumns,
    StarFusionTableColumns,
    MetadataTableColumns,
    RNASeqQCTableColumns,
    WGLaneLevelDataColumns,
    WGCallReadyDataColumns,
    RNASeqQCMergedTableColumns,
)
from plot import Plot

class Table:
    base_db_path = "/scratch2/groups/gsi/staging/qcetl_v1/"
    title: str
    blurb = ""
    headings: Dict[str,str]
    columns: Dict[str, int]
    data: Dict
    source_table: str
    source_db: str
    process: str
    plots = {}
    data = {}
    pipeline_step: str
    glossary: Dict[str,str]
    pct_stats = set()

    def __init__(self, input_file):
        with open(input_file) as f:
            Table.data = json.load(f)
        Table.sample_ids = list(Table.data.keys())
        Table.sample_ids.sort()

    def get_select(self):
        indices = {}
        select_block = ""
        for index, (key, value) in enumerate(self.columns.items()):
            indices[key] = index
            select_block = select_block + f"{value}, "
        return select_block[:-2], indices

    def add_plot_data(self, plot, val, id):
        if plot in self.plots.keys():
            self.plots[plot].add_data(val, id)
    
    def generate_plots(self):
        for plot in self.plots.values():
            plot.generate_plot(self.process)
    
    def get_sample_id(self, cur, case, swid, pk, table_index=0):
        row = cur.execute(
            f"""
            select "Tissue Type", "Tissue Origin", "Library Design", "Group ID"
            from {self.source_table[table_index]}
            where "{pk}" like '%{swid}%';
            """
        ).fetchall()[0]
        return f"{case}_{row[0]}_{row[1]}_{row[2]}_{row[3]}"

    
    def get_data(self):
        con = sqlite3.connect(self.base_db_path + self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for case in Table.sample_ids:
            wfr = list(Table.data[case]["analysis"][self.pipeline_step][self.process].keys())[0]
            entry = {}
            row = cur.execute(
                f"""
                select {select_block} from {self.source_table[0]} where "Workflow Run SWID" like '%{wfr}%';
                """
            ).fetchall()[0]
            entry["sample_id"] = self.get_sample_id(cur, case, wfr, "Workflow Run SWID")
            for column in self.columns.keys():
                entry[column] = (
                    row[indices[column]] * 100
                    if column in self.pct_stats
                    else row[indices[column]]
                )
                entry[column] = (
                    entry[column]
                    if isinstance(entry[column], str)
                    else round(entry[column], 2)
                )

                self.add_plot_data(column, entry[column], entry["sample_id"])
            
            # for key, value in entry.items():
            #     if not isinstance(value, str):
            #         entry[key] = f"{value:,}"
            data.append(entry)
        
        cur.close()
        con.close()
        return data

    def load_context(self):
        context = {
            "title": self.title,
            "columns": self.columns,
            "headings": self.headings,
            "data": self.get_data(),
            "blurb": self.blurb,
            "glossary": self.glossary,
        }
        return context

class WGLaneLevelData(Table):
    def __init__(self, title, lib_type, tiss_origin, blurb=""):
        self.title = title
        self.blurb = blurb
        self.lib_type = lib_type
        self.tiss_origin = tiss_origin
        self.headings = {
            WGLaneLevelDataColumns.SampleID: "Sample ID",
            WGLaneLevelDataColumns.Lane: "Lane",
            WGLaneLevelDataColumns.CoverageDedup: "Coverage Depth",
            WGLaneLevelDataColumns.InsertSizeAvg: "Insert Size",
            WGLaneLevelDataColumns.MarkDupPctDup: "Duplication",
            WGLaneLevelDataColumns.TotalClusters: "Read Pairs",
            WGLaneLevelDataColumns.MappedReads: "Mapped Reads",
        }
        self.glossary = {
            WGLaneLevelDataColumns.CoverageDedup: "Mean depth of coverage corrected for duplication",
            WGLaneLevelDataColumns.InsertSizeAvg: "Mean size of the sequenced insert",
            WGLaneLevelDataColumns.MarkDupPctDup: "percent of reads marked as duplicates",
            WGLaneLevelDataColumns.TotalClusters: "Number of read pairs generated",
            WGLaneLevelDataColumns.MappedReads: "Percent of reads mapping to the genomic reference",

        }
        self.columns = {
            WGLaneLevelDataColumns.CoverageDedup: "\"coverage deduplicated\"",
            WGLaneLevelDataColumns.InsertSizeAvg: "\"insert size average\"",
            WGLaneLevelDataColumns.MarkDupPctDup: "\"mark duplicates_PERCENT_DUPLICATION\"",
            WGLaneLevelDataColumns.TotalClusters: "\"total clusters\"",
            WGLaneLevelDataColumns.MappedReads: "(1 - (CAST(\"unmapped reads meta\" as FLOAT) / CAST(\"total input reads meta\" as FLOAT)))",
        }
        self.pipeline_step = "alignments_WG.lanelevel"
        self.source_table = ["dnaseqqc_dnaseqqc_5", "bamqc4_bamqc4_5"]
        self.source_db = ["dnaseqqc", "bamqc4"]
        self.process = "bwaMem"
        self.pct_stats = set(
            [
                WGLaneLevelDataColumns.MarkDupPctDup,
                WGLaneLevelDataColumns.MappedReads,
            ]
        )
        self.plots = {
            WGLaneLevelDataColumns.CoverageDedup: Plot(
                title="Coverage Depth",
                x_axis="Sample ID",
                y_axis="Coverage Depth",
            ),
            WGLaneLevelDataColumns.InsertSizeAvg: Plot(
                title="Insert Size",
                x_axis="Sample ID",
                y_axis="Insert Size",
            ),
            WGLaneLevelDataColumns.MarkDupPctDup: Plot(
                title="Duplication",
                x_axis="Sample ID",
                y_axis="Duplication",
                is_pct=True,
            ),
            WGLaneLevelDataColumns.TotalClusters: Plot(
                title="Read Pairs",
                x_axis="Sample ID",
                y_axis="Read Pairs",
            ),
            WGLaneLevelDataColumns.MappedReads: Plot(
                title="Mapped Reads",
                x_axis="Sample ID",
                y_axis="Mapped Reads",
                is_pct=True,
            ),

        }


    def get_sample_id(self, case, swid):
        lims_lane_mapping = Table.data[case][self.lib_type][self.tiss_origin]
        for id, value in lims_lane_mapping.items():
            if swid in value.keys():
                return id, value[swid]["run"]
        raise Exception("There is no Sample ID associated with the limkey")

    def get_data(self):
        dnaseqqc_con = sqlite3.connect(self.base_db_path + self.source_db[0] + "/latest")
        dnaseqqc_cur = dnaseqqc_con.cursor()
        bamqc4_con = sqlite3.connect(self.base_db_path+ self.source_db[1] + "/latest")
        bamqc4_cur = bamqc4_con.cursor()

        select_block, indices = self.get_select()
        data = []

        print(select_block)

        for id in Table.sample_ids:
            lims_keys = []
            for key in Table.data[id][self.lib_type][self.tiss_origin].keys():
                lims_keys = lims_keys + list(Table.data[id][self.lib_type][self.tiss_origin][key].keys())

            for lims in lims_keys:
                target_cur = dnaseqqc_cur
                src_table_index = 0
                entry = {}
                rows = target_cur.execute(
                    f"""
                    select {select_block} from {self.source_table[0]} where "Pinery Lims ID" like '%{lims}%';
                    """).fetchall()
                if 0 == len(rows):
                    target_cur = bamqc4_cur
                    src_table_index = 1
                    rows = target_cur.execute(
                        f"""
                        select {select_block} from {self.source_table[1]} where "Pinery Lims ID" like '%{lims}%';
                        """).fetchall()
                for row in rows:
                    entry[WGLaneLevelDataColumns.SampleID], entry[WGLaneLevelDataColumns.Lane] = self.get_sample_id(id, lims)
                    for column in self.columns.keys():
                        entry[column] = (
                            row[indices[column]] * 100
                            if column in self.pct_stats
                            else row[indices[column]]
                        )
                        entry[column] = (
                            entry[column]
                            if isinstance(entry[column], str)
                            else round(entry[column], 2)
                        )
                        self.add_plot_data(column, entry[column], entry[WGLaneLevelDataColumns.SampleID])                
                
                # for key, value in entry.items():
                #     if not isinstance(value, str):
                #         entry[key] = f"{value:,}"
                data.append(entry)

        dnaseqqc_cur.close()
        dnaseqqc_con.close()
        bamqc4_cur.close()
        bamqc4_con.close()

        return data

class WGCallReadyData(Table):
    def __init__(self, title, lib_type, tiss_origin, blurb=""):
        self.title = title
        self.blurb = blurb
        self.lib_type = lib_type
        self.tiss_origin = tiss_origin
        self.headings = {
            WGCallReadyDataColumns.SampleID: "Sample ID",
            WGCallReadyDataColumns.CoverageDedup: "Coverage Depth",
            WGCallReadyDataColumns.MarkDupPctDup: "Duplication",
            WGCallReadyDataColumns.TotalClusters: "Read Pairs",
            WGCallReadyDataColumns.MappedReads: "Mapped Reads",
            WGCallReadyDataColumns.NumLimsKeys: "Lanes Sequenced",
        }
        self.glossary = {
            WGCallReadyDataColumns.CoverageDedup: "Mean depth of coverage corrected for duplication",
            WGCallReadyDataColumns.MarkDupPctDup: "Percent of reads marked as duplicates",
            WGCallReadyDataColumns.TotalClusters: "Number of read pairs generated",
            WGCallReadyDataColumns.MappedReads: "Percent of reads mapping to the genomic reference",
            WGCallReadyDataColumns.NumLimsKeys: "Number of lanes of sequencing merged to call ready",
        }
        self.columns = {
            WGCallReadyDataColumns.CoverageDedup: "\"coverage deduplicated\"",
            WGCallReadyDataColumns.MarkDupPctDup: "\"mark duplicates_PERCENT_DUPLICATION\"",
            WGCallReadyDataColumns.TotalClusters: "\"total clusters\"",
            WGCallReadyDataColumns.MappedReads: """
            (1 - CAST("unmapped reads meta" as FLOAT)
            /CAST("total input reads meta" as FLOAT))
            """,
        }
        self.source_table = ["bamqc4merged_bamqc4merged_5"]
        self.source_db = "bamqc4merged"
        self.process = "bamMergePreprocessing_by_tumor_group"
        self.pct_stats = set(
            [
                WGCallReadyDataColumns.MarkDupPctDup,
                WGCallReadyDataColumns.MappedReads,
            ]
        )
        self.plots = {
            WGCallReadyDataColumns.CoverageDedup: Plot(
                "Coverage Depth",
                "Sample ID",
                "Coverage Depth",
            ),
            WGCallReadyDataColumns.MarkDupPctDup: Plot(
                "Duplication",
                "Sample ID",
                "Duplication",
                is_pct=True,
            ),
            WGCallReadyDataColumns.TotalClusters: Plot(
                "Read Pairs",
                "Sample ID",
                "Read Pairs",
            ),
            WGCallReadyDataColumns.MappedReads: Plot(
                "Mapped Reads",
                "Sample ID",
                "Mapped Reads",
                is_pct=True,
            )
        }

    def get_data(self):
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            lims_keys = []
            for key in Table.data[id][self.lib_type][self.tiss_origin].keys():
                lims_keys = lims_keys + list(Table.data[id][self.lib_type][self.tiss_origin][key].keys())
            lims_keys.sort()
            lims = "[\"" + '\", \"'.join(lims_keys) + "\"]"
            entry = {}
            # print(id, lims)
            row = cur.execute(
                f"""
                select {select_block} from {self.source_table[0]} where "Merged Pinery Lims ID"  like '%{lims}%';
                """
            ).fetchall()[0]

            entry["sample_id"] = self.get_sample_id(cur, id, lims, "Merged Pinery Lims ID")
            
            for column in self.columns.keys():
                entry[column] = (
                    row[indices[column]] * 100
                    if column in self.pct_stats
                    else row[indices[column]]
                )
                entry[column] = (
                    entry[column]
                    if isinstance(entry[column], str)
                    else round(entry[column], 2)
                )
                self.add_plot_data(column, entry[column], entry["sample_id"])
                
            entry[WGCallReadyDataColumns.NumLimsKeys] = len(lims_keys)

            # for key, value in entry.items():
            #     if not isinstance(value, str):
            #         entry[key] = f"{value:,}"   
            data.append(entry)
        
        cur.close()
        con.close()
        return data

class DellyTable(Table):
    def __init__(self):
        self.title = "Genomic Structural Variants"
        self.headings = {
            DellyTableColumns.SampleID: "Sample ID",
            DellyTableColumns.NumCalls: "SV Calls",
            DellyTableColumns.NumPASS: "SV PASS Calls",
            DellyTableColumns.NumBND: "Translocations",
            DellyTableColumns.NumDEL: "Deletions",
            DellyTableColumns.NumDUP: "Duplications",
            DellyTableColumns.NumINS: "Insertions",
            DellyTableColumns.NumINV: "Inversions",
        }
        self.columns = {
            DellyTableColumns.NumCalls: "\"num_calls\"",
            DellyTableColumns.NumPASS: "\"num_PASS\"",
            DellyTableColumns.NumBND: "\"num_BND\"",
            DellyTableColumns.NumDEL: "\"num_DEL\"",
            DellyTableColumns.NumDUP: "\"num_DUP\"",
            DellyTableColumns.NumINS: "\"num_INS\"",
            DellyTableColumns.NumINV: "\"num_INV\"",
        }
        self.pipeline_step = "calls.structuralvariants"
        self.source_table = ["analysis_delly_analysis_delly_1"]
        self.source_db = "analysis_delly"
        self.process = "delly_matched_by_tumor_group"
        self.plots = {
            DellyTableColumns.NumPASS: Plot(
                title="SV Calls",
                x_axis="Sample IDs",
                y_axis="SV Calls"
            ),
        }
        self.glossary = {
            DellyTableColumns.NumCalls: "The number of somatic structural variant calls identified by delly",
            DellyTableColumns.NumPASS: "The number of structural variant calls marked as PASS",
            DellyTableColumns.NumBND: "The number of PASS translocation calls",
            DellyTableColumns.NumDEL: "The number of PASS deletion calls",
            DellyTableColumns.NumDUP: "The number of PASS duplication calls",
            DellyTableColumns.NumINS: "The number of PASS insertions calls",
            DellyTableColumns.NumINV: "The number of PASS insertions calls",
        }

class MetadataTable(Table):
    def __init__(self):
        self.title = "Cases"
        self.blurb = ""
        self.headings = {
            MetadataTableColumns.Case: "Case",
            MetadataTableColumns.LibraryDesign: "Library Type",
            MetadataTableColumns.TissueType : "Sample Type",
            MetadataTableColumns.TissueOrigin : "Tissue Origin",
            MetadataTableColumns.SampleID: "Sample ID",
            MetadataTableColumns.ExternalID: "External ID",
        }
        self.columns = {}
        self.glossary ={}
        self.pipeline_step = "calls.fusions"
        self.source_table = ["analysis_starfusion_analysis_starfusion_1"]
        self.source_db = "analysis_starfusion"
        self.process = "starfusion"

    def get_data(self):
        data = []
        
        for case in Table.data.values():
            ids = (
                list(case["WT"]["Tumour"].keys())
                + list(case["WG"]["Tumour"].keys())
                + list(case["WG"]["Normal"].keys())
            )

            for id in ids:
                metadata = id.split("_")
                data.append(
                    {
                        MetadataTableColumns.Case: metadata[0] + "_" + metadata[1],
                        MetadataTableColumns.TissueType: metadata[2],
                        MetadataTableColumns.TissueOrigin: metadata[3],
                        MetadataTableColumns.LibraryDesign: metadata[4],
                        MetadataTableColumns.ExternalID: case["external_id"],
                        MetadataTableColumns.SampleID: id,
                    }
                )
        return data
        
class Mutect2Table(Table):
    def __init__(self):
        self.title = "Mutations"
        self.blurb = ""
        self.headings = {
            Mutect2TableColumns.SampleID: "Sample ID",
            Mutect2TableColumns.NumCalls: "Mutation Calls",
            Mutect2TableColumns.NumPASS: "Mutation PASS Calls",
            Mutect2TableColumns.NumSNPs: "snvs",
            Mutect2TableColumns.NumIndels: "indels",
            Mutect2TableColumns.TITVRatio: "Ti/Tv",
        }
        self.glossary = {
            Mutect2TableColumns.NumCalls: "The number of somatic mutation calls (snvs + indels)  identified by mutect2",
            Mutect2TableColumns.NumPASS: "The number of somatic mutation calls marked as PASS",
            Mutect2TableColumns.NumSNPs: "The number of somatic mutation calls marked as PASS",
            Mutect2TableColumns.NumIndels: "The number of structural variant calls marked as PASS",
            Mutect2TableColumns.TITVRatio: "Transition to Transversion ratio",
        }
        self.columns = {
            Mutect2TableColumns.NumCalls: "num_calls",
            Mutect2TableColumns.NumPASS: "num_PASS",
            Mutect2TableColumns.NumSNPs: "num_SNPs",
            Mutect2TableColumns.NumIndels: "num_indels",
            Mutect2TableColumns.TITVRatio: "titv_ratio",
        }
        self.pipeline_step = "calls.mutations"
        self.source_table = ["analysis_mutect2_analysis_mutect2_1"]
        self.source_db = "analysis_mutect2"
        self.process = "mutect2_matched_by_tumor_group"
        self.plots = {
            Mutect2TableColumns.NumPASS: Plot(
                title="Mutation Calls",
                x_axis="Sample IDs",
                y_axis="Mutation Calls"
            ),
            Mutect2TableColumns.TITVRatio: Plot(
                title="Ti/Tv",
                x_axis="Sample IDs",
                y_axis="Ti/Tv"
            ),
        }

class RNASeqQCTable(Table):
    def __init__(self):
        self.title = "WT Raw Sequence Data"
        self.blurb = ""
        self.headings = {
            RNASeqQCTableColumns.SampleID: "Sample ID",
            RNASeqQCTableColumns.Lane: "Lane",
            RNASeqQCTableColumns.PctCodingBases: "Percent Coding",
            RNASeqQCTableColumns.TotalClusters: "Read Pairs",
            RNASeqQCTableColumns.MappedReads: "Mapped Reads",
            RNASeqQCTableColumns.RRNAContamination: "rRNA Contamination",

        }
        self.glossary = {
            RNASeqQCTableColumns.PctCodingBases: "Percentage of bases mapping to the coding regions of the genome",
            RNASeqQCTableColumns.TotalClusters: "Number of read pairs generated",
            RNASeqQCTableColumns.MappedReads: "Percentage of reads mapping to the genomic reference",
            RNASeqQCTableColumns.RRNAContamination: "Pecentage of reads mapping to ribosomal RNA",

        }
        self.columns = {
            RNASeqQCTableColumns.PctCodingBases: "\"PCT_CODING_BASES\"",
            RNASeqQCTableColumns.TotalClusters: "\"total clusters\"",
            RNASeqQCTableColumns.MappedReads: """
                (1 - CAST(\"unmapped reads\" as FLOAT)
                /CAST(\"total reads\" as FLOAT))
            """,
            RNASeqQCTableColumns.RRNAContamination: """
                (CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT))
            """,
        }
        self.source_table = ["rnaseqqc2_rnaseqqc2_2"]
        self.source_db = "rnaseqqc2"
        self.process = "star_lane_level"
        self.pipeline_step = "alignments_WT.lanelevel"
        self.pct_stats = set(
            [
                RNASeqQCTableColumns.MappedReads,
                RNASeqQCTableColumns.RRNAContamination,
            ]
        )
        self.plots = {
            RNASeqQCTableColumns.PctCodingBases: Plot(
                title="Percent Coding",
                x_axis="Sample ID",
                y_axis="Percent Coding",
                is_pct=True,
            ),
            RNASeqQCTableColumns.TotalClusters: Plot(
                title="Reads Pair",
                x_axis="Sample ID",
                y_axis="Reads Pair",
            ),
            RNASeqQCTableColumns.MappedReads: Plot(
                title="Mapped Reads",
                x_axis="Sample ID",
                y_axis="Mapped Reads",
                is_pct=True,
            ),
            RNASeqQCTableColumns.RRNAContamination: Plot(
                title="rRNA Contamination",
                x_axis="Sample ID",
                y_axis="rRNA Contamination",
                is_pct=True,
            ),

        }
    
    def get_sample_id(self, case, swid):
        lims_lane_mapping = Table.data[case]["WT"]["Tumour"]
        for id, value in lims_lane_mapping.items():
            if swid in value.keys():
                return id, value[swid]["run"]
        raise Exception("There is no Sample ID associated with the limkey")

    def get_data(self):
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        print(select_block)

        for id in Table.sample_ids:
            lims_keys = []
            for key in Table.data[id]["WT"]["Tumour"].keys():
                lims_keys = lims_keys + list(Table.data[id]["WT"]["Tumour"][key].keys())

            for lims in lims_keys:
                entry = {}
                rows = cur.execute(
                    f"""
                    select {select_block} from {self.source_table[0]} where "Pinery Lims ID" = '{lims}';
                    """).fetchall()
                for row in rows:
                    entry[RNASeqQCTableColumns.SampleID], entry[RNASeqQCTableColumns.Lane] = self.get_sample_id(id, lims)
                    for column in self.columns.keys():
                        entry[column] = (
                            row[indices[column]] * 100
                            if column in self.pct_stats
                            else row[indices[column]]
                        )
                        entry[column] = (
                            entry[column]
                            if isinstance(entry[column], str)
                            else round(entry[column], 2)
                        )
                        self.add_plot_data(column, entry[column], entry["sample_id"])

                # for key, value in entry.items():
                #     if not isinstance(value, str):
                #         entry[key] = f"{value:,}"
                data.append(entry)

        cur.close()
        con.close()
        return data

class RNASeqQCMergedTable(Table):
    def __init__(self):
        self.title = "WT Call Ready Alignments"
        self.blurb = ""
        self.headings = {
            RNASeqQCMergedTableColumns.SampleID: "Sample ID",
            RNASeqQCMergedTableColumns.PctCodingBases: "Percent Coding",
            RNASeqQCMergedTableColumns.TotalClusters: "Read Pairs",
            RNASeqQCMergedTableColumns.MappedReads: "Mapped Reads",
            RNASeqQCMergedTableColumns.RRNAContamination: "rRNA Contamination",
            RNASeqQCMergedTableColumns.NumLimsKeys: "Lanes Sequenced",
        }
        self.glossary = {
            RNASeqQCMergedTableColumns.PctCodingBases: "Percentage of bases mapping to the coding regions of the genome",
            RNASeqQCMergedTableColumns.TotalClusters: "	Number of read pairs generated",
            RNASeqQCMergedTableColumns.MappedReads: "Percentage of reads mapping to the genomic reference",
            RNASeqQCMergedTableColumns.RRNAContamination: "Pecentage of reads mapping to ribosomal RNA",
            RNASeqQCMergedTableColumns.NumLimsKeys: "Number of lanes of sequencing merged to call ready",
        }
        self.columns = {
            RNASeqQCMergedTableColumns.PctCodingBases: "\"PCT_CODING_BASES\"",
            RNASeqQCMergedTableColumns.TotalClusters: "\"total clusters\"",
            RNASeqQCMergedTableColumns.MappedReads: """
                (1 - CAST("unmapped reads" as FLOAT)
                /CAST("total reads" as FLOAT))
            """,
            RNASeqQCMergedTableColumns.RRNAContamination: """
                (CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT))
            """,
        }
        self.pct_stats = set(
            [
                RNASeqQCMergedTableColumns.MappedReads,
                RNASeqQCMergedTableColumns.RRNAContamination,
            ]
        )
        self.pipeline_step = "alignments_WT.callready"
        self.source_table = ["rnaseqqc2merged_rnaseqqc2merged_2"]
        self.source_db = "rnaseqqc2merged"
        self.process = "star_call_ready"
        self.plots = {
            RNASeqQCMergedTableColumns.PctCodingBases: Plot(
                title="Percent Coding",
                x_axis="Sample ID",
                y_axis="Percent Coding",
                is_pct=True,
            ),
            RNASeqQCMergedTableColumns.TotalClusters: Plot(
                title="Read Pairs",
                x_axis="Sample ID",
                y_axis="Read Pairs",
            ),
            RNASeqQCMergedTableColumns.MappedReads: Plot(
                title="Mapped Reads",
                x_axis="Sample ID",
                y_axis="Mapped Reads",
                is_pct=True,
            ),
            RNASeqQCMergedTableColumns.RRNAContamination: Plot(
                title="rRNA Contamination",
                x_axis="Sample ID",
                y_axis="rRNA Contamination",
                is_pct=True,
            ),
        }

    def get_data(self):
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            limskeys = list(Table.data[id]["analysis"][self.pipeline_step][self.process].values())[0]["limkeys"].split(':')
            limskeys.sort()
            lims = "[\"" + "\", \"".join(limskeys) + "\"]"

            entry = {}
            for row in cur.execute(
                f"""
                select {select_block} from {self.source_table[0]} where "Merged Pinery Lims ID"  = '{lims}';
                """):
                entry["sample_id"] = self.get_sample_id(cur, id, lims, "Merged Pinery Lims ID")
                for column in self.columns.keys():
                    entry[column] = (
                        row[indices[column]] * 100
                        if column in self.pct_stats
                        else row[indices[column]]
                    )
                    entry[column] = (
                        entry[column]
                        if isinstance(entry[column], str)
                        else round(entry[column], 2)
                    )
                    self.add_plot_data(column, entry[column], entry[RNASeqQCMergedTableColumns.SampleID]) 

            entry[RNASeqQCMergedTableColumns.NumLimsKeys] = len(limskeys)
            
            # for key, value in entry.items():
            #     if not isinstance(value, str):
            #         entry[key] = f"{value:,}"
            data.append(entry)
        
        cur.close()
        con.close()
        return data

class RSEMTable(Table):
    def __init__(self):
        self.title = "Gene Expression"
        self.blurb = ""
        self.headings = {
            RSEMTableColumns.SampleID: "Sample ID",
            RSEMTableColumns.Total: "Total Genes",
            RSEMTableColumns.PctNonZero: "Percent Expressed",
            RSEMTableColumns.Q0_05: "TPM, 5th percentile",
            RSEMTableColumns.Q0_5: "Median TPM",
            RSEMTableColumns.Q0_95: "TPM, 95th percentile",
        }
        self.glossary = {
            RSEMTableColumns.Total: "The total number of genes in the gene model",
            RSEMTableColumns.PctNonZero: "The percentage of genes expressed (non-zero)",
            RSEMTableColumns.Q0_05: "5th percentile of the Transcripts per Million scores",
            RSEMTableColumns.Q0_5: "Median Transcripts per Million score",
            RSEMTableColumns.Q0_95: "95th percentile of the Transcripts per Million score",
        }
        self.columns = {
            RSEMTableColumns.Total: "\"total\"",
            RSEMTableColumns.PctNonZero: "\"pct_non_zero\"",
            RSEMTableColumns.Q0_05: "\"Q0.05\"",
            RSEMTableColumns.Q0_5: "\"Q0.5\"",
            RSEMTableColumns.Q0_95: "\"Q0.95\"",
        }
        self.pct_stats = set(
            [
                RSEMTableColumns.PctNonZero,
            ]
        )
        self.pipeline_step = "calls.expression"
        self.source_table = ["analysis_rsem_analysis_rsem_1"]
        self.source_db = "analysis_rsem"
        self.process = "rsem"
        self.plots = {
            "pct_non_zero": Plot(
                title="Percent Expressed",
                x_axis="Sample ID",
                y_axis="Percent Expressed",
                is_pct=True,
            ),
            "Q0.5": Plot(
                title="Median TMP",
                x_axis="Sample ID",
                y_axis="Median TMP",
            ),
        }

class SequenzaTable(Table):
    def __init__(self):
        self.title = "Copy Number Alterations"
        self.blurb = ""
        self.headings = {
            SequenzaTableColumns.SampleID: "Sample ID",
            SequenzaTableColumns.Cellularity: "Cellularity",
            SequenzaTableColumns.Ploidy: "Ploidy",
            SequenzaTableColumns.FGA: "FGA"
        }
        self.glossary = {
            SequenzaTableColumns.Cellularity: "cellularity estimate (gamma=500)",
            SequenzaTableColumns.Ploidy: "ploidy estimate (gamma = 500)",
            SequenzaTableColumns.FGA: "Fraction of the genome altered (gamma = 500)"
        }
        self.columns = {
            SequenzaTableColumns.Cellularity: "\"cellularity\"",
            SequenzaTableColumns.Ploidy: "\"ploidy\"",
        }
        self.pipeline_step = "calls.copynumber"
        self.source_table = ["analysis_sequenza_analysis_sequenza_alternative_solutions_1", "analysis_sequenza_analysis_sequenza_gamma_500_fga_1"]
        self.source_db = "analysis_sequenza"
        self.process = "sequenza_by_tumor_group"
        self.pct_stats = set(
            [
               SequenzaTableColumns.FGA, 
            ]
        )
        self.plots = {
            SequenzaTableColumns.Cellularity: Plot(
                title="Cellularity",
                x_axis="Sample ID",
                y_axis="Cellularity",
            ),
            SequenzaTableColumns.Ploidy: Plot(
                title="Ploidy",
                x_axis="Sample ID",
                y_axis="Ploidy",
            ),
            SequenzaTableColumns.FGA: Plot(
                title="FGA",
                x_axis="Sample ID",
                y_axis="FGA",
            ),
        }
    
    def get_data(self):
        con = sqlite3.connect(self.base_db_path + self.source_db + "/latest")
        cur = con.cursor()
        select_block, indices = self.get_select()
        data = []
        for id in Table.sample_ids:
            wfr = list(Table.data[id]["analysis"][self.pipeline_step][self.process].keys())[0]

            entry = {}
            row = cur.execute(
                f"""
                select {select_block} from {self.source_table[0]} where "Workflow Run SWID" like '%{wfr}%' and gamma = 500;
                """
            ).fetchall()[0]

            entry["sample_id"] = self.get_sample_id(cur, id, wfr, "Workflow Run SWID")
            for column in self.columns.keys():
                entry[column] = (
                    row[indices[column]] * 100
                    if column in self.pct_stats
                    else row[indices[column]]
                )
                entry[column] = (
                    entry[column]
                    if isinstance(entry[column], str)
                    else round(entry[column], 2)
                )
                self.add_plot_data(column, entry[column], entry["sample_id"])

            #get FGA
            row = cur.execute(
                f"""
                select fga from {self.source_table[1]} where "Workflow Run SWID" like '%{wfr}%';
                """
            ).fetchall()[0]
            entry[SequenzaTableColumns.FGA] = round(row[0],2)
            self.add_plot_data(SequenzaTableColumns.FGA, row[0], entry["sample_id"])

            # for key, value in entry.items():
            #     if not isinstance(value, str):
            #         entry[key] = f"{value:,}"
            data.append(entry)
        
        cur.close()
        con.close()
        return data

class StarFusionTable(Table):
    def __init__(self):
        self.title = "Gene Fusions"
        self.blurb = ""
        self.headings = {
            StarFusionTableColumns.SampleID: "Sample ID",
            StarFusionTableColumns.NumRecords: "Fusion Calls",
        }
        self.glossary = {
            StarFusionTableColumns.NumRecords: "Number of gene fusions identified by StarFusion",
        }
        self.columns = {
            StarFusionTableColumns.NumRecords: "\"num_records\"",
        }
        self.pipeline_step = "calls.fusions"
        self.source_table = ["analysis_starfusion_analysis_starfusion_1"]
        self.source_db = "analysis_starfusion"
        self.process = "starfusion"
        self.plots = {
            StarFusionTableColumns.NumRecords: Plot(
                title="Fusion Calls",
                x_axis="Sample ID",
                y_axis="Fusion Calls"
            )
        }
