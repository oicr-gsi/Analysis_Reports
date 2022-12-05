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
    BAMQC4TableColumns,
    BAMQC4MergedTableColumns,
    RNASeqQCMergedTableColumns,
)
from plot import Plot

class Table:
    base_db_path = "/scratch2/groups/gsi/staging/qcetl_v1/"
    title = ""
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
            select_block = select_block + f"\"{value}\", "
        return select_block[:-2], indices

    def add_plot_data(self, plot, val, id):
        if plot in self.plots.keys():
            self.plots[plot].add_data(val, id)
    
    def generate_plots(self):
        for plot in self.plots.values():
            plot.generate_plot(self.process)
    


    
    def get_data(self):
        con = sqlite3.connect(self.base_db_path + self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            wfr = list(Table.data[id]["analysis"][self.pipeline_step][self.process].keys())[0]
            entry = {}
            row = cur.execute( # forloop automatically ignores empty SELECTs
                f"""
                select {select_block} from {self.source_table} where "Workflow Run SWID" like '%{wfr}%';
                """
            ).fetchall()[0]
            entry["sample_id"] = id
            for column in self.columns.keys():
                entry[column] = (
                    row[indices[column]]
                    if isinstance(row[indices[column]], str)
                    else round(row[indices[column]], 2)
                )
                self.add_plot_data(column, entry[column], id)
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
        }
        return context

class BAMQC4Table(Table):
    def __init__(self):
        self.title = "BAMQC4"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            BAMQC4TableColumns.SampleID: "Sample ID",
            BAMQC4TableColumns.CoverageDedup: "Coverage Depth",
            BAMQC4TableColumns.InsertSizeAvg: "Insert Size",
            BAMQC4TableColumns.MarkDupPctDup: "Duplication",
            BAMQC4TableColumns.TotalClusters: "Reads Pairs",
            BAMQC4TableColumns.MappedReads: "Mapped Reads",
        }
        self.columns = {
            BAMQC4TableColumns.CoverageDedup: "coverage deduplicated",
            BAMQC4TableColumns.InsertSizeAvg: "insert size average",
            BAMQC4TableColumns.MarkDupPctDup: "mark duplicates_PERCENT_DUPLICATION",
            BAMQC4TableColumns.TotalClusters: "total clusters",
        }
        self.pipeline_step = "alignments_WG.lanelevel"
        self.source_table = ["dnaseqqc_dnaseqqc_5", "bamqc4_bamqc4_5"]
        self.source_db = ["dnaseqqc", "bamqc4"]
        self.process = "bwaMem"
        self.plots = {}

    def get_data(self):
        dnaseqqc_con = sqlite3.connect(self.base_db_path + self.source_db[0] + "/latest")
        dnaseqqc_cur = dnaseqqc_con.cursor()
        bamqc4_con = sqlite3.connect(self.base_db_path+ self.source_db[1] + "/latest")
        bamqc4_cur = bamqc4_con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            lims_keys = [
                Table.data[id]["analysis"][self.pipeline_step][self.process][key]["limkeys"]
                for key in Table.data[id]["analysis"][self.pipeline_step][self.process].keys()
            ]

            for lims in lims_keys:
                target_cur = dnaseqqc_cur
                src_table_index = 0
                entry = {}
                rows = target_cur.execute(
                    f"""
                    select {select_block} from {self.source_table[0]} where "Pinery Lims ID" = '{lims}';
                    """).fetchall()
                if 0 == len(rows):
                    target_cur = bamqc4_cur
                    src_table_index = 1
                    rows = target_cur.execute(
                        f"""
                        select {select_block} from {self.source_table[1]} where "Pinery Lims ID" = '{lims}';
                        """).fetchall()
                for row in rows:

                    entry["sample_id"] = id
                    for column in self.columns.keys():
                        entry[column] = (
                            row[indices[column]]
                            if isinstance(row[indices[column]], str)
                            else round(row[indices[column]], 2)
                        )
                        self.add_plot_data(column, entry[column], id)                   
                for row in target_cur.execute(
                f"""
                    select (CAST("mapped reads" as FLOAT) /CAST("total reads" as FLOAT)) from {self.source_table[src_table_index]} where "Pinery Lims ID" = '{lims}';
                    """
                ):
                    entry[BAMQC4TableColumns.MappedReads] = row[0]
                    self.add_plot_data(BAMQC4TableColumns.MappedReads, row[0], id)
                data.append(entry)

        dnaseqqc_cur.close()
        dnaseqqc_con.close()
        bamqc4_cur.close()
        bamqc4_con.close()

        return data

class BAMQC4MergedTable(Table):
    def __init__(self):
        self.title = "BAMQC4Merged"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            BAMQC4MergedTableColumns.SampleID: "Sample ID",
            BAMQC4MergedTableColumns.CoverageDedup: "Coverage Depth",
            BAMQC4MergedTableColumns.MarkDupPctDup: "Duplication",
            BAMQC4MergedTableColumns.TotalClusters: "Reads Pairs",
            # BAMQC4MergedTableColumns.MappedReads: "Mapped Reads",
        }
        self.columns = {
            BAMQC4MergedTableColumns.CoverageDedup: "coverage deduplicated",
            BAMQC4MergedTableColumns.MarkDupPctDup: "mark duplicates_PERCENT_DUPLICATION",
            BAMQC4MergedTableColumns.TotalClusters: "total clusters",
        }
        self.source_table = "bamqc4merged_bamqc4merged_5"
        self.source_db = "bamqc4merged"
        self.process = "bamMergePreprocessing_by_tumor_group"

    def get_data(self):
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            lims = "[\"" + '\", \"'.join(Table.data[id][self.process]["limskeys"]) + "\"]"

            for row in cur.execute( # forloop automatically ignores empty SELECTs
                f"""
                select {select_block} from {self.source_table} where "Merged Pinery Lims ID"  = '{lims}';
                """):
                entry = {}
                entry["sample_id"] = id
                for column in self.columns.keys():
                    entry[column] = (
                        row[indices[column]]
                        if isinstance(row[indices[column]], str)
                        else round(row[indices[column]], 2)
                    )
                    self.add_plot_data(column, entry[column], id)  
                data.append(entry)
        
        cur.close()
        con.close()
        return data

class DellyTable(Table):
    def __init__(self):
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
            DellyTableColumns.NumCalls: "num_calls",
            DellyTableColumns.NumPASS: "num_PASS",
            DellyTableColumns.NumBND: "num_BND",
            DellyTableColumns.NumDEL: "num_DEL",
            DellyTableColumns.NumDUP: "num_DUP",
            DellyTableColumns.NumINS: "num_INS",
            DellyTableColumns.NumINV: "num_INV",
        }
        self.pipeline_step = "calls.structuralvariants"
        self.source_table = "analysis_delly_analysis_delly_1"
        self.source_db = "analysis_delly"
        self.process = "delly_matched_by_tumor_group"
        self.plots = {
            DellyTableColumns.NumPASS: Plot(
                title="SV Calls",
                x_axis="Sample IDs",
                y_axis="SV Calls"
            ),
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
        self.columns = {
            MetadataTableColumns.TissueOrigin: "Tissue Origin",
            MetadataTableColumns.TissueType: "Tissue Type",
            MetadataTableColumns.LibraryDesign: "Library Design",
        }
        self.pipeline_step = "calls.fusions"
        self.source_table = "analysis_starfusion_analysis_starfusion_1"
        self.source_db = "analysis_starfusion"
        self.process = "starfusion"

    def get_data(self):
        return super().get_data()
        

class Mutect2Table(Table):
    def __init__(self):
        self.headings = {
            Mutect2TableColumns.SampleID: "Sample ID",
            Mutect2TableColumns.NumCalls: "Mutation Calls",
            Mutect2TableColumns.NumPASS: "Mutation PASS Calls",
            Mutect2TableColumns.NumSNPs: "snvs",
            Mutect2TableColumns.NumIndels: "indels",
            Mutect2TableColumns.TITVRatio: "Ti/Tv",
        }
        self.columns = {
            Mutect2TableColumns.NumCalls: "num_calls",
            Mutect2TableColumns.NumPASS: "num_PASS",
            Mutect2TableColumns.NumSNPs: "num_SNPs",
            Mutect2TableColumns.NumIndels: "num_indels",
            Mutect2TableColumns.TITVRatio: "titv_ratio",
        }
        self.pipeline_step = "calls.mutations"
        self.source_table = "analysis_mutect2_analysis_mutect2_1"
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
        self.title = ""
        self.blurb = ""
        self.headings = {
            RNASeqQCTableColumns.SampleID: "Sample ID",
            RNASeqQCTableColumns.PctCodingBases: "Percent Coding",
            RNASeqQCTableColumns.TotalClusters: "Reads Pairs",
            RNASeqQCTableColumns.MappedReads: "Mapped Reads",
            RNASeqQCTableColumns.RRNAContamination: "rRNA Contamination",
        }
        self.columns = {
            RNASeqQCTableColumns.PctCodingBases: "PCT_CODING_BASES",
            RNASeqQCTableColumns.TotalClusters: "total clusters",
        }
        self.source_table = "rnaseqqc2_rnaseqqc2_2"
        self.source_db = "rnaseqqc2"
        self.process = "star_lane_level"
        self.pipeline_step = "alignments_WT.lanelevel"
        self.plots = {
            RNASeqQCTableColumns.PctCodingBases: Plot(
                title="Percent Coding",
                x_axis="Sample ID",
                y_axis="Percent Coding",
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
            ),
            RNASeqQCTableColumns.RRNAContamination: Plot(
                title="rRna Contamination",
                x_axis="Sample ID",
                y_axis="rRna Contamination",
            ),

        }
    
    def get_data(self):
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            lims_keys = [
                Table.data[id]["analysis"][self.pipeline_step][self.process][key]["limkeys"]
                for key in Table.data[id]["analysis"][self.pipeline_step][self.process].keys()
            ]

            for lims in lims_keys:
                entry = {}
                rows = cur.execute(
                    f"""
                    select {select_block} from {self.source_table} where "Pinery Lims ID" = '{lims}';
                    """).fetchall()
                for row in rows:
                    entry["sample_id"] = id
                    for column in self.columns.keys():
                        entry[column] = (
                            row[indices[column]]
                            if isinstance(row[indices[column]], str)
                            else round(row[indices[column]], 2)
                        )
                        self.add_plot_data(column, entry[column], id)                     


                for row in cur.execute(
                    f"""
                        select (1 - CAST("unmapped reads" as FLOAT) /CAST("total reads" as FLOAT)) * 100 from {self.source_table} where "Pinery Lims ID" like '%{lims}%';
                        """
                    ):
                    entry[RNASeqQCTableColumns.MappedReads] = round(row[0], 2)
                    self.add_plot_data(RNASeqQCTableColumns.MappedReads, row[0], id)  

                for row in cur.execute(
                f"""
                    select
                    (CAST("rrna contamination properly paired" as FLOAT)
                    /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT))
                    from {self.source_table}
                    where "Pinery Lims ID" like '%{lims}%';
                    """
                ):
                    entry[RNASeqQCTableColumns.RRNAContamination] = round(row[0], 2)
                    self.add_plot_data(RNASeqQCTableColumns.RRNAContamination, row[0], id)  

                data.append(entry)

        cur.close()
        con.close()
        return data

class RNASeqQCMergedTable(Table):
    def __init__(self):
        self.title = ""
        self.blurb = ""
        self.headings = {
            RNASeqQCMergedTableColumns.SampleID: "Sample ID",
            RNASeqQCMergedTableColumns.PctCodingBases: "Percent Coding",
            RNASeqQCMergedTableColumns.TotalClusters: "Total Clusters",
            RNASeqQCMergedTableColumns.MappedReads: "Mapped Reads",
            RNASeqQCMergedTableColumns.RRNAContamination: "rRNA Contamination",
            RNASeqQCMergedTableColumns.NumLimsKeys: "Lanes Sequenced",
        }
        self.columns = {
            RNASeqQCMergedTableColumns.PctCodingBases: "PCT_CODING_BASES",
            RNASeqQCMergedTableColumns.TotalClusters: "total clusters",

        }
        self.pipeline_step = "alignments_WT.callready"
        self.source_table = "rnaseqqc2merged_rnaseqqc2merged_2"
        self.source_db = "rnaseqqc2merged"
        self.process = "star_call_ready"

    def get_data(self):
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            limskeys = list(Table.data[id]["analysis"][self.pipeline_step][self.process].values())[0]["limkeys"].split(':')
            lims = "[\"" + "\", \"".join(limskeys) + "\"]"

            print(lims)
            entry = {}
            for row in cur.execute( # forloop automatically ignores empty SELECTs
                f"""
                select {select_block} from {self.source_table} where "Merged Pinery Lims ID"  = '{lims}';
                """):
                entry["sample_id"] = id
                for column in self.columns.keys():
                    entry[column] = (
                        row[indices[column]]
                        if isinstance(row[indices[column]], str)
                        else round(row[indices[column]], 2)
                    )

            for row in cur.execute(
                f"""
                    select (CAST("mapped reads" as FLOAT) /CAST("total reads" as FLOAT)) from {self.source_table} where "Merged Pinery Lims ID" = '{lims}';
                    """
                ):
                entry[RNASeqQCMergedTableColumns.MappedReads] = round(row[0], 2)

            for row in cur.execute(
               f"""
                select
                (CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT))
                from {self.source_table}
                where "Merged Pinery Lims ID" = '{lims}';
                """
            ):
                entry[RNASeqQCMergedTableColumns.RRNAContamination] = round(row[0], 2)
            
            entry[RNASeqQCMergedTableColumns.NumLimsKeys] = len(limskeys)
            data.append(entry)
        
        cur.close()
        con.close()
        return data

class RSEMTable(Table):
    def __init__(self):
        self.title = "Call Summary"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            RSEMTableColumns.SampleID: "Sample ID",
            RSEMTableColumns.Total: "Total Genes",
            RSEMTableColumns.PctNonZero: "Percent Expressed",
            RSEMTableColumns.Q0_05: "TPM, 5th percentile",
            RSEMTableColumns.Q0_5: "Median TPM",
            RSEMTableColumns.Q0_95: "TPM, 95th percentile",
        }
        self.columns = {
            RSEMTableColumns.Total: "total",
            RSEMTableColumns.PctNonZero: "pct_non_zero",
            RSEMTableColumns.Q0_05: "Q0.05",
            RSEMTableColumns.Q0_5: "Q0.5",
            RSEMTableColumns.Q0_95: "Q0.95",
        }
        self.pipeline_step = "calls.expression"
        self.source_table = "analysis_rsem_analysis_rsem_1"
        self.source_db = "analysis_rsem"
        self.process = "rsem"

class SequenzaTable(Table):
    def __init__(self):
        self.title = "Alternative Solutions"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            SequenzaTableColumns.SampleID: "Sample ID",
            SequenzaTableColumns.Cellularity: "Cellularity",
            SequenzaTableColumns.Ploidy: "Ploidy",
            SequenzaTableColumns.SLPP: "SLPP",
            SequenzaTableColumns.FGA: "FGA"
        }
        self.columns = {
            SequenzaTableColumns.Cellularity: "cellularity",
            SequenzaTableColumns.Ploidy: "ploidy",
            SequenzaTableColumns.SLPP: "SLPP",
            # SequenzaTableColumns.FGA: "fga" this column is added manually
        }
        self.pipeline_step = "calls.copynumber"
        self.source_table = ["analysis_sequenza_analysis_sequenza_alternative_solutions_1", "analysis_sequenza_analysis_sequenza_gamma_500_fga_1"]
        self.source_db = "analysis_sequenza"
        self.process = "sequenza_by_tumor_group"
    
    def get_data(self):
        con = sqlite3.connect(self.base_db_path + self.source_db + "/latest")
        cur = con.cursor()
        select_block, indices = self.get_select()
        data = []
        for id in Table.sample_ids:
            wfr = list(Table.data[id]["analysis"][self.pipeline_step][self.process].keys())[0]

            entry = {}
            row = cur.execute( # forloop automatically ignores empty SELECTs
                f"""
                select {select_block} from {self.source_table[0]} where "Workflow Run SWID" like '%{wfr}%' and gamma = 500;
                """
            ).fetchall()[0]

            entry["sample_id"] = id
            for column in self.columns.keys():
                entry[column] = (
                    row[indices[column]]
                    if isinstance(row[indices[column]], str)
                    else round(row[indices[column]], 2)
                )
                self.add_plot_data(column, entry[column], id)
            
            row = cur.execute(
                f"""
                select fga from {self.source_table[1]} where "Workflow Run SWID" like '%{wfr}%';
                """
            ).fetchall()[0]
            entry[SequenzaTableColumns.FGA] = row[0]
            self.add_plot_data(SequenzaTableColumns.FGA, row[0], id)
            data.append(entry)
        
        cur.close()
        con.close()
        return data

class StarFusionTable(Table):
    def __init__(self):
        self.title = ""
        self.blurb = ""
        self.headings = {
            StarFusionTableColumns.SampleID: "Sample ID",
            StarFusionTableColumns.NumRecords: "Fusion Calls",
        }
        self.columns = {
            StarFusionTableColumns.NumRecords: "num_records",
        }
        self.pipeline_step = "calls.fusions"
        self.source_table = "analysis_starfusion_analysis_starfusion_1"
        self.source_db = "analysis_starfusion"
        self.process = "starfusion"
        self.plots = {
            "num_records": Plot(
                title="Fusion Calls",
                x_axis="Sample ID",
                y_axis="Fusion Calls"
            )
        }
