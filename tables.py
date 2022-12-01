from typing import Dict, Type, List, Callable, Union, Tuple, Set, Any
import sqlite3
import json
from table_columns import (
    SequenzaFGATableColumns,
    RSEMTableColumns,
    SequenzaAltSolnTableColumns,
    DellyTableColumns,
    Mutect2TableColumns,
    StarFusionTableColumns,
    MetadataTableColumns,
    Mutect2TITVStatsTableColumns,
    RNASeqQCTableColumns,
    BAMQC4TableColumns,
    BAMQC4MergedTableColumns,
    RNASeqQCMergedTableColumns,
)
from plot import Plot

class Table:
    base_db_path = "/.mounts/labs/gsiprojects/gsi/gsiusers/jqian/gsi-qc-etl/gsiqcetl/"
    title: str
    blurb: str
    headings: Dict[str,str]
    columns: Dict[str, int]
    data: Dict
    source_table: str
    source_db: str
    process: str
    plots: Dict[str, Plot]

    def __init__(self, input_file):
        self.data = {}
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

    def get_plot_data(self, plot, val, id):
        if plot in self.plots.keys():
            self.plots[plot].data["x"].append(id)
            self.plots[plot].data["y"].append(val)

    def get_data(self):
        con = sqlite3.connect(self.base_db_path + self.source_db + "/temp/"+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            wfr = Table.data[id][self.process]["wfrun"]

            for row in cur.execute( # forloop automatically ignores empty SELECTs
                f"""
                select {select_block} from {self.source_table} where "Workflow Run SWID" like '%{wfr}%';
                """):
                entry = {}
                entry["sample_id"] = id
                for column in self.columns.keys():
                    entry[column] = (
                        row[indices[column]]
                        if isinstance(row[indices[column]], str)
                        else round(row[indices[column]], 2)
                    )
                    self.get_plot_data(column, entry[column], id)
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
        self.source_table = ["dnaseqqc_dnaseqqc_5", "bamqc4_bamqc4_5"]
        self.source_db = ["dnaseqqc", "bamqc4"]
        self.process = "bwaMem"

    def get_data(self):
        dnaseqqc_con = sqlite3.connect("/scratch2/groups/gsi/production/qcetl_v1/"+ self.source_db[0] + "/latest")
        dnaseqqc_cur = dnaseqqc_con.cursor()
        bamqc4_con = sqlite3.connect("/scratch2/groups/gsi/production/qcetl_v1/"+ self.source_db[1] + "/latest")
        bamqc4_cur = bamqc4_con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            target_cur = dnaseqqc_cur
            src_table_index = 0

            lims = Table.data[id][self.process]["limskeys"][0]
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

                entry["sample_id"] = id
                for column in self.columns.keys():
                    entry[column] = (
                        row[indices[column]]
                        if isinstance(row[indices[column]], str)
                        else round(row[indices[column]], 2)
                    )
                
            for row in target_cur.execute(
               f"""
                select (CAST("mapped reads" as FLOAT) /CAST("total reads" as FLOAT)) from {self.source_table[src_table_index]} where "Pinery Lims ID" like '%{lims}%';
                """
            ):
                entry[BAMQC4TableColumns.MappedReads] = row[0]
            data.append(entry)
        dnaseqqc_cur.close()
        dnaseqqc_con.close()
        bamqc4_cur.close()
        bamqc4_con.close()

        return data

class BAMQC4MergedTable(Table):
    def __init__(self):
        self.base_db_path = "/scratch2/groups/gsi/production/qcetl_v1/"
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
                data.append(entry)
        
        cur.close()
        con.close()
        return data

class DellyTable(Table):
    def __init__(self):
        self.title = "Delly"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            DellyTableColumns.SampleID: "Sample ID",
            DellyTableColumns.NumCalls: "Number of Calls",
            DellyTableColumns.NumPASS: "Number of Pass Calls",
            DellyTableColumns.NumBND: "Number of BND",
            DellyTableColumns.NumDEL: "Number of DEL",
            DellyTableColumns.NumDUP: "Number of DUP",
            DellyTableColumns.NumINS: "Number of INS",
            DellyTableColumns.NumINV: "Number of INV",
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
        self.source_table = "analysis_delly_analysis_delly_1"
        self.source_db = "analysis_delly"
        self.process = "delly_matched_by_tumor_group"

class MetadataTable(Table):
    def __init__(self):
        self.title = "Metadata"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            MetadataTableColumns.SampleID: "Sample ID",
            MetadataTableColumns.GroupID: "Group ID",
            MetadataTableColumns.TissueType : "Tissue Type",
            MetadataTableColumns.TissueOrigin : "Tissue Origin",
            MetadataTableColumns.LibraryDesign: "Library Design",
        }
        self.columns = {
            MetadataTableColumns.GroupID: "Group ID",
            MetadataTableColumns.TissueOrigin: "Tissue Origin",
            MetadataTableColumns.TissueType: "Tissue Type",
            MetadataTableColumns.LibraryDesign: "Library Design",
        }
        self.source_table = "analysis_starfusion_analysis_starfusion_1"
        self.source_db = "analysis_starfusion"
        self.process = "starfusion"

class Mutect2Table(Table):
    def __init__(self):
        self.title = "Mutect2"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            Mutect2TableColumns.SampleID: "Sample ID",
            Mutect2TableColumns.NumCalls: "Number of Calls",
            Mutect2TableColumns.NumPASS: "Number of PASS calls",
            Mutect2TableColumns.NumSNPs: "Number of SNPs",
            Mutect2TableColumns.NumMultiSNPs: "Number of multi-SNPs",
            Mutect2TableColumns.NumIndels: "Number of indels",
            Mutect2TableColumns.NumMNPs: "Number of MNPs",
        }
        self.columns = {
            Mutect2TableColumns.NumCalls: "num_calls",
            Mutect2TableColumns.NumPASS: "num_PASS",
            Mutect2TableColumns.NumSNPs: "num_SNPs",
            Mutect2TableColumns.NumMultiSNPs: "num_multi_SNPs",
            Mutect2TableColumns.NumIndels: "num_indels",
            Mutect2TableColumns.NumMNPs: "num_MNPs",
        }
        self.data = {}
        self.source_table = "analysis_mutect2_analysis_mutect2_1"
        self.source_db = "analysis_mutect2"
        self.process = "mutect2_matched_by_tumor_group"

class Mutect2TITVStatsTable(Table):
    def __init__(self):
        self.title = "Mutect2 TI/TV Summary"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            Mutect2TITVStatsTableColumns.SampleID: "Sample ID",
            Mutect2TITVStatsTableColumns.TITVRatio: "TI/TV ratio",
            Mutect2TITVStatsTableColumns.PctTI: "Percentage of TI",
            Mutect2TITVStatsTableColumns.PctTV: "Percentage of TV",
        }
        self.columns = {
            Mutect2TITVStatsTableColumns.TITVRatio: "titv_ratio",
            Mutect2TITVStatsTableColumns.PctTI: "pct_ti",
            Mutect2TITVStatsTableColumns.PctTV: "pct_tv",
        }
        self.data = {}
        self.source_table = "analysis_mutect2_analysis_mutect2_1"
        self.source_db = "analysis_mutect2"
        self.process = "mutect2_matched_by_tumor_group"

class RNASeqQCTable(Table):
    def __init__(self):
        self.base_db_path = "/scratch2/groups/gsi/production/qcetl_v1/"
        self.title = "RNASeqQC"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
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
    
    def get_data(self):
        con = sqlite3.connect("/scratch2/groups/gsi/production/qcetl_v1/"+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            # print(f"ID: {id}: ")
            lims = Table.data[id][self.process]["limskeys"][0]
            entry = {}
            for row in cur.execute( # forloop automatically ignores empty SELECTs
                f"""
                select {select_block} from {self.source_table} where "Pinery Lims ID" like '%{lims}%';
                """):
                # print(row)
                # print("\n")
 
                entry["sample_id"] = id
                for column in self.columns.keys():
                    entry[column] = (
                        row[indices[column]]
                        if isinstance(row[indices[column]], str)
                        else round(row[indices[column]], 2)
                    )

            for row in cur.execute(
                f"""
                    select (CAST("mapped reads" as FLOAT) /CAST("total reads" as FLOAT)) from {self.source_table} where "Pinery Lims ID" like '%{lims}%';
                    """
                ):
                entry[RNASeqQCTableColumns.MappedReads] = row[0]

            for row in cur.execute(
               f"""
                select
                (CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT))
                from {self.source_table}
                where "Pinery Lims ID" like '%{lims}%';
                """
            ):
                entry[RNASeqQCTableColumns.RRNAContamination] = row[0]

            data.append(entry)
            # print("+++++++++++++++++++++++++++++++++")

        cur.close()
        con.close()
        return data

class RNASeqQCMergedTable(Table):
    def __init__(self):
        self.base_db_path = "/scratch2/groups/gsi/production/qcetl_v1/"
        self.title = "RNASeqQCMerged"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
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
        self.data = {}
        self.source_table = "rnaseqqc2merged_rnaseqqc2merged_2"
        self.source_db = "rnaseqqc2merged"
        self.process = "star_call_ready"

    def get_data(self):
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for id in Table.sample_ids:
            lims = "[\"" + '\", \"'.join(Table.data[id][self.process]["limskeys"]) + "\"]"

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
                entry[RNASeqQCMergedTableColumns.MappedReads] = row[0]

            for row in cur.execute(
               f"""
                select
                (CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT))
                from {self.source_table}
                where "Merged Pinery Lims ID" = '{lims}';
                """
            ):
                entry[RNASeqQCMergedTableColumns.RRNAContamination] = row[0]
            
            entry[RNASeqQCMergedTableColumns.NumLimsKeys] = len(Table.data[id][self.process]["limskeys"])
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
            RSEMTableColumns.Total: "Total",
            RSEMTableColumns.PctNonZero: "Percent of Non-Zero Records",
            RSEMTableColumns.Q0: "Q0",
            RSEMTableColumns.Q0_05: "Q0.05",
            RSEMTableColumns.Q0_1: "Q0.1",
            RSEMTableColumns.Q0_25: "Q0.25",
            RSEMTableColumns.Q0_5: "Q0.5",
            RSEMTableColumns.Q0_75: "Q0.75",
            RSEMTableColumns.Q0_9: "Q0.9",
            RSEMTableColumns.Q0_95: "Q0.95",
            RSEMTableColumns.Q1: "Q1", 
        }
        self.columns = {
            RSEMTableColumns.Total: "total",
            RSEMTableColumns.PctNonZero: "pct_non_zero",
            RSEMTableColumns.Q0: "Q0",
            RSEMTableColumns.Q0_05: "Q0.05",
            RSEMTableColumns.Q0_1: "Q0.1",
            RSEMTableColumns.Q0_25: "Q0.25",
            RSEMTableColumns.Q0_5: "Q0.5",
            RSEMTableColumns.Q0_75: "Q0.75",
            RSEMTableColumns.Q0_9: "Q0.9",
            RSEMTableColumns.Q0_95: "Q0.95",
            RSEMTableColumns.Q1: "Q1",
        }
        self.source_table = "analysis_rsem_analysis_rsem_1"
        self.source_db = "analysis_rsem"
        self.process = "rsem"

class SequenzaAltSolnTable(Table):
    def __init__(self):
        self.title = "Alternative Solutions"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            SequenzaAltSolnTableColumns.SampleID: "Sample ID",
            SequenzaAltSolnTableColumns.Index: "Index",
            SequenzaAltSolnTableColumns.Cellularity: "Cellularity",
            SequenzaAltSolnTableColumns.Ploidy: "Ploidy",
            SequenzaAltSolnTableColumns.SLPP: "SLPP",
            SequenzaAltSolnTableColumns.Gamma: "Gamma",
        }
        self.columns = {
            SequenzaAltSolnTableColumns.Index: "index",
            SequenzaAltSolnTableColumns.Cellularity: "cellularity",
            SequenzaAltSolnTableColumns.Ploidy: "ploidy",
            SequenzaAltSolnTableColumns.SLPP: "SLPP",
            SequenzaAltSolnTableColumns.Gamma: "gamma",
        }
        self.source_table = "analysis_sequenza_analysis_sequenza_alternative_solutions_1"
        self.source_db = "analysis_sequenza"
        self.process = "sequenza_by_tumor_group"

class SequenzaFGATable(Table):
    def __init__(self):
        self.title = "FGA with Gamma = 500"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
                SequenzaFGATableColumns.SampleID: "Sample ID",
                SequenzaFGATableColumns.FGA: "Fraction Genome Altered",
            }
        self.columns = {
            SequenzaFGATableColumns.FGA: "fga",
            }
        self.source_table = "analysis_sequenza_analysis_sequenza_gamma_500_fga_1"
        self.source_db = "analysis_sequenza"
        self.process = "sequenza_by_tumor_group"

class StarFusionTable(Table):
    def __init__(self):
        self.title = "Call Summary"
        self.blurb = """
            Lorem ipsum dolor sit amet, consectetur adipiscing elit.
            Curabitur sit amet commodo lorem.
            Vivamus mollis mauris urna, eget volutpat ex pharetra vel.
            Suspendisse ultrices commodo feugiat.
            """
        self.headings = {
            StarFusionTableColumns.SampleID: "Sample ID",
            StarFusionTableColumns.NumRecords: "Number of Records",
        }
        self.columns = {
            StarFusionTableColumns.NumRecords: "num_records",
        }
        self.source_table = "analysis_starfusion_analysis_starfusion_1"
        self.source_db = "analysis_starfusion"
        self.process = "starfusion"
