import json
import sys
import os
import pandas as pd
import sqlite3
from table_columns import (
    CommonColumns,
    DellyTableColumns,
    Mutect2TableColumns,
    RSEMTableColumns,
    StarFusionTableColumns,
)

# The Table class defines each table that is generated
class Table:
    base_db_path: str       # base path to databases that are queried
    title: str              # title of table
    blurb = ""              # descriptive blurb of table
    headings: dict          # headings for each column to be displayed on table
    columns: dict           # columns we want from sql table. Must match EXACTLY
    source_table: str       # table we query from
    source_db: str          # database we query from
    glossary: dict          # Dict[name of column, definition]

    def load_context(self, workflow_ids, base_db_path):
        '''
        None -> dict[str, Any]
        
        Returns a dict used for loading the html in jinja2 templating
        '''
        context = {
            "title": self.title,
            "blurb": self.blurb,
            "headings": self.headings,
            "columns": self.columns,
            "glossary": self.glossary,
        }

        # Call to extract metrics and get data
        table_data = extract_metrics(self.__class__, workflow_ids, base_db_path) 
        context["data"] = table_data.to_dict(orient='records')

        return context

# DellyTable class defines a table for the delly workflow
class DellyTable(Table):
    def __init__(self):
        self.title = "Genomic Structural Variants"
        self.headings = {
            DellyTableColumns.Case: "Case",
            DellyTableColumns.NumCalls: "SV Calls",
            DellyTableColumns.NumPASS: "SV PASS Calls",
            DellyTableColumns.NumBND: "Translocations",
            DellyTableColumns.NumDEL: "Deletions",
            DellyTableColumns.NumDUP: "Duplications",
            DellyTableColumns.NumINS: "Insertions",
            DellyTableColumns.NumINV: "Inversions",
        }
        self.columns = {
            DellyTableColumns.Case: "\"Donor\"",
            DellyTableColumns.NumCalls: "\"num_calls\"",
            DellyTableColumns.NumPASS: "\"num_PASS\"",
            DellyTableColumns.NumBND: "\"num_BND\"",
            DellyTableColumns.NumDEL: "\"num_DEL\"",
            DellyTableColumns.NumDUP: "\"num_DUP\"",
            DellyTableColumns.NumINS: "\"num_INS\"",
            DellyTableColumns.NumINV: "\"num_INV\"",
        }
        self.source_table = ["analysis_delly_analysis_delly_1"]
        self.source_db = "analysis_delly"
        self.glossary = {
            DellyTableColumns.NumCalls: "The number of somatic structural variant calls identified by delly",
            DellyTableColumns.NumPASS: "The number of structural variant calls marked as PASS",
            DellyTableColumns.NumBND: "The number of PASS translocation calls",
            DellyTableColumns.NumDEL: "The number of PASS deletion calls",
            DellyTableColumns.NumDUP: "The number of PASS duplication calls",
            DellyTableColumns.NumINS: "The number of PASS insertions calls",
            DellyTableColumns.NumINV: "The number of PASS inversions calls",
        }

# Mutect2Table class defines a table for the Mutect2 workflow
class Mutect2Table(Table):
    def __init__(self):
        self.title = "Somatic Mutations"
        self.headings = {
            Mutect2TableColumns.Case: "Case",
            Mutect2TableColumns.NumCalls: "Calls",
            Mutect2TableColumns.NumPASS: "PASS Calls",
            Mutect2TableColumns.NumSNPs: "SNPs",
            Mutect2TableColumns.NumIndels: "Indels",
            Mutect2TableColumns.TITVRatio: "Ti/Tv Ratio",
        }
        self.columns = {
            Mutect2TableColumns.Case: "\"Donor\"",
            Mutect2TableColumns.NumCalls: "\"num_calls\"",
            Mutect2TableColumns.NumPASS: "\"num_PASS\"",
            Mutect2TableColumns.NumSNPs: "\"num_SNPs\"",
            Mutect2TableColumns.NumIndels: "\"num_indels\"",
            Mutect2TableColumns.TITVRatio: "\"titv_ratio\"",
        }
        self.source_table = ["analysis_mutect2_analysis_mutect2_1"]
        self.source_db = "analysis_mutect2"
        self.glossary = {
            Mutect2TableColumns.NumCalls: "Total number of calls identified by Mutect2",
            Mutect2TableColumns.NumPASS: "Total number of PASS calls",
            Mutect2TableColumns.NumSNPs: "Number of SNP calls",
            Mutect2TableColumns.NumIndels: "Number of insertion/deletion calls",
            Mutect2TableColumns.TITVRatio: "Ratio of transition vs transversion mutations",
        }

# RSEMTable class defines a table for the rsem workflow
class RSEMTable(Table):
    def __init__(self):
        self.title = "Gene Expression"
        self.headings = {
            RSEMTableColumns.Case: "Case",
            RSEMTableColumns.Total: "Total Reads",
            RSEMTableColumns.PctNonZero: "Percent Non-zero",
            RSEMTableColumns.Q0_05: "0.05 Quantile",
            RSEMTableColumns.Q0_5: "0.5 Quantile",
            RSEMTableColumns.Q0_95: "0.95 Quantile",
        }
        self.columns = {
            RSEMTableColumns.Case: "\"Donor\"",
            RSEMTableColumns.Total: "\"total\"",
            RSEMTableColumns.PctNonZero: "\"pct_non_zero\"",
            RSEMTableColumns.Q0_05: "\"Q0.05\"",
            RSEMTableColumns.Q0_5: "\"Q0.5\"",
            RSEMTableColumns.Q0_95: "\"Q0.95\"",
        }
        self.source_table = ["analysis_rsem_analysis_rsem_1"]
        self.source_db = "analysis_rsem"
        self.glossary = {
            RSEMTableColumns.Total: "Total number of reads assigned to a gene",
            RSEMTableColumns.PctNonZero: "Percentage of non-zero read counts",
            RSEMTableColumns.Q0_05: "Expression at the 0.05 quantile",
            RSEMTableColumns.Q0_5: "Expression at the 0.5 quantile",
            RSEMTableColumns.Q0_95: "Expression at the 0.95 quantile",
        }

# StarFusionTable class defines a table for the star fusion workflow
class StarFusionTable(Table):
    def __init__(self):
        self.title = "Gene Fusions"
        self.headings = {
            StarFusionTableColumns.Case: "Case",
            StarFusionTableColumns.NumRecords: "Fusion Calls",
        }
        self.columns = {
            StarFusionTableColumns.Case: "\"Donor\"",
            StarFusionTableColumns.NumRecords: "\"num_records\"",
        }
        self.source_table = ["analysis_starfusion_analysis_starfusion_1"]
        self.source_db = "analysis_starfusion"
        self.glossary = {
            StarFusionTableColumns.NumRecords: "Number of gene fusions identified by StarFusion",
        }  


def extract_metrics(table_class, workflow_ids, base_db_path):
    '''
    Fetch data from the database based on the provided table class, workflow_ids, and base_db_path.
    '''
    table_obj = table_class()

    con = sqlite3.connect(base_db_path + table_obj.source_db + "/latest")
    cur = con.cursor()

    query = f'''
    SELECT {', '.join(table_obj.columns.values())}
    FROM {table_obj.source_table[0]}
    WHERE "Workflow Run SWID" LIKE "vidarr:%/run/%" AND "Workflow Run SWID" LIKE ?;
    '''
    
    extracted_metrics = []

    for workflow_id in workflow_ids:
        cur.execute(query, (f"%{workflow_id}",))
        rows = cur.fetchall()
        extracted_metrics.extend(rows)

    columns = [desc[0] for desc in cur.description]
    res = pd.DataFrame(extracted_metrics, columns=columns)
    res.columns = [table_obj.headings.get(col, col) for col in columns]

    cur.close()
    con.close()

    return res
