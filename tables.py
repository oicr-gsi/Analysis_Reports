import json
import sys
import os
import gzip
import re
import pandas as pd
import sqlite3
from table_columns import (
    CommonColumns,
    CasesTableColumns,
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

# CasesTable class defines the Cases table
class CasesTable(Table):
    def __init__(self):
        self.title = "Donors"
        self.headings = {
            CasesTableColumns.Case: "Donor",
            CasesTableColumns.GroupID: "Group ID",
            CasesTableColumns.LibraryType: "Library Type",
            CasesTableColumns.TissueType: "Tissue Type",
            CasesTableColumns.TissueOrigin: "Tissue Origin",
            CasesTableColumns.TissuePreparation: "Tissue Preparation",
            CasesTableColumns.ExternalID: "External ID",
            CasesTableColumns.SampleID: "Sample ID",
        }
        self.columns = {
            CasesTableColumns.Case: "\"Donor\"",
            CasesTableColumns.GroupID: "\"Group ID\"",
            CasesTableColumns.LibraryType: "\"Library Type\"",
            CasesTableColumns.TissueType: "\"Tissue Type\"",
            CasesTableColumns.TissueOrigin: "\"Tissue Origin\"",
            CasesTableColumns.TissuePreparation: "\"Tissue Preparation\"",
            CasesTableColumns.ExternalID: "\"External ID\"",
            CasesTableColumns.SampleID: "\"Sample ID\"",
        }
        self.source_file = "/scratch2/groups/gsi/production/vidarr/vidarr_files_report_latest.tsv.gz" 
        self.glossary = {}
        # Glossary templates
        self.tissue_types = {
            'X': 'Xenograft derived from some tumour.',
            'U': 'Unspecified', 
            'T': 'Unclassified tumour', 
            'S': 'Serum from blood where clotting proteins have been removed',
            'R': 'Reference or non-tumour, non-diseased tissue sample.',
            'P': 'Primary tumour', 
            'O': 'Organoid', 
            'n': 'Unknown', 
            'M': 'Metastatic tumour',
            'F': 'Fibroblast cells', 
            'E': 'Endothelial cells', 
            'C': 'Cell line derived from a tumour',
            'B': 'Benign tumour', 
            'A': 'Cells taken from Ascites fluid'
        }
        self.tissue_origin = {
            'Ab': 'Abdomen', 'Ad': 'Adipose', 'Ae': 'Adnexa', 'Ag': 'Adrenal', 'An': 'Anus',
            'Ao': 'Anorectal', 'Ap': 'Appendix', 'As': 'Ascites', 'At': 'Astrocytoma', 'Av': 'Ampulla',
            'Ax': 'Axillary', 'Ba': 'Back', 'Bd': 'Bile', 'Bi': 'Biliary', 'Bl': 'Bladder',
            'Bm': 'Bone', 'Bn': 'Brain', 'Bo': 'Bone', 'Br': 'Breast', 'Bu': 'Buccal',
            'Bw': 'Bowel', 'Cb': 'Cord', 'Cc': 'Cecum', 'Ce': 'Cervix', 'Cf': 'Cell-Free', 'Ch': 'Chest',
            'Cj': 'Conjunctiva', 'Ck': 'Cheek', 'Cn': 'Central', 'Co': 'Colon', 'Cr': 'Colorectal',
            'Cs': 'Cul-de-sac', 'Ct': 'Circulating', 'Di': 'Diaphragm', 'Du': 'Duodenum',
            'En': 'Endometrial', 'Ep': 'Epidural', 'Es': 'Esophagus', 'Ey': 'Eye', 'Fa': 'Fallopian',
            'Fb': 'Fibroid', 'Fs': 'Foreskin', 'Ft': 'Foot', 'Ga': 'Gastric', 'Gb': 'Gallbladder',
            'Ge': 'Gastroesophageal', 'Gi': 'Gastrointestinal', 'Gj': 'Gastrojejunal', 'Gn': 'Gingiva',
            'Gt': 'Genital', 'Hp': 'Hypopharynx', 'Hr': 'Heart', 'Ic': 'Ileocecum', 'Il': 'Ileum',
            'Ki': 'Kidney', 'La': 'Lacrimal', 'Lb': 'Limb', 'Le': 'Leukocyte', 'Lg': 'Leg',
            'Li': 'Large', 'Ln': 'Lymph', 'Lp': 'Lymphoblast', 'Lu': 'Lung', 'Lv': 'Liver', 
            'Lx': 'Larynx', 'Ly': 'Lymphocyte', 'Md': 'Mediastinum', 'Me': 'Mesenchyme', 'Mn': 'Mandible',
            'Mo': 'Mouth', 'Ms': 'Mesentary', 'Mu': 'Muscle', 'Mx': 'Maxilla', 'Nk': 'Neck',
            'nn': 'Unknown', 'No': 'Nose', 'Np': 'Nasopharynx', 'Oc': 'Oral', 'Om': 'Omentum',
            'Or': 'Orbit', 'Ov': 'Ovary', 'Pa': 'Pancreas', 'Pb': 'Peripheral', 'Pc': 'Pancreatobiliary',
            'Pd': 'Parathyroid', 'Pe': 'Pelvic', 'Pg': 'Parotid', 'Ph': 'Paratracheal', 'Pi': 'Penis',
            'Pl': 'Plasma', 'Pm': 'Peritoneum', 'Pn': 'Peripheral', 'Po': 'Peri-aorta', 'Pr': 'Prostate',
            'Pt': 'Palate', 'Pu': 'Pleura', 'Py': 'Periampullary', 'Ra': 'Right', 'Rc': 'Rectosigmoid',
            'Re': 'Rectum', 'Ri': 'Rib', 'Rp': 'Retroperitoneum', 'Sa': 'Saliva', 'Sb': 'Small',
            'Sc': 'Scalp', 'Se': 'Serum', 'Sg': 'Salivary', 'Si': 'Small', 'Sk': 'Skin', 'Sm': 'Skeletal',
            'Sn': 'Spine', 'So': 'Soft', 'Sp': 'Spleen', 'Sr': 'Serosa', 'Ss': 'Sinus', 'St': 'Stomach',
            'Su': 'Sternum', 'Ta': 'Tail', 'Te': 'Testes', 'Tg': 'Thymic', 'Th': 'Thymus',
            'Tn': 'Tonsil', 'To': 'Throat', 'Tr': 'Trachea', 'Tu': 'Tongue', 'Ty': 'Thyroid',
            'Uc': 'Urachus', 'Ue': 'Ureter', 'Um': 'Umbilical', 'Up': 'Urine', 'Ur': 'Urethra',
            'Us': 'Urine', 'Ut': 'Uterus', 'Uw': 'Urine', 'Vg': 'Vagina', 'Vu': 'Vulva', 'Wm': 'Worm'
        }
        self.library_design = {
            'WT': 'Whole Transcriptome', 'WG': 'Whole Genome', 'TS': 'Targeted Sequencing',
            'TR': 'Total RNA', 'SW': 'Shallow Whole Genome', 'SM': 'smRNA', 'SC': 'Single Cell',
            'NN': 'Unknown', 'MR': 'mRNA', 'EX': 'Exome', 'CT': 'ctDNA', 'CM': 'cfMEDIP',
            'CH': 'ChIP-Seq', 'BS': 'Bisulphite Sequencing', 'AS': 'ATAC-Seq'
        }

    def load_context(self, workflow_ids, base_db_path):
        '''
        Load the context for the Cases Table by querying the TSV file.
        '''
        context = {
            "title": self.title,
            "headings": self.headings,
            "columns": self.columns,
            "glossary": self.glossary,
        }

        # Query the TSV file and get the cases data
        cases_data = query_provenance_file(self.source_file, workflow_ids)
        cases_data['Sample ID'] = cases_data.apply(
            lambda row: f"{row['Donor']}_{row['Tissue Origin']}_{row['Tissue Type']}_{row['Library Type']}_{row['Group ID']}",
            axis=1
        )

        # Get the unique values present in the dataset for each category
        ttypes = cases_data['Tissue Type'].unique()
        torigins = cases_data['Tissue Origin'].unique()
        ltypes = cases_data['Library Type'].unique()

        # For each type, get the full form description
        ttype = {key: self.tissue_types.get(key, 'Unknown') for key in ttypes}
        torigin = {key: self.tissue_origin.get(key, 'Unknown') for key in torigins}
        ltype = {key: self.library_design.get(key, 'Unknown') for key in ltypes}

        # Update glossary with full descriptions based on the data in the table
        self.glossary[CasesTableColumns.TissueType] = "\n".join([f"{k}: {v}" for k, v in ttype.items()])
        self.glossary[CasesTableColumns.TissueOrigin] = "\n".join([f"{k}: {v}" for k, v in torigin.items()])
        self.glossary[CasesTableColumns.LibraryType] = "\n".join([f"{k}: {v}" for k, v in ltype.items()])
        
        context["data"] = cases_data.to_dict(orient='records')
        return context

# DellyTable class defines a table for the delly workflow
class DellyTable(Table):
    def __init__(self):
        self.title = "Genomic Structural Variants"
        self.headings = {
            DellyTableColumns.Case: "Donor",
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
            Mutect2TableColumns.Case: "Donor",
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
            RSEMTableColumns.Case: "Donor",
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
            StarFusionTableColumns.Case: "Donor",
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

def query_provenance_file(file_provenance_path, workflow_ids):
    workflow = re.compile('|'.join(workflow_ids))
    case = []
    with gzip.open(file_provenance_path, 'rt') as f:
        header = f.readline().strip().split('\t')
        col_indices = {
            'Workflow Run SWID': header.index('Workflow Run SWID'),
            'Root Sample Name': header.index('Root Sample Name'),
            'Sample Attributes': header.index('Sample Attributes')
        }

        for line in f:
            row = line.strip().split('\t')
            if any(workflow_id in row[col_indices['Workflow Run SWID']] for workflow_id in workflow_ids):
                workflow_swid = row[col_indices['Workflow Run SWID']]
                donor = row[col_indices['Root Sample Name']]
                sample_attributes = row[col_indices['Sample Attributes']]

                metadata_dict = {}
                for item in sample_attributes.split(';'):
                    if '=' in item:
                        key, value = item.split('=', 1)
                        metadata_dict[key] = value

                case.append({
                    'Donor': donor,
                    'Group ID': metadata_dict.get('geo_group_id'),
                    'Library Type': metadata_dict.get('geo_library_source_template_type'),
                    'Tissue Type': metadata_dict.get('geo_tissue_type'),
                    'Tissue Origin': metadata_dict.get('geo_tissue_origin'),
                    'Tissue Preparation': metadata_dict.get('geo_tissue_preparation'),
                    'External ID': metadata_dict.get('geo_external_name'),
                })

    # Convert the list of cases to a DataFrame
    cases = pd.DataFrame(case).drop_duplicates()
    return cases
