import pandas as pd
import sqlite3
from typing import List
from table_columns import (
    CommonColumns,
    CasesTableColumns,
    DellyTableColumns,
    Mutect2TableColumns,
    RSEMTableColumns,
    StarFusionTableColumns,
    WGCallReadyTableColumns,
    WGLaneLevelTableColumns,
    WTCallReadyTableColumns,
    WTLaneLevelTableColumns,
)
from plot import(
    Plot,
)

NUM_DP = 2
# The Table class defines each table that is generated
class Table:
    base_db_path: str       # base path to databases that are queried
    title: str              # title of table
    blurb = ""              # descriptive blurb of table
    headings: dict          # headings for each column to be displayed on table
    columns: dict           # columns we want from sql table. Must match EXACTLY
    source_table: str       # table we query from
    source_db: str          # database we query from
    process: List[str]      # workflow names
    plots = {}
    glossary: dict          # Dict[name of column, definition]
    pct_stats = set()
    pipeline_step: str

    def load_context(self, workflow_ids, base_db_path, cases_data):
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
        table_data = get_metrics(self.__class__, workflow_ids, base_db_path) 
        table_data.columns = table_data.columns.str.strip('"')
        table_data = table_data.drop(columns=['SampleID'])

        cases_data = cases_data.drop(columns=['Workflow Run SWID', 'LIMS ID']).drop_duplicates()

        # Get SampleIDs from FPR
        table_data = table_data.merge(cases_data[['Donor', 'SampleID']], on='Donor', how='left')
        column_order = ['Donor', 'SampleID'] + [col for col in table_data.columns if col not in ['Donor', 'SampleID']]
        table_data = table_data[column_order]

        if table_data.empty:
            context["data"] = []
        else:  
            context["data"] = table_data.to_dict(orient='records')

        return context

# CasesTable class defines the Cases table
class CasesTable(Table):
    def __init__(self):
        self.title = "Donors"
        self.headings = {
            CasesTableColumns.Case: "Donor",
            CasesTableColumns.SampleID: "SampleID",
            CasesTableColumns.GroupID: "Group ID",
            CasesTableColumns.LibraryType: "Library Type",
            CasesTableColumns.TissueType: "Tissue Type",
            CasesTableColumns.TissueOrigin: "Tissue Origin",
            CasesTableColumns.TissuePreparation: "Tissue Preparation",
            CasesTableColumns.ExternalID: "External ID",
        }
        self.columns = {
            CasesTableColumns.Case: "\"Donor\"",
            CasesTableColumns.SampleID: "\"SampleID\"",
            CasesTableColumns.GroupID: "\"Group ID\"",
            CasesTableColumns.LibraryType: "\"Library Type\"",
            CasesTableColumns.TissueType: "\"Tissue Type\"",
            CasesTableColumns.TissueOrigin: "\"Tissue Origin\"",
            CasesTableColumns.TissuePreparation: "\"Tissue Preparation\"",
            CasesTableColumns.ExternalID: "\"External ID\"",
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

    def load_context(self, workflow_ids, base_db_path, cases_data):
        '''
        Load the context for the Cases Table. 
        '''
        context = {
            "title": self.title,
            "headings": self.headings,
            "columns": self.columns,
            "glossary": self.glossary,
        }

        cases_data = cases_data.drop(columns=['Workflow Run SWID', 'LIMS ID']).drop_duplicates()
        column_order = ['Donor', 'SampleID'] + [col for col in cases_data.columns if col not in ['Donor', 'SampleID']]
        cases_data = cases_data[column_order]

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
        
        if cases_data.empty:
            context["data"] = []
        else:  
            context["data"] = cases_data.to_dict(orient='records')

        return context

# DellyTable class defines a table for the delly workflow
class DellyTable(Table):
    def __init__(self):
        self.title = "Genomic Structural Variants"
        self.headings = {
            DellyTableColumns.Case: "Donor",
            DellyTableColumns.SampleID: "SampleID",
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
            DellyTableColumns.SampleID: "\"SampleID\"",
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
        self.process = ["delly_matched_by_tumor_group", "delly"]
        self.plots = {
            DellyTableColumns.NumPASS: Plot(
                title="SV PASS Calls",
                x_axis="SampleIDs",
                y_axis="SV PASS Calls"
            ),
        }
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
            Mutect2TableColumns.SampleID: "SampleID",
            Mutect2TableColumns.NumCalls: "Calls",
            Mutect2TableColumns.NumPASS: "PASS Calls",
            Mutect2TableColumns.NumSNPs: "SNPs",
            Mutect2TableColumns.NumIndels: "Indels",
            Mutect2TableColumns.TITVRatio: "Ti/Tv Ratio",
        }
        self.columns = {
            Mutect2TableColumns.Case: "\"Donor\"",
            Mutect2TableColumns.SampleID: "\"SampleID\"",
            Mutect2TableColumns.NumCalls: "\"num_calls\"",
            Mutect2TableColumns.NumPASS: "\"num_PASS\"",
            Mutect2TableColumns.NumSNPs: "\"num_SNPs\"",
            Mutect2TableColumns.NumIndels: "\"num_indels\"",
            Mutect2TableColumns.TITVRatio: "\"titv_ratio\"",
        }
        self.pipeline_step = "calls.mutations"
        self.source_table = ["analysis_mutect2_analysis_mutect2_1"]
        self.source_db = "analysis_mutect2"
        self.process = ["mutect2_matched_by_tumor_group", "mutect2"]
        self.plots = {
            Mutect2TableColumns.NumPASS: Plot(
                title="Mutation Calls",
                x_axis="SampleIDs",
                y_axis="Mutation Calls"
            ),
            Mutect2TableColumns.TITVRatio: Plot(
                title="Ti/Tv",
                x_axis="SampleIDs",
                y_axis="Ti/Tv",
                lo=0,
            ),
        }
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
            RSEMTableColumns.SampleID: "SampleID",
            RSEMTableColumns.Total: "Total Reads",
            RSEMTableColumns.PctNonZero: "Percent Non-zero",
            RSEMTableColumns.Q0_05: "0.05 Quantile",
            RSEMTableColumns.Q0_5: "0.5 Quantile",
            RSEMTableColumns.Q0_95: "0.95 Quantile",
        }
        self.columns = {
            RSEMTableColumns.Case: "\"Donor\"",
            RSEMTableColumns.SampleID: "\"SampleID\"",
            RSEMTableColumns.Total: "\"total\"",
            RSEMTableColumns.PctNonZero: "\"pct_non_zero\"",
            RSEMTableColumns.Q0_05: "\"Q0.05\"",
            RSEMTableColumns.Q0_5: "\"Q0.5\"",
            RSEMTableColumns.Q0_95: "\"Q0.95\"",
        }
        self.pipeline_step = "calls.expression"
        self.source_table = ["analysis_rsem_analysis_rsem_1"]
        self.source_db = "analysis_rsem"
        self.process = ["rsem"]
        self.pct_stats = set(
            [
                RSEMTableColumns.PctNonZero,
            ]
        )
        self.plots = {
            "pct_non_zero": Plot(
                title="Percent Expressed",
                x_axis="SampleID",
                y_axis="Percent Expressed (%)",
                hi=100,
                lo=0,
            ),
            "Q0.5": Plot(
                title="Median TPM",
                x_axis="SampleID",
                y_axis="Median TPM",
            ),
        }
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
            StarFusionTableColumns.SampleID: "SampleID",
            StarFusionTableColumns.NumRecords: "Fusion Calls",
        }
        self.columns = {
            StarFusionTableColumns.Case: "\"Donor\"",
            StarFusionTableColumns.SampleID: "\"SampleID\"",
            StarFusionTableColumns.NumRecords: "\"num_records\"",
        }
        self.pipeline_step = "calls.fusions"
        self.source_table = ["analysis_starfusion_analysis_starfusion_1"]
        self.source_db = "analysis_starfusion"
        self.process = ["starfusion", "starFusion"]
        self.plots = {
            StarFusionTableColumns.NumRecords: Plot(
                title="Fusion Calls",
                x_axis="SampleID",
                y_axis="Fusion Calls"
            )
        }
        self.glossary = {
            StarFusionTableColumns.NumRecords: "Number of gene fusions identified by StarFusion",
        }  

#WGCallReadyTable class defines a table for the bamqc4merged workflow
class WGCallReadyTable(Table):
    def __init__(self):
        self.title = "Whole Genome Libraries, tumour and matched normal"
        self.headings = {
            WGCallReadyTableColumns.Case: "Donor",
            WGCallReadyTableColumns.SampleID: "SampleID",
            WGCallReadyTableColumns.CoverageDedup: "Coverage Depth",
            WGCallReadyTableColumns.MarkDupPctDup: "Duplication (%)",
            WGCallReadyTableColumns.TotalClusters: "Read Pairs",
            WGCallReadyTableColumns.MappedReads: "Mapped Reads (%)",
        }
        self.columns = {
            WGCallReadyTableColumns.Case: "\"Donor\"",
            WGCallReadyTableColumns.SampleID: "\"SampleID\"",
            WGCallReadyTableColumns.CoverageDedup: "\"coverage deduplicated\"",
            WGCallReadyTableColumns.MarkDupPctDup: "\"mark duplicates_PERCENT_DUPLICATION\"",
            WGCallReadyTableColumns.TotalClusters: "\"total clusters\"",
            WGCallReadyTableColumns.MappedReads: """
            ROUND((1 - CAST("unmapped reads meta" as FLOAT) / CAST("total input reads meta" as FLOAT)) * 100, 2)
            """,
        }
        self.source_table = ["bamqc4merged_bamqc4merged_5"]
        self.source_db = "bamqc4merged"
        self.glossary = {
            WGCallReadyTableColumns.CoverageDedup: "Mean depth of coverage corrected for duplication",
            WGCallReadyTableColumns.MarkDupPctDup: "Percent of reads marked as duplicates",
            WGCallReadyTableColumns.TotalClusters: "Number of read pairs generated",
            WGCallReadyTableColumns.MappedReads: "Percent of reads mapping to the genomic reference",
        }
    
    def load_context(self, workflow_ids, base_db_path, cases_data):
        '''
        Load the context for the WG Call Ready Table. 
        '''
        context = {
            "title": self.title,
            "headings": self.headings,
            "columns": self.columns,
            "glossary": self.glossary,
        }
        table_data = CallReady_metrics(self.__class__, cases_data, base_db_path) 
        
        if table_data.empty:
            context["data"] = []
        else:  
            context["data"] = table_data.to_dict(orient='records')

        return context

#WGLaneLevelTable class defines a table for the bamqc4 workflow
class WGLaneLevelTable(Table):
    def __init__(self):
        self.title = "Whole Genome Libraries, tumour and matched normal"
        self.headings = {
            WGLaneLevelTableColumns.Case: "Donor",
            WGLaneLevelTableColumns.SampleID: "SampleID",
            WGLaneLevelTableColumns.CoverageDedup: "Coverage Depth",
            WGLaneLevelTableColumns.InsertSizeAvg: "Insert Size",
            WGLaneLevelTableColumns.MarkDupPctDup: "Duplication (%)",
            WGLaneLevelTableColumns.TotalClusters: "Read Pairs",
            WGLaneLevelTableColumns.MappedReads: "Mapped Reads (%)",
        }
        self.columns = {
            WGLaneLevelTableColumns.Case: "\"Donor\"",
            WGLaneLevelTableColumns.SampleID: "\"sample\"",
            WGLaneLevelTableColumns.CoverageDedup: "\"coverage deduplicated\"",
            WGLaneLevelTableColumns.InsertSizeAvg: "\"insert size average\"",
            WGLaneLevelTableColumns.MarkDupPctDup: "\"mark duplicates_PERCENT_DUPLICATION\"",
            WGLaneLevelTableColumns.TotalClusters: "\"total clusters\"",
            WGLaneLevelTableColumns.MappedReads: """
            ROUND((1 - CAST("unmapped reads meta" as FLOAT) / CAST("total input reads meta" as FLOAT)) * 100, 2)
            """,
        }
        self.source_table = ["bamqc4_bamqc4_5"]
        self.source_db = "bamqc4"
        self.glossary = {
            WGLaneLevelTableColumns.CoverageDedup: "Mean depth of coverage corrected for duplication",
            WGLaneLevelTableColumns.InsertSizeAvg: "Mean size of the sequenced insert",
            WGLaneLevelTableColumns.MarkDupPctDup: "Percent of reads marked as duplicates",
            WGLaneLevelTableColumns.TotalClusters: "Number of read pairs generated",
            WGLaneLevelTableColumns.MappedReads: "Percent of reads mapping to the genomic reference",

        }
    
    def get_data(self, cases_data, base_db_path):
        '''
        Fetch metrics data for WG libraries and map Donor using SampleID.
        '''
        sample_ids = cases_data['SampleID'].tolist()
        con = sqlite3.connect(base_db_path + self.source_db + "/latest")
        cur = con.cursor()

        query = f'''
        SELECT {', '.join(self.columns.values())}
        FROM {self.source_table[0]}
        WHERE "sample" = ?
        '''

        extracted_metrics = []
        for sample_id in sample_ids:
            cur.execute(query, (sample_id.strip(),))
            rows = cur.fetchall()
            extracted_metrics.extend(rows)

        columns = [desc[0] for desc in cur.description]
        res = pd.DataFrame(extracted_metrics, columns=columns)
        res.columns = [self.headings.get(col, col) for col in columns]

        cur.close()
        con.close()

        res.columns = res.columns.str.strip('"')
        res = res.drop(columns=['Donor'], errors='ignore')
        res.rename(columns={res.columns[0]: 'SampleID'}, inplace=True)

        # Prepare donor mapping and merge
        cases_data = cases_data.drop(columns=['Workflow Run SWID', 'LIMS ID'], errors='ignore').drop_duplicates()
        donor_map = cases_data.set_index('SampleID')['Donor']
        res['Donor'] = res['SampleID'].map(donor_map)

        # Reorder columns
        col_order = ['Donor', 'SampleID'] + [col for col in res.columns if col not in ['Donor', 'SampleID']]
        res = res[col_order].drop_duplicates()

        return res

    def load_context(self, workflow_ids, base_db_path, cases_data):
        '''
        Load the context for the WG Lane Level Table. 
        '''
        context = {
            "title": self.title,
            "headings": self.headings,
            "columns": self.columns,
            "glossary": self.glossary,
        }

        table_data = self.get_data(cases_data, base_db_path)
        
        if table_data.empty:
            context["data"] = []
        else:  
            context["data"] = table_data.to_dict(orient='records')

        return context


#WTCallReadyTable class defines a table for the rnaseqqc2merged workflow
class WTCallReadyTable(Table):
    def __init__(self):
        self.title = "Whole Transcriptome Libraries, tumour only"
        self.headings = {
            WTCallReadyTableColumns.Case: "Donor",
            WTCallReadyTableColumns.SampleID: "SampleID",
            WTCallReadyTableColumns.PctCodingBases: "Percent Coding (%)",
            WTCallReadyTableColumns.TotalClusters: "Read Pairs",
            WTCallReadyTableColumns.MappedReads: "Mapped Reads (%)",
            WTCallReadyTableColumns.RRNAContamination: "rRNA Contamination (%)",
        }
        self.columns = {
            WTCallReadyTableColumns.Case: "\"Donor\"",
            WTCallReadyTableColumns.SampleID: "\"SampleID\"",
            WTCallReadyTableColumns.PctCodingBases: "\"PCT_CODING_BASES\"",
            WTCallReadyTableColumns.TotalClusters: "\"total clusters\"",
            WTCallReadyTableColumns.MappedReads: """
            ROUND((1 - CAST("unmapped reads" as FLOAT)/CAST("total reads" as FLOAT)) * 100,2)
            """,
            WTCallReadyTableColumns.RRNAContamination: """
                ROUND((CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT)), 2)
            """,
        }
        self.source_table = ["rnaseqqc2merged_rnaseqqc2merged_3"]
        self.source_db = "rnaseqqc2merged"
        self.glossary = {
            WTCallReadyTableColumns.PctCodingBases: "Percentage of bases mapping to the coding regions of the genome",
            WTCallReadyTableColumns.TotalClusters: "Number of read pairs generated",
            WTCallReadyTableColumns.MappedReads: "Percentage of reads mapping to the genomic reference",
            WTCallReadyTableColumns.RRNAContamination: "Pecentage of reads mapping to ribosomal RNA",
        }

    def load_context(self, workflow_ids, base_db_path, cases_data):
        '''
        Load the context for the WT Call Ready Table. 
        '''
        context = {
            "title": self.title,
            "headings": self.headings,
            "columns": self.columns,
            "glossary": self.glossary,
        }
        table_data = CallReady_metrics(self.__class__, cases_data, base_db_path) 
        
        if table_data.empty:
            context["data"] = []
        else:  
            context["data"] = table_data.to_dict(orient='records')

        return context

#WTLaneLevelTable class defines a table for the rnaseqqc2 workflow
class WTLaneLevelTable(Table):
    def __init__(self):
        self.title = "Whole Transcriptome Libraries, tumour only"
        self.blurb = ""
        self.headings = {
            WTLaneLevelTableColumns.Case: "Donor",
            WTLaneLevelTableColumns.SampleID: "SampleID",
            WTLaneLevelTableColumns.PctCodingBases: "Percent Coding (%)",
            WTLaneLevelTableColumns.TotalClusters: "Read Pairs",
            WTLaneLevelTableColumns.MappedReads: "Mapped Reads (%)",
            WTLaneLevelTableColumns.RRNAContamination: "rRNA Contamination (%)",

        }
        self.columns = {
            WTLaneLevelTableColumns.Case: "\"Donor\"",
            WTLaneLevelTableColumns.SampleID: "\"sample\"",
            WTLaneLevelTableColumns.PctCodingBases: "\"PCT_CODING_BASES\"",
            WTLaneLevelTableColumns.TotalClusters: "\"total clusters\"",
            WTLaneLevelTableColumns.MappedReads: """
            ROUND((1 - CAST("unmapped reads" as FLOAT)/CAST("total reads" as FLOAT)) * 100,2)
            """,
            WTLaneLevelTableColumns.RRNAContamination: """
                ROUND((CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT)), 2)
            """,
        }
        self.source_table = ["rnaseqqc2_rnaseqqc2_3"]
        self.source_db = "rnaseqqc2"
        self.glossary = {
            WTLaneLevelTableColumns.PctCodingBases: "Percentage of bases mapping to the coding regions of the genome",
            WTLaneLevelTableColumns.TotalClusters: "Number of read pairs generated",
            WTLaneLevelTableColumns.MappedReads: "Percentage of reads mapping to the genomic reference",
            WTLaneLevelTableColumns.RRNAContamination: "Pecentage of reads mapping to ribosomal RNA",

        }
    
    def get_data(self, cases_data, base_db_path):
        '''
        Fetch metrics data for WT libraries and map Donor and SampleID using LIMS ID.
        '''
        lims_ids = cases_data['LIMS ID']

        con = sqlite3.connect(base_db_path + self.source_db + "/latest")
        cur = con.cursor()

        query = f'''
        SELECT {', '.join(self.columns.values())}
        FROM {self.source_table[0]}
        WHERE "Pinery Lims ID" = ?
        '''

        extracted_metrics = []
        lims_id_tracker = []

        for lims_id_list in lims_ids:
            lims_id_values = lims_id_list.split(',')

            for lims_id in lims_id_values:
                lims_id = lims_id.strip()
                cur.execute(query, (lims_id,))
                rows = cur.fetchall()
                extracted_metrics.extend(rows)
                lims_id_tracker.extend([lims_id] * len(rows))  

        if not extracted_metrics:
            cur.close()
            con.close()
            return pd.DataFrame()  

        columns = [desc[0] for desc in cur.description]
        res = pd.DataFrame(extracted_metrics, columns=columns)
        res.columns = [self.headings.get(col, col) for col in columns]

        res.columns = res.columns.str.strip('"') 
        res = res.drop(columns=['Donor', 'sample'])
        res['LIMS ID'] = lims_id_tracker  

        cur.close()
        con.close()

        # Explode cases_data so each LIMS ID gets its own row
        cases = cases_data.copy()
        cases['LIMS ID'] = cases['LIMS ID'].str.split(',')
        cases = cases.explode('LIMS ID')
        cases['LIMS ID'] = cases['LIMS ID'].str.strip()

        # Merge Donor and SampleID from exploded cases_data using LIMS ID
        cases = cases[['LIMS ID', 'Donor', 'SampleID']].drop_duplicates()
        res = res.merge(cases, on='LIMS ID', how='left')
        res.drop(columns=['LIMS ID'], inplace=True)

        # Reorder columns
        column_order = ['Donor', 'SampleID'] + [col for col in res.columns if col not in ['Donor', 'SampleID', 'LIMS ID']]
        res = res[column_order].drop_duplicates()

        return res

    def load_context(self, workflow_ids, base_db_path, cases_data):
        '''
        Load the context for the WT Call Ready Table. 
        '''
        context = {
            "title": self.title,
            "headings": self.headings,
            "columns": self.columns,
            "glossary": self.glossary,
        }
        table_data = self.get_data(cases_data, base_db_path)
        
        if table_data.empty:
            context["data"] = []
        else:  
            context["data"] = table_data.to_dict(orient='records')

        return context



def get_metrics(table_class, workflow_ids, base_db_path):
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


def CallReady_metrics(table_class, cases_data, base_db_path):
    '''
    Fetch data from the WG/WT caches based on the provided table class, workflow_ids, and base_db_path.
    '''
    table_obj = table_class()
    cases = []
    cases = cases_data['Donor']

    con = sqlite3.connect(base_db_path + table_obj.source_db + "/latest")
    cur = con.cursor()

    query = f'''
    SELECT {', '.join(table_obj.columns.values())}
    FROM {table_obj.source_table[0]}
    WHERE "Donor" = ?
    '''
    
    extracted_metrics = []

    for case in cases:
        cur.execute(query, (case,))
        rows = cur.fetchall()
        extracted_metrics.extend(rows)

    columns = [desc[0] for desc in cur.description]
    res = pd.DataFrame(extracted_metrics, columns=columns)
    res.columns = [table_obj.headings.get(col, col) for col in columns]

    cur.close()
    con.close()

    res.columns = res.columns.str.strip('"')
    res = res.drop(columns=['SampleID'])

    # Get SampleIDs from FPR
    res = res.merge(cases_data[['Donor', 'SampleID']], on='Donor', how='left').drop_duplicates()
    column_order = ['Donor', 'SampleID'] + [col for col in res.columns if col not in ['Donor', 'SampleID']]
    res = res[column_order]

    return res

