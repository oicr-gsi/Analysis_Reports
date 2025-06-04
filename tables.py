from typing import List, Dict, Any
import pandas as pd
import os
from gsiqcetl import QCETLMultiCache
import gsiqcetl.column
import logging
import sqlite3
from datetime import date
from table_columns import (
    CommonColumns,
    CasesTableColumns,
    DellyTableColumns,
    PurpleTableColumns,
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
    title: str              # title of table
    blurb = ""              # descriptive blurb of table
    headings: dict          # headings for each column to be displayed on table
    columns: dict           # columns we want from sql table. Must match EXACTLY
    gsiqcetl_dirs: str      # qc-etl cache that we query from
    data = {}               # data to be displayed in table
    plots = {}              # plots to be generated for the table data
    glossary: dict          # Dict[name of column, definition]
    pipeline_step: str
    
    def get_data(self, cache, cases_data):
        cache = cache.copy()
        cache['SWID'] = cache[self.col.WorkflowRunSWID].str.extract(r'([^/]+$)')

        cases_data.rename(columns={'Workflow Run ID': 'SWID'}, inplace=True)
        swid = cases_data['SWID'].unique()
        data = cache[cache['SWID'].isin(swid)].copy()

        for col in data.select_dtypes(include='object').columns:
            if data[col].apply(lambda x: isinstance(x, list)).any():
                data.drop(columns=[col], inplace=True)
        
        merge_col = ['SWID', 'Tissue Type']
        data = data.merge(cases_data[['SWID', 'Tissue Type', 'SampleID']], on=merge_col, how='left')
        data = data[[col.strip('"') for col in self.columns.values() if col.strip('"') in data.columns]].copy()
        col = ['Donor', 'SampleID'] + [col for col in data.columns if col not in ['Donor', 'SampleID']]
        data = data[col].drop_duplicates()
        data[data.select_dtypes(include='float').columns] = data.select_dtypes(include='float').round(2)

        rename_map = {}
        for key, o_col in self.columns.items():
            o_col = o_col.strip('"')
            if o_col in data.columns and key in self.headings:
                rename_map[o_col] = self.headings[key]

        data.rename(columns=rename_map, inplace=True)

        return data
    
    def add_plot_data(self, data):
        if self.plots and not data.empty:
            plots = {}
            for col_key, plot in self.plots.items():
                if plot.y_axis in data.columns:
                    try:
                        plot_data = data[[plot.x_axis, plot.y_axis]].dropna()
                        if not plot_data.empty:
                            # Create a unique filename for the plot
                            plot_prefix = self.pipeline_step.replace('.', '_')
                            plot_filename = f"{plot_prefix}_{col_key.replace(' ', '_')}.{date.today().strftime('%Y-%m-%d')}._plot.png"
                            plot_path = plot.generate_plots(plot_data, plot_filename)
                            plots[col_key] = {
                                "fig_path": plot_path,
                                "title": plot.title
                            }

                    except Exception as e:
                        print(f"Plot generation failed for {col_key}: {e}")
                
            return plots

    def add_Seqplot_data(self, data):
        if self.plots and not data.empty:
            plots = {}
            for col_key, plot in self.plots.items():
                if plot.y_axis in data.columns:
                    try:
                        data[plot.y_axis] = pd.to_numeric(data[plot.y_axis], errors='coerce')
                        plot_data = data[[plot.x_axis, plot.y_axis, 'Sample Type']].dropna()
                        if not plot_data.empty:
                            # Create a unique filename for the plot
                            plot_prefix = self.pipeline_step.replace('.', '_')
                            plot_filename = f"{plot_prefix}_{col_key.replace(' ', '_')}.{date.today().strftime('%Y-%m-%d')}._plot.png"
                            plot_path = plot.generate_Seqplots(plot_data, plot_filename)
                            plots[col_key] = {
                                "fig_path": plot_path,
                                "title": plot.title
                            }

                    except Exception as e:
                        print(f"Plot generation failed for {col_key}: {e}")
                
            return plots
    
    def get_context(self, cases_data):
        return {
            "title": self.title,
            "blurb": self.blurb,
            "headings": self.headings,
            "columns": self.columns,
            "glossary": self.glossary,
        }


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

    def case_data(self, cases_data):
        cases_data = cases_data.drop(columns=['Workflow Run ID', 'LIMS ID', 'Sample Name']).drop_duplicates()
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
        
        return cases_data

    def load_context(self, cases_data):
        '''
        Load the context for the Cases Table. 
        '''
        context = self.get_context(cases_data)
        data = self.case_data(cases_data)

        if data.empty:
            context["data"] = []
        else:
            data = data.sort_values(by=['Donor', 'Library Type'])  
            context["data"] = data.to_dict(orient='records')

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
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/staging/qcetl_v1', '/.mounts/labs/gsi/gsiqcetl_archival/staging/ro']
        self.col = gsiqcetl.column.AnalysisDellyColumn
        self.plots = {
            DellyTableColumns.NumPASS: Plot(
                title="SV PASS Calls",
                x_axis="SampleID",
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
    
    def load_context(self, cases_data):
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        delly = load_cache(etl_caches, 'analysis_delly', 'analysis_delly',
            gsiqcetl.column.AnalysisDellyColumn.MergedPineryLimsID, True)

        data = self.get_data(delly, cases_data) 
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:  
            data = data[data['SampleID'].str.contains('_WG_')].drop_duplicates()
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_plot_data(data)

        return context

# PurpleTable class defines a table for the purple workflow
class PurpleTable(Table):
    def __init__(self):
        self.title = "Purity/Ploidy Assessment"
        self.headings = {
            PurpleTableColumns.Case: "Donor",
            PurpleTableColumns.SampleID: "SampleID",
            PurpleTableColumns.Purity: "Purity",
            PurpleTableColumns.Ploidy: "Ploidy",
            PurpleTableColumns.Pga: "PGA",
        }
        self.columns = {
            PurpleTableColumns.Case: "\"Donor\"",
            PurpleTableColumns.SampleID: "\"SampleID\"",
            PurpleTableColumns.Purity: "\"purity\"",
            PurpleTableColumns.Ploidy: "\"ploidy\"",
            PurpleTableColumns.Pga: "\"PGA\"",
        }
        self.pipeline_step = "calls.purityploidyestimates"
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/staging/qcetl_v1']
        self.col = gsiqcetl.column.AnalysisPurpleColumn
        self.plots = {
            PurpleTableColumns.Purity: Plot(
                title="Purity",
                x_axis="SampleID",
                y_axis="Purity"
            ),
            PurpleTableColumns.Ploidy: Plot(
                title="Ploidy",
                x_axis="SampleID",
                y_axis="Ploidy"
            ),
            PurpleTableColumns.Pga: Plot(
                title="Percent Genome Altered",
                x_axis="SampleID",
                y_axis="PGA"
            ),
        }
        self.glossary = {
            PurpleTableColumns.Purity: "Purity of tumor in the sample",
            PurpleTableColumns.Ploidy: "Average ploidy of the tumor sample after adjusting for purity",
            PurpleTableColumns.Pga: "Percent of the genome that is altered in the tumor sample",
        }
    
    def load_context(self, cases_data):
        print(self.gsiqcetl_dirs)
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        purple = load_cache(etl_caches, 'analysis_purple', 'analysis_purple',
            gsiqcetl.column.AnalysisPurpleColumn.MergedPineryLimsID, True)

        data = self.get_data(purple, cases_data) 
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:  
            data = data[data['SampleID'].str.contains('_WG_')].drop_duplicates()
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_plot_data(data)

        return context

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
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/staging/qcetl_v1', '/.mounts/labs/gsi/gsiqcetl_archival/staging/ro']
        self.col = gsiqcetl.column.AnalysisMutect2Column
        self.plots = {
            Mutect2TableColumns.NumPASS: Plot(
                title="Mutation Calls",
                x_axis="SampleID",
                y_axis="PASS Calls"
            ),
            Mutect2TableColumns.TITVRatio: Plot(
                title="Ti/Tv",
                x_axis="SampleID",
                y_axis="Ti/Tv Ratio",
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

    def load_context(self, cases_data):
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        mutect2 = load_cache(etl_caches, 'analysis_mutect2', 'analysis_mutect2',
            gsiqcetl.column.AnalysisMutect2Column.MergedPineryLimsID, True)

        data = self.get_data(mutect2, cases_data) 
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:  
            data = data[data['SampleID'].str.contains('_WG_')].drop_duplicates()
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_plot_data(data)

        return context

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
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/staging/qcetl_v1', '/.mounts/labs/gsi/gsiqcetl_archival/staging/ro']
        self.col = gsiqcetl.column.AnalysisRSEMColumn
        self.plots = {
            RSEMTableColumns.PctNonZero: Plot(
                title="Percent Expressed",
                x_axis="SampleID",
                y_axis="Percent Non-zero",
                hi=100,
                lo=0,
            ),
            RSEMTableColumns.Q0_5: Plot(
                title="Median TPM",
                x_axis="SampleID",
                y_axis="0.5 Quantile",
            ),
        }
        self.glossary = {
            RSEMTableColumns.Total: "Total number of reads assigned to a gene",
            RSEMTableColumns.PctNonZero: "Percentage of non-zero read counts",
            RSEMTableColumns.Q0_05: "Expression at the 0.05 quantile",
            RSEMTableColumns.Q0_5: "Expression at the 0.5 quantile",
            RSEMTableColumns.Q0_95: "Expression at the 0.95 quantile",
        }
    
    def load_context(self, cases_data):
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        rsem = load_cache(etl_caches, 'analysis_rsem', 'analysis_rsem',
            gsiqcetl.column.AnalysisRSEMColumn.MergedPineryLimsID, True)

        data = self.get_data(rsem, cases_data) 
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:  
            data = data[data['SampleID'].str.contains('_WT_')].drop_duplicates()
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_plot_data(data)

        return context

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
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/staging/qcetl_v1', '/.mounts/labs/gsi/gsiqcetl_archival/staging/ro']
        self.col = gsiqcetl.column.AnalysisStarFusionColumn
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
    
    def load_context(self, cases_data):
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        starfusion = load_cache(etl_caches, 'analysis_starfusion', 'analysis_starfusion',
            gsiqcetl.column.AnalysisStarFusionColumn.MergedPineryLimsID, True)

        data = self.get_data(starfusion, cases_data) 
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:  
            data = data[data['SampleID'].str.contains('_WT_')].drop_duplicates()
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_plot_data(data)

        return context

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
            WGCallReadyTableColumns.SampleType: "Sample Type",
        }
        self.columns = {
            WGCallReadyTableColumns.Case: "\"Donor\"",
            WGCallReadyTableColumns.SampleID: "\"SampleID\"",
            WGCallReadyTableColumns.CoverageDedup: "\"coverage deduplicated\"",
            WGCallReadyTableColumns.MarkDupPctDup: "\"mark duplicates_PERCENT_DUPLICATION\"",
            WGCallReadyTableColumns.TotalClusters: "\"total clusters\"",
            WGCallReadyTableColumns.MappedReads: "\"MappedReads\"",
            WGCallReadyTableColumns.SampleType: "\"Sample Type\"",
        }
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/production/qcetl_v1', '/.mounts/labs/gsi/gsiqcetl_archival/production/ro']
        self.pipeline_step = "alignments_WG.CallReady"
        self.plots = {
            WGCallReadyTableColumns.CoverageDedup: Plot(
                "Coverage Depth",
                "SampleID",
                "Coverage Depth",
            ),
            WGCallReadyTableColumns.MarkDupPctDup: Plot(
                "Duplication",
                "SampleID",
                "Duplication (%)",
            ),
            WGCallReadyTableColumns.TotalClusters: Plot(
                "Read Pairs",
                "SampleID",
                "Read Pairs",
            ),
            WGCallReadyTableColumns.MappedReads: Plot(
                "Mapped Reads",
                "SampleID",
                "Mapped Reads (%)",
            )
        }
        self.glossary = {
            WGCallReadyTableColumns.CoverageDedup: "Mean depth of coverage corrected for duplication",
            WGCallReadyTableColumns.MarkDupPctDup: "Percent of reads marked as duplicates",
            WGCallReadyTableColumns.TotalClusters: "Number of read pairs generated",
            WGCallReadyTableColumns.MappedReads: "Percent of reads mapping to the genomic reference",
        }
    
    def get_data(self, bamqc4merged, cases_data):
        def derive(data):
            data['MappedReads'] = (
                    (1 - data["unmapped reads meta"].astype(float) /
                    data["total input reads meta"].astype(float)) * 100
                    ).round(2)
            return data
        
        return get_seq_metrics(
            cache=bamqc4merged,
            cases_data=cases_data,
            column=self.columns,
            derived_col=derive,
            join_col='Sample Name',
            rename_col=self.headings
        )
    
    def load_context(self, cases_data):
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        bamqc4merged_columns = gsiqcetl.column.BamQc4MergedColumn
        bamqc4merged = load_cache(etl_caches, 'bamqc4merged', 'bamqc4merged',
            bamqc4merged_columns.Donor, True)

        data = self.get_data(bamqc4merged, cases_data) 
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:  
            data = data[data['SampleID'].str.contains('_WG_')]
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_Seqplot_data(data)

        return context

#WGLaneLevelTable class defines a table for the bamqc4 workflow
class WGLaneLevelTable(Table):
    def __init__(self):
        self.title = "Whole Genome Libraries, tumour and matched normal"
        self.headings = {
            WGLaneLevelTableColumns.Case: "Donor",
            WGLaneLevelTableColumns.SampleID: "SampleID",
            WGLaneLevelTableColumns.Lane: "Sequencing Run",
            WGLaneLevelTableColumns.CoverageDedup: "Coverage Depth",
            WGLaneLevelTableColumns.InsertSizeAvg: "Insert Size",
            WGLaneLevelTableColumns.MarkDupPctDup: "Duplication (%)",
            WGLaneLevelTableColumns.TotalClusters: "Read Pairs",
            WGLaneLevelTableColumns.MappedReads: "Mapped Reads (%)",
            WGLaneLevelTableColumns.SampleType: "Sample Type",
        }
        self.columns = {
            WGLaneLevelTableColumns.Case: "\"Donor\"",
            WGLaneLevelTableColumns.SampleID: "\"SampleID\"",
            WGLaneLevelTableColumns.Lane: "\"Lane\"",
            WGLaneLevelTableColumns.CoverageDedup: "\"coverage deduplicated\"",
            WGLaneLevelTableColumns.InsertSizeAvg: "\"insert size average\"",
            WGLaneLevelTableColumns.MarkDupPctDup: "\"mark duplicates_PERCENT_DUPLICATION\"",
            WGLaneLevelTableColumns.TotalClusters: "\"total clusters\"",
            WGLaneLevelTableColumns.MappedReads: "\"MappedReads\"",
            WGLaneLevelTableColumns.SampleType: "\"Sample Type\"",
        }
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/production/qcetl_v1', '/.mounts/labs/gsi/gsiqcetl_archival/production/ro']
        self.pipeline_step = "alignments_WG.lanelevel"
        self.plots = {
            WGLaneLevelTableColumns.CoverageDedup: Plot(
                title="Coverage Depth",
                x_axis="SampleID",
                y_axis="Coverage Depth",
            ),
            WGLaneLevelTableColumns.InsertSizeAvg: Plot(
                title="Insert Size",
                x_axis="SampleID",
                y_axis="Insert Size",
            ),
            WGLaneLevelTableColumns.MarkDupPctDup: Plot(
                title="Duplication",
                x_axis="SampleID",
                y_axis="Duplication (%)",
            ),
            WGLaneLevelTableColumns.TotalClusters: Plot(
                title="Read Pairs",
                x_axis="SampleID",
                y_axis="Read Pairs",
            ),
            WGLaneLevelTableColumns.MappedReads: Plot(
                title="Mapped Reads",
                x_axis="SampleID",
                y_axis="Mapped Reads (%)",
            ),
        }
        self.glossary = {
            WGLaneLevelTableColumns.CoverageDedup: "Mean depth of coverage corrected for duplication",
            WGLaneLevelTableColumns.InsertSizeAvg: "Mean size of the sequenced insert",
            WGLaneLevelTableColumns.MarkDupPctDup: "Percent of reads marked as duplicates",
            WGLaneLevelTableColumns.TotalClusters: "Number of read pairs generated",
            WGLaneLevelTableColumns.MappedReads: "Percent of reads mapping to the genomic reference",

        }
    
    def get_data(self, bamqc4, cases_data):
        def derive(data):
            data['MappedReads'] = (
                    (1 - data["unmapped reads meta"].astype(float) /
                    data["total input reads meta"].astype(float)) * 100
                    ).round(2)
            return data
        
        return get_seq_metrics(
            cache=bamqc4,
            cases_data=cases_data,
            column=self.columns,
            derived_col=derive,
            join_col='Sample Name',
            add_lane=True,
            lane_col_params={
                'run_alias': 'Run Alias',
                'lane_number': 'Lane Number'
            },
            rename_col=self.headings
        )

    def load_context(self, cases_data):
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        bamqc4_columns = gsiqcetl.column.BamQc4Column
        bamqc4 = load_cache(etl_caches, 'bamqc4', 'bamqc4', bamqc4_columns.Barcodes) 

        data = self.get_data(bamqc4, cases_data)
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:  
            data = data[data['SampleID'].str.contains('_WG_')]
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_Seqplot_data(data)

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
            WTCallReadyTableColumns.SampleType: "Sample Type",
        }
        self.columns = {
            WTCallReadyTableColumns.Case: "\"Donor\"",
            WTCallReadyTableColumns.SampleID: "\"SampleID\"",
            WTCallReadyTableColumns.PctCodingBases: "\"PCT_CODING_BASES\"",
            WTCallReadyTableColumns.TotalClusters: "\"total clusters\"",
            WTCallReadyTableColumns.MappedReads: "\"MappedReads\"",
            WTCallReadyTableColumns.RRNAContamination: "\"rrnacontaminationpercent\"",
            WTCallReadyTableColumns.SampleType: "\"Sample Type\"",
        }
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/production/qcetl_v1', '/.mounts/labs/gsi/gsiqcetl_archival/production/ro']
        self.pipeline_step = "alignments_WT.callready"
        self.plots = {
            WTCallReadyTableColumns.PctCodingBases: Plot(
                title="Percent Coding",
                x_axis="SampleID",
                y_axis="Percent Coding (%)",
            ),
            WTCallReadyTableColumns.TotalClusters: Plot(
                title="Read Pairs",
                x_axis="SampleID",
                y_axis="Read Pairs",
            ),
            WTCallReadyTableColumns.MappedReads: Plot(
                title="Mapped Reads",
                x_axis="SampleID",
                y_axis="Mapped Reads (%)",
            ),
            WTCallReadyTableColumns.RRNAContamination: Plot(
                title="rRNA Contamination",
                x_axis="SampleID",
                y_axis="rRNA Contamination (%)",
            ),
        }
        self.glossary = {
            WTCallReadyTableColumns.PctCodingBases: "Percentage of bases mapping to the coding regions of the genome",
            WTCallReadyTableColumns.TotalClusters: "Number of read pairs generated",
            WTCallReadyTableColumns.MappedReads: "Percentage of reads mapping to the genomic reference",
            WTCallReadyTableColumns.RRNAContamination: "Pecentage of reads mapping to ribosomal RNA",
        }
    def get_data(self, rnaseqqc2merged, cases_data):
        def derive(data):
            data['MappedReads'] = (
                    (1 - data["unmapped reads"].astype(float) /
                    data["total reads"].astype(float)) * 100
                    ).round(2)
            data['rrnacontaminationpercent'] = ((
                    data["rrna contamination properly paired"].astype(float) /
                    data["rrna contamination in total (QC-passed reads + QC-failed reads)"].astype(float)
                ) * 100).round(2)
            data['PCT_CODING_BASES'] = ((data["PCT_CODING_BASES"].astype(float)) * 100).round(2)
            return data
        
        return get_seq_metrics(
            cache=rnaseqqc2merged,
            cases_data=cases_data,
            column=self.columns,
            derived_col=derive,
            join_col='LIMS ID',
            rename_col=self.headings
        )

    def load_context(self, cases_data):
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        rnaseqqc2merged_columns = gsiqcetl.column.RnaSeqQc2MergedColumn
        rnaseqqc2merged = load_cache(etl_caches, 'rnaseqqc2merged', 'rnaseqqc2merged',
            rnaseqqc2merged_columns.Donor, True)

        data = self.get_data(rnaseqqc2merged, cases_data)
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:  
            data = data[data['SampleID'].str.contains('_WT_')]
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_Seqplot_data(data)

        return context

#WTLaneLevelTable class defines a table for the rnaseqqc2 workflow
class WTLaneLevelTable(Table):
    def __init__(self):
        self.title = "Whole Transcriptome Libraries, tumour only"
        self.blurb = ""
        self.headings = {
            WTLaneLevelTableColumns.Case: "Donor",
            WTLaneLevelTableColumns.SampleID: "SampleID",
            WTLaneLevelTableColumns.Lane: "Sequencing Run",
            WTLaneLevelTableColumns.PctCodingBases: "Percent Coding (%)",
            WTLaneLevelTableColumns.TotalClusters: "Read Pairs",
            WTLaneLevelTableColumns.MappedReads: "Mapped Reads (%)",
            WTLaneLevelTableColumns.RRNAContamination: "rRNA Contamination (%)",
            WTLaneLevelTableColumns.SampleType: "Sample Type",
        }
        self.columns = {
            WTLaneLevelTableColumns.Case: "\"Donor\"",
            WTLaneLevelTableColumns.SampleID: "\"SampleID\"",
            WTLaneLevelTableColumns.Lane: "\"Lane\"",
            WTLaneLevelTableColumns.PctCodingBases: "\"PCT_CODING_BASES\"",
            WTLaneLevelTableColumns.TotalClusters: "\"total clusters\"",
            WTLaneLevelTableColumns.MappedReads: "\"MappedReads\"",
            WTLaneLevelTableColumns.RRNAContamination: "\"rrnacontaminationpercent\"",
            WTLaneLevelTableColumns.SampleType: "\"Sample Type\"",
        }
        self.gsiqcetl_dirs = ['/scratch2/groups/gsi/production/qcetl_v1', '/.mounts/labs/gsi/gsiqcetl_archival/production/ro']
        self.pipeline_step = "alignments_WT.lanelevel"
        self.plots = {
            WTLaneLevelTableColumns.PctCodingBases: Plot(
                title="Percent Coding",
                x_axis="SampleID",
                y_axis="Percent Coding (%)",
            ),
            WTLaneLevelTableColumns.TotalClusters: Plot(
                title="Reads Pairs",
                x_axis="SampleID",
                y_axis="Read Pairs",
            ),
            WTLaneLevelTableColumns.MappedReads: Plot(
                title="Mapped Reads",
                x_axis="SampleID",
                y_axis="Mapped Reads (%)",
            ),
            WTLaneLevelTableColumns.RRNAContamination: Plot(
                title="rRNA Contamination",
                x_axis="SampleID",
                y_axis="rRNA Contamination (%)",
            ),
        }
        self.glossary = {
            WTLaneLevelTableColumns.PctCodingBases: "Percentage of bases mapping to the coding regions of the genome",
            WTLaneLevelTableColumns.TotalClusters: "Number of read pairs generated",
            WTLaneLevelTableColumns.MappedReads: "Percentage of reads mapping to the genomic reference",
            WTLaneLevelTableColumns.RRNAContamination: "Pecentage of reads mapping to ribosomal RNA",

        }
    
    def get_data(self, rnaseqqc2, cases_data):
        def derive(data):
            data['MappedReads'] = (
                    (1 - data["unmapped reads"].astype(float) /
                    data["total reads"].astype(float)) * 100
                    ).round(2)
            data['rrnacontaminationpercent'] = ((
                    data["rrna contamination properly paired"].astype(float) /
                    data["rrna contamination in total (QC-passed reads + QC-failed reads)"].astype(float)
                ) * 100).round(2)
            data['PCT_CODING_BASES'] = ((data["PCT_CODING_BASES"].astype(float)) * 100).round(2)
            return data
        
        return get_seq_metrics(
            cache=rnaseqqc2,
            cases_data=cases_data,
            column=self.columns,
            derived_col=derive,
            join_col='LIMS ID',
            add_lane=True,
            lane_col_params={
                'run_alias': 'Run Alias',
                'lane_number': 'Lane Number'
            },
            rename_col=self.headings
        )

    def load_context(self, cases_data):
        context = self.get_context(cases_data)
        etl_caches = QCETLMultiCache(self.gsiqcetl_dirs)
        rnaseqqc2_columns = gsiqcetl.column.RnaSeqQc2Column
        rnaseqqc2 = load_cache(etl_caches, 'rnaseqqc2', 'rnaseqqc2',
            rnaseqqc2_columns.Barcodes)

        data = self.get_data(rnaseqqc2, cases_data)
        
        if data.empty:
            context["data"] = []
            context["plots"] = {}
        else:
            data = data[data['SampleID'].str.contains('_WT_')]
            data = data.sort_values(by=['Donor', 'SampleID'])
            context["data"] = data.to_dict(orient='records')
            context["plots"] = self.add_Seqplot_data(data)

        return context

def get_seq_metrics(cache, cases_data, column, derived_col, join_col, add_lane=False, lane_col_params=None, rename_col=None):
    cases = cases_data.copy()
    if join_col == 'Sample Name':
        query = cases['Sample Name'].unique()
        data = cache[cache['library'].isin(query)].copy()

    elif join_col == 'LIMS ID':
        cases['LIMS ID'] = cases['LIMS ID'].str.split(',')
        cases = cases.explode('LIMS ID').reset_index(drop=True)
        cases['LIMS ID'] = cases['LIMS ID'].str.strip()
        lims_col = None
        if 'Merged Pinery Lims ID' in cache.columns:
            lims_col = 'Merged Pinery Lims ID'
        elif 'Pinery Lims ID' in cache.columns:
            lims_col = 'Pinery Lims ID'
        else:
            raise ValueError("No valid LIMS ID column found in cache")

        query = cases['LIMS ID'].unique()
        data = cache[cache[lims_col].isin(query)].copy()
    
    data.drop(columns=[col for col in ['Donor', 'SampleID', 'Sample Type'] if col in data.columns], inplace=True)
    
    if derived_col:
        data = derived_col(data)

    if add_lane and lane_col_params:
        data['Lane'] = data[lane_col_params['run_alias']] + "_lane_" + data[lane_col_params['lane_number']].astype(str)
    
    if join_col == 'LIMS ID':
        data = data.merge(cases[['LIMS ID', 'Donor', 'SampleID']], left_on=lims_col, right_on='LIMS ID', how='left')
    elif join_col == 'Sample Name':
        data = data.merge(cases[['Sample Name', 'Donor', 'SampleID']], left_on='library', right_on='Sample Name', how='left')
    
    data['Sample Type'] = data['SampleID'].apply(lambda x: 'Matched Normal' if '_R_' in str(x) else 'Tumor').str.strip()
    
    data = data[[col.strip('"') for col in column.values() if col.strip('"') in data.columns]].copy()
    col = ['Donor', 'SampleID'] + [col for col in data.columns if col not in ['Donor', 'SampleID']]
    data = data[col].drop_duplicates()
    
    data[data.select_dtypes(include='float').columns] = data.select_dtypes(include='float').round(2)

    if rename_col:
        rename_map = {}

        for k, v in column.items():
            stripped_v = v.strip('"')
            if stripped_v in data.columns and k in rename_col:
                rename_map[stripped_v] = rename_col[k]

        for k in rename_col:
            if k in data.columns:
                rename_map[k] = rename_col[k]

        data.rename(columns=rename_map, inplace=True)

    return data


def load_cache(etl_caches, cache_version: str, cache_name: str, id_column, merged: bool = False):
    try:
        version = etl_caches.load_same_version(cache_version).remove_missing(cache_name)
        cache = version.unique(cache_name)
        if id_column not in cache:
            logging.warning(f"'{id_column}' column not found in cache: {cache_name}.{cache_version}")
            return pd.DataFrame()
        else:
            if merged:
                single_id_column = "Pinery Lims ID"
                cache[single_id_column] = cache[id_column]
                cache = cache.explode(single_id_column)
            else:
                single_id_column = id_column
            cache.set_index(single_id_column, inplace=True, drop=False)
            cache.sort_index(inplace=True)
            return cache
    except Exception:
            logging.exception(f'Error loading cache: {cache_version}.{cache_name}')
            return pd.DataFrame()

