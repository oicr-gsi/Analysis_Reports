from typing import Dict, List
import sqlite3
import json
from table_columns import (
    CommonColumns,
    RSEMTableColumns,
    SequenzaTableColumns,
    DellyTableColumns,
    Mutect2TableColumns,
    StarFusionTableColumns,
    CasesTableColumns,
    WTLaneLevelTableColumns,
    WGLaneLevelTableColumns,
    WGCallReadyTableColumns,
    WTCallReadyTableColumns,
)
from plot import (
    Plot,
    SeqPlot,
)

NUM_DP = 2 # number of decimal points

# The Table class defines each table that is generated
class Table:
    base_db_path: str       # base path to databases that are queried
    title: str              # title of table
    blurb = ""              # descriptive blurb of table
    headings: Dict[str,str] # headings for each column to be displayed on table
    columns: Dict[str, int] # columns we want from sql table. Must match EXACTLY
    data: Dict              # data from input file
    source_table: str       # table we query from
    source_db: str          # database we query from
    process: List[str]      # workflow names
    plots = {}              # Plots to generate for this table
    data = {}
    pipeline_step: str
    glossary: Dict[str,str] # Dict[name of column, definition]
    pct_stats = set()       # set of columns where the data is a percentage (i.e. 0 < data < 1) WHEN IT IS PULLED FROM database
                            # some stats are already multipled by 100 and should NOT be added

    def __init__(self, input_file, use_stage):
        env = "staging" if use_stage else "production"

        with open(input_file) as f:
            Table.data = json.load(f)
            Table.project = Table.data["project"]
            Table.release = Table.data["release"]
            Table.data = Table.data["cases"]
        Table.cases = list(Table.data.keys())
        Table.cases.sort()
        Table.base_db_path = f"/scratch2/groups/gsi/{env}/qcetl_v1/" 

    def get_select(self):
        """
        None -> str, dict[str, int]
        
        Creates an index dictionary that maps what column of the table corresponds to
        the i-th index of the select block.

        Ex. the select block "num1, num2, num3, num4" will produce an index dictionary
            {
                num1: 0,
                num2: 1,
                num3: 2, 
                num4: 3,
            }
    
        """
        indices = {}
        select_block = ""
        for index, (key, value) in enumerate(self.columns.items()):
            indices[key] = index
            select_block = select_block + f"{value}, "
        return select_block[:-2], indices

    def add_plot_data(self, plot, val, id):
        """
        (str, Any, str) -> None
        
        Checks if plot is a graph we want to generate for this table and then
        adds the data-point (val, id) to the correct plot data.

        Parameters
        ----------
        - plot (str): name of the column we want to plot
        - val  (any): y-value being plotted
        - id (str): x-value being plotted
    
        """
        if plot in self.plots.keys():
            self.plots[plot].add_data(val, id)
    
    def get_sample_id(self, cur, case, swid, pk, table_index=0):
        """
        (SQLCursor, str, str, str, int) -> str
        
        Gets the sample id for case where the pk (primary key) equals swid

        Parameters
        ----------
        - cur (SQLCursor): SQL cursor connected to the correct database
        - case (str): the case of the sample being queried
        - swid (str): the swid of the sample being queried
        - pk (str): the primary key that swid refers to
        - swid (str): the value of the primary key of the sample being queried
    
        """
        try:
            row = cur.execute(
                f"""
                select "Tissue Type", "Tissue Origin", "Library Design", "Group ID"
                from {self.source_table[table_index]}
                where "{pk}" like '%{swid}%';
                """
            ).fetchall()[0]
            return f"{case}_{row[0]}_{row[1]}_{row[2]}_{row[3]}"
        except:
            print(f"No data found for {case}, where {pk} = {swid}")
            return "nd"

    def get_row_data(self, indices, row, table_cols, entry):
        """
        (dict, tuple, ColumnObject, dict) -> dict
        
        Returns the a dict of the an entire row of values for the table

        Parameters
        ----------
        - indices (dict): dictionary that maps the column with what index it is in the SQL query
        - row (tuple): the row that is returned from the sql query
        - table_cols (ColumnObject): table-specific column object, found in table_columns.py
        - entry (dict): dictionary that maps the column with its value for the current row
    
        """
        for column in self.columns.keys():
            entry[column] = (
                #multiply by 100 if the current column is a percentage that is in decimal form
                row[indices[column]] * 100
                if column in self.pct_stats
                else row[indices[column]]
            )
            entry[column] = (
                #if entry is not a string, assume its an int/float and round to NUM_DP decimal places
                entry[column]
                if isinstance(entry[column], str)
                else round(entry[column], NUM_DP)
            )
            #add each data point to the correct graph
            self.add_plot_data(column, entry[column], entry[table_cols.SampleID])
        return entry

    def get_nd_entry(self, indices, case, source_table):
        """
        (dict, str, str) -> dict
        
        Creates a row that has nd ("no data") for each column. This is only
        run when the SQL query fails (i.e. row[0] throws an error)

        Parameters
        ----------
        - indices (dict): dictionary that maps the column with what index it
                          is in the SQL query
        - case (str): the case that is being queried
        - source_table (str): the name of the SQL table being queried    
        """
        print(f"No data found for {case} from {source_table}")
        entry = {}
        for column in indices.keys():
            entry[column] = "nd"
        return entry

    def get_data(self):
        """
        None -> list[dict]
        
        Queries self.table for the values defined in self.columns and 
        returns a list of dict with the column mapped to its value. Each 
        dict represents a row in the table and each kvp in the dict
        represent a column in the specific row
        """
        #connect to the correct database
        con = sqlite3.connect(self.base_db_path + self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for case in Table.cases:
            context  = {
                case: []
            }
            #get wfr, used as primary key to query SQL table
            wfr = None
            for limkey, run_info in Table.data[case]["analysis"][self.pipeline_step].items():
                if run_info["wf"] in self.process:
                    wfr = limkey
            if not wfr:
                raise Exception(f"No limkey found for {case} -- {self.process}")

            try:
                row = cur.execute(
                    f"""
                    select {select_block}
                    from {self.source_table[0]}
                    where "Workflow Run SWID" like '%{wfr}%';
                    """
                ).fetchall()[0]
                entry = {}
                entry[CommonColumns.Case] = case
                entry[CommonColumns.SampleID] = self.get_sample_id(cur, case, wfr, "Workflow Run SWID")
                context[case].append(self.get_row_data(indices, row, CommonColumns, entry))
            except:
                entry = self.get_nd_entry(indices, case ,self.source_table[0])
                entry[CommonColumns.Case] = case
                entry[CommonColumns.SampleID] = self.get_sample_id(cur, case, wfr, "Workflow Run SWID")
                context[case].append(entry)
            data.append(context)
        
        cur.close()
        con.close()
        return data

    def load_context(self):
        """
        None -> dict[str, Any]
        
        Returns a dict used for loading the html in jinja2 templating
        """
        context = {
            "title": self.title,
            "headings": self.headings,
            "data": self.get_data(),
            "blurb": self.blurb,
            "glossary": self.glossary,
        }
        return context

# SeqTable class defines a table for Raw Sequence data (WG and WT) and
# Call Ready Alignment Data (WG and WT)
# Main addition from Table class is that SeqTable can specify in the
# plots/tables if the sample is Matched Normal or Tumour
class SeqTable(Table):
    def add_plot_data(self, plot, sample_type, val, id):
        """
        (str, str, Any, str) -> None
        
        Checks if plot is a graph we want to generate for this table and then
        adds the data-point (val, id) to the correct plot data.

        Parameters
        ----------
        - plot (str): name of the column we want to plot
        - sample_type (str): specifies what type the sample is (Normal, Tumour)
        - val  (any): y-value being plotted
        - id (str): x-value being plotted
    
        """
        if plot in self.plots.keys():
            self.plots[plot].add_data(sample_type, val, id)

    def get_row_data(self, indices, row, table_cols, entry, sample_type):
        """
        (dict, tuple, ColumnObject, dict, str) -> dict
        
        Returns the a dict of the an entire row of values for the table

        Parameters
        ----------
        - indices (dict): dictionary that maps the column with what index it is in the SQL query
        - row (tuple): the row that is returned from the sql query
        - table_cols (ColumnObject): table-specific column object, found in table_columns.py
        - entry (dict): dictionary that maps the column with its value for the current row
        - sample_type (str): specifies what type the sample is (Normal, Tumour)

        """
        for column in self.columns.keys():
            entry[column] = (
                row[indices[column]] * 100
                if column in self.pct_stats
                else row[indices[column]]
            )
            entry[column] = (
                entry[column]
                if isinstance(entry[column], str)
                else round(entry[column],NUM_DP)
            )
            self.add_plot_data(column, sample_type, entry[column], entry[table_cols.Case])
        return entry

# CasesTable class defines the Cases table
class CasesTable(Table):
    def __init__(self):
        self.title = "Cases"
        self.blurb = ""
        self.headings = {
            CasesTableColumns.Case: "Case",
            CasesTableColumns.LibraryDesign: "Library Type",
            CasesTableColumns.TissueType : "Sample Type",
            CasesTableColumns.TissueOrigin : "Tissue Origin",
            CasesTableColumns.SampleID: "Sample ID",
            CasesTableColumns.ExternalID: "External ID",
        }
        self.columns = {}
        self.pipeline_step = ""
        self.source_table = []
        self.source_db = ""
        self.process = ""
        self.glossary = {}

    def get_data(self):
        """
        None -> list[dict]
        
        Queries self.table for the values defined in self.columns and returns a list of dict
        with the column mapped to its value. Each dict represents a row in the table and each 
        kvp in the dict represent a column in the specific row
        """
        #takes the glossary and tranforms it into a string for easier display in templating
        def _get_glossary_string(glossary):
            keys = list(glossary.keys())
            keys.sort()
            entries = [f"{k}: {glossary[k]}" for k in keys]
            return entries

        metadata_indices = {
            "proj": 0,
            "proj_num": 1,
            "tiss_origin": 2,
            "tiss_type": 3,
            "lib": 4,
        }

        # definitions are from the Configuration tab in MISO
        tissue_types = {
            'X': 'Xenograft derived from some tumour. Note: may not necessarily be a mouse xenograft',
            'U': 'Unspecified', 'T': 'Unclassifed tumour', 'S': 'Serum from blood where clotting proteins have been removed',
            'R': 'Reference or non-tumour, non-diseased tissue sample. Typically used as a donor-specific comparison to a diseased tissue, usually a cancer',
            'P': 'Primary tumour', 'O': 'Organoid', 'n': 'Unknown', 'M': 'Metastatic tumour',
            'F': 'Fibroblast cells', 'E': 'Endothelial cells', 'C': 'Cell line derived from a tumour',
            'B': 'Benign tumour', 'A': 'Cells taken from Ascites fluid'
        }
        tissue_origin = {
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
            'Gt': 'Genital', 'Hp': 'Hypopharynx', 'Hr': 'Heart', 'Ic': 'ileocecum', 'Il': 'Ileum',
            'Ki': 'Kidney', 'La': 'Lacrimal', 'Lb': 'Limb', 'Le': 'Leukocyte', 'Lg': 'Leg',
            'Li': 'Large', 'Ln': 'Lymph', 'Lp': 'Lymphoblast', 'Lu': 'Lung', 'Lv': 'Liver', 
            'Lx': 'Larynx', 'Ly': 'Lymphocyte', 'Md': 'Mediastinum', 'Me': 'Mesenchyme', 'Mn': 'Mandible',
            'Mo': 'Mouth', 'Ms': 'Mesentary', 'Mu': 'Muscle', 'Mx': 'Maxilla', 'Nk': 'Neck',
            'nn': 'Unknown', 'No': 'Nose', 'Np': 'Nasopharynx', 'Oc': 'Oral', 'Om': 'Omentum',
            'Or': 'Orbit', 'Ov': 'Ovary', 'Pa': 'Pancreas', 'Pb': 'Peripheral', 'Pc': 'Pancreatobiliary',
            'Pd': 'Parathyroid', 'Pe': 'Pelvic', 'Pg': 'Parotid', 'Ph': 'Paratracheal', 'Pi': 'Penis',
            'Pl': 'Plasma', 'Pm': 'Peritoneum', 'Pn': 'Peripheral', 'Po': 'Peri-aorta', 'Pr': 'Prostate',
            'Pt': 'Palate', 'Pu': 'Pleura', 'Py': 'periampullary', 'Ra': 'Right', 'Rc': 'Rectosigmoid',
            'Re': 'Rectum', 'Ri': 'Rib', 'Rp': 'Retroperitoneum', 'Sa': 'Saliva', 'Sb': 'Small',
            'Sc': 'Scalp', 'Se': 'Serum', 'Sg': 'Salivary', 'Si': 'Small', 'Sk': 'Skin', 'Sm': 'Skeletal',
            'Sn': 'Spine', 'So': 'Soft', 'Sp': 'Spleen', 'Sr': 'Serosa', 'Ss': 'Sinus', 'St': 'Stomach',
            'Su': 'Sternum', 'Ta': 'Tail', 'Te': 'Testes', 'Tg': 'Thymic', 'Th': 'Thymus',
            'Tn': 'Tonsil', 'To': 'Throat', 'Tr': 'Trachea', 'Tu': 'Tongue', 'Ty': 'Thyroid',
            'Uc': 'Urachus', 'Ue': 'Ureter', 'Um': 'Umbilical', 'Up': 'Urine', 'Ur': 'Urethra',
            'Us': 'Urine', 'Ut': 'Uterus', 'Uw': 'Urine', 'Vg': 'Vagina', 'Vu': 'Vulva', 'Wm': 'Worm'
        }
        library_design = {
            'WT': 'Whole Transcriptome', 'WG': 'Whole Genome', 'TS': 'Targeted Sequencing',
            'TR': 'Total RNA', 'SW': 'Shallow Whole Genome', 'SM': 'smRNA', 'SC': 'Single Cell',
            'NN': 'Unknown', 'MR': 'mRNA', 'EX': 'Exome', 'CT': 'ctDNA', 'CM': 'cfMEDIP',
            'CH': 'ChIP-Seq', 'BS': 'Bisulphite Sequencing', 'AS': 'ATAC-Seq'
            }

        data = []
        glossary = {
            CasesTableColumns.TissueType: {},
            CasesTableColumns.TissueOrigin: {},
            CasesTableColumns.LibraryDesign: {},
        }

        for case_id, case in Table.data.items():
            context = {
                case_id: []
            }
            ids = (
                list(case["WG"]["Normal"].keys())
                + list(case["WG"]["Tumour"].keys())
                + list(case["WT"]["Tumour"].keys())
            )

            for id in ids:
                metadata = id.split("_")
                context[case_id].append(
                    {
                        CasesTableColumns.Case: (
                            metadata[metadata_indices["proj"]]
                            + "_" +
                            metadata[metadata_indices["proj_num"]]
                        ),
                        CasesTableColumns.TissueType: metadata[
                            metadata_indices["tiss_type"]
                        ],
                        CasesTableColumns.TissueOrigin: metadata[
                            metadata_indices["tiss_origin"]
                        ],
                        CasesTableColumns.LibraryDesign: metadata[
                            metadata_indices["lib"]
                        ],
                        CasesTableColumns.ExternalID: case["external_id"],
                        CasesTableColumns.SampleID: id,
                    }
                )            
                glossary[CasesTableColumns.TissueType][metadata[metadata_indices["tiss_type"]]] = tissue_types[metadata[metadata_indices["tiss_type"]]]
                glossary[CasesTableColumns.TissueOrigin][metadata[metadata_indices["tiss_origin"]]] = tissue_origin[metadata[metadata_indices["tiss_origin"]]]
                glossary[CasesTableColumns.LibraryDesign][metadata[metadata_indices["lib"]]] = library_design[metadata[metadata_indices["lib"]]]
            data.append(context)

            for key, value in glossary.items():
                self.glossary[key] = ", ".join(_get_glossary_string(value)) + "."

        data = sorted(data, key=lambda d: list(d.keys()))
        return data
    
    def load_context(self):
        """
        None -> dict[str, Any]
        
        Returns a dict used for loading the html in jinja2 templating
        """
        context = {
            "title": self.title,
            "headings": self.headings,
            "data": self.get_data(),
            "blurb": self.blurb,
            "glossary": self.glossary,
        }
        return context

#DellyTable class defines a table for the delly workflow
class DellyTable(Table):
    def __init__(self):
        self.title = "Genomic Structural Variants"
        self.headings = {
            DellyTableColumns.Case: "Case",
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
        self.process = ["delly_matched_by_tumor_group", "delly"]
        self.plots = {
            DellyTableColumns.NumPASS: Plot(
                title="SV PASS Calls",
                x_axis="Sample IDs",
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

#Mutect2Table class defines a table for the mutect2 workflow   
class Mutect2Table(Table):
    def __init__(self):
        self.title = "Mutations"
        self.blurb = ""
        self.headings = {
            Mutect2TableColumns.Case: "Case",
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
        self.process = ["mutect2_matched_by_tumor_group", "mutect2"]
        self.plots = {
            Mutect2TableColumns.NumPASS: Plot(
                title="Mutation Calls",
                x_axis="Sample IDs",
                y_axis="Mutation Calls"
            ),
            Mutect2TableColumns.TITVRatio: Plot(
                title="Ti/Tv",
                x_axis="Sample IDs",
                y_axis="Ti/Tv",
                lo=0,
            ),
        }

#RSEMTable defines a table for the RSEM workflow
class RSEMTable(Table):
    def __init__(self):
        self.title = "Gene Expression"
        self.blurb = ""
        self.headings = {
            RSEMTableColumns.Case: "Case",
            RSEMTableColumns.SampleID: "Sample ID",
            RSEMTableColumns.Total: "Total Genes",
            RSEMTableColumns.PctNonZero: "Percent Expressed (%)",
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
        self.process = ["rsem"]
        self.plots = {
            "pct_non_zero": Plot(
                title="Percent Expressed",
                x_axis="Sample ID",
                y_axis="Percent Expressed (%)",
                hi=100,
                lo=0,
            ),
            "Q0.5": Plot(
                title="Median TPM",
                x_axis="Sample ID",
                y_axis="Median TPM",
            ),
        }

#SequenzaTable class defines a table for the Sequenza workflow
# SequenzaTable pulls data from two SQL databases for one table
class SequenzaTable(Table):
    def __init__(self):
        self.title = "Copy Number Alterations"
        self.blurb = ""
        self.headings = {
            SequenzaTableColumns.Case: "Case",
            SequenzaTableColumns.SampleID: "Sample ID",
            SequenzaTableColumns.Cellularity: "Cellularity",
            SequenzaTableColumns.Ploidy: "Ploidy",
            SequenzaTableColumns.FGA: "FGA (%)"
        }
        self.glossary = {
            SequenzaTableColumns.Cellularity: "Cellularity estimate (gamma=500)",
            SequenzaTableColumns.Ploidy: "Ploidy estimate (gamma = 500)",
            SequenzaTableColumns.FGA: "Fraction of the genome altered (gamma = 500)"
        }
        self.columns = {
            SequenzaTableColumns.Cellularity: "\"cellularity\"",
            SequenzaTableColumns.Ploidy: "\"ploidy\"",
        }
        self.pipeline_step = "calls.copynumber"
        self.source_table = [
            "analysis_sequenza_analysis_sequenza_alternative_solutions_1",
            "analysis_sequenza_analysis_sequenza_gamma_500_fga_1"
        ]
        self.source_db = "analysis_sequenza"
        self.process = ["sequenza_by_tumor_group", "sequenza"]
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
                y_axis="FGA (%)",
                hi=100,
                lo=0,
            ),
        }
    
    def get_data(self):
        """
        None -> list[dict]
        
        Queries self.table for the values defined in self.columns and returns a list of dict
        with the column mapped to its value. Each dict represents a row in the table and each 
        kvp in the dict represent a column in the specific row
        """
        #connect to database
        con = sqlite3.connect(self.base_db_path + self.source_db + "/latest")
        cur = con.cursor()
        select_block, indices = self.get_select()
        data = []
    
        for case in Table.cases:
            context = {
                case: []
            }
            wfr = None
            for limkey, run_info in Table.data[case]["analysis"][self.pipeline_step].items():
                if run_info["wf"] in self.process:
                    wfr = limkey
            if not wfr:
                raise Exception(f"No limkey found for {case} -- {self.process}")
            entry = {}
            entry[SequenzaTableColumns.Case] = case
            entry[SequenzaTableColumns.SampleID] = self.get_sample_id(
                cur, case, wfr, "Workflow Run SWID"
            )

            #get FGA
            try:
                row = cur.execute(
                    f"""
                    select fga from {self.source_table[1]} where "Workflow Run SWID" like '%{wfr}%';
                    """
                ).fetchall()[0]
                entry[SequenzaTableColumns.FGA] = round(row[0] * 100, NUM_DP)            
                self.add_plot_data(
                    SequenzaTableColumns.FGA,
                    entry[SequenzaTableColumns.FGA],
                    entry["sample_id"]
                )
            except:
                print(f"No data found for {case} in {self.source_table[1]}")
                entry[SequenzaTableColumns.FGA] = "nd"
            
            #get rest of column values
            try:
                row = cur.execute(
                    f"""
                    select {select_block} from {self.source_table[0]} where "Workflow Run SWID" like '%{wfr}%' and gamma = 500;
                    """
                ).fetchall()[0]
                context[case].append(
                    self.get_row_data(indices, row, SequenzaTableColumns, entry)
                )
            except:
                for column in indices.keys():
                    print(f"No data found for {case} in {self.source_table[0]}")
                    entry[column] = "nd"
                context[case].append(entry)
            data.append(context)
        
        cur.close()
        con.close()
        return data

#StarFusion class defines a table for the starfusion workflow
class StarFusionTable(Table):
    def __init__(self):
        self.title = "Gene Fusions"
        self.blurb = ""
        self.headings = {
            StarFusionTableColumns.Case: "Case",
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
        self.process = ["starfusion", "starFusion"]
        self.plots = {
            StarFusionTableColumns.NumRecords: Plot(
                title="Fusion Calls",
                x_axis="Sample ID",
                y_axis="Fusion Calls"
            )
        }

#WGCallReadyTable class defines a table for the bamqc4merged workflow
class WGCallReadyTable(SeqTable):
    def __init__(self):
        self.title = "Whole Genome Libraries, tumour and matched normal"
        self.blurb = ""
        self.headings = {
            WGCallReadyTableColumns.Case: "Case",
            WGCallReadyTableColumns.SampleID: "Sample ID",
            WGCallReadyTableColumns.SampleType: "Sample Type",
            WGCallReadyTableColumns.CoverageDedup: "Coverage Depth",
            WGCallReadyTableColumns.MarkDupPctDup: "Duplication (%)",
            WGCallReadyTableColumns.TotalClusters: "Read Pairs",
            WGCallReadyTableColumns.MappedReads: "Mapped Reads (%)",
            WGCallReadyTableColumns.NumLimsKeys: "Lanes Sequenced",
        }
        self.glossary = {
            WGCallReadyTableColumns.CoverageDedup: "Mean depth of coverage corrected for duplication",
            WGCallReadyTableColumns.MarkDupPctDup: "Percent of reads marked as duplicates",
            WGCallReadyTableColumns.TotalClusters: "Number of read pairs generated",
            WGCallReadyTableColumns.MappedReads: "Percent of reads mapping to the genomic reference",
            WGCallReadyTableColumns.NumLimsKeys: "Number of lanes of sequencing merged to call ready",
        }
        self.columns = {
            WGCallReadyTableColumns.CoverageDedup: "\"coverage deduplicated\"",
            WGCallReadyTableColumns.MarkDupPctDup: "\"mark duplicates_PERCENT_DUPLICATION\"",
            WGCallReadyTableColumns.TotalClusters: "\"total clusters\"",
            WGCallReadyTableColumns.MappedReads: """
            (1 - CAST("unmapped reads meta" as FLOAT)
            /CAST("total input reads meta" as FLOAT))
            """,
        }
        self.source_table = ["bamqc4merged_bamqc4merged_5"]
        self.source_db = "bamqc4merged"
        self.process = ["bamMergePreprocessing_by_tumor_group"]
        self.pct_stats = set(
            [
                WGCallReadyTableColumns.MarkDupPctDup,
                WGCallReadyTableColumns.MappedReads,
            ]
        )
        self.plots = {
            WGCallReadyTableColumns.CoverageDedup: SeqPlot(
                "Coverage Depth",
                "Sample ID",
                "Coverage Depth",
            ),
            WGCallReadyTableColumns.MarkDupPctDup: SeqPlot(
                "Duplication",
                "Sample ID",
                "Duplication (%)",
                hi=100,
                lo=0,
            ),
            WGCallReadyTableColumns.TotalClusters: SeqPlot(
                "Read Pairs",
                "Sample ID",
                "Read Pairs",
            ),
            WGCallReadyTableColumns.MappedReads: SeqPlot(
                "Mapped Reads",
                "Sample ID",
                "Mapped Reads (%)",
                hi=100,
                lo=0,
            )
        }
        self.sample_types = {
            "Normal": "Matched Normal",
            "Tumour": "Tumour",
        }

    def get_data(self):
        """
        None -> list[dict]
        
        Queries self.table for the values defined in self.columns and
        returns a list of dict with the column mapped to its value.
        Each dict represents a row in the table and each kvp in the 
        dict represent a column in the specific row
        """
        #connect to database
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for case in Table.cases:
            context = {
                case: []
            }
            for stype, display_type in self.sample_types.items():
                #primary key is a string of list of limkeys
                lim_keys = []
                for key in Table.data[case]["WG"][stype].keys():
                    lim_keys = lim_keys + list(Table.data[case]["WG"][stype][key].keys())
                lim_keys.sort()
                lims = "[\"" + '\", \"'.join(lim_keys) + "\"]"

                try:
                    row = cur.execute(
                        f"""
                        select {select_block} from {self.source_table[0]} where "Merged Pinery Lims ID"  like '%{lims}%';
                        """
                    ).fetchall()[0]
                    entry = {}
                    entry[WGCallReadyTableColumns.Case] = case
                    entry[WGCallReadyTableColumns.SampleType] = display_type
                    entry[WGCallReadyTableColumns.SampleID] = self.get_sample_id(
                        cur, case, lims, "Merged Pinery Lims ID"
                    )                
                    entry[WGCallReadyTableColumns.NumLimsKeys] = len(lim_keys)
                    context[case].append(
                        self.get_row_data(
                            indices,
                            row,
                            WGCallReadyTableColumns,
                            entry,
                            stype
                        )
                    )
                except:
                    entry = self.get_nd_entry(
                        indices,
                        case,
                        self.source_table[0]
                    )
                    entry[WGCallReadyTableColumns.Case] = case
                    entry[WGCallReadyTableColumns.SampleType] = display_type
                    entry[WGCallReadyTableColumns.SampleID] = self.get_sample_id(
                        cur, case, lims, "Merged Pinery Lims ID"
                    )   
                    context[case].append(entry)
            data.append(context)
        cur.close()
        con.close()
        data = sorted(data, key=lambda d: list(d.keys()))
        return data

#WGLaneLevelTable class defines a table for the bamqc4 workflow
class WGLaneLevelTable(SeqTable):
    def __init__(self):
        self.title = "Whole Genome Libraries, tumour and matched normal"
        self.blurb = ""
        self.headings = {
            WGLaneLevelTableColumns.Case: "Case",
            WGLaneLevelTableColumns.SampleID: "Sample ID",
            WGCallReadyTableColumns.SampleType: "Sample Type",
            WGLaneLevelTableColumns.Lane: "Sequencing Run",
            WGLaneLevelTableColumns.CoverageDedup: "Coverage Depth",
            WGLaneLevelTableColumns.InsertSizeAvg: "Insert Size",
            WGLaneLevelTableColumns.MarkDupPctDup: "Duplication (%)",
            WGLaneLevelTableColumns.TotalClusters: "Read Pairs",
            WGLaneLevelTableColumns.MappedReads: "Mapped Reads (%)",
        }
        self.glossary = {
            WGLaneLevelTableColumns.CoverageDedup: "Mean depth of coverage corrected for duplication",
            WGLaneLevelTableColumns.InsertSizeAvg: "Mean size of the sequenced insert",
            WGLaneLevelTableColumns.MarkDupPctDup: "Percent of reads marked as duplicates",
            WGLaneLevelTableColumns.TotalClusters: "Number of read pairs generated",
            WGLaneLevelTableColumns.MappedReads: "Percent of reads mapping to the genomic reference",

        }
        self.columns = {
            WGLaneLevelTableColumns.CoverageDedup: "\"coverage deduplicated\"",
            WGLaneLevelTableColumns.InsertSizeAvg: "\"insert size average\"",
            WGLaneLevelTableColumns.MarkDupPctDup: "\"mark duplicates_PERCENT_DUPLICATION\"",
            WGLaneLevelTableColumns.TotalClusters: "\"total clusters\"",
            WGLaneLevelTableColumns.MappedReads: "(1 - (CAST(\"unmapped reads meta\" as FLOAT) / CAST(\"total input reads meta\" as FLOAT)))",
        }
        self.pipeline_step = "alignments_WG.lanelevel"
        self.source_table = ["dnaseqqc_dnaseqqc_5", "bamqc4_bamqc4_5"]
        self.source_db = ["dnaseqqc", "bamqc4"]
        self.process = ["bwaMem"]
        self.pct_stats = set(
            [
                WGLaneLevelTableColumns.MarkDupPctDup,
                WGLaneLevelTableColumns.MappedReads,
            ]
        )
        self.plots = {
            WGLaneLevelTableColumns.CoverageDedup: SeqPlot(
                title="Coverage Depth",
                x_axis="Sample ID",
                y_axis="Coverage Depth",
            ),
            WGLaneLevelTableColumns.InsertSizeAvg: SeqPlot(
                title="Insert Size",
                x_axis="Sample ID",
                y_axis="Insert Size",
            ),
            WGLaneLevelTableColumns.MarkDupPctDup: SeqPlot(
                title="Duplication",
                x_axis="Sample ID",
                y_axis="Duplication (%)",
                hi=100,
                lo=0,
            ),
            WGLaneLevelTableColumns.TotalClusters: SeqPlot(
                title="Read Pairs",
                x_axis="Sample ID",
                y_axis="Read Pairs",
            ),
            WGLaneLevelTableColumns.MappedReads: SeqPlot(
                title="Mapped Reads",
                x_axis="Sample ID",
                y_axis="Mapped Reads (%)",
                hi=100,
                lo=0,
            ),
        }
        self.sample_types = {
            "Normal": "Matched Normal",
            "Tumour": "Tumour",
        }

    def get_sample_id(self, stype, case, swid):
        """
        (str, str, str) -> str, str
        
        Gets the sample ID and lane run data of the sample of type stype that matches
        case and swid

        Parameters
        ----------
        - stype (str): sample type (Normal or Tumour)
        - case (str): the case of the sample being queried
        - swid (str): the swid of the sample being queried
    
        """
        lims_lane_mapping = Table.data[case]["WG"][stype]
        for id, value in lims_lane_mapping.items():
            if swid in value.keys():
                return id, value[swid]["run"]
        raise Exception("There is no Sample ID associated with the limkey")

    def get_row(self, target_cur, src_table_index, select_block, lims):
        """
        (SQLCursor, int, str, str) -> List(tuple)
        
        Gets the columns in select_block from the src_table_index-th table in self.table
        that matches lims

        Parameters
        ----------
        - target_cur (SQLCursor): SQL cursor connected to the correct database
        - src_table_index (int): the index of the table being queried
        - select_block (str): the columns being selected from the SQL table
        - lims (str): the primary key of the sample
    
        """
        rows = target_cur.execute(
            f"""
            select {select_block}
            from {self.source_table[src_table_index]}
            where "Pinery Lims ID" like '%{lims}%';
            """).fetchall()
        return rows

    def get_data(self):
        """
        None -> list[dict]
        
        Queries self.table for the values defined in self.columns and returns a list of dict
        with the column mapped to its value. Each dict represents a row in the table and each 
        kvp in the dict represent a column in the specific row

        Note: get_data needs to check two tables: dnaseqqc and bamqc4. If data is not found in dnaseqqc, 
                then query bamqc4
        """
        dnaseqqc_con = sqlite3.connect(
            self.base_db_path
            + self.source_db[0]
            + "/latest"
        )
        dnaseqqc_cur = dnaseqqc_con.cursor()

        bamqc4_con = sqlite3.connect(
            self.base_db_path
            + self.source_db[1]
            + "/latest"
        )
        bamqc4_cur = bamqc4_con.cursor()

        select_block, indices = self.get_select()
        data = []

        for case in Table.cases:
            context = {
                case: []
            }
            for stype, display_type in self.sample_types.items():
                lims_keys = []
                for key in Table.data[case]["WG"][stype].keys():
                    lims_keys = lims_keys + list(
                        Table.data[case]["WG"][stype][key].keys()
                    )

                for lims in lims_keys:
                    target_cur = dnaseqqc_cur
                    src_table_index = 0

                    rows = self.get_row(
                        target_cur,
                        src_table_index,
                        select_block,
                        lims
                    )
                    if 0 == len(rows):
                        target_cur = bamqc4_cur
                        src_table_index = 1
                        rows = self.get_row(
                            target_cur,
                            src_table_index,
                            select_block,
                            lims
                        )

                    if len(rows) > 1:
                        raise Exception(
                            f"Multiple rows returned for limkeys {lims}"
                        )

                    try:
                        entry = {}
                        (
                            entry[WGLaneLevelTableColumns.SampleID],
                            entry[WGLaneLevelTableColumns.Lane]
                        ) = self.get_sample_id(stype, case, lims)
                        entry[WGLaneLevelTableColumns.Case] = case
                        entry[WGLaneLevelTableColumns.SampleType] = display_type              
                        context[case].append(
                            self.get_row_data(
                                indices,
                                rows[0],
                                WGLaneLevelTableColumns,
                                entry,
                                stype,
                            )
                        )
                    except:
                        entry = self.get_nd_entry(indices, case, self.source_table[src_table_index])
                        entry[WGLaneLevelTableColumns.Case] = case
                        entry[WGLaneLevelTableColumns.SampleType] = display_type
                        (
                            entry[WGLaneLevelTableColumns.SampleID],
                            entry[WGLaneLevelTableColumns.Lane]
                        ) = self.get_sample_id(stype, case, lims)
                        context[case].append(entry)    

            data.append(context)
        dnaseqqc_cur.close()
        dnaseqqc_con.close()
        bamqc4_cur.close()
        bamqc4_con.close()
        data = sorted(data, key=lambda d: list(d.keys()))
        return data

#WTCallReadyTable class defines a table for the rnaseqqc2merged workflow
class WTCallReadyTable(SeqTable):
    def __init__(self):
        self.title = "Whole Transcriptome Libraries, tumour only"
        self.blurb = ""
        self.headings = {
            WTCallReadyTableColumns.Case: "Case",
            WTCallReadyTableColumns.SampleID: "Sample ID",
            WTCallReadyTableColumns.PctCodingBases: "Percent Coding (%)",
            WTCallReadyTableColumns.TotalClusters: "Read Pairs",
            WTCallReadyTableColumns.MappedReads: "Mapped Reads (%)",
            WTCallReadyTableColumns.RRNAContamination: "rRNA Contamination (%)",
            WTCallReadyTableColumns.NumLimsKeys: "Lanes Sequenced",
        }
        self.glossary = {
            WTCallReadyTableColumns.PctCodingBases: "Percentage of bases mapping to the coding regions of the genome",
            WTCallReadyTableColumns.TotalClusters: "Number of read pairs generated",
            WTCallReadyTableColumns.MappedReads: "Percentage of reads mapping to the genomic reference",
            WTCallReadyTableColumns.RRNAContamination: "Pecentage of reads mapping to ribosomal RNA",
            WTCallReadyTableColumns.NumLimsKeys: "Number of lanes of sequencing merged to call ready",
        }
        self.columns = {
            WTCallReadyTableColumns.PctCodingBases: "\"PCT_CODING_BASES\"",
            WTCallReadyTableColumns.TotalClusters: "\"total clusters\"",
            WTCallReadyTableColumns.MappedReads: """
                (1 - CAST("unmapped reads" as FLOAT)
                /CAST("total reads" as FLOAT))
            """,
            WTCallReadyTableColumns.RRNAContamination: """
                (CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT))
            """,
        }
        self.pct_stats = set(
            [
                WTCallReadyTableColumns.MappedReads,
                WTCallReadyTableColumns.RRNAContamination,
            ]
        )
        self.pipeline_step = "alignments_WT.callready"
        self.source_table = ["rnaseqqc2merged_rnaseqqc2merged_2"]
        self.source_db = "rnaseqqc2merged"
        self.process = ["star_call_ready", "STAR"]
        self.plots = {
            WTCallReadyTableColumns.PctCodingBases: SeqPlot(
                title="Percent Coding",
                x_axis="Sample ID",
                y_axis="Percent Coding (%)",
                hi=100,
                lo=0,
            ),
            WTCallReadyTableColumns.TotalClusters: SeqPlot(
                title="Read Pairs",
                x_axis="Sample ID",
                y_axis="Read Pairs",
            ),
            WTCallReadyTableColumns.MappedReads: SeqPlot(
                title="Mapped Reads",
                x_axis="Sample ID",
                y_axis="Mapped Reads (%)",
                hi=100,
                lo=0,
            ),
            WTCallReadyTableColumns.RRNAContamination: SeqPlot(
                title="rRNA Contamination",
                x_axis="Sample ID",
                y_axis="rRNA Contamination (%)",
                hi=100,
                lo=0,
            ),
        }
        self.sample_types = {
            "Tumour": "Tumour",
        }

    def get_data(self):
        """
        None -> list[dict]
        
        Queries self.table for the values defined in self.columns and returns a list of dict
        with the column mapped to its value. Each dict represents a row in the table and each 
        kvp in the dict represent a column in the specific row

        Note: get_data needs to check two tables: dnaseqqc and bamqc4. If data is not found in dnaseqqc, 
                then query bamqc4
        """
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for case in Table.cases:
            context = {
                case: []
            }
            for stype in self.sample_types.keys():
                #primary key is a string of list of limkeys
                limkeys = list(
                    Table.data[case]["analysis"][self.pipeline_step].values()
                )[0]["limkeys"].split(':')
                limkeys.sort()
                lims = "[\"" + "\", \"".join(limkeys) + "\"]"

                try:
                    row = cur.execute(
                        f"""
                        select {select_block} from {self.source_table[0]} where "Merged Pinery Lims ID"  = '{lims}';
                        """
                    ).fetchall()[0]
                    entry = {}
                    entry[WTCallReadyTableColumns.Case] = case
                    entry[WTCallReadyTableColumns.SampleID] = self.get_sample_id(
                        cur, case, lims, "Merged Pinery Lims ID"
                    )
                    entry[WTCallReadyTableColumns.NumLimsKeys] = len(limkeys)
                    context[case].append(self.get_row_data(indices, row, WTCallReadyTableColumns, entry, stype))
                except:
                    entry = self.get_nd_entry(indices, case, self.source_table[0])
                    entry[WTCallReadyTableColumns.Case] = case
                    entry[WTCallReadyTableColumns.SampleID] = self.get_sample_id(
                        cur, case, lims, "Merged Pinery Lims ID"
                    )
                    entry[WTCallReadyTableColumns.NumLimsKeys] = len(limkeys) 
                    context[case].append(entry)
            data.append(context)

        cur.close()
        con.close()
        return data

#WTLaneLevelTable class defines a table for the rnaseqqc2 workflow
class WTLaneLevelTable(SeqTable):
    def __init__(self):
        self.title = "Whole Transcriptome Libraries, tumour only"
        self.blurb = ""
        self.headings = {
            WTLaneLevelTableColumns.Case: "Case",
            WTLaneLevelTableColumns.SampleID: "Sample ID",
            WTLaneLevelTableColumns.Lane: "Sequencing Run",
            WTLaneLevelTableColumns.PctCodingBases: "Percent Coding (%)",
            WTLaneLevelTableColumns.TotalClusters: "Read Pairs",
            WTLaneLevelTableColumns.MappedReads: "Mapped Reads (%)",
            WTLaneLevelTableColumns.RRNAContamination: "rRNA Contamination (%)",

        }
        self.glossary = {
            WTLaneLevelTableColumns.PctCodingBases: "Percentage of bases mapping to the coding regions of the genome",
            WTLaneLevelTableColumns.TotalClusters: "Number of read pairs generated",
            WTLaneLevelTableColumns.MappedReads: "Percentage of reads mapping to the genomic reference",
            WTLaneLevelTableColumns.RRNAContamination: "Pecentage of reads mapping to ribosomal RNA",

        }
        self.columns = {
            WTLaneLevelTableColumns.PctCodingBases: "\"PCT_CODING_BASES\"",
            WTLaneLevelTableColumns.TotalClusters: "\"total clusters\"",
            WTLaneLevelTableColumns.MappedReads: """
                (1 - CAST(\"unmapped reads\" as FLOAT)
                /CAST(\"total reads\" as FLOAT))
            """,
            WTLaneLevelTableColumns.RRNAContamination: """
                (CAST("rrna contamination properly paired" as FLOAT)
                /CAST("rrna contamination in total (QC-passed reads + QC-failed reads)" as FLOAT))
            """,
        }
        self.source_table = ["rnaseqqc2_rnaseqqc2_2"]
        self.source_db = "rnaseqqc2"
        self.process = ["star_lane_level", "STAR"]
        self.pipeline_step = "alignments_WT.lanelevel"
        self.pct_stats = set(
            [
                WTLaneLevelTableColumns.MappedReads,
                WTLaneLevelTableColumns.RRNAContamination,
            ]
        )
        self.plots = {
            WTLaneLevelTableColumns.PctCodingBases: SeqPlot(
                title="Percent Coding",
                x_axis="Sample ID",
                y_axis="Percent Coding (%)",
                hi=100,
                lo=0,
            ),
            WTLaneLevelTableColumns.TotalClusters: SeqPlot(
                title="Reads Pair",
                x_axis="Sample ID",
                y_axis="Reads Pair",
            ),
            WTLaneLevelTableColumns.MappedReads: SeqPlot(
                title="Mapped Reads",
                x_axis="Sample ID",
                y_axis="Mapped Reads (%)",
                hi=100,
                lo=0,
            ),
            WTLaneLevelTableColumns.RRNAContamination: SeqPlot(
                title="rRNA Contamination",
                x_axis="Sample ID",
                y_axis="rRNA Contamination (%)",
                hi=100,
                lo=0,
            ),
        }
        self.sample_types = {
            "Tumour": "Tumour",
        }
    
    def get_sample_id(self, case, swid):
        """
        (str, str) -> str, str
        
        Gets the sample ID and lane run data of the sample of type stype that matches
        case and swid

        Parameters
        ----------
        - case (str): the case of the sample being queried
        - swid (str): the swid of the sample being queried
    
        """
        lims_lane_mapping = Table.data[case]["WT"]["Tumour"]
        for id, value in lims_lane_mapping.items():
            if swid in value.keys():
                return id, value[swid]["run"]
        raise Exception("There is no Sample ID associated with the limkey")

    def get_data(self):
        """
        None -> list[dict]
        
        Queries self.table for the values defined in self.columns and
        returns a list of dict with the column mapped to its value. Each
        dict represents a row in the table and each kvp in the dict
        represent a column in the specific row

        Note: get_data needs to check two tables: dnaseqqc and bamqc4.
        If data is not found in dnaseqqc, then query bamqc4
        """
        con = sqlite3.connect(self.base_db_path+ self.source_db + "/latest")
        cur = con.cursor()

        select_block, indices = self.get_select()
        data = []

        for case in Table.cases:
            context = {
                case: []
            }
            for stype in self.sample_types.keys():
                lims_keys = []
                for key in Table.data[case]["WT"]["Tumour"].keys():
                    lims_keys = lims_keys + list(
                        Table.data[case]["WT"]["Tumour"][key].keys()
                    )

                for lims in lims_keys:
                    rows = cur.execute(
                        f"""
                        select {select_block} from {self.source_table[0]} where "Pinery Lims ID" = '{lims}';
                        """).fetchall()
                    if len(rows) > 1:
                        raise Exception(f"Multiple rows returned for limkeys {lims}")
                    
                    try:
                        entry = {}
                        entry[WTLaneLevelTableColumns.Case] = case
                        (entry[WTLaneLevelTableColumns.SampleID],
                        entry[WTLaneLevelTableColumns.Lane]
                        ) = self.get_sample_id(case, lims)
                        context[case].append(
                            self.get_row_data(
                            indices, rows[0], WTLaneLevelTableColumns, entry, stype
                            )
                        )
                    except:
                        entry = self.get_nd_entry(indices, case, self.source_table[0])
                        entry[WTLaneLevelTableColumns.Case] = case
                        (
                            entry[WTLaneLevelTableColumns.SampleID],
                            entry[WTLaneLevelTableColumns.Lane]
                        ) = self.get_sample_id(case, lims)
                        context[case].append(entry)
            data.append(context)
        cur.close()
        con.close()
        return data
