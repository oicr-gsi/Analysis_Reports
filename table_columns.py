'''
A new column object is defined for each table.
These objects provide a unqiue identifier for each column of the table and 
are used to map the heading, sql column name, and order of the columns together.
These column names appear in the jinja2 context and should follow naming conventions
for json keys if possible. 
'''

# add columns that are used in multiple tables here for re-usability 
class CommonColumns:
    Case = "Donor"

class CasesTableColumns:
    Case = CommonColumns.Case
    GroupID = "Group ID"
    LibraryType = "Library Type"
    TissueType = "Tissue Type"
    TissueOrigin = "Tissue Origin"
    TissuePreparation = "Tissue Preparation"
    ExternalID = "External ID"
    SampleID = "Sample ID"

class DellyTableColumns:
    Case = CommonColumns.Case
    NumCalls = "num_calls"
    NumPASS = "num_PASS"
    NumBND = "num_BND"
    NumDEL = "num_DEL"
    NumDUP = "num_DUP"
    NumINS = "num_INS"
    NumINV = "num_INV"

class Mutect2TableColumns:
    Case = CommonColumns.Case
    NumCalls = "num_calls"
    NumPASS = "num_PASS"
    NumSNPs = "num_SNPs"
    NumIndels = "num_indels"
    TITVRatio = "titv_ratio"

class RSEMTableColumns:
    Case = CommonColumns.Case
    Total = "total"
    PctNonZero = "pct_non_zero"
    Q0 = "Q0"
    Q0_05 = "Q0.05"
    Q0_1 = "Q0.1"
    Q0_25 = "Q0.25"
    Q0_5 = "Q0.5"
    Q0_75 = "Q0.75"
    Q0_9 = "Q0.9"
    Q0_95 = "Q0.95"
    Q1 = "Q1"

class StarFusionTableColumns:
    Case = CommonColumns.Case
    NumRecords = "num_records"