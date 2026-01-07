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
    SampleID = "sample"    

class CasesTableColumns:
    Case = CommonColumns.Case
    ExternalID = "External ID"
    SampleID = "SampleID"
    LibraryType = "Library Type"
    TissueOrigin = "Tissue Origin"
    TissueType = "Tissue Type"
    GroupID = "Group ID"
    TissuePreparation = "Tissue Preparation"
    

class WGLaneLevelTableColumns:
    SampleID = CommonColumns.SampleID
    Lane = "Lane"
    CoverageDedup = "coverage deduplicated"
    InsertSizeAvg = "insert size average"
    MarkDupPctDup = "mark duplicates_PERCENT_DUPLICATION"
    TotalClusters = "total clusters"
    MappedReads = "mapped reads"
    SampleType = "Sample Type"

class WTLaneLevelTableColumns:
    SampleID = CommonColumns.SampleID
    Lane = "Lane"
    PctCodingBases = "PCT_CODING_BASES"
    TotalClusters = "total clusters"
    MappedReads = "mapped reads"
    RRNAContamination = "rrna_contam"
    SampleType = "Sample Type"

class WGCallReadyTableColumns:
    SampleID = "SampleID"
    CoverageDedup = "coverage deduplicated"
    MarkDupPctDup = "mark duplicates_PERCENT_DUPLICATION"
    TotalClusters = "total clusters"
    MappedReads = "mapped reads"
    SampleType = "Sample Type"

class WTCallReadyTableColumns:
    SampleID = "SampleID"
    PctCodingBases = "PCT_CODING_BASES"
    TotalClusters = "total clusters"
    MappedReads = "mapped reads"
    RRNAContamination = "rrna_contam"

class Mutect2TableColumns:
    SampleID = "SampleID"
    NumCalls = "num_calls"
    NumPASS = "num_PASS"
    NumSNPs = "num_SNPs"
    NumIndels = "num_indels"
    TITVRatio = "titv_ratio"

class DellyTableColumns:
    SampleID = "SampleID"
    NumCalls = "num_calls"
    NumPASS = "num_PASS"
    NumBND = "num_BND"
    NumDEL = "num_DEL"
    NumDUP = "num_DUP"
    NumINS = "num_INS"
    NumINV = "num_INV"

class PurpleTableColumns:
    SampleID = "SampleID"
    Purity = "purity"
    Ploidy = "ploidy"
    Pga = "PGA"

class MrdTableColumns:
    SampleID = "SampleID"
    SitesDetected = "sites_detected"
    CancerDetected = "cancer_detected"
    CandidateSNPs = "sample_candidate_SNPs"
    SampleCoverage = "sample_coverage"

class RSEMTableColumns:
    SampleID = "SampleID"
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
    SampleID = "SampleID"
    NumRecords = "num_records"