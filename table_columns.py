"""
A new column object is defined for each table.
These objects provide a unqiue identifier for each column of the table and 
are used to map the heading, sql column name, and order of the columns together.
These column names appear in the jinja2 context and should follow naming conventions
for json keys if possible. 
"""

# add columns that are used in multiple tables here for re-usability 
class CommonColumns:
    SampleID = "sample_id"
    TotalClusters = "total_clusters"
    MappedReads = "mapped_reads"
    Case = "case"

class CasesTableColumns:
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    ExternalID = "external_id"
    TissueType = "tissue_type"
    TissueOrigin = "tissue_origin"
    LibraryDesign = "library_design"

class DellyTableColumns:
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    NumCalls = "num_calls"
    NumPASS = "num_PASS"
    NumBND = "num_BND"
    NumDEL = "num_DEL"
    NumDUP = "num_DUP"
    NumINS = "num_INS"
    NumINV = "num_INV"

class Mutect2TableColumns:
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    NumCalls = "num_calls"
    NumPASS = "num_PASS"
    NumSNPs = "num_SNPs"
    NumIndels = "num_indels"
    TITVRatio = "titv_ratio"

class RSEMTableColumns:
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
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

class SequenzaTableColumns:
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    Index = "index"
    Cellularity = "cellularity"
    Ploidy = "ploidy"
    SLPP = "slpp"
    Gamma = "gamma"
    FGA = "fga"

class StarFusionTableColumns:
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    NumRecords = "num_records"

class WGLaneLevelTableColumns:
    SampleType = "sample_type"
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    CoverageDedup = "coverage_dedup"
    InsertSizeAvg = "insert_size_avg"
    MarkDupPctDup = "mark_dup_pct_dup"
    TotalClusters = CommonColumns.TotalClusters
    MappedReads = CommonColumns.MappedReads
    Lane = "lane"
    
class WGCallReadyTableColumns:
    SampleType = WGLaneLevelTableColumns.SampleType
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    CoverageDedup = WGLaneLevelTableColumns.CoverageDedup
    MarkDupPctDup = WGLaneLevelTableColumns.MarkDupPctDup
    TotalClusters = CommonColumns.TotalClusters
    MappedReads = WGLaneLevelTableColumns.MappedReads
    NumLimsKeys = "num_limskeys"


class WTLaneLevelTableColumns:
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    PctCodingBases = "pct_coding_bases"
    TotalClusters = CommonColumns.TotalClusters
    MappedReads = CommonColumns.MappedReads
    RRNAContamination = "rrna_contam"
    Lane = "lane"

class WTCallReadyTableColumns:
    Case = CommonColumns.Case
    SampleID = CommonColumns.SampleID
    PctCodingBases = WTLaneLevelTableColumns.PctCodingBases
    TotalClusters = WTLaneLevelTableColumns.TotalClusters
    MappedReads = WTLaneLevelTableColumns.MappedReads
    RRNAContamination = WTLaneLevelTableColumns.RRNAContamination
    NumLimsKeys = "num_limskeys"
