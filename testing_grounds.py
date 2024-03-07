# Standard library imports

# Related third party imports
from Bio.SeqFeature import CompoundLocation
import pyranges as pr
import polars as pl

from BarCodeLibrary import BarCodeLibrary
from BowtieRunner import BowtieRunner
from GenBankParser import GenBankParser
from Logger import Logger
from PySamParser import PySamParser
from GenBankParser import PAMFinder

### Example usage ###
logger = Logger()
barcodes = BarCodeLibrary("Example_Libraries/CN-32-zmo.tsv", column="spacer")
genbank = GenBankParser("GCA_003054575.1.gb")
pam = PAMFinder(genbank.records, "NGNC", "downstream")

with BowtieRunner() as bowtie:
    bowtie.make_fasta(genbank.records)
    bowtie.make_fastq(barcodes.barcodes)
    bowtie.create_index()
    bowtie.align(num_mismatches=1, num_threads=12)
    sam = PySamParser(bowtie.sam_path)
    targets = sam.ranges.join(genbank.ranges)

targets.PAM = targets.df.apply(lambda row: pam.get_pam_seq(row), axis=1)

targets.Targeting = targets.PAM.apply(lambda x: pam.pam_matches(x))

# "Source" targets are the targets that are in the genome, agnostic to genes or features

# Immediately reject Barcodes that:
# 1. Do not have a PAM
# 2. Are not targeting
# 3. Are not mapped
# 4. Are not unique (i.e. have multiple targets)

source_unique_targets = (
    targets.df[
        (targets.df["Type"] == "source")
        & (targets.df["Targeting"] == True)
        & (targets.df["Mapped"] == True)
    ]
    .loc[lambda df: ~df.duplicated(subset=["Barcode"])]
    .reset_index(drop=True)
)


# "Feature" targets are the targets that are in the genome and are annotated to a feature

feature_targets = (
    targets.df[
        (targets.df["Type"] != "source")
        & (targets.df["Targeting"] == True)
        & (targets.df["Mapped"] == True)
    ]
    .assign(
        Offset=lambda df: df.apply(
            lambda row: {"+": row.Start - row.Start_b, "-": row.End_b - row.End}.get(
                row.Strand_b, None
            ),
            axis=1,
        ),
        Overlap=lambda df: df.apply(
            lambda row: max(min(row.End, row.End_b) - max(row.Start, row.Start_b), 0),
            axis=1,
        ),
    )
    .reset_index(drop=True)
)

# Features that are not in multiple sites, but might have multiple overlapping annotations

feature_unique_targets = (
    feature_targets[feature_targets["Barcode"].isin(source_unique_targets.Barcode)]
    .sort_values(["Chromosome", "Start", "End"])
    .reset_index(drop=True)
)

# Features that are not in multiple sites and have unique annotations

feature_unambiguous_targets = feature_unique_targets[
    ~feature_unique_targets.duplicated(subset=["Barcode"]).reset_index(drop=True)
]
