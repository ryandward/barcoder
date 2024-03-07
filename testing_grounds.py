# Standard library imports

# Related third party imports
from Bio.SeqFeature import CompoundLocation
import pyranges as pr
import polars as pl

from BarCodeLibrary import BarCodeLibrary
from BowtieRunner import BowtieRunner
from CRISPRiLibrary import CRISPRiLibrary
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

crispr_guides = CRISPRiLibrary(targets.df, pam)
source_unique_targets = crispr_guides.get_source_unique_targets()
feature_targets = crispr_guides.get_feature_targets()
feature_unique_targets = crispr_guides.get_feature_unique_targets(source_unique_targets)
feature_unambiguous_targets = crispr_guides.get_feature_unambiguous_targets(
    feature_unique_targets
)
