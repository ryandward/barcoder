# Standard library imports
# (No standard library imports in this selection)

# Related third party imports
from Bio.SeqFeature import CompoundLocation
from BarCodeLibrary import BarCodeLibrary
from BowtieRunner import BowtieRunner
from CRISPRiLibrary import CRISPRiLibrary
from GenBankParser import GenBankParser
from PAMProcessor import PAMProcessor, GuideFinder, PAMFinder
from Logger import Logger
from PySamParser import PySamParser
import pyranges as pr

### Example usage ###
logger = Logger()

# barcodes = BarCodeLibrary("Example_Libraries/CN-32-zmo.tsv", column="spacer")

genbank = GenBankParser("GCA_003054575.1.gb")

finder = GuideFinder(genbank.records, "GGGGGGG", "downstream", 20)

guide_sequences = finder.find_guides_from_pam()

barcodes = BarCodeLibrary()

barcodes.load_from_list(guide_sequences)

pam = PAMFinder(genbank.records, "GGGGGGG", "downstream")

with BowtieRunner() as bowtie:
    bowtie.make_fasta(genbank.records)
    bowtie.make_fastq(barcodes.barcodes)  # Convert barcodes to FASTQ format.
    bowtie.create_index()  # Create an index for the Bowtie alignment.
    bowtie.align(num_mismatches=1, num_threads=12)
    sam = PySamParser(bowtie.sam_path)
    targets = sam.ranges.join(genbank.ranges)

guides = CRISPRiLibrary(targets.df, pam)

print(guides.unique_targets)
print(guides.unambiguous_targets)
