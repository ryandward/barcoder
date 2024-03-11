# Standard library imports
# (No standard library imports in this selection)

# Related third party imports
from Bio.SeqFeature import CompoundLocation
from BarCodeLibrary import BarCodeLibrary
from BowtieRunner import BowtieRunner
from CRISPRiLibrary import CRISPRiLibrary
from GenBankParser import GenBankParser
from PAMProcessor import PAMProcessor, GuideFinderByPAM, PAMRetriever
from Logger import Logger
from PySamParser import PySamParser
import pyranges as pr

### Example usage ###
# Initialize a Logger instance to handle logging throughout the script.
logger = Logger()

# Create a BarCodeLibrary instance by reading a TSV file containing barcodes.
# The 'column' parameter specifies which column in the TSV to use for barcodes.
# barcodes = BarCodeLibrary("Example_Libraries/CN-32-zmo.tsv", column="spacer")

# Parse a GenBank file to extract genetic information using GenBankParser.
genbank = GenBankParser("GCA_003054575.1.gb")

finder = GuideFinderByPAM(genbank.records, "GGGGGG", "downstream", 32)

pam_sequences = finder.find_pam_sequences()

print(pam_sequences[1:20])

# # Find Protospacer Adjacent Motifs (PAMs) using the PAMFinder class.
# # It searches for the specified PAM sequence ('NGNC') downstream of the genetic records.
# pam = PAMRetriever(genbank.records, "NGNC", "downstream")

# # Use BowtieRunner as a context manager to handle the lifecycle of the Bowtie alignment tool.
# # This block will create necessary FASTA and FASTQ files, build an index, and perform alignment.
# with BowtieRunner() as bowtie:
#     # Convert GenBank records to FASTA format.
#     bowtie.make_fasta(genbank.records)
#     bowtie.make_fastq(barcodes.barcodes)  # Convert barcodes to FASTQ format.
#     bowtie.create_index()  # Create an index for the Bowtie alignment.
#     # Perform the alignment with specified parameters.
#     bowtie.align(num_mismatches=1, num_threads=12)
#     # Parse the resulting SAM file to extract alignment data.
#     sam = PySamParser(bowtie.sam_path)
#     # Join the alignment data with GenBank ranges.
#     targets = sam.ranges.join(genbank.ranges)

# # Create a CRISPRiLibrary instance to manage CRISPR interference (CRISPRi) guides.
# guides = CRISPRiLibrary(targets.df, pam)

# # Print the unique targets for the CRISPRi library.
# print(guides.unique_targets)

# # Print unambiguous targets for the CRISPRi library.
# print(guides.unambiguous_targets)
