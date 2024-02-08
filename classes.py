from math import e
import select
from tracemalloc import start
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
import Bio.SeqFeature
from functools import lru_cache
import json
from itertools import islice
import pandas as pd
from rich.console import Console
from rich.logging import RichHandler
import logging
from functools import wraps
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import tempfile
import os
import pysam
import subprocess

# Configure logging
logging.basicConfig(
    level="NOTSET", format="%(message)s", datefmt="[%X]", handlers=[RichHandler()]
)


# Create a base class for logging
class LoggingBase:
    def __init__(self):
        self.console = Console()

    @staticmethod
    def logged_property(func):
        @property
        @wraps(func)
        def wrapper(self):
            self.console.log(
                f"{self.__class__.__name__} is creating '{func.__name__}'..."
            )
            return func(self)

        return wrapper


class GenBankReader:
    def __init__(self, filename):
        self.filename = filename

    @property
    @lru_cache(maxsize=None)
    def records(self):
        return SeqIO.index(self.filename, "genbank")


# Class to process GenBank files
class GenBankProcessor:
    def __init__(self, records):
        self.records = records
        self.console = Console(stderr=True, style="bold blue")

    @property
    @lru_cache(maxsize=None)
    def organisms(self):
        return {
            id: record.annotations.get("organism", None)
            for id, record in self.records.items()
        }

    @property
    @lru_cache(maxsize=None)
    def seq_lens(self):
        return {id: len(record.seq) for id, record in self.records.items()}

    @property
    @lru_cache(maxsize=None)
    def topologies(self):
        return {
            id: record.annotations.get("topology", None)
            for id, record in self.records.items()
        }

    @property
    @lru_cache(maxsize=None)
    def num_genes(self):
        return {
            id: len([f for f in record.features if f.type == "gene"])
            for id, record in self.records.items()
        }

    @property
    @lru_cache(maxsize=None)
    def overhangs(self):
        return {
            id: (100_000 if self.topologies[id] == "circular" else 0)
            for id, _ in self.records.items()
        }

    @property
    @lru_cache(maxsize=None)
    def locus_tree(self):

        locus_tree = IntervalTree()

        for (
            id,
            record,
        ) in self.records.items():
            for feature in record.features:
                # Check if the feature location is a CompoundLocation
                if isinstance(feature.location, CompoundLocation):
                    # If so, use the parts of the location
                    locations = feature.location.parts
                else:
                    # If not, use the location itself
                    locations = [feature.location]

                # Iterate over the locations
                for location in locations:
                    start = int(location.start)
                    end = int(location.end)
                    interval = Interval(
                        start,
                        end,
                        {
                            "id": id,
                            "feature": feature,  # Store the entire SeqFeature object
                            "sequence": record.seq[start:end],
                            # "strand": feature.location.strand,
                            # "feature_type": feature.type,  # Include the feature type
                        },
                    )
                    locus_tree.add(interval)

        return locus_tree

    def write_fasta(self, filename):
        # Write the records to a FASTA file
        with open(filename, "w") as fasta_file:
            SeqIO.write(self.records.values(), fasta_file, "fasta")


from Bio import SeqIO
import csv


class BarCode:
    def __init__(self, sequence: Seq):
        if not isinstance(sequence, Seq):
            raise ValueError("sequence must be a Seq object")
        self.sequence = sequence


import csv
from Bio import SeqIO

import csv
from Bio import SeqIO


class BarCodeLibraryReader:
    def __init__(self, filename, barcode_column=None):
        self.filename = filename
        self.barcode_column = barcode_column

    def read_barcodes(self):
        if self.filename.endswith(".fasta"):
            return self._read_fasta()
        elif self.filename.endswith(".tsv"):
            return self._read_tsv()
        else:
            raise ValueError(f"Unsupported file format: {self.filename}")

    def _read_fasta(self):
        barcodes = []
        for record in SeqIO.parse(self.filename, "fasta"):
            barcodes.append(str(record.seq))
        return barcodes

    def _read_tsv(self):
        if self.barcode_column is None:
            raise ValueError("A barcode column must be specified for TSV files")
        barcodes = []
        with open(self.filename, "r") as file:
            reader = csv.reader(file, delimiter="\t")
            header = next(reader)
            if self.barcode_column not in header:
                raise ValueError(f"Column '{self.barcode_column}' not found in file")
            barcode_index = header.index(self.barcode_column)
            for row in reader:
                barcodes.append(row[barcode_index])
        return barcodes


class BarCodeLibrary:
    def __init__(self, filename=None):
        self._barcodes = set()
        if filename is not None:
            self.load(filename)

    @property
    @lru_cache(maxsize=None)
    def barcodes(self):
        return self._barcodes

    def add(self, sequence):
        barcode = sequence
        self._barcodes.add(barcode)

    def remove(self, sequence):
        barcode = sequence
        self._barcodes.remove(barcode)

    def load(self, filename):
        reader = BarCodeLibraryReader(filename)
        for sequence in reader.read_barcodes():
            self.add(sequence)

    @property
    @lru_cache(maxsize=None)
    def size(self):
        return len(self._barcodes)


class BowtieRunner(LoggingBase):
    def __init__(self):
        super().__init__()
        self.temp_dir = tempfile.TemporaryDirectory()
        self._fasta_path = None
        self._index_path = None
        self._fastq_path = None
        self._sam_path = None

    @property
    def index_path(self):
        if self._index_path is None:
            self._index_path = tempfile.mktemp(dir=self.temp_dir.name)
        return self._index_path

    @property
    def fasta_path(self):
        if self._fasta_path is None:
            self._fasta_path = self.index_path + ".fasta"
        return self._fasta_path

    @property
    def fastq_path(self):
        if self._fastq_path is None:
            self._fastq_path = self.index_path + ".fastq"
        return self._fastq_path
    
    @property
    def sam_path(self):
        if self._sam_path is None:
            self._sam_path = self.index_path + ".sam"
        return self._sam_path    

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.temp_dir.cleanup()

    def write_fasta(self, records):
        # Write the records to a FASTA file
        self.console.log(f"Writing FASTA file {self.fasta_path}...")
        try:
            with open(self.fasta_path, "w") as fasta_file:
                SeqIO.write(records.values(), fasta_file, "fasta")
        except Exception as e:
            raise BowtieError("Failed to write FASTA file") from e

    def write_fastq(self, barcodes):
        self.console.log(f"Writing FASTQ file {self.fastq_path}...")
        try:
            with open(self.fastq_path, "a") as output_handle:
                for barcode in barcodes:
                    # Create a SeqRecord for the barcode
                    record = SeqRecord(Seq(barcode))
                    # Assign a quality score of 40 to each base in the barcode
                    record.letter_annotations["phred_quality"] = [40] * len(barcode)
                    # Write the SeqRecord to the FASTQ file
                    SeqIO.write([record], output_handle, "fastq")
        except Exception as e:
            raise BowtieError("Failed to write FASTQ file") from e

    def create_index(self):
        if not os.path.exists(self.fasta_path) or os.path.getsize(self.fasta_path) == 0:
            raise BowtieError("BowtieRunner.fasta_path does not exist or is an empty")

        bowtie_build_command = ["bowtie-build", self.fasta_path, self.index_path]

        self.console.log(f"Creating Bowtie index {bowtie_build_command}...")

        try:
            subprocess.run(
                bowtie_build_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            self.console.log("Bowtie index created successfully.")
        except subprocess.CalledProcessError as e:
            raise BowtieError("Bowtie failed to index.") from e

    def align(self, num_mismatches=0, num_threads=1):
        bowtie_align_command = [
            "bowtie",
            "-S",
            "-k 100",
            "--nomaqround",
            "-p",
            str(num_threads),
            "--tryhard",
            "-v",
            str(num_mismatches),
            "-x",
            self.index_path,
            self.fastq_path,
            "-S",
            self.sam_path
        ]
        try:
            subprocess.run(
                bowtie_align_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            self.console.log("Bowtie alignment complete.")
        except subprocess.CalledProcessError as e:
            raise BowtieError("Bowtie failed to align.") from e


class BowtieError(Exception):
    """Exception raised for errors in the BowtieRunner.
    Attributes: message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


bcl = BarCodeLibrary("A_E_coli.fasta")
gbr = GenBankReader("GCA_000005845.2.gb")
gbp = GenBankProcessor(gbr.records)

with BowtieRunner() as bowtie:
    bowtie.write_fasta(gbr.records)
    bowtie.write_fastq(bcl.barcodes)
    bowtie.create_index()
    bowtie.align(3, 12)


# def run_bowtie_and_parse(
#     sgrna_fastq_file_name,
#     topological_fasta_file_name,
#     locus_map,
#     num_mismatches,
#     num_threads,
# ):
#     results = []
#     # Create a temporary directory for the bowtie index files
#     with tempfile.TemporaryDirectory() as temp_dir:
#         # Create a temporary name for the genome index
#         genome_index_temp_name = tempfile.NamedTemporaryFile(
#             dir=temp_dir, delete=False
#         ).name
#         index_prefix = os.path.join(temp_dir, genome_index_temp_name)

#         with open(os.devnull, "w") as devnull:

#             bowtie_build_command = [
#                 "bowtie-build",
#                 topological_fasta_file_name,
#                 index_prefix,
#             ]

#             bowtie_build_process = subprocess.Popen(
#                 bowtie_build_command,
#                 stdout=devnull,
#                 stderr=devnull
#             )
#             bowtie_build_process.wait()

#             # Create a temporary file for the bowtie output
#             with tempfile.NamedTemporaryFile(delete=False) as bowtie_output_temp_file:

#                 bowtie_command = [
#                     "bowtie",
#                     "-S",
#                     "-k 100",
#                     "--nomaqround",
#                     "-p",
#                     str(num_threads),
#                     "--tryhard",
#                     "-v",
#                     str(num_mismatches),
#                     "-x",
#                     index_prefix,
#                     sgrna_fastq_file_name,
#                 ]

#                 bowtie_process = subprocess.Popen(
#                     bowtie_command,
#                     stdout=subprocess.PIPE,
#                     stderr=devnull,
#                     universal_newlines=True
#                 )

#         if bowtie_process.stdout is None:
#             raise RuntimeError("Bowtie was unable to start. Check your installation.")

#         if bowtie_process.stdout is not None:
#                 with pysam.AlignmentFile(bowtie_process.stdout, "r") as samfile:
#                     results = parse_sam_output(
#                         samfile,
#                         locus_map,
#                         topological_fasta_file_name,
#                         args.genome_file,
#                         args.pam,
#                         args.pam_direction,
#                     )

#                 bowtie_process.wait()

#                 # Delete the temporary file after we're done with it
#                 os.remove(bowtie_output_temp_file.name)

#     # Return the results of parse_sam_output
#     if results is []:
#         raise RuntimeError("No results were returned from the Bowtie process. Check your input files and parameters.")
#     else:
#         return results


# read example_ranges.tsv as a pandas dataframe and convert to a list of tuples

# ranges = pd.read_csv("example_ranges.tsv", sep="\t")

# print(ranges)

# ranges = [tuple(x) for x in ranges.values]

# for start, end in ranges:
#     overlapping_intervals = gbp.locus_tree.overlap(start, end)
#     for interval in overlapping_intervals:
#         feature = interval.data["feature"]
#         if isinstance(feature.location, Bio.SeqFeature.CompoundLocation):
#             print(f"Query range: {start}-{end}")
#             print(f"Overlapping interval: {interval}")
#             print(f"Feature type: {interval.data['feature'].type}")
#             print(f"Qualifiers: {interval.data['feature'].qualifiers}")
#             print()

# Assuming gbp.fastas is a dictionary of sequences

# Run bowtie and parse the output using parse_sam_output and temp files
