# Standard library imports
from collections import defaultdict
from functools import lru_cache, wraps
from itertools import islice
from math import e
import csv
import json
import logging
import os
import re
import subprocess
import sys
import tempfile
from textwrap import indent
from tracemalloc import start

# Related third party imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord
import Bio.SeqFeature
from intervaltree import Interval, IntervalTree
import pandas as pd
import pysam
from rich.console import Console
from rich.logging import RichHandler

os.environ["TMPDIR"] = "/tmp/ramdisk"


class LoggingBase:
    def __init__(self):
        logging.basicConfig(
            level=logging.NOTSET,
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(console=Console(stderr=True))],
        )
        self.logger = logging.getLogger("rich")

    def format_numbers(self, message):
        if isinstance(message, str):
            # Find all standalone numbers in the string
            numbers = re.findall(r"\b\d+\b", message)
            # Replace each number with its formatted version
            for number in numbers:
                formatted_number = f"{int(number):,}"
                message = message.replace(number, formatted_number)
        elif isinstance(message, int):
            # Format the number with commas as thousand separators
            message = f"{message:,}"
        return message

    def info(self, message):
        message = self.format_numbers(message)
        self.logger.info(message)

    def warn(self, message):
        message = self.format_numbers(message)
        self.logger.warning(message)


class GenBankReader:
    def __init__(self, filename):
        self.filename = filename

    @property
    @lru_cache(maxsize=None)
    def records(self):
        return SeqIO.index(self.filename, "genbank")


class ReadInterval(Interval):
    @property
    def read(self):
        return self.data


class FeatureInterval(Interval):
    @property
    def feature(self):
        return self.data


# Class to process GenBank files
class GenBankParser:
    def __init__(self, filename):
        self.reader = GenBankReader(filename)
        self.records = self.reader.records

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
    def trees(self):
        trees = defaultdict(IntervalTree)

        for id, record in self.records.items():
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
                    # Create an Interval object for the feature
                    # Pack the entire feature object into the interval's data attribute
                    interval = FeatureInterval(start, end, data=feature)

                    trees[id].add(interval)

        return trees

    def write_fasta(self, filename):
        # Write the records to a FASTA file
        with open(filename, "w") as fasta_file:
            SeqIO.write(self.records.values(), fasta_file, "fasta")


class BarCode:
    def __init__(self, sequence: Seq):
        # if not isinstance(sequence, Seq):
        #     raise ValueError("sequence must be a Seq object")
        self.sequence = sequence


class BarCodeLibraryReader:
    def __init__(self, filename, column=None):
        self.filename = filename
        self.column = column

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
        if self.column is None:
            raise ValueError("A barcode column must be specified for TSV files")
        barcodes = []
        with open(self.filename, "r") as file:
            reader = csv.reader(file, delimiter="\t")
            header = next(reader)
            if self.column not in header:
                raise ValueError(f"Column '{self.column}' not found in file")
            barcode_index = header.index(self.column)
            for row in reader:
                barcodes.append(row[barcode_index])
        return barcodes


class BarCodeLibrary(LoggingBase):
    def __init__(self, filename=None, **kwargs):
        super().__init__()
        self._barcodes = set()
        self.kwargs = kwargs
        if filename is not None:
            self.reader = BarCodeLibraryReader(filename, **self.kwargs)
            self.load()

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

    def load(self):
        try:
            for sequence in self.reader.read_barcodes():
                self.add(sequence)
            self.info(
                f"Loaded {self.size} barcodes from {os.path.abspath(self.reader.filename)} ..."
            )
        except Exception as e:
            raise BarCodeLibraryError("Failed to load barcodes") from e

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
        self.info(f"Writing FASTA file {self.fasta_path} ...")
        try:
            with open(self.fasta_path, "w") as fasta_file:
                SeqIO.write(records.values(), fasta_file, "fasta")
        except Exception as e:
            raise BowtieError("Failed to write FASTA file") from e

    def write_fastq(self, barcodes):
        self.info(f"Writing FASTQ file {self.fastq_path} ...")
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

        self.info(f"Creating Bowtie index {bowtie_build_command} ...")

        try:
            subprocess.run(
                bowtie_build_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            self.info("Bowtie index created successfully.")
        except subprocess.CalledProcessError as e:
            raise BowtieError("Bowtie failed to index.") from e

    def align(self, num_mismatches=0, num_threads=os.cpu_count()):
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
            self.sam_path,
        ]
        try:
            subprocess.run(
                bowtie_align_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=True,
            )
            self.info("Bowtie alignment complete.")
        except subprocess.CalledProcessError as e:
            raise BowtieError("Bowtie failed to align.") from e


class PySamParser:
    def __init__(self, filename):
        self.filename = filename

    def read_sam(self):
        with pysam.AlignmentFile(self.filename, "r") as samfile:
            for read in samfile:
                yield read

    @property
    @lru_cache(maxsize=None)
    def trees(self):
        trees = defaultdict(IntervalTree)
        for read in self.read_sam():
            # Get the IntervalTree for this read's reference_name
            tree = trees[read.reference_name]
            # Add an interval to the tree with the read's start and end positions
            # and pack the entire read object into the interval's data attribute
            tree.add(Interval(read.reference_start, read.reference_end, data=read))
        return trees


class BowtieError(Exception):
    """Exception raised for errors in the BowtieRunner.
    Attributes: message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class BarCodeLibraryError(Exception):
    """Exception raised for errors in the BowtieRunner.
    Attributes: message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


def find_overlaps(PySamParser, GenBankParser):
    overlaps = []
    for chromosome, feature_tree in GenBankParser.trees.items():
        read_tree = PySamParser.trees.get(chromosome)
        if read_tree is None:
            continue
        for feature_interval in feature_tree:
            overlapping_intervals = read_tree[
                feature_interval.begin : feature_interval.end
            ]
            for overlap in overlapping_intervals:
                feature = (
                    feature_interval.data
                )  # This is the feature from GenBankParser
                read = overlap.data
                gene_name = (
                    feature.qualifiers.get("gene")[0]
                    if "gene" in feature.qualifiers
                    else None
                )
                locus_tag = (
                    feature.qualifiers.get("locus_tag")[0]
                    if "locus_tag" in feature.qualifiers
                    else None
                )
                overlaps.append(
                    {
                        "read": {
                            "chromosome": chromosome,
                            "start": read.reference_start,
                            "end": read.reference_end,
                            "sequence": read.query_sequence,
                            "is_reverse": read.is_reverse,
                            "mismatches": (
                                read.get_tag("NM") if read.has_tag("NM") else "0"
                            ),
                        },
                        "feature": {
                            "chromosome": chromosome,
                            "start": int(feature.location.start),
                            "end": int(feature.location.end),
                            "locus_tag": locus_tag,
                            "gene_name": gene_name,
                            "strand": feature.location.strand,
                        },
                    }
                )
    return overlaps


barcodes = BarCodeLibrary(
    "Example_Libraries/unambiguous-20-NGG-eco.tsv", column="spacer"
)
genbank = GenBankParser("GCA_000005845.2.gb")

with BowtieRunner() as bowtie:
    bowtie.write_fasta(genbank.records)

    bowtie.write_fastq(barcodes.barcodes)

    bowtie.create_index()

    bowtie.align(num_mismatches=1, num_threads=12)

    sam = PySamParser(bowtie.sam_path)

    overlaps = find_overlaps(sam, genbank)

    print(json.dumps(overlaps, indent=4))
