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
    def fastas(self):
        return {id: record.seq for id, record in self.records.items()}

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
                strand = feature.location.strand

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
                            "strand": strand,
                            "feature_type": feature.type,  # Include the feature type
                        },
                    )
                    locus_tree.add(interval)

        return locus_tree


from Bio import SeqIO
import csv


class Barcode:
    def __init__(self, sequence: Seq):
        if not isinstance(sequence, Seq):
            raise ValueError("sequence must be a Seq object")
        self.sequence = sequence


import csv
from Bio import SeqIO

import csv
from Bio import SeqIO


class BarcodeLibraryReader:
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


class BarcodeLibrary:
    def __init__(self):
        self.barcodes = set()

    def add(self, sequence):
        barcode = sequence
        self.barcodes.add(barcode)

    def remove(self, sequence):
        barcode = sequence
        self.barcodes.remove(barcode)

    def load_from_file(self, filename):
        reader = BarcodeLibraryReader(filename)
        for sequence in reader.read_barcodes():
            self.add(sequence)

    @property
    @lru_cache(maxsize=None)
    def size(self):
        return len(self.barcodes)


bl = BarcodeLibrary()
bl.load_from_file("A_E_coli.fasta")
# print(bl.size)

gbr = GenBankReader("GCA_000005845.2.gb")

gbp = GenBankProcessor(gbr.records)

# # read example_ranges.tsv as a pandas dataframe and convert to a list of tuples

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
#             print(f"Feature type: {interval.data['feature_type']}")
#             print(f"Qualifiers: {interval.data['feature'].qualifiers}")
#             print()

# Assuming gbp.fastas is a dictionary of sequences
