from Logger import Logger
from Bio import SeqIO
from Bio.Seq import Seq
import csv
import os
from functools import lru_cache


class BarCode:
    def __init__(self, sequence: Seq):
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


class BarCodeLibrary(Logger):
    def __init__(self, filename=None, barcodes=None, **kwargs):
        super().__init__()
        self._barcodes = set()
        self.kwargs = kwargs
        if filename is not None:
            self.reader = BarCodeLibraryReader(filename, **self.kwargs)
            self.load()
        if barcodes is not None:
            self.load_from_list(barcodes)

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
        self.warn(
            "Genome circularity is not yet implemented. Barcodes spanning the origin will be missed!"
        )
        try:
            for sequence in self.reader.read_barcodes():
                self.add(sequence)
            self.info(
                f"Loaded {self.size} barcodes from {os.path.abspath(self.reader.filename)} ..."
            )
        except Exception as e:
            raise BarCodeLibraryError("Failed to load barcodes") from e

    def load_from_list(self, barcodes):
        for sequence in barcodes:
            self.add(sequence)
        self.info(f"Loaded {self.size} barcodes from list ...")

    @property
    @lru_cache(maxsize=None)
    def size(self):
        return len(self._barcodes)


class BarCodeLibraryError(Exception):
    """Exception raised for errors in the BarCodeLibrary.
    Attributes: message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
