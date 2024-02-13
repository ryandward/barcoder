# Standard library imports
from calendar import c
from functools import lru_cache, wraps
import csv
import locale
import json
import logging
import os
import re
import shutil
import subprocess
import tempfile

# Related third party imports
from babel.numbers import format_decimal
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord
from numpy import std
import pandas as pd
import pyranges as pr
import pysam
from rich.console import Console
from rich.highlighter import JSONHighlighter
from rich.logging import RichHandler
from rich.theme import Theme


user_locale = locale.getdefaultlocale()[0]  # Get the user's locale


class Logger:

    SUBPROC = 25  # Between INFO (20) and WARNING (30)

    def __init__(self):
        console = Console(
            stderr=True, theme=Theme({"logging.level.subproc": "bold blue"})
        )

        logging.basicConfig(
            level=logging.NOTSET,
            format="%(message)s",
            datefmt="[%X]",
            handlers=[RichHandler(console=console)],
        )
        self.logger = logging.getLogger("rich")
        logging.addLevelName(self.SUBPROC, "SUBPROC")

    def format_numbers(self, message):

        if isinstance(message, str):
            lines = message.splitlines()
            for i, line in enumerate(lines):
                words = line.split()
                for j, word in enumerate(words):
                    try:
                        # Try to convert the word to a float
                        num = float(word)
                        # If successful, format the number and replace the word
                        words[j] = format_decimal(num, locale=user_locale)
                    except ValueError:
                        # If the word can't be converted to a float, ignore it
                        pass
                # Join the words back into a single string
                lines[i] = " ".join(words)
            # Join the lines back into a single string
            message = "\n".join(lines)
        elif isinstance(message, int):
            message = format_decimal(message, locale=user_locale)
        return message

    def info(self, message):
        message = self.format_numbers(message)
        self.logger.info(message)

    def warn(self, message):
        message = self.format_numbers(message)
        self.logger.warning(message)

    def subproc(self, message, *args, **kwargs):
        message = self.format_numbers(message)
        if not message:  # Check if the message is empty
            message = "No errors reported"
        if self.logger.isEnabledFor(self.SUBPROC):
            self.logger._log(self.SUBPROC, message, args, **kwargs)

    def json(self, data):
        json_str = json.dumps(data, indent=4)
        self.logger.info(json_str, extra={"highlighter": JSONHighlighter()})


class GenBankReader:
    def __init__(self, filename):
        self.filename = filename

    @property
    @lru_cache(maxsize=None)
    def records(self):
        with open(self.filename, "r") as handle:
            return SeqIO.to_dict(SeqIO.parse(handle, "genbank"))


class GenBankParser(Logger):
    def __init__(self, filename):
        super().__init__()
        self.reader = GenBankReader(filename)
        self.records = self.reader.records
        self.info(f"Found the following records:")
        self.json(self.organisms)

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
    def ranges(self):
        data = []

        for id, record in self.records.items():
            for feature in record.features:
                if feature.type not in ["source", "gene"]:
                    continue
                # Check if the feature location is a CompoundLocation
                if isinstance(feature.location, CompoundLocation):
                    # If so, use the parts of the location
                    locations = feature.location.parts
                else:
                    # If not, use the location itself
                    locations = [feature.location]

                for location in locations:
                    start = int(location.start)
                    end = int(location.end)
                    strand = location.strand

                    interval_dict = {
                        "Chromosome": id,
                        "Start": start,
                        "End": end,
                        "Strand": "+" if strand == 1 else "-" if strand == -1 else ".",
                        "Locus_Tag": feature.qualifiers.get("locus_tag", [None])[0],
                        "Gene": feature.qualifiers.get("gene", [None])[0],
                        "Type": feature.type,
                    }

                    data.append(interval_dict)

        df = pd.DataFrame(data)
        ranges = pr.PyRanges(df)
        return ranges

    def make_fasta(self, filename):
        # Write the records to a FASTA file
        with open(filename, "w") as fasta_file:
            SeqIO.write(self.records.values(), fasta_file, "fasta")

    def get_pam_sequence(self, row, pam_length, direction):
        # Fetch the sequence for the range
        sequence = self.records[row.Chromosome].seq[row.Start : row.End]

        # If the strand is "-", get the reverse complement of the sequence
        if row.Strand == "-":
            sequence = sequence.reverse_complement()

        # Get the PAM sequence
        if direction == "upstream":
            if row.Strand == "+":
                pam_sequence = self.records[row.Chromosome].seq[
                    row.Start - pam_length : row.Start
                ]
            else:
                pam_sequence = self.records[row.Chromosome].seq[
                    row.End : row.End + pam_length
                ]
        elif direction == "downstream":
            if row.Strand == "+":
                pam_sequence = self.records[row.Chromosome].seq[
                    row.End : row.End + pam_length
                ]
            else:
                pam_sequence = self.records[row.Chromosome].seq[
                    row.Start - pam_length : row.Start
                ]
        else:
            raise ValueError("direction must be 'upstream' or 'downstream'")

        # If the strand is "-", get the reverse complement of the PAM sequence
        if row.Strand == "-":
            pam_sequence = pam_sequence.reverse_complement()

        return str(pam_sequence)


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

    @property
    @lru_cache(maxsize=None)
    def size(self):
        return len(self._barcodes)


class BowtieRunner(Logger):
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
            temp_file = tempfile.NamedTemporaryFile(
                dir=self.temp_dir.name, delete=False
            )
            self._index_path = temp_file.name
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

    def make_fasta(self, records):
        # Write the records to a FASTA file
        self.info(f"Writing FASTA file {self.fasta_path} ...")
        try:
            with open(self.fasta_path, "w") as fasta_file:
                SeqIO.write(records.values(), fasta_file, "fasta")
        except Exception as e:
            raise BowtieError("Failed to write FASTA file") from e

    def make_fastq(self, barcodes):
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

        try:
            bowtie_build_path = shutil.which("bowtie-build")
        except subprocess.CalledProcessError as e:
            raise BowtieError("Bowtie not found in PATH") from e

        bowtie_build_command = ["bowtie-build", self.fasta_path, self.index_path]
        self.info(f"Creating index using {bowtie_build_path} ...")
        self.json(bowtie_build_command)

        try:
            bowtie_build_complete = subprocess.run(
                bowtie_build_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                check=True,
            )
            self.subproc(f"{bowtie_build_path} reporting:")
            self.subproc(bowtie_build_complete.stderr.decode("utf-8"))

        except subprocess.CalledProcessError as e:
            raise BowtieError("Failed to index") from e

    def align(self, num_mismatches=0, num_threads=os.cpu_count()):

        try:
            bowtie_path = shutil.which("bowtie")
        except subprocess.CalledProcessError as e:
            raise BowtieError("Bowtie not found in PATH") from e

        bowtie_align_command = [
            "bowtie",
            "-S",
            "-a",
            "--nomaqround",
            "--mm",
            "--tryhard",
            "--quiet",
            "--best",
            f"-p{str(num_threads)}",
            f"-v{str(num_mismatches)}",
            f"-x{self.index_path}",
            self.fastq_path,
            self.sam_path,
        ]

        self.info(f"Performing alignment using {bowtie_path} ...")
        self.json(bowtie_align_command)

        try:
            bowtie_complete = subprocess.run(
                bowtie_align_command,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                check=True,
            )
            self.subproc(f"{bowtie_path} reporting:")
            self.subproc(bowtie_complete.stderr.decode("utf-8"))

        except subprocess.CalledProcessError as e:
            raise BowtieError(f"{bowtie_path} failed") from e


class PySamParser:
    def __init__(self, filename):
        self.filename = filename

    def read_sam(self):
        with pysam.AlignmentFile(self.filename, "r") as samfile:
            for read in samfile:
                yield read

    @property
    @lru_cache(maxsize=None)
    def ranges(self):
        data = []

        for read in self.read_sam():

            interval_dict = {
                "Chromosome": read.reference_name,
                "Start": read.reference_start,
                "End": read.reference_end,
                "Mapped": True if not read.is_unmapped else False,
                "Strand": (
                    "+"
                    if read.is_reverse == 0
                    else "-" if read.is_reverse == 1 else "."
                ),
                "Barcode": read.query_sequence,
                "Mismatches": (read.get_tag("NM") if read.has_tag("NM") else "0"),
            }

            data.append(interval_dict)

        df = pd.DataFrame(data)
        ranges = pr.PyRanges(df)
        return ranges


class BowtieError(Exception):
    """Exception raised for errors in the BowtieRunner.
    Attributes: message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


class BarCodeLibraryError(Exception):
    """Exception raised for errors in the BarCodeLibrary.
    Attributes: message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message


### Example usage ###
logger = Logger()
barcodes = BarCodeLibrary("Example_Libraries/CN-32-zmo.tsv", column="spacer")
genbank = GenBankParser("GCA_003054575.1.gb")


with BowtieRunner() as bowtie:
    bowtie.make_fasta(genbank.records)
    bowtie.make_fastq(barcodes.barcodes)
    bowtie.create_index()
    bowtie.align(num_mismatches=0, num_threads=12)

    sam = PySamParser(bowtie.sam_path)

    # TODO: This is a work in progress. I am currently exploring how to create a class.
    pam = "NGNC"
    direction = "downstream"

    targets = sam.ranges.join(genbank.ranges)

    ## compile stats ##

    # barcodes with genome targets, here the genome is called "source"
    source_targets = targets[targets.Type == "source"]
    
    # Create a DataFrame with unique combinations of 'Barcode', 'Chromosome', 'Start', 'End'
    unique_genomic_sites = source_targets.df.drop_duplicates(subset=['Barcode', 'Chromosome', 'Start', 'End'])

    # Find barcodes that appear more than once in the unique combinations
    barcode_occurrences_in_unique_sites = unique_genomic_sites['Barcode'].value_counts()
    barcodes_in_multiple_sites = barcode_occurrences_in_unique_sites[barcode_occurrences_in_unique_sites > 1].index

    multiple_site_barcodes_df = source_targets[source_targets.Barcode.isin(barcodes_in_multiple_sites)].df

    # Add the 'Sites' column
    multiple_site_barcodes_df['Sites'] = multiple_site_barcodes_df['Barcode'].map(barcode_occurrences_in_unique_sites)

    # Convert DataFrame back to PyRanges
    multiple_site_barcodes = pr.PyRanges(multiple_site_barcodes_df)
    
    

    def convert_pam(pam):
        return pam.replace("N", "[GACT]")

    pam_search = convert_pam(pam)

    targets.PAM = targets.df.apply(
        lambda row: genbank.get_pam_sequence(row, len(pam), direction), axis=1
    )

    targets.Offset = targets.df.apply(
        lambda row: (
            row.Start - row.Start_b if row.Strand_b == "+" else row.End_b - row.End
        ),
        axis=1,
    )

    targets.Overlap = targets.df.apply(
        lambda row: min(row.End, row.End_b) - max(row.Start, row.Start_b), axis=1
    )

    targets_with_pam = targets[targets.PAM.str.contains(pam_search, regex=True)]

    if targets_with_pam.empty:
        logger.warn(f"No targets with matching PAM: '{pam}'")

    else:
        # here we're leveraging that source contains the whole genome to get the barcode occurrences
        pam_barcode_occurrences = targets_with_pam[
            targets_with_pam.Type == "source"
        ].Barcode.value_counts()

        gene_targets_with_pam = targets_with_pam[targets_with_pam.Type == "gene"]

        logger.info(f"Interpolated {len(sam.ranges)} barcodes from aligner ...")
        logger.info(
            f"Found {len(targets[targets.Type == 'source'])} barcodes in the genome ..."
        )
        logger.info(
            f"Found {len(pam_barcode_occurrences)} targets with matching PAM: '{pam}'"
        )
        # logger.info(f"Marking {len(gene_targets_with_pam)} gene targets with matching PAM: '{pam}'")
