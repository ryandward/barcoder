from Logger import Logger

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


import os
import shutil
import subprocess
import tempfile


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
        

class BowtieError(Exception):
    """Exception raised for errors in the BowtieRunner.
    Attributes: message -- explanation of the error
    """

    def __init__(self, message):
        self.message = message
