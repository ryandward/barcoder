import re
from Bio import SeqIO
import pyranges as pr
import pandas as pd
from Logger import Logger
from Bio.SeqFeature import CompoundLocation
from functools import lru_cache


class GenBankReader:
    def __init__(self, filename):
        self.filename = filename

    @property
    @lru_cache(maxsize=None)
    def records(self):
        with open(self.filename, "r") as handle:
            return SeqIO.to_dict(SeqIO.parse(handle, "genbank"))


class PAMFinder:
    def __init__(self, records, pam, direction):
        self.records = records
        self.pam = pam
        self.pam_length = len(pam)
        self.pam_pattern = self.pam.replace("N", "[ATCG]")

        self.direction = direction

    def get_pam_seq(self, row):
        # Fetch the sequence for the range
        sequence = self.records[row.Chromosome].seq[row.Start : row.End]

        # If the strand is "-", get the reverse complement of the sequence
        if row.Strand == "-":
            sequence = sequence.reverse_complement()

        # Get the PAM sequence
        if self.direction == "upstream":

            if row.Strand == "+":
                pam_sequence = self.records[row.Chromosome].seq[
                    row.Start - self.pam_length : row.Start
                ]
            else:
                pam_sequence = self.records[row.Chromosome].seq[
                    row.End : row.End + self.pam_length
                ]
        elif self.direction == "downstream":
            if row.Strand == "+":
                pam_sequence = self.records[row.Chromosome].seq[
                    row.End : row.End + self.pam_length
                ]
            else:
                pam_sequence = self.records[row.Chromosome].seq[
                    row.Start - self.pam_length : row.Start
                ]
        else:
            raise ValueError("direction must be 'upstream' or 'downstream'")

        # If the strand is "-", get the reverse complement of the PAM sequence
        if row.Strand == "-":
            pam_sequence = pam_sequence.reverse_complement()

        return str(pam_sequence)

    def pam_matches(self, sequence):
        # check if the sequence matches the PAM pattern
        return bool(re.search(self.pam_pattern, sequence))


class GenBankParser(Logger):
    def __init__(self, filename):
        super().__init__()
        self.reader = GenBankReader(filename)
        self.records = self.reader.records
        # self.pam_finder = PAMFinder(self.records)
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

    def find_gene_name_for_locus(self, locus_tag):
        # Iterate through all records in the GenBank file
        for record_id, record in self.records.items():
            for feature in record.features:
                if feature.type == "gene":
                    # Check if this gene feature has a locus tag that matches the locus_tag argument
                    if (
                        "locus_tag" in feature.qualifiers
                        and feature.qualifiers["locus_tag"][0] == locus_tag
                    ):
                        # If a gene name is available, return it, otherwise return the locus tag
                        return feature.qualifiers.get("gene", [locus_tag])[0]
        # Return None or locus_tag if not found; depends on how you want to handle not found cases
        return None
