import pandas as pd
import pyranges as pr
import pysam
from functools import lru_cache
from Bio.Seq import Seq


def rev_comp(seq):
    return str(Seq(seq).reverse_complement())


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

            if read.is_reverse == 0:
                strand = "+"
                sequence = read.query_sequence
            elif read.is_reverse == 1:
                strand = "-"
                sequence = rev_comp(read.query_sequence)
            else:
                strand = "."
                sequence = None

            interval_dict = {
                "Chromosome": read.reference_name,
                "Start": read.reference_start,
                "End": read.reference_end,
                "Mapped": True if not read.is_unmapped else False,
                "Strand": strand,
                "Barcode": sequence,
                "Mismatches": (read.get_tag("NM") if read.has_tag("NM") else "0"),
            }

            data.append(interval_dict)

        df = pd.DataFrame(data)
        ranges = pr.PyRanges(df)
        return ranges
