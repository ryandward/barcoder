from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation
from functools import lru_cache
import json
from itertools import islice
from rich.console import Console
from rich.logging import RichHandler
import logging
from functools import wraps
from collections import defaultdict

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


# Class to process GenBank files
class GenBankProcessor:
    def __init__(self, filename):
        self.filename = filename
        self.console = Console(stderr=True, style="bold blue")

    @LoggingBase.logged_property
    @lru_cache(maxsize=None)
    def records(self):
        return SeqIO.index(self.filename, "genbank")

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

    @LoggingBase.logged_property
    @lru_cache(maxsize=None)
    def locus_map(self):

        def update_locus_map(
            locus_map, id, position, locus_tag, gene_name, start, end, feature
        ):
            if position not in locus_map[id]:
                locus_map[id][position] = []
            locus_map[id][position].append(
                (locus_tag, gene_name, int(start), int(end), feature.location.strand)
            )
            return locus_map

        overhang_continue = defaultdict(int)
        locus_map = defaultdict(dict)

        for id, record in self.records.items():
            locus_map[id] = {}
            overhang_continue[id] = 0

            for feature in record.features:
                if feature.type == "gene":
                    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                    gene_name = feature.qualifiers.get("gene", [None])[0]

                    if isinstance(feature.location, CompoundLocation) and any(
                        part.start == 0 or part.end == len(record.seq)
                        for part in feature.location.parts
                    ):

                        # Find segments that wrap around the genome
                        genome_end_segment = next(
                            part
                            for part in feature.location.parts
                            if part.end == len(record.seq)
                        )
                        genome_start_segment = next(
                            part for part in feature.location.parts if part.start == 0
                        )

                        # Adjust positions
                        adj_start = int(genome_end_segment.start)
                        adj_end = int(genome_start_segment.end + len(record.seq))
                        overhang_continue[id] = int(genome_start_segment.end)

                        # Iterate over the adjusted range
                        for position in range(adj_start, adj_end):
                            update_locus_map(
                                locus_map,
                                id,
                                position,
                                locus_tag,
                                gene_name,
                                int(adj_start),
                                int(adj_end),
                                feature,
                            )

                    else:
                        # normal genes
                        # Check the type of feature.location before the loop
                        if isinstance(feature.location, CompoundLocation):
                            locations = feature.location.parts
                        else:
                            locations = [feature.location]

                        # Now you can iterate over locations without creating a new list in each iteration
                        for part_location in locations:
                            # rest of your code
                            for position in range(
                                int(part_location.start), int(part_location.end)
                            ):
                                update_locus_map(
                                    locus_map,
                                    id,
                                    position,
                                    locus_tag,
                                    gene_name,
                                    int(part_location.start),
                                    int(part_location.end),
                                    feature,
                                )

                                if (
                                    overhang_continue[id]
                                    <= position
                                    < self.overhangs[id]
                                ):
                                    update_locus_map(
                                        locus_map,
                                        id,
                                        position + len(record.seq),
                                        locus_tag,
                                        gene_name,
                                        int(part_location.start) + len(record.seq),
                                        int(part_location.end) + len(record.seq),
                                        feature,
                                    )

        return locus_map


# Create a GenBankProcessor object
gbp = GenBankProcessor("GCA_003054575.1.gb")

