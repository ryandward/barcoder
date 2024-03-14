import re


class PAMProcessor:
    def __init__(self, records, pam, direction):
        self.records = records
        self.pam = pam.replace("N", "[ATCG]")
        self.direction = direction

    def get_sequence(self, row):
        sequence = self.records[row.Chromosome].seq[row.Start : row.End]
        if row.Strand == "-":
            sequence = sequence.reverse_complement()
        return sequence

    def get_strand(self, strand_symbol):
        # Normalize the input to handle common variations
        normalized_strand = str(strand_symbol).lower().strip()
        if normalized_strand in ("+", "1", "+1", "fwd", "forward"):
            return 1
        elif normalized_strand in ("-", "-1", "rev", "reverse"):
            return -1
        else:
            raise ValueError(f"Unrecognized strand symbol: {strand_symbol}")


class GuideFinder:
    def __init__(self, records, pam, direction, length):
        self.records = records
        self.pam = pam.replace("N", "[ATCG]")
        self.direction = direction
        self.length = length

    def find_guides_from_pam(self):
        guide_sequences = []

        for record_id, record in self.records.items():
            for sequence in [record.seq, record.seq.reverse_complement()]:
                sequence_str = str(sequence)
                for match in re.finditer(self.pam, sequence_str):
                    start = match.start()
                    end = match.end()

                    if self.direction == "downstream":
                        guide_sequence = sequence_str[
                            max(0, start - self.length) : start
                        ]
                    elif self.direction == "upstream":
                        guide_sequence = sequence_str[
                            end : min(end + self.length, len(sequence_str))
                        ]
                    else:
                        raise ValueError("Direction must be 'upstream' or 'downstream'")

                    guide_sequences.append(guide_sequence)

        return guide_sequences


class PAMFinder(PAMProcessor):
    def __init__(self, records, pam, direction):
        super().__init__(records, pam, direction)
        self.pam_length = len(pam)

    def get_pam_seq(self, row):
        sequence = self.get_sequence(row)
        strand = self.get_strand(row.Strand)

        if self.direction == "upstream":
            pam_sequence = (
                self.records[row.Chromosome].seq[row.End : row.End + self.pam_length]
                if strand == 1
                else self.records[row.Chromosome].seq[
                    row.Start - self.pam_length : row.Start
                ]
                # .reverse_complement()
            )

        elif self.direction == "downstream":
            pam_sequence = (
                self.records[row.Chromosome].seq[row.End : row.End + self.pam_length]
                if strand == 1
                else self.records[row.Chromosome].seq[
                    row.Start - self.pam_length : row.Start
                ]
                # .reverse_complement()
            )
        else:
            raise ValueError("direction must be 'upstream' or 'downstream'")

        if strand == -1:
            pam_sequence = pam_sequence.reverse_complement()

        return str(pam_sequence)

    def pam_matches(self, sequence):
        return bool(re.search(self.pam, sequence))
