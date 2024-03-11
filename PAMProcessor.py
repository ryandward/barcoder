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


class GuideFinderByPAM(PAMProcessor):
    def __init__(self, records, pam, direction, length):
        super().__init__(records, pam, direction)
        self.length = length

    def find_pam_sequences(self):
        pam_sequences = []

        for record_id, record in self.records.items():
            for strand, sequence in [
                (+1, record.seq),
                (-1, record.seq.reverse_complement()),
            ]:
                sequence_str = str(sequence)
                for match in re.finditer(self.pam, sequence_str):
                    start = match.start()
                    end = match.end()

                    if (
                        self.direction == "downstream"
                        and strand == 1
                        or self.direction == "upstream"
                        and strand == -1
                    ):
                        if start >= self.length:
                            guide_sequence = sequence_str[start - self.length : start]
                            pam_sequences.append(guide_sequence)
                    elif (
                        self.direction == "upstream"
                        and strand == 1
                        or self.direction == "downstream"
                        and strand == -1
                    ):
                        if end + self.length <= len(sequence_str):
                            guide_sequence = sequence_str[end : end + self.length]
                            pam_sequences.append(guide_sequence)

        return pam_sequences


class PAMRetriever(PAMProcessor):
    def __init__(self, records, pam, direction):
        super().__init__(records, pam, direction)
        self.pam_length = len(pam)

    def get_pam_seq(self, row):
        sequence = self.get_sequence(row)

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

        if row.Strand == "-":
            pam_sequence = pam_sequence.reverse_complement()

        return str(pam_sequence)

    def pam_matches(self, sequence):
        return bool(re.search(self.pam, sequence))
