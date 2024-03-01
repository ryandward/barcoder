import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def create_topological_fasta(
    genbank_file_name, topological_fasta_file_name, overhang_length=0
):
    topological_records = []
    open_func = gzip.open if genbank_file_name.endswith(".gz") else open

    with open_func(genbank_file_name, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            if record.annotations.get("topology", None) == "circular":
                overhang_length = 100_000

            new_seq = record.seq + record.seq[:overhang_length]
            topological_record = SeqRecord(
                Seq(str(new_seq)), id=record.id, description=record.description
            )
            topological_records.append(topological_record)

    with open(topological_fasta_file_name, "w") as output_handle:
        SeqIO.write(topological_records, output_handle, "fasta")


def create_fake_topological_fastq(topological_fasta_file_name, fastq_file):
    if topological_fasta_file_name.endswith(".gz"):
        with gzip.open(topological_fasta_file_name, "rt") as input_handle, open(
            fastq_file, "w"
        ) as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")
    else:
        with open(topological_fasta_file_name, "rt") as input_handle, open(
            fastq_file, "w"
        ) as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")
