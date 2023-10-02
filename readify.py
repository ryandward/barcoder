import gzip
import random
import sys

def readify(file1, file2):
    paired_reads = []
    
    with gzip.open(file1, 'rt') as f1, gzip.open(file2, 'rt') as f2:
        while True:
            f1.readline()  # Skip header
            seq1 = f1.readline().strip()
            f1.readline()  # Skip separator
            f1.readline()  # Skip quality
            
            f2.readline()  # Skip header
            seq2 = f2.readline().strip()
            f2.readline()  # Skip separator
            f2.readline()  # Skip quality
            
            if not seq1:
                break
            
            if all(seq1.count(base) / len(seq1) < 0.5 for base in "ATGC") and \
               all(seq2.count(base) / len(seq2) < 0.5 for base in "ATGC"):
                paired_reads.append((seq1, seq2))
    
    sampled_paired_reads = random.sample(paired_reads, 1000)
    remaining_paired_reads = [x for x in paired_reads if x not in sampled_paired_reads]
    
    remaining_paired_reads.sort()
    
    final_paired_reads = sampled_paired_reads + remaining_paired_reads

    output_file1 = file1.replace(".fastq.gz", ".reads")
    output_file2 = file2.replace(".fastq.gz", ".reads")
    
    with open(output_file1, "w") as out1, open(output_file2, "w") as out2:
        for seq1, seq2 in final_paired_reads:
            out1.write(f"{seq1}\n")
            out2.write(f"{seq2}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_to_reads.py file1.fastq.gz file2.fastq.gz")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    readify(file1, file2)


# zstd compression to implement later
# import os
# import glob
# import gzip
# import pyzstd

# def convert_fastq_gz_to_reads_zst(fastq_gz_file):
#     base_name = os.path.splitext(os.path.splitext(fastq_gz_file)[0])[0]
#     reads_file = base_name + '.reads'
    
#     with gzip.open(fastq_gz_file, 'rt') as f_in, open(reads_file, 'w') as f_out:
#         while True:
#             try:
#                 next(f_in)  # Skip @SEQUENCE_ID
#                 seq = next(f_in).strip()  # Capture SEQUENCE
#                 next(f_in)  # Skip +
#                 next(f_in)  # Skip QUALITY
#                 f_out.write(seq + '\n')
#             except StopIteration:
#                 break

#     # Compress the reads file into .reads.zst
#     with open(reads_file, 'rb') as f_in, pyzstd.open(reads_file + '.zst', 'wb') as f_out:
#         f_out.write(f_in.read())
#     os.remove(reads_file)

# if __name__ == "__main__":
#     for fastq_gz_file in glob.glob('*.fastq.gz'):
#         convert_fastq_gz_to_reads_zst(fastq_gz_file)

