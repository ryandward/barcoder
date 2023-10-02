import os
import glob
import gzip
import pyzstd

def convert_fastq_gz_to_reads_zst(fastq_gz_file):
    base_name = os.path.splitext(os.path.splitext(fastq_gz_file)[0])[0]
    reads_file = base_name + '.reads'
    
    with gzip.open(fastq_gz_file, 'rt') as f_in, open(reads_file, 'w') as f_out:
        while True:
            try:
                next(f_in)  # Skip @SEQUENCE_ID
                seq = next(f_in).strip()  # Capture SEQUENCE
                next(f_in)  # Skip +
                next(f_in)  # Skip QUALITY
                f_out.write(seq + '\n')
            except StopIteration:
                break

    # Compress the reads file into .reads.zst
    with open(reads_file, 'rb') as f_in, pyzstd.open(reads_file + '.zst', 'wb') as f_out:
        f_out.write(f_in.read())
    os.remove(reads_file)

if __name__ == "__main__":
    for fastq_gz_file in glob.glob('*.fastq.gz'):
        convert_fastq_gz_to_reads_zst(fastq_gz_file)
