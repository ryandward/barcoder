import gzip
import sys
from concurrent.futures import ThreadPoolExecutor
import pyzstd

def parallel_sort(paired_reads):
    chunk_size = len(paired_reads) // 4
    chunks = [paired_reads[i:i + chunk_size] for i in range(0, len(paired_reads), chunk_size)]
    
    with ThreadPoolExecutor() as executor:
        sorted_chunks = list(executor.map(sorted, chunks))
        
    return sorted_chunks

def merge_sorted_chunks(sorted_chunks):
    from heapq import merge
    return list(merge(*sorted_chunks))

def read_and_process(file1, file2):
    sampled_paired_reads = []
    remaining_paired_reads = []
    
    with gzip.open(file1, 'rt') as f1, gzip.open(file2, 'rt') as f2:
        for i in range(1000):
            f1.readline()
            seq1 = f1.readline().strip()
            f1.readline()
            f1.readline()

            f2.readline()
            seq2 = f2.readline().strip()
            f2.readline()
            f2.readline()
            
            sampled_paired_reads.append((seq1, seq2))
        
        while True:
            f1.readline()
            seq1 = f1.readline().strip()
            f1.readline()
            f1.readline()

            f2.readline()
            seq2 = f2.readline().strip()
            f2.readline()
            f2.readline()
            
            if not seq1:
                break
                
            remaining_paired_reads.append((seq1, seq2))
    
    sorted_chunks = parallel_sort(remaining_paired_reads)
    sorted_remaining_paired_reads = merge_sorted_chunks(sorted_chunks)
    
    final_paired_reads = sampled_paired_reads + sorted_remaining_paired_reads
    
    output_file1 = file1.replace(".fastq.gz", ".reads")
    output_file2 = file2.replace(".fastq.gz", ".reads")
    
    with open(output_file1, "wb") as out1, open(output_file2, "wb") as out2:
        for seq1, seq2 in final_paired_reads:
            out1.write(f"{seq1}\n".encode())
            out2.write(f"{seq2}\n".encode())

    with open(output_file1, 'rb') as f_in, pyzstd.open(output_file1 + '.zst', 'wb') as f_out:
        f_out.write(f_in.read())

    with open(output_file2, 'rb') as f_in, pyzstd.open(output_file2 + '.zst', 'wb') as f_out:
        f_out.write(f_in.read())

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_to_reads.py file1.fastq.gz file2.fastq.gz")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]

    read_and_process(file1, file2)
