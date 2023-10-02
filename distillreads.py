import gzip
import sys
from concurrent.futures import ThreadPoolExecutor
import pyzstd


def compress_and_write(output_file):
    with open(output_file, 'rb') as f_in, pyzstd.open(output_file + '.zst', 'wb') as f_out:
        f_out.write(f_in.read())


def parallel_sort(paired_reads):
    chunk_size = len(paired_reads) // 4
    chunks = [paired_reads[i:i + chunk_size] for i in range(0, len(paired_reads), chunk_size)]
    
    with ThreadPoolExecutor() as executor:
        sorted_chunks = list(executor.map(sorted, chunks))
        
    return sorted_chunks


def merge_sorted_chunks(sorted_chunks):
    from heapq import merge
    return list(merge(*sorted_chunks))


def read_and_process(*files):
    sampled_paired_reads = []
    remaining_paired_reads = []

    file_handles = [gzip.open(file, 'rt') for file in files]
    
    for _ in range(1000):
        read_tuple = []
        for f in file_handles:
            f.readline()
            seq = f.readline().strip()
            f.readline()
            f.readline()
            read_tuple.append(seq)
        
        sampled_paired_reads.append(tuple(read_tuple))
    
    while True:
        read_tuple = []
        for f in file_handles:
            f.readline()
            seq = f.readline().strip()
            f.readline()
            f.readline()
            read_tuple.append(seq)

        if not all(read_tuple):
            break

        remaining_paired_reads.append(tuple(read_tuple))

    for f in file_handles:
        f.close()
        
    sorted_chunks = parallel_sort(remaining_paired_reads)
    sorted_remaining_paired_reads = merge_sorted_chunks(sorted_chunks)

    final_paired_reads = sampled_paired_reads + sorted_remaining_paired_reads

    output_files = [file.replace(".fastq.gz", ".reads") for file in files]
    
    with open(output_files[0], "wb") as out1:
        for read_tuple in final_paired_reads:
            out1.write(f"{read_tuple[0]}\n".encode())
            
    for i, output_file in enumerate(output_files[1:], start=1):
        with open(output_file, "wb") as out:
            for read_tuple in final_paired_reads:
                out.write(f"{read_tuple[i]}\n".encode())
                
        compress_and_write(output_file)
        
    compress_and_write(output_files[0])


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python convert_to_reads.py file1.fastq.gz [file2.fastq.gz ...]")
        sys.exit(1)

    files = sys.argv[1:]
    read_and_process(*files)
