import argparse
import gzip
import os
import platform
import sys
from collections import Counter, defaultdict
from contextlib import nullcontext
from datetime import datetime
from multiprocessing import Pool, cpu_count
from typing import List, Tuple
from collections import Counter
from typing import List, Set, Tuple

from collections import Counter
from typing import List, Tuple, Set

import zstandard as zstd
from Bio.Seq import Seq
from rich.console import Console
from rich.table import Table

import rich

def safe_len(s):
    return 0 if s is None else len(s)

def read_fasta(fasta_file) -> Set[str]:
    barcodes = set()
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    with open_func(fasta_file, 'rt') as f:
        for line in f:
            if not line.startswith(">"):
                barcodes.add(line.strip())
    return barcodes

def open_file(file_path, mode):
    if file_path.endswith('.gz'):
        return gzip.open(file_path, mode)
    elif file_path.endswith('.zst'):
        return zstd.open(file_path, mode)
    elif file_path.endswith('.fastq'):
        return open(file_path, mode)
    elif file_path.endswith('.reads'):
        return open(file_path, mode)
    else:
        raise ValueError("Unsupported file type.")

def read_in_chunks(file1, file2=None, chunk_size=2**16) -> Tuple[List[str], List[str]]:

    reads1, reads2 = [], []
    
    # Remove compression extension if present
    stripped_file1 = file1
    if file1.endswith('.gz') or file1.endswith('.zst'):
        stripped_file1 = os.path.splitext(file1)[0]
        
    # Determine the file type based on the stripped extension
    if stripped_file1.endswith('.fastq'):
        file_type = 'fastq'
    elif stripped_file1.endswith('.reads'):
        file_type = 'reads'
    else:
        raise ValueError("Unsupported file type. Must be '.fastq' or '.reads'.")
    
    with open_file(file1, 'rt') as f1, (open_file(file2, 'rt') if file2 else nullcontext()) as f2:
        iters = [iter(f1), iter(f2) if f2 else iter([])]
        
        while True:
            try:
                if file_type == 'fastq':
                    next(iters[0])  # Skip @SEQUENCE_ID
                    reads1.append(next(iters[0]).strip())  # Capture SEQUENCE
                    next(iters[0])  # Skip +
                    next(iters[0])  # Skip QUALITY

                    if file2:
                        next(iters[1])  # Skip @SEQUENCE_ID
                        reads2.append(next(iters[1]).strip())  # Capture SEQUENCE
                        next(iters[1])  # Skip +
                        next(iters[1])  # Skip QUALITY

                elif file_type == 'reads':
                    reads1.append(next(iters[0]).strip())  # Capture SEQUENCE

                    if file2:
                        reads2.append(next(iters[1]).strip())  # Capture SEQUENCE
                
               
                if len(reads1) >= chunk_size:
                    yield (reads1[:chunk_size], reads2[:chunk_size] if reads2 else None)
                    reads1, reads2 = reads1[chunk_size:], reads2[chunk_size:] if reads2 else []
                    
            except StopIteration:
                break
                
        if reads1:
            yield (reads1, reads2 if reads2 else None)


def contains_bc(seq, barcodes, bc_len) -> bool:
    if seq is None:
        return False
    for i in range(0, len(seq) - bc_len + 1):
        kmer = seq[i:i + bc_len]
        if kmer in barcodes:
            return True
    return False


def sample_data(file1, file2, barcodes, is_paired) -> Tuple[int, bool, Set[str], Set[str]]:
    bc_len = len(next(iter(barcodes)))
    valid_samples1 = set()
    valid_samples2 = set()
    total_reads_sampled = 0
    new_reads_sampled = 0

    chunk_generator = read_in_chunks(file1, file2 if is_paired else None, chunk_size=len(barcodes))

    for read1_chunk, read2_chunk in chunk_generator:            
        for read1, read2 in zip(read1_chunk, read2_chunk if read2_chunk else [None]*len(read1_chunk)):
            total_reads_sampled += 1
            
            if read2 and (read1 in valid_samples1 or (read2 in valid_samples2)):
                continue
            elif read1 in valid_samples1:
                continue

            new_reads_sampled += 1

            bc_F_in_read1 = contains_bc(read1, barcodes, bc_len) if read1 else None
            bc_F_in_read2 = contains_bc(read2, barcodes, bc_len) if read2 else None

            bc_R_in_read1 = contains_bc(read1[::-1 or 0].translate(str.maketrans("ATCGN", "TAGCN")), barcodes, bc_len) if read1 else None
            bc_R_in_read2 = contains_bc(read2[::-1 or 0].translate(str.maketrans("ATCGN", "TAGCN")), barcodes, bc_len) if read2 else None
            
            if is_paired:
                if bc_F_in_read1 and bc_R_in_read2:
                    valid_samples1.add(read1)
                    valid_samples2.add(read2)
                if bc_R_in_read1 and bc_F_in_read2:
                    valid_samples1.add(read1)
                    valid_samples2.add(read2)
            else:
                if bc_F_in_read1:
                    valid_samples1.add(read1)
                    valid_samples2 = None # Placeholder until we can determine the orientation of read1
                if bc_R_in_read1:
                    valid_samples1.add(read1)
                    valid_samples2 = None # Placeholder until we can determine the orientation of read1

            if len(valid_samples1) > len(barcodes):
                break

    total_count1 = sum(1 for read in valid_samples1 if read and contains_bc(read, barcodes, bc_len)) if valid_samples1 else 0
    total_count2 = sum(1 for read in valid_samples2 if read and contains_bc(read, barcodes, bc_len)) if valid_samples2 else 0

    need_swap = False
    if is_paired:
        need_swap = total_count2 > total_count1
    else:  # Single-end
        # Count barcodes in reverse complement of read1
        rc_count = sum(1 for read in valid_samples1 if contains_bc(read[::-1].translate(str.maketrans("ATCGN", "TAGCN")), barcodes, bc_len))
        need_swap = rc_count > total_count1

    return new_reads_sampled, need_swap, valid_samples1, valid_samples2

def find_start_positions(reads, barcodes, bc_len, rev=False):

    if reads is None:
        return None
    
    offset_counts = Counter()
    for read in reads:

        for i in range(len(read) - bc_len + 1):
            kmer = read[i:i+bc_len]
            if rev:
                kmer = kmer[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            if kmer in barcodes:
                offset_counts[i] += 1
    return offset_counts.most_common(1)[0][0] if offset_counts else None


def generate_kmers(seq: str, k: int) -> Set[str]:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def find_ends(reads: List[str], barcodes: Set[str], start: int, bc_len: int, max_flank: int = 4, rev: bool = False) -> Tuple[str, str]:
    # Initialize maximum lengths for L_flank and R_flank flanks; reverse if the 'rev' flag is set.
    L_max, R_max = (max_flank, max_flank) if not rev else (max_flank, max_flank)[::-1]

    # Counters to store occurrences of L_flank and R_flank flanking sequences.
    L_flanks, R_flanks = Counter(), Counter()

    def update_flanks(side: str, seq: str, max_len: int):
        """Update counters with sub-sequences of varying lengths."""
        counts = L_flanks if side == "L_flank" else R_flanks
        for i in range(max_len, 0, -1):
            truncated = seq[-i:] if side == "L_flank" else seq[:i]
            counts[truncated] += 1

    for read in reads:
        # Generate kmers for the current read
        kmers = generate_kmers(read, bc_len)
        # Identify overlapping barcodes with the read's kmers
        overlapping_bcs = barcodes & kmers
        if not overlapping_bcs:
            continue

        # Find the starting position of the barcode in the read
        common_start = read.find(next(iter(overlapping_bcs)))
        # Update maximum possible L_flank and R_flank flanks based on the found position
        L_max = min(L_max, common_start)
        R_max = min(R_max, len(read) - (common_start + bc_len))

        # Extract potential flanking sequences from the read
        L_flank = read[common_start - L_max: common_start]
        R_flank = read[common_start + bc_len: common_start + bc_len + R_max]

        # Update flank counters with sequences of varying lengths
        update_flanks("L_flank", L_flank, L_max)
        update_flanks("R_flank", R_flank, R_max)

    def extract_best_flank(counts: Counter) -> str:
        """Identify the most suitable flanking sequence. 
        To justify using a shorter flank, it must be at least twice as common as the next longer one."""
        most_common_next = None
        for bc_len in range(max_flank, 0, -1):
            # Find sequences of the current length from the counter
            potential_seqs = [seq for seq in counts if len(seq) == bc_len]
            if not potential_seqs:
                continue
            # Identify the most common sequence of the current length
            most_common = max(potential_seqs, key=lambda x: counts[x])
            # If it's the maximum length or meets the "twice as common" criteria, return it
            if bc_len == max_flank or counts[most_common] >= 2 * counts[most_common_next]:
                return most_common
            most_common_next = most_common
        return None

    # Extract the most appropriate L_flank and R_flank flanks using the criteria defined above
    L_most_common = extract_best_flank(L_flanks)
    R_most_common = extract_best_flank(R_flanks)

    return L_most_common, R_most_common


def process_chunk(chunk, bcs_with_flanks_fwd, bcs_with_flanks_rev, L_fwd_start, L_rev_start, bc_len, L_fwd, R_fwd, L_rev, R_rev, need_swap) -> Tuple[Counter, int]:
    
    counts = Counter()

    if need_swap: reads2, reads1 = chunk
    else: reads1, reads2 = chunk

    def validate_read(seq_with_flanks, L_flank, R_flank, rev=False):
        if rev: in_bcs_with_flanks = seq_with_flanks in bcs_with_flanks_rev
        else: in_bcs_with_flanks = seq_with_flanks in bcs_with_flanks_fwd

        seq = seq_with_flanks[safe_len(L_flank):safe_len(seq_with_flanks) - safe_len(R_flank)]
        has_flanks = seq_with_flanks.startswith(L_flank or '') and seq_with_flanks.endswith(R_flank or '')
        
        return in_bcs_with_flanks, has_flanks, seq

    L_fwd_len = safe_len(L_fwd)
    R_fwd_len = safe_len(R_fwd)
    L_rev_len = safe_len(L_rev)
    R_rev_len = safe_len(R_rev)

    def process_paired_end(reads1, reads2, L_fwd_start, L_fwd_len, R_fwd_len, L_rev_start, L_rev_len, R_rev_len, L_fwd, R_fwd, L_rev, R_rev):
        if len(reads1) != len(reads2):
            raise ValueError("Length of reads1 and reads2 should be the same for paired-end data.")
        
        for record_fwd, record_rev in zip(reads1, reads2):
            if 'N' in record_fwd:
                continue
            
            seq_with_flanks_fwd = record_fwd[L_fwd_start:L_fwd_start + L_fwd_len + bc_len + R_fwd_len]
            seq_with_flanks_rev = record_rev[L_rev_start:L_rev_start + L_rev_len + bc_len + R_rev_len]
            
            in_bcs_with_flanks_fwd, has_flanks_fwd, seq1 = validate_read(seq_with_flanks_fwd, L_fwd, R_fwd)
            in_bcs_with_flanks_rev, has_flanks_rev, seq2 = validate_read(seq_with_flanks_rev, L_rev, R_rev, rev=True)

            if seq1 != seq2[::-1].translate(str.maketrans("ATCGN", "TAGCN")): 
                continue
            
            if in_bcs_with_flanks_fwd and in_bcs_with_flanks_rev and has_flanks_fwd and has_flanks_rev:
                counts[seq1] += 1
            
            elif has_flanks_fwd and has_flanks_rev:
                counts[seq1 + "*"] += 1
    
    def process_single_end(reads, start_pos, len_L, len_R, L_flank, R_flank, reverse=False):
        for record in reads:
            if 'N' in record:
                continue
            
            seq_with_flanks = record[start_pos:start_pos + len_L + bc_len + len_R]
            in_bcs_with_flanks, has_flanks, seq = validate_read(seq_with_flanks, L_flank, R_flank, rev=reverse)
            
            if in_bcs_with_flanks and has_flanks:
                counts[seq] += 1
            elif has_flanks:
                counts[seq + "*"] += 1

    if reads1 and reads2:  # Paired-end processing
        process_paired_end(reads1, reads2, L_fwd_start, L_fwd_len, R_fwd_len, L_rev_start, L_rev_len, R_rev_len, L_fwd, R_fwd, L_rev, R_rev)
    elif reads1:  # Single-end processing, fwd orientation
        process_single_end(reads1, L_fwd_start, L_fwd_len, R_fwd_len, L_fwd, R_fwd)
    elif reads2:  # Single-end processing, reverse complement
        process_single_end(reads2, L_rev_start, L_rev_len, R_rev_len, L_rev, R_rev, reverse=True)

    return (counts, len(reads1) if reads1 else len(reads2))


def main(args):
    console = Console(stderr=True, highlight=True)
    console.log("[bold red]Initializing heuristic barcode counting[/bold red]...")

    reads1 = args.file1
    reads2 = args.file2

    num_threads = cpu_count()//2

    # Reading FASTA File
    console.log("Reading barcodes...")
    barcodes = read_fasta(args.fasta_file)
    bcs_rev = {barcode[::-1].translate(str.maketrans("ATCGN", "TAGCN")) for barcode in barcodes}
    
    bc_len = len(next(iter(barcodes)))

    is_paired = bool(args.file2)
        
    # Reading reads Files
    console.log("Sampling initial reads for orientation...")

    # Initialize
    console.log("Determining orientation of reads...")
    new_reads_sampled, need_swap, sample1, sample2 = sample_data(reads1, reads2, barcodes, is_paired)
    console.log(f"Sampled {new_reads_sampled:,} unique reads and found {safe_len(sample1):,} File_1 and {safe_len(sample2):,} File_2 valid matches...")
    console.log(f"Swapping reads..." if need_swap else f"Proceeding with reads as is...")                

    # Apply the swap logic
    if need_swap:
        sample1, sample2 = sample2, sample1
        reads1, reads2 = reads2, reads1

    # Initialize to None
    bc_start1 = None
    bc_start2 = None        

    # Skip finding barcode starts for single-end that needed a swap
    if sample1 is not None:
        console.log("Finding fwd coordinates...")
        bc_start1 = find_start_positions(sample1, barcodes, bc_len)
        if bc_start1 is None:
            console.log("[bold red]No barcodes found in sample 1. Exiting.[/bold red]")
            sys.exit(1)

    # For paired-end or single-end that needed a swap
    if sample2:
        console.log("Finding reverse coordinates...")
        bc_start2 = find_start_positions(sample2, barcodes, bc_len, rev=True)
        if bc_start2 is None:
            console.log("[bold red]No barcodes found in sample 2. Exiting.[/bold red]")
            sys.exit(1)
    else:
        # Set bc_start2 to an empty list to be consistent with bc_start1
        bc_start2 = None

    # Find flanking sequences
    if sample1 is not None:
        console.log("Identifying forward read flanks...")
        L_fwd, R_fwd = find_ends(sample1, barcodes, bc_start1, bc_len)
        L_fwd_start = bc_start1 - len(L_fwd) if L_fwd else 0
    else:
        L_fwd, R_fwd = None, None
        L_fwd_start = None

    if sample2 is not None:
        console.log("Identifying reverse read flanks...")
        L_rev, R_rev = find_ends(sample2, bcs_rev, bc_start2, bc_len, rev=True)
        L_rev_start =  bc_start2 - len(L_rev) if L_rev else 0
    else:
        L_rev, R_rev = None, None
        L_rev_start = None


   # Calculate reverse complements
    L_rev_rev = L_rev[::-1].translate(str.maketrans("ATCGN", "TAGCN")) if L_rev else None
    R_rev_rev = R_rev[::-1].translate(str.maketrans("ATCGN", "TAGCN")) if R_rev else None

    # Check if the fwd and reverse flanking sequences are reverse complements
    L_min_len = R_min_len = 0

    if L_fwd and R_rev_rev:
        L_min_len = min(len(L_fwd), len(R_rev_rev))
        if L_fwd[-L_min_len:] != R_rev_rev[:L_min_len]:
            console.log("[bold red]Error: Forward and reverse L_flank flanking sequences are not reverse complements.[/bold red]")
            sys.exit(1)
        
    if R_fwd and L_rev_rev:
        R_min_len = min(len(R_fwd), len(L_rev_rev))
        if R_fwd[:R_min_len] != L_rev_rev[:R_min_len]:
            console.log("[bold red]Error: Forward and reverse R_flank flanking sequences are not reverse complements.[/bold red]")
            sys.exit(1)

    def add_flank(barcodes, L_flank=None, R_flank=None, rev=False):
        if rev:
            barcodes = {barcode[::-1].translate(str.maketrans("ATCGN", "TAGCN")) for barcode in barcodes}
        L_flank, R_flank = (L_flank or ""), (R_flank or "")
        return {L_flank + barcode + R_flank for barcode in barcodes}
    
    # Add flanks to the ends of the barcodes, if present

    bcs_with_flanks_fwd = add_flank(barcodes, L_fwd, R_fwd)
    bcs_with_flanks_rev = add_flank(bcs_rev, L_rev, R_rev)

    console.log("Executing high-throughput read analysis...")
    chunk_generator = read_in_chunks(args.file1, args.file2 if is_paired else None)

    # Create argument generator
    args_generator = ((chunk, bcs_with_flanks_fwd, bcs_with_flanks_rev, L_fwd_start, L_rev_start, bc_len, L_fwd, R_fwd, L_rev, R_rev, need_swap) for chunk in chunk_generator)

    # Execute multiprocessing
    with Pool(num_threads) as pool:
        results = pool.starmap(process_chunk, args_generator)
    # Collating Results
    console.log("[bold red]Collating results[/bold red]...")
    
    doc_bcs = Counter()
    undoc_bcs = Counter()

    total_reads = 0

    for result, chunk_size in results:
        total_reads += chunk_size
        for barcode, count in result.items():
            if barcode.endswith('*'):
                undoc_bcs[barcode] += count
            else:
                doc_bcs[barcode] += count

    if is_paired:
        file1_filename = os.path.basename(args.file1) if not need_swap else os.path.basename(args.file2)
        file2_filename = os.path.basename(args.file2) if not need_swap else os.path.basename(args.file1)
    else:
        file1_filename = os.path.basename(args.file1) if not need_swap else None
        file2_filename = None if not need_swap else os.path.basename(args.file1)

    # Table

    # Create a single table with enhanced styles
    combined_table = Table(
        # title="Summary",
        box=rich.table.box.SIMPLE_HEAVY,
        caption="Finished at [u]{}[/u]".format(datetime.now()),
        title_style="bold bright_white",
        caption_style="bold white",
        header_style="bold bright_white",
        border_style="bold bright_white",
        show_header=True
    )
    # Define columns with justifications
    combined_table.add_column(os.path.basename(sys.argv[0]), justify="right", style="white", min_width=30)
    combined_table.add_column("Summary", justify="right", style="bold bright_white", min_width=20)

    # Input & Configuration Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_magenta]Input & Config[/bold bright_magenta]", "")

    if args.fasta_file: combined_table.add_row("Barcodes", f"[bold]{os.path.basename(args.fasta_file)}[/bold]")
    if file1_filename: combined_table.add_row("Forward Reads", f"[bold]{file1_filename}[/bold]")
    if file2_filename: combined_table.add_row("Reverse Reads", f"[bold]{file2_filename}[/bold]")
    if num_threads: combined_table.add_row("Threads", f"[bold]{num_threads}[/bold]")
    if platform.system: combined_table.add_row("Operating System", f"[bold]{platform.system()}[/bold]")

    # Heuristic Statistics Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_blue]Heuristics[/bold bright_blue]", "")

    if bc_len: combined_table.add_row("Barcode Length", f"[bold]{bc_len}[/bold]")
    if bc_start1:combined_table.add_row("Forward Offset", f"[bold]{bc_start1}[/bold]")
    if bc_start2:combined_table.add_row("Reverse Offset", f"[bold]{bc_start2}[/bold]")
    if L_fwd or R_fwd: combined_table.add_row("Forward Flanks", f"[bold]{L_fwd}...{R_fwd}[/bold]")
    if L_rev or R_rev: combined_table.add_row("Reverse Flanks", f"[bold]{L_rev}...{R_rev}[/bold]")

    # Numeric Statistics Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_green]Barcode Alignment Stats[/bold bright_green]", "")

    combined_table.add_row("Declared Barcodes", f"[bold]{len(barcodes):,}[/bold]")
    combined_table.add_row("Seen Documented Barcodes", f"[bold]{len(doc_bcs):,}[/bold]")
    combined_table.add_row("Unseen Documented Barcodes", f"[bold]{(len(barcodes) - len(doc_bcs)):,}[/bold]")
    combined_table.add_row("Found Undocumented Barcodes", f"[bold]{len(undoc_bcs):,}[/bold]")
    combined_table.add_row("Total Reads", f"[bold]{total_reads:,}[/bold]")
    combined_table.add_row("Documented Barcode Reads", f"[bold]{sum(doc_bcs.values()):,}[/bold]")
    combined_table.add_row("Undocumented Barcode Reads", f"[bold]{sum(undoc_bcs.values()):,}[/bold]")
    combined_table.add_row("Documented Fraction", f"[bold]{(sum(doc_bcs.values()) / total_reads if total_reads != 0 else 0):.3f}[/bold]")
    combined_table.add_row("Undocumented Fraction", f"[bold]{(sum(undoc_bcs.values()) / total_reads if total_reads != 0 else 0):.3f}[/bold]", end_section=True)

    # Sequence Information Sub-heading  
    combined_table.add_section()
    top_doc_bcs = min(5, len(doc_bcs))
    combined_table.add_row(f"[bold bright_yellow]Top {top_doc_bcs} Documented Barcodes[/bold bright_yellow]", "")
    for idx, (barcode, count) in enumerate(doc_bcs.most_common(top_doc_bcs)):
        end_section = idx == (top_doc_bcs - 1)
        combined_table.add_row(barcode, f"{count:,}", end_section=end_section)

    combined_table.add_section()
    top_undoc_bcs = min(5, len(undoc_bcs))
    combined_table.add_row(f"[bold bright_red]Top {top_undoc_bcs} Undocumented Barcodes[/bold bright_red]", "")
    for idx, (barcode, count) in enumerate(undoc_bcs.most_common(top_undoc_bcs)):
        end_section = idx == (top_undoc_bcs - 1)
        combined_table.add_row(barcode, f"{count:,}", end_section=end_section)

    # Print the combined table
    console.log(combined_table)

    for barcode, count in doc_bcs.items():
        print("\t".join([barcode, str(count)]))

    for barcode, count in undoc_bcs.items():
        print("\t".join([barcode, str(count)]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Barcodes.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file.')
    parser.add_argument('file1', type=str, help='First reads file: FASTQ or raw reads.')
    parser.add_argument('file2', type=str, nargs='?', default=None, help='Second reads file: FASTQ or raw reads (optional).')
    args = parser.parse_args()
    main(args)