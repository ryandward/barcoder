import argparse
import gzip
import os
import platform
import sys
from collections import Counter, defaultdict
from contextlib import nullcontext
from datetime import datetime
from multiprocessing import Pool, cpu_count

import rich
import zstandard as zstd
from Bio.Seq import Seq
from rich.console import Console
from rich.table import Table


def read_fasta(fasta_file):
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

def read_in_chunks(file1, file2=None, chunk_size=4096):

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

def count_barcodes(seq, barcode_length, barcodes):
    # Generate k-mers and count matches with barcodes
    kmers = {seq[i:i + barcode_length] for i in range(len(seq) - barcode_length + 1)}
    return len(barcodes & kmers)

def contains_barcode(seq, barcodes, barcode_length):
    """
    Slides a window of size barcode_length across the sequence and 
    checks for the presence of any matching barcode.
    """
    if seq is None:
        return False
    for i in range(0, len(seq) - barcode_length + 1):
        kmer = seq[i:i + barcode_length]
        if kmer in barcodes:
            return True
    return False

def count_barcodes_in_sequence(seq, barcodes, barcode_length):
    """
    Counts occurrences of barcodes both in forward and reverse complement of the sequence.
    """
    if seq is None:
        return 0
    count = 0
    if contains_barcode(seq, barcodes, barcode_length):
        count += 1

    seq_rc = seq[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
    if contains_barcode(seq_rc, barcodes, barcode_length):
        count += 1

    return count

def determine_forward_read_with_generator(file1, file2, barcodes, is_paired_end):
    barcode_length = len(next(iter(barcodes)))
    valid_samples1 = []
    valid_samples2 = []
    reads_sampled = 0

    chunk_generator = read_in_chunks(file1, file2 if is_paired_end else None, chunk_size=2**16)

    for s1, s2 in chunk_generator:            
        for seq1, seq2 in zip(s1, s2 if s2 else [None]*len(s1)):
            reads_sampled += 1

            valid_seq1_fwd = contains_barcode(seq1, barcodes, barcode_length) if seq1 else None
            valid_seq2_fwd = contains_barcode(seq2, barcodes, barcode_length) if seq2 else None

            valid_seq1_rev = contains_barcode(seq1[::-1 or 0].translate(str.maketrans("ATCGN", "TAGCN")), barcodes, barcode_length) if seq1 else None
            valid_seq2_rev = contains_barcode(seq2[::-1 or 0].translate(str.maketrans("ATCGN", "TAGCN")), barcodes, barcode_length) if seq2 else None
            
            if is_paired_end:
                if valid_seq1_fwd and valid_seq2_rev:
                    valid_samples1.append(seq1)
                    valid_samples2.append(seq2)
                elif valid_seq1_rev and valid_seq2_fwd:
                    valid_samples1.append(seq1)
                    valid_samples2.append(seq2)
            else:
                if valid_seq1_fwd:
                    valid_samples1.append(seq1)
                    valid_samples2 = None # Placeholder until we can determine the orientation of sample1
                elif valid_seq1_rev:
                    valid_samples1.append(seq1)
                    valid_samples2 = None # Placeholder until we can determine the orientation of sample1

            if len(valid_samples1) >= 2**16 and len(valid_samples1) >= reads_sampled * 0.5:
                break

    total_count1 = sum(1 for seq in valid_samples1 if contains_barcode(seq, barcodes, barcode_length)) if valid_samples1 else 0
    total_count2 = sum(1 for seq in valid_samples2 if seq and contains_barcode(seq, barcodes, barcode_length)) if valid_samples2 else 0

    needs_swap = False
    if is_paired_end:
        needs_swap = total_count2 > total_count1
    else:  # Single-end
        # Count barcodes in reverse complement of sample1
        rc_count = sum(1 for seq in valid_samples1 if contains_barcode(seq[::-1].translate(str.maketrans("ATCGN", "TAGCN")), barcodes, barcode_length))
        needs_swap = rc_count > total_count1

    return reads_sampled, needs_swap, valid_samples1, valid_samples2


def find_ends(reads, start, length, reverse_strand=False):

    # Initialize max values; reverse if needed
    max_left, max_right = (4, 4) if not reverse_strand else (4, 4)[::-1]

    # Initialize Counters for left and right flanking regions
    lefts, rights = Counter(), Counter()

    for read in reads:
        max_left = min(max_left, start)
        max_right = min(max_right, len(read) - (start + length))

        # Extract flanking sequences
        left = read[start - max_left: start]
        right = read[start + length: start + length + max_right]

        # Update Counters
        if left:
            lefts[left] += 1
        if right:
            rights[right] += 1

    # Find the most common flanking sequence
    left_most_common = lefts.most_common(1)[0][0] if lefts else None
    right_most_common = rights.most_common(1)[0][0] if rights else None

    return left_most_common, right_most_common

def process_chunk(chunk, barcodes, barcode_start1, barcode_start2, barcode_length, left1, right1, left2, right2, need_swap):

    reads1, reads2 = chunk

    if need_swap:
        reads1, reads2 = reads2, reads1

    counts = defaultdict(int)
    left1_len = len(left1) if left1 else 0
    right1_len = len(right1) if right1 else 0
    left2_len = len(left2) if left2 else 0
    right2_len = len(right2) if right2 else 0

    if reads1 and reads2:  # Paired-end processing
        for rec1, rec2 in zip(reads1, reads2):
            candidate1 = rec1[barcode_start1:barcode_start1 + barcode_length]
            candidate2 = rec2[barcode_start2:barcode_start2 + barcode_length]
            candidate2_rc = candidate2[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            if candidate2_rc == candidate1:
                if candidate1 in barcodes:
                    counts[candidate1[left1_len:-right1_len or None]] += 1
                elif candidate1.startswith(left1 or '') and candidate1.endswith(right1 or ''):
                    barcode = candidate1[left1_len:-right1_len or None] + "*"  
                    counts[barcode] += 1

    elif reads1:  # Single-end processing, forward orientation
        for rec1 in reads1:
            candidate1 = rec1[barcode_start1:barcode_start1 + barcode_length]
            if candidate1 in barcodes:
                counts[candidate1[left1_len:-right1_len or None]] += 1
            elif candidate1.startswith(left1 or '') and candidate1.endswith(right1 or ''):
                barcode = candidate1[left1_len:-right1_len or None] + "*"  
                counts[barcode] += 1

    elif reads2:  # Single-end processing, reverse complement
        for rec2 in reads2:
            candidate2 = rec2[barcode_start2:barcode_start2 + barcode_length]
            candidate2_rc = candidate2[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            if candidate2 in barcodes:
                counts[candidate2_rc[right2_len:-left2_len or None]] += 1
            elif candidate2.startswith(left2 or '') and candidate2.endswith(right2 or ''):
                barcode = candidate2_rc[right2_len:-left2_len or None] + "*"  
                counts[barcode] += 1
    
    if reads1 and reads2 and len(reads1) != len(reads2):
        raise ValueError("Length of reads1 and reads2 should be the same for paired-end data.")

    return (counts, len(reads1) if reads1 else len(reads2))

def find_start_positions(reads, barcodes, barcode_length, reverse_strand=False):

    if reads is None:
        return None
    
    offset_counts = Counter()
    for read in reads:

        for i in range(len(read) - barcode_length + 1):
            kmer = read[i:i+barcode_length]
            if reverse_strand:
                kmer = kmer[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            if kmer in barcodes:
                offset_counts[i] += 1
    return offset_counts.most_common(1)[0][0] if offset_counts else None

def main(args):
    console = Console(stderr=True, highlight=True)
    console.log("[bold red]Initializing heuristic barcode counting[/bold red]...")

    num_threads = cpu_count()

    # Reading FASTA File
    console.log("Reading barcodes...")
    barcodes = read_fasta(args.fasta_file)
    barcode_length = len(next(iter(barcodes)))

    is_paired_end = bool(args.fastq2)
    
    # Reading FASTQ Files
# Using the function:
    is_paired_end = bool(args.fastq2)
        
    # Reading FASTQ Files
    console.log("Sampling initial reads for orientation...")

    # Initialize
    need_swap = False 

    # Your existing logging and swapping logic
    console.log("Determining orientation of reads...")
    reads_sampled, need_swap, sample1, sample2 = determine_forward_read_with_generator(args.fastq1, args.fastq2, barcodes, is_paired_end)
    
    console.log(f"Sampled {reads_sampled} reads and found {len(sample1)} valid matches...")
    console.log(f"Swapping reads..." if need_swap else f"Proceeding with reads as is...")                

    # Apply the swap logic
    if need_swap:
        sample1, sample2 = sample2, sample1

    # Initialize to None
    barcode_start1 = None
    barcode_start2 = None        

    # Skip finding barcode starts for single-end that needed a swap
    if sample1 is not None:
        console.log("Finding forward coordinates...")
        barcode_start1 = find_start_positions(sample1, barcodes, barcode_length)
        if barcode_start1 is None:
            console.log("[bold red]No barcodes found in sample 1. Exiting.[/bold red]")
            sys.exit(1)

    # For paired-end or single-end that needed a swap
    if sample2:
        console.log("Finding reverse coordinates...")
        barcode_start2 = find_start_positions(sample2, barcodes, barcode_length, reverse_strand=True)
        if barcode_start2 is None:
            console.log("[bold red]No barcodes found in sample 2. Exiting.[/bold red]")
            sys.exit(1)
    else:
        barcode_start2 = None

    # Find flanking sequences
    if sample1 is not None:
        console.log("Identifying forward read junctions...")
        left1, right1 = find_ends(sample1, barcode_start1, barcode_length)
        barcode_start1 -= len(left1) if left1 else 0
    else:
        left1, right1 = None, None

    if sample2 is not None:
        console.log("Identifying reverse read junctions...")
        left2, right2 = find_ends(sample2, barcode_start2, barcode_length, reverse_strand=True)
        barcode_start2 -= len(left2) if left2 else 0
    else:
        left2, right2 = None, None

   # Calculate reverse complements
    left2_rc = left2[::-1].translate(str.maketrans("ATCGN", "TAGCN")) if left2 else None
    right2_rc = right2[::-1].translate(str.maketrans("ATCGN", "TAGCN")) if right2 else None

    # Check if the forward and reverse flanking sequences are reverse complements
    min_len_left = min_len_right = 0

    if left1 and right2_rc:
        min_len_left = min(len(left1), len(right2_rc))
        # console.log(f"Comparing {left1[-min_len_left:]} and {right2_rc[:min_len_left]}")
        if left1[-min_len_left:] != right2_rc[:min_len_left]:
            console.log("[bold red]Error: Forward and reverse left flanking sequences are not reverse complements.[/bold red]")
            sys.exit(1)
        
    if right1 and left2_rc:
        min_len_right = min(len(right1), len(left2_rc))
        # console.log(f"Comparing {right1[:min_len_right]} and {left2_rc[:min_len_right]}")
        if right1[:min_len_right] != left2_rc[:min_len_right]:
            console.log("[bold red]Error: Forward and reverse right flanking sequences are not reverse complements.[/bold red]")
            sys.exit(1)

    # Update barcodes
    console.log("Associating barcodes with read junctions...")
    if left1 and right1:
        barcodes = {left1 + barcode + right1 for barcode in barcodes}
    elif left2 and right2:
        barcodes = {left2 + barcode[::-1].translate(str.maketrans("ATCGN", "TAGCN")) + right2 for barcode in barcodes}
    barcode_length = len(next(iter(barcodes)))

    console.log("Executing high-throughput read analysis...")
    chunk_generator = read_in_chunks(args.fastq1, args.fastq2 if is_paired_end else None)

    # Create argument generator
    args_generator = ((chunk, barcodes, barcode_start1, barcode_start2, barcode_length, left1, right1, left2, right2, need_swap) for chunk in chunk_generator)

    # Execute multiprocessing
    with Pool(num_threads) as pool:
        results = pool.starmap(process_chunk, args_generator)
    # Collating Results
    console.log("[bold red]Collating results[/bold red]...")
    
    documented_barcodes = Counter()
    undocumented_barcodes = Counter()

    total_reads = 0

    for result, chunk_size in results:
        total_reads += chunk_size
        for barcode, count in result.items():
            if barcode.endswith('*'):
                undocumented_barcodes[barcode] += count
            else:
                documented_barcodes[barcode] += count

    if is_paired_end:
        fastq1_filename = os.path.basename(args.fastq1) if not need_swap else os.path.basename(args.fastq2)
        fastq2_filename = os.path.basename(args.fastq2) if not need_swap else os.path.basename(args.fastq1)
    else:
        fastq1_filename = os.path.basename(args.fastq1) if not need_swap else None
        fastq2_filename = None if not need_swap else os.path.basename(args.fastq1)

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

    # Rows for Input & Configuration
    if args.fasta_file:
        combined_table.add_row("Barcodes", f"[bold]{os.path.basename(args.fasta_file)}[/bold]")
    if fastq1_filename:
        combined_table.add_row("Forward Reads", f"[bold]{fastq1_filename}[/bold]")
    if fastq2_filename:
        combined_table.add_row("Reverse Reads", f"[bold]{fastq2_filename}[/bold]")
    combined_table.add_row("Threads", f"[bold]{num_threads}[/bold]")
    combined_table.add_row("Operating System", f"[bold]{platform.system()}[/bold]")

    # Heuristic Statistics Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_blue]Heuristics[/bold bright_blue]", "")

    if barcode_length:
        len_left = len(left1) if left1 else len(right2) if right2 else 0
        len_right = len(right1) if right1 else len(left2) if left2 else None
        combined_table.add_row("Barcode Length", f"[bold]{barcode_length - len_left - (len_right if len_right is not None else 0)}[/bold]")
    if barcode_start1:
        combined_table.add_row("Forward Offset", f"[bold]{barcode_start1 + len(left1)}[/bold]")
    if barcode_start2:
        combined_table.add_row("Reverse Offset", f"[bold]{barcode_start2 + len(right2)}[/bold]")
    if left1 and right1:
        combined_table.add_row("Forward Junction", f"[bold]{left1}...{right1}[/bold]")
    if left2 and right2:
        combined_table.add_row("Reverse Junction", f"[bold]{left2}...{right2}[/bold]")

    # Numeric Statistics Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_green]Barcode Alignment Stats[/bold bright_green]", "")

    # Rows for Numeric Statistics
    combined_table.add_row("Barcodes Declared", f"[bold]{len(barcodes)}[/bold]")
    combined_table.add_row("Documented Barcodes Found", f"[bold]{len(documented_barcodes)}[/bold]")
    combined_table.add_row("Undocumented Barcodes Found", f"[bold]{len(undocumented_barcodes)}[/bold]")
    combined_table.add_row("Total Reads", f"[bold]{total_reads}[/bold]")
    combined_table.add_row("Documented Barcode Reads", f"[bold]{sum(documented_barcodes.values())}[/bold]")
    combined_table.add_row("Undocumented Barcode Reads", f"[bold]{sum(undocumented_barcodes.values())}[/bold]")
    combined_table.add_row("Documented Fraction", f"[bold]{(sum(documented_barcodes.values()) / total_reads if total_reads != 0 else 0):.4f}[/bold]")
    combined_table.add_row("Undocumented Fraction", f"[bold]{(sum(undocumented_barcodes.values()) / total_reads if total_reads != 0 else 0):.4f}[/bold]", end_section=True)

    # Sequence Information Sub-heading  
    combined_table.add_section()

    # Add documented_barcodes to the main table
    combined_table.add_section()
    top_documented_barcodes = min(5, len(documented_barcodes))
    combined_table.add_row(f"[bold bright_yellow]Top {top_documented_barcodes} Documented Barcodes[/bold bright_yellow]", "")
    for idx, (barcode, count) in enumerate(documented_barcodes.most_common(top_documented_barcodes)):
        end_section = idx == (top_documented_barcodes - 1)
        combined_table.add_row(barcode, str(count), end_section=end_section)

    # Add undocumented_barcodes to the main table
    combined_table.add_section()
    top_undocumented_barcodes = min(5, len(undocumented_barcodes))
    combined_table.add_row(f"[bold bright_red]Top {top_undocumented_barcodes} Undocumented Barcodes[/bold bright_red]", "")
    for idx, (barcode, count) in enumerate(undocumented_barcodes.most_common(top_undocumented_barcodes)):
        end_section = idx == (top_undocumented_barcodes - 1)
        combined_table.add_row(barcode, str(count), end_section=end_section)

    # Print the combined table
    console.log(combined_table)

    for barcode, count in documented_barcodes.items():
        print("\t".join([barcode, str(count)]))

    for barcode, count in undocumented_barcodes.items():
        print("\t".join([barcode, str(count)]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Barcodes.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file.')
    parser.add_argument('fastq1', type=str, help='First FASTQ file.')
    parser.add_argument('fastq2', type=str, nargs='?', default=None, help='Second FASTQ file (optional).')
    args = parser.parse_args()
    main(args)