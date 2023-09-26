import argparse
import gzip
import os
import platform
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime
from multiprocessing import Pool, cpu_count

import rich
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

def fastq_reader(proc_stdout, should_decode):
    while True:
        next(proc_stdout, None)
        seq_line = next(proc_stdout, None)
        next(proc_stdout, None)
        next(proc_stdout, None)
        if seq_line is None:
            break
        yield seq_line.decode().strip() if should_decode else seq_line.strip()

def read_fastq(fastq1_file, fastq2_file, num_threads):
    if fastq2_file:

        paired_reads1, paired_reads2 = [], []
        proc1 = proc2 = None

        if fastq1_file.endswith('.gz'):
            try:
                proc1 = subprocess.Popen(["pigz", "-dc", fastq1_file], stdout=subprocess.PIPE)
                proc2 = subprocess.Popen(["pigz", "-dc", fastq2_file], stdout=subprocess.PIPE)
                f1, f2 = proc1.stdout, proc2.stdout
            except FileNotFoundError:
                f1 = gzip.open(fastq1_file, 'rt')
                f2 = gzip.open(fastq2_file, 'rt')
        else:
            f1 = open(fastq1_file, 'rt')
            f2 = open(fastq2_file, 'rt')

        should_decode = bool(proc1 and proc2)

        for seq1, seq2 in zip(fastq_reader(f1, should_decode), fastq_reader(f2, should_decode)):
            paired_reads1.append(seq1)
            paired_reads2.append(seq2)

        if proc1 and proc2:
            proc1.wait()
            proc2.wait()

        chunk_size = len(paired_reads1) // num_threads
        return [(paired_reads1[i:i+chunk_size], paired_reads2[i:i+chunk_size]) for i in range(0, len(paired_reads1), chunk_size)]
    
    else:
        reads1 = []
        if fastq1_file.endswith('.gz'):
            f1 = gzip.open(fastq1_file, 'rt')
        else:
            f1 = open(fastq1_file, 'rt')
        for seq1 in fastq_reader(f1, False):
            reads1.append(seq1)
        chunk_size = len(reads1) // num_threads
        return [(reads1[i:i+chunk_size], None) for i in range(0, len(reads1), chunk_size)]

def determine_forward_read(sample1, sample2, barcodes):
    count1, count2 = 0, 0
    barcode_length = len(next(iter(barcodes)))

    if sample2 is None:  # For single-end reads
        for seq in sample1[:250]:
            kmers = {seq[i:i + barcode_length] for i in range(len(seq) - barcode_length + 1)}
            if barcodes & kmers:
                count1 += 1
            
            # Check the reverse complement
            seq_rc = seq[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            kmers_rc = {seq_rc[i:i + barcode_length] for i in range(len(seq_rc) - barcode_length + 1)}
            if barcodes & kmers_rc:
                count2 += 1

        return count2 > count1  # True if more barcodes found in the reverse complement

    else:  # For paired-end reads
        for seq1, seq2 in zip(sample1[:250], sample2[:250]):
            kmers1 = {seq1[i:i + barcode_length] for i in range(len(seq1) - barcode_length + 1)}
            kmers2 = {seq2[i:i + barcode_length] for i in range(len(seq2) - barcode_length + 1)}

            if barcodes & kmers1:
                count1 += 1
            if barcodes & kmers2:
                count2 += 1

        return count2 > count1  # True if more barcodes found in sample2, implying need for swap

from collections import Counter


def find_ends(reads, start, length):
    lefts, rights = [], []

    for read in reads[:250]:
        left = read[start - 4: start]
        right = read[start + length: start + length + 4]
        lefts.append(left)
        rights.append(right)

    left = Counter(lefts).most_common(1)[0][0]
    right = Counter(rights).most_common(1)[0][0]
    return left, right

from collections import defaultdict


def process_chunk(chunk, barcodes, barcode_start1=None, barcode_start2=None, barcode_length=None):
    reads1, reads2 = chunk
    counts = defaultdict(int)
    unexpected_sequences = defaultdict(int)

    if reads1 and reads2:  # Paired-end processing
        for rec1, rec2 in zip(reads1, reads2):
            candidate1 = rec1[barcode_start1:barcode_start1 + barcode_length]
            candidate2 = rec2[barcode_start2:barcode_start2 + barcode_length]
            candidate2_rc = candidate2[::-1].translate(str.maketrans("ATCGN", "TAGCN"))

            if candidate1 in barcodes and candidate2_rc == candidate1:
                counts[candidate1[4:-4]] += 1
            elif candidate1 not in barcodes and candidate2_rc == candidate1:
                unexpected_sequences[candidate1] += 1

    elif reads1:  # Single-end processing, forward orientation
        for rec1 in reads1:
            candidate1 = rec1[barcode_start1:barcode_start1 + barcode_length]
            if candidate1 in barcodes:
                counts[candidate1[4:-4]] += 1
            else:
                unexpected_sequences[candidate1] += 1

    elif reads2:  # Single-end processing, reverse complement
        for rec2 in reads2:
            candidate2 = rec2[barcode_start2:barcode_start2 + barcode_length]
            candidate2_rc = candidate2[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            if candidate2 in barcodes:
                counts[candidate2_rc[4:-4]] += 1
            else:
                unexpected_sequences[candidate2_rc] += 1

    return counts, unexpected_sequences

def find_start_positions(reads, barcodes, barcode_length, is_read2=False):
    if reads is None:
        return None
    
    offset_counts = Counter()
    for read in reads[:250]:

        for i in range(len(read) - barcode_length + 1):
            kmer = read[i:i+barcode_length]
            if is_read2:
                kmer = str(Seq(kmer).reverse_complement())
            if kmer in barcodes:
                offset_counts[i] += 1
    return offset_counts.most_common(1)[0][0] if offset_counts else None

def main(args):
    console = Console(stderr=True, highlight=True)
    console.log("[bold red]Initializing heuristic barcode counting")

    num_threads = cpu_count()

    # Reading FASTA File
    console.log("Reading barcodes...")
    barcodes = read_fasta(args.fasta_file)
    barcode_length = len(next(iter(barcodes)))

    is_paired_end = bool(args.fastq2)
    
    # Reading FASTQ Files
    console.log("Reading results from FASTQ...")
    chunks = read_fastq(args.fastq1, args.fastq2 if is_paired_end else None, num_threads)
    sample1, sample2 = chunks[0]

    # Initialize
    need_swap = False  

    # Determine if a swap is needed
    if is_paired_end:
        console.log("Finding orientation of paired-end reads...")
        need_swap = determine_forward_read(sample1, sample2, barcodes)
    else:
        console.log("Determining orientation of single-end reads...")
        need_swap = determine_forward_read(sample1, None, barcodes)

    # Apply the swap logic
    if need_swap:
        if is_paired_end:
            chunks = [(reads2, reads1) for reads1, reads2 in chunks]
        else:
            # For single-end reads, make them the new 'sample2' by leaving 'sample1' as None
            chunks = [(None, reads1) for reads1, _ in chunks]

    # Always get sample1 and sample2 from the first chunk
    sample1, sample2 = chunks[0]

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
        console.log("Finding reverse coordinates")
        barcode_start2 = find_start_positions(sample2, barcodes, barcode_length, is_read2=True)
        if barcode_start2 is None:
            console.log("[bold red]No barcodes found in sample 2. Exiting.[/bold red]")
            sys.exit(1)
    else:
        barcode_start2 = None

    
    # Find flanking sequences
    if sample1 is not None:
        console.log("Finding forward junctions...")
        left1, right1 = find_ends(sample1, barcode_start1, barcode_length)
        barcode_start1 -= len(left1)
    else:
        left1, right1 = None, None

    if sample2:
        console.log("Finding reverse junctions...")
        left2, right2 = find_ends(sample2, barcode_start2, barcode_length)
        barcode_start2 -= len(left2)
    else:
        left2, right2 = None, None

    # Update barcodes
    console.log("Mapping junctions to barcodes...")
    if left1 and right1:
        barcodes = {left1 + barcode + right1 for barcode in barcodes}
    elif left2 and right2:
        barcodes = {left2 + barcode[::-1].translate(str.maketrans("ATCGN", "TAGCN")) + right2 for barcode in barcodes}
    barcode_length = len(next(iter(barcodes)))

    # Processing Chunks
    console.log("Counting barcodes at junctions...")
    with Pool(num_threads) as pool:
        results = pool.starmap(process_chunk, [(chunk, barcodes, barcode_start1, barcode_start2, barcode_length) for chunk in chunks])

    # Collating Results
    console.log("Preparing results...")
    final_counts = Counter()
    final_unexpected_sequences = Counter()
    for counts, unexpected_sequences in results:
        final_counts.update(counts)
        final_unexpected_sequences.update(unexpected_sequences)

    # Filtering and Creating New Barcodes
    console.log("Heuristic filtering and identifying new barcodes...")
    new_barcodes = Counter()
    left, right = (left1, right1) if left1 is not None else (left2, right2)
    if left1 is None and left2 is not None:
        left, right = right[::-1].translate(str.maketrans("ATCGN", "TAGCN")), left[::-1].translate(str.maketrans("ATCGN", "TAGCN"))


    # Initialize new_barcodes
    new_barcodes = Counter()

    if left and right:  # Check if either pair is not None
        for seq, count in final_unexpected_sequences.items():
            if seq.startswith(left) and seq.endswith(right):
                barcode = seq[4:-4] + "*"                 
                new_barcodes[barcode] = count

    # final_counts = final_counts + new_barcodes

    total_reads = sum(len(chunk[0]) if chunk[0] is not None else len(chunk[1]) for chunk in chunks)

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
        combined_table.add_row("Barcode Length", f"[bold]{barcode_length-8}[/bold]")
    if barcode_start1:
        combined_table.add_row("Forward Offset", f"[bold]{barcode_start1 + 4}[/bold]")
    if barcode_start2:
        combined_table.add_row("Reverse Offset", f"[bold]{barcode_start2 + 4}[/bold]")
    if left1 and right1:
        combined_table.add_row("Forward Junction", f"[bold]{left1}...{right1}[/bold]")
    if left2 and right2:
        combined_table.add_row("Reverse Junction", f"[bold]{left2}...{right2}[/bold]")

    # Numeric Statistics Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_green]Barcode Alignment Stats[/bold bright_green]", "")

    # Rows for Numeric Statistics
    combined_table.add_row("Barcodes Declared", f"[bold]{len(barcodes)}[/bold]")
    combined_table.add_row("Documented Barcodes Found", f"[bold]{len(final_counts)}[/bold]")
    combined_table.add_row("Undocumented Barcodes Found", f"[bold]{len(new_barcodes)}[/bold]")
    combined_table.add_row("Total Reads", f"[bold]{total_reads}[/bold]")
    combined_table.add_row("Documented Barcode Reads", f"[bold]{sum(final_counts.values())}[/bold]")
    combined_table.add_row("Undocumented Barcode Reads", f"[bold]{sum(new_barcodes.values())}[/bold]")
    combined_table.add_row("Documented Fraction", f"[bold]{(sum(final_counts.values()) / total_reads if total_reads != 0 else 0):.4f}[/bold]")
    combined_table.add_row("Undocumented Fraction", f"[bold]{(sum(new_barcodes.values()) / total_reads if total_reads != 0 else 0):.4f}[/bold]", end_section=True)

    # Sequence Information Sub-heading  
    combined_table.add_section()

    # Add final_counts to the main table
    combined_table.add_section()
    top_final_counts = min(5, len(final_counts))
    combined_table.add_row(f"[bold bright_yellow]Top {top_final_counts} Documented Barcodes[/bold bright_yellow]", "")
    for idx, (barcode, count) in enumerate(final_counts.most_common(top_final_counts)):
        end_section = idx == (top_final_counts - 1)
        combined_table.add_row(barcode, str(count), end_section=end_section)

    # Add new_barcodes to the main table
    combined_table.add_section()
    top_new_barcodes = min(5, len(new_barcodes))
    combined_table.add_row(f"[bold bright_red]Top {top_new_barcodes} Undocumented Barcodes[/bold bright_red]", "")
    for idx, (barcode, count) in enumerate(new_barcodes.most_common(top_new_barcodes)):
        end_section = idx == (top_new_barcodes - 1)
        combined_table.add_row(barcode, str(count), end_section=end_section)

    # Print the combined table
    console.log(combined_table)

    for barcode, count in final_counts.items():
        print("\t".join([barcode, str(count)]))

    for barcode, count in new_barcodes.items():
        print("\t".join([barcode, str(count)]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process Barcodes.')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file.')
    parser.add_argument('fastq1', type=str, help='First FASTQ file.')
    parser.add_argument('fastq2', type=str, nargs='?', default=None, help='Second FASTQ file (optional).')
    args = parser.parse_args()
    main(args)

