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

def count_barcodes(seq, barcode_length, barcodes):
    # Generate k-mers and count matches with barcodes
    kmers = {seq[i:i + barcode_length] for i in range(len(seq) - barcode_length + 1)}
    return len(barcodes & kmers)

def determine_forward_read(sample1, sample2, barcodes):
    barcode_length = len(next(iter(barcodes)))
    count1, count2 = 0, 0

    for i in range(min(250, len(sample1))):
        seq1 = sample1[i]
        seq2 = sample2[i] if sample2 else None  # seq2 is None for single-end

        # Count barcode matches in seq1
        count1 += count_barcodes(seq1, barcode_length, barcodes)
        
        # Early exit if count1 is overwhelmingly large
        if count1 > 125:
            return False
        
        # For paired-end or reverse complement of single-end
        if seq2:
            # Count barcode matches in seq2
            count2 += count_barcodes(seq2, barcode_length, barcodes)
            
            # Early exit if count2 is overwhelmingly large
            if count2 > 125:
                return True
        else:
            # Count barcode matches in reverse complement of seq1
            seq1_rc = seq1[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            count2 += count_barcodes(seq1_rc, barcode_length, barcodes)
            
            # Early exit if count2 is overwhelmingly large
            if count2 > 125:
                return True

    # True if more barcodes found in seq2 or reverse complement of seq1
    return count2 > count1


def find_ends(reads, start, length, reverse_strand=False):
    # Initialize max values; reverse if needed
    max_left, max_right = (4, 4) if not reverse_strand else (4, 4)[::-1]

    # Initialize Counters for left and right flanking regions
    lefts, rights = Counter(), Counter()

    for read in reads[:250]:
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

def process_chunk(chunk, barcodes, barcode_start1=None, barcode_start2=None, barcode_length=None, left1=None, right1=None, left2=None, right2=None):
    reads1, reads2 = chunk
    counts = defaultdict(int)

    if reads1 and reads2:  # Paired-end processing
        for rec1, rec2 in zip(reads1, reads2):
            candidate1 = rec1[barcode_start1:barcode_start1 + barcode_length]
            candidate2 = rec2[barcode_start2:barcode_start2 + barcode_length]
            candidate2_rc = candidate2[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            if candidate2_rc == candidate1:
                if candidate1 in barcodes:
                    counts[candidate1[len(left1):-len(right1) or None]] += 1
                else:
                    if candidate1.startswith(left1) and candidate1.endswith(right1):
                        barcode = candidate1[len(left1):-len(right1) or None] + "*"  
                        counts[barcode] += 1

    elif reads1:  # Single-end processing, forward orientation
        for rec1 in reads1:
            candidate1 = rec1[barcode_start1:barcode_start1 + barcode_length]
            if candidate1 in barcodes:
                counts[candidate1[len(left1):-len(right1) or None]] += 1
            else:
                if candidate1.startswith(left1) and candidate1.endswith(right1):
                    barcode = candidate1[len(left1):-len(right1) or None] + "*"  
                    counts[barcode] += 1

    elif reads2:  # Single-end processing, reverse complement
        for rec2 in reads2:
            candidate2 = rec2[barcode_start2:barcode_start2 + barcode_length]
            candidate2_rc = candidate2[::-1].translate(str.maketrans("ATCGN", "TAGCN"))
            if candidate2 in barcodes:
                counts[candidate2_rc[len(right2):-len(left2) or None]] += 1
            else:
                # print(candidate2, left2, right2, file=sys.stderr)
                if candidate2.startswith(left2) and candidate2.endswith(right2):
                    barcode = candidate2_rc[len(right2):-len(left2) or None] + "*"  
                    counts[barcode] += 1
    
    return counts


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
        barcode_start1 -= len(left1) if left1 else 0
    else:
        left1, right1 = None, None

    if sample2 is not None:
        console.log("Finding reverse junctions...")
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
        console.log(f"Comparing {left1[-min_len_left:]} and {right2_rc[:min_len_left]}")
        if left1[-min_len_left:] != right2_rc[:min_len_left]:

            console.log("[bold red]Error: Forward and reverse left flanking sequences are not reverse complements.[/bold red]")
            sys.exit(1)
        

    if right1 and left2_rc:
        min_len_right = min(len(right1), len(left2_rc))
        console.log(f"Comparing {right1[:min_len_right]} and {left2_rc[:min_len_right]}")
        if right1[:min_len_right] != left2_rc[:min_len_right]:

            console.log("[bold red]Error: Forward and reverse right flanking sequences are not reverse complements.[/bold red]")
            sys.exit(1)

        min_junction_len = min(min_len_left, min_len_right) 
        left1, right1 = left1[-min_junction_len:], right1[:min_junction_len] if left1 and right1 else None
        left2, right2 = left2[-min_junction_len:], right2[:min_junction_len] if left2 and right2 else None

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
        results = pool.starmap(process_chunk, [(chunk, barcodes, barcode_start1, barcode_start2, barcode_length, left1, right1, left2, right2) for chunk in chunks])


    # Collating Results
    console.log("Preparing results...")
    
    documented_barcodes = Counter()
    undocumented_barcodes = Counter()

    for result in results:
        for barcode, count in result.items():
            if barcode.endswith('*'):
                undocumented_barcodes[barcode] += count
            else:
                documented_barcodes[barcode] += count


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

