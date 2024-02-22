import argparse
import gzip
import os
import platform
import sys
from collections import Counter, defaultdict
from contextlib import nullcontext
from datetime import datetime
from multiprocessing import Pool, cpu_count
from typing import List, Set, Tuple

import rich
import zstandard as zstd
from Bio.Seq import Seq
from rich.console import Console
from rich.table import Table
import rich.table
from Bio import SeqIO


def rev_comp(sequence: str) -> str:
    return sequence[::-1].translate(str.maketrans("ATCGN", "TAGCN"))


def safe_len(s):
    return 0 if s is None else len(s)


def generate_kmers(seq: str, k: int) -> Set[str]:
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def read_fasta(fasta_file) -> Set[str]:
    barcodes = set()
    open_func = gzip.open if fasta_file.endswith(".gz") else open
    with open_func(fasta_file, "rt") as f:
        for line in f:
            if not line.startswith(">"):
                barcodes.add(line.strip())
    return barcodes


def validate_barcodes(barcodes):
    # If the input is a string, treat it as a filename
    if isinstance(barcodes, str):
        sequences = list(SeqIO.parse(barcodes, "fasta"))
    # If the input is a list, treat it as a list of sequences
    elif isinstance(barcodes, list):
        sequences = barcodes
    else:
        raise ValueError("The input should be a filename or a list of barcodes.")

    # Check if there are at least 10 sequences
    if len(sequences) < 10:
        raise ValueError(
            "The input contains fewer than 10 sequences. Please provide at least 10 short barcodes."
        )

    # Check if the sequences are short (e.g., no more than 1,000 bases)
    for seq in sequences:
        if len(seq) > 1000:
            raise ValueError(
                f'The sequence "{seq.id}" is longer than 1,000 bases. Provide a list or fasta file of short barcodes.'
            )


def open_reads_file(file_path, mode):
    if file_path.endswith(".gz"):
        return gzip.open(file_path, mode)
    elif file_path.endswith(".zst"):
        return zstd.open(file_path, mode)
    elif file_path.endswith(".fastq") or file_path.endswith(".fq"):
        return open(file_path, mode)
    elif file_path.endswith(".reads"):
        return open(file_path, mode)
    else:
        raise ValueError(
            f'"{file_path}" does not appear to be a supported reads file: .fastq or .reads.'
        )


def read_in_chunks(file1, file2=None, chunk_size=2**16) -> Tuple[List[str], List[str]]:
    reads1, reads2 = [], []

    # Remove compression extension if present
    stripped_file1 = file1
    if file1.endswith(".gz") or file1.endswith(".zst"):
        stripped_file1 = os.path.splitext(file1)[0]

    # Determine the file type based on the stripped extension
    if stripped_file1.endswith(".fastq"):
        file_type = "fastq"
    elif stripped_file1.endswith(".reads"):
        file_type = "reads"
    else:
        raise ValueError("Unsupported file type. Must be '.fastq' or '.reads'.")

    with open_reads_file(file1, "rt") as f1, (
        open_reads_file(file2, "rt") if file2 else nullcontext()
    ) as f2:
        iters = [iter(f1), iter(f2) if f2 else iter([])]

        while True:
            try:
                if file_type == "fastq":
                    next(iters[0])  # Skip @SEQUENCE_ID
                    reads1.append(next(iters[0]).strip())  # Capture SEQUENCE
                    next(iters[0])  # Skip +
                    next(iters[0])  # Skip QUALITY

                    if file2:
                        next(iters[1])  # Skip @SEQUENCE_ID
                        reads2.append(next(iters[1]).strip())  # Capture SEQUENCE
                        next(iters[1])  # Skip +
                        next(iters[1])  # Skip QUALITY

                elif file_type == "reads":
                    reads1.append(next(iters[0]).strip())  # Capture SEQUENCE

                    if file2:
                        reads2.append(next(iters[1]).strip())  # Capture SEQUENCE

                if len(reads1) >= chunk_size:
                    yield (reads1[:chunk_size], reads2[:chunk_size] if reads2 else None)
                    reads1, reads2 = reads1[chunk_size:], (
                        reads2[chunk_size:] if reads2 else []
                    )

            except StopIteration:
                break

        if reads1:
            yield (reads1, reads2 if reads2 else None)


def sample_data(file1, file2, barcodes, is_paired):

    rev_barcodes = set(rev_comp(bc) for bc in barcodes)
    bc_len = len(next(iter(barcodes)))
    chunk_generator = read_in_chunks(
        file1, file2 if is_paired else None, chunk_size=10 * len(barcodes)
    )
    processed_count1 = 0
    processed_count2 = 0
    valid_count1 = 0
    valid_count2 = 0

    # Overall Counters
    global_read1_orientations = Counter()
    global_read1_offsets = Counter()
    global_read2_orientations = Counter()
    global_read2_offsets = Counter()

    # Sets to determine sampling depth
    valid_reads1 = set()
    valid_reads2 = set()

    seen_reads = set()

    # Keep track of the unique barcodes seen in any orientation
    unique_barcodes = set()

    for read1_chunk, read2_chunk in chunk_generator:

        # Lists to collect orientations and offsets for this chunk
        read1_orientations = []
        read1_offsets = []

        read2_orientations = []
        read2_offsets = []
        valid_kmers = set()

        for read1, read2 in zip(
            read1_chunk, read2_chunk if read2_chunk else [None] * len(read1_chunk)
        ):

            # Bypass pairs of reads if the read1 or read2 has been seen before
            if read1 in seen_reads or (read2 and read2 in seen_reads):
                continue

            seen_reads.add(read1)
            if read2:
                seen_reads.add(read2)

            processed_count1 += 1
            if read2:
                processed_count2 += 1

            for i in range(len(read1) - bc_len + 1):
                kmer = read1[i : i + bc_len]
                if kmer in barcodes:

                    if kmer in valid_kmers:
                        processed_count1 -= 1
                        continue
                    valid_count1 += 1
                    valid_kmers.add(kmer)
                    unique_barcodes.add(kmer)
                    read1_orientations.append("forward")
                    read1_offsets.append(i)
                    valid_reads1.add(read1)

                elif kmer in rev_barcodes:

                    if kmer in valid_kmers:
                        processed_count1 -= 1
                        continue
                    valid_count1 += 1
                    valid_kmers.add(kmer)
                    unique_barcodes.add(kmer)
                    read1_orientations.append("reverse")
                    read1_offsets.append(i)
                    valid_reads1.add(read1)

                if is_paired:
                    kmer2 = read2[i : i + bc_len]
                    if kmer2 in barcodes:
                        if kmer2 in valid_kmers:
                            processed_count2 -= 1
                            continue
                        valid_count2 += 1
                        valid_kmers.add(kmer2)
                        unique_barcodes.add(kmer)
                        read2_orientations.append("forward")
                        read2_offsets.append(i)
                        valid_reads2.add(read2)

                    elif kmer2 in rev_barcodes:
                        if kmer2 in valid_kmers:
                            processed_count2 -= 1
                            continue
                        valid_count2 += 1
                        valid_kmers.add(kmer2)
                        unique_barcodes.add(kmer)
                        read2_orientations.append("reverse")
                        read2_offsets.append(i)
                        valid_reads2.add(read2)

        # print(len(seen_reads), len(valid_kmers), len(valid_reads1), len(valid_reads2), valid_count1, valid_count2, file=sys.stderr)

        # Update the global counters after processing each chunk
        global_read1_orientations.update(read1_orientations)
        global_read1_offsets.update(read1_offsets)

        global_read2_orientations.update(read2_orientations)
        global_read2_offsets.update(read2_offsets)

        read1_offsets_common = Counter(global_read1_offsets).most_common(2)
        read2_offsets_common = Counter(global_read2_offsets).most_common(2)

        # Check conditions for two reads
        # if ((read2 and valid_count1 >= 10 * len(barcodes) and valid_count2 >= 10 * len(barcodes)) or (len(seen_reads) >= 10 * len(valid_kmers) and len(valid_kmers) != 0)):
        if (
            len(valid_reads1) >= len(barcodes)
            and len(valid_reads2) >= len(barcodes)
            and (
                valid_count1 >= 10 * len(barcodes)
                and valid_count2 >= 10 * len(barcodes)
            )
            or (len(seen_reads) >= 10 * len(valid_kmers) and len(valid_kmers) != 0)
        ):
            # Break if there's a single dominant offset or the most common is 5 times more frequent than the second
            # print("TRUE", file=sys.stderr)
            if (len(read1_offsets_common) == 1 and len(read2_offsets_common) == 1) or (
                len(read1_offsets_common) > 1
                and len(read2_offsets_common) > 1
                and read1_offsets_common[0][1] >= 2 * read1_offsets_common[1][1]
                and read2_offsets_common[0][1] >= 2 * read2_offsets_common[1][1]
            ):
                break

        # Check conditions for a single read
        # elif (not read2 and valid_count1 >= 250 * len(barcodes) or (len(seen_reads) >= 10 * len(valid_kmers) and len(valid_kmers) != 0)):
        elif (not read2 and len(valid_reads1) >= len(barcodes)) and (
            valid_count1 >= 10 * len(barcodes)
            or (len(seen_reads) >= 10 * len(valid_kmers) and len(valid_kmers) != 0)
        ):
            # Break if there's a single dominant offset or the most common is 5 times more frequent than the second
            if len(read1_offsets_common) == 1 or (
                len(read1_offsets_common) > 1
                and read1_offsets_common[0][1] >= 2 * read1_offsets_common[1][1]
            ):
                break

    read1_orientation = (
        Counter(global_read1_orientations).most_common(1)[0][0] if read1 else None
    )
    read1_offset = Counter(global_read1_offsets).most_common(1)[0][0] if read1 else None

    read2_orientation = (
        Counter(global_read2_orientations).most_common(1)[0][0] if read2 else None
    )
    read2_offset = Counter(global_read2_offsets).most_common(1)[0][0] if read2 else None

    if read1_orientation == "forward" or read2_orientation == "reverse":
        need_swap = False
        return (
            processed_count1 + processed_count2,
            read1_offset,
            read2_offset,
            valid_reads1,
            valid_reads2,
            unique_barcodes,
            need_swap,
        )

    elif read1_orientation == "reverse" or read2_orientation == "forward":
        need_swap = True
        return (
            processed_count1 + processed_count2,
            read2_offset,
            read1_offset,
            valid_reads2,
            valid_reads1,
            unique_barcodes,
            need_swap,
        )

    else:
        raise ValueError(
            "Unable to determine orientation of reads. Please check the input files."
        )


def find_flanks(
    reads: List[str], start: int, bc_len: int, max_flank: int = 10
) -> Tuple[str, str]:
    L_flanks, R_flanks = Counter(), Counter()

    def update_flanks(side: str, seq: str, max_len: int):
        """Update counters with sub-sequences of varying lengths."""
        counts = L_flanks if side == "L_flank" else R_flanks
        for i in range(max_len, 0, -1):
            truncated = seq[-i:] if side == "L_flank" else seq[:i]
            counts[truncated] += 1

    for read in reads:
        # Extract potential flanking sequences from the read based on start and bc_len
        L_flank = (
            read[start - max_flank : start] if start - max_flank >= 0 else read[0:start]
        )
        R_flank = read[start + bc_len : start + bc_len + max_flank]

        # Update flank counters with sequences of varying lengths
        update_flanks("L_flank", L_flank, len(L_flank))
        update_flanks("R_flank", R_flank, len(R_flank))

    def extract_best_flank(counts: Counter) -> str:
        most_common_prev = None
        for fl_len in range(max_flank, 0, -1):
            # print(f"Potential flank sequences: {counts}", file=sys.stderr)

            potential_seqs = [seq for seq in counts if len(seq) == fl_len]
            if not potential_seqs:
                continue
            most_common = max(potential_seqs, key=lambda x: counts[x])
            if fl_len == max_flank:
                if (
                    most_common_prev is None
                    or counts[most_common] > 3 * counts[most_common_prev]
                ):
                    return most_common
            else:
                if (
                    most_common_prev is not None
                    and counts[most_common] > 3 * counts[most_common_prev]
                ):
                    return most_common
                if most_common_prev is None or (
                    counts[most_common] * 3 < counts[most_common_prev]
                ):
                    most_common_prev = most_common
        return None

    L_most_common = extract_best_flank(L_flanks)
    R_most_common = extract_best_flank(R_flanks)
    return L_most_common, R_most_common


def process_chunk(
    chunk,
    bcs_with_flanks_fwd,
    bcs_with_flanks_rev,
    L_fwd_start,
    L_rev_start,
    bc_len,
    L_fwd,
    R_fwd,
    L_rev,
    R_rev,
    need_swap,
) -> Tuple[Counter, int]:

    counts = Counter()

    if need_swap:
        reads2, reads1 = chunk
    else:
        reads1, reads2 = chunk

    def validate_read(seq_with_flanks, L_flank, R_flank, rev=False):
        if rev:
            in_bcs_with_flanks = seq_with_flanks in bcs_with_flanks_rev
        else:
            in_bcs_with_flanks = seq_with_flanks in bcs_with_flanks_fwd

        seq = seq_with_flanks[
            safe_len(L_flank) : safe_len(seq_with_flanks) - safe_len(R_flank)
        ]
        has_flanks = seq_with_flanks.startswith(
            L_flank or ""
        ) and seq_with_flanks.endswith(R_flank or "")

        return in_bcs_with_flanks, has_flanks, seq

    L_fwd_len = safe_len(L_fwd)
    R_fwd_len = safe_len(R_fwd)
    L_rev_len = safe_len(L_rev)
    R_rev_len = safe_len(R_rev)

    def process_paired_end(
        reads1,
        reads2,
        L_fwd_start,
        L_fwd_len,
        R_fwd_len,
        L_rev_start,
        L_rev_len,
        R_rev_len,
        L_fwd,
        R_fwd,
        L_rev,
        R_rev,
    ):
        if len(reads1) != len(reads2):
            raise ValueError(
                "Length of reads1 and reads2 must be the same for paired-end data."
            )

        for record_fwd, record_rev in zip(reads1, reads2):
            if "N" in record_fwd:
                continue

            seq_with_flanks_fwd = record_fwd[
                L_fwd_start : L_fwd_start + L_fwd_len + bc_len + R_fwd_len
            ]
            seq_with_flanks_rev = record_rev[
                L_rev_start : L_rev_start + L_rev_len + bc_len + R_rev_len
            ]

            in_bcs_with_flanks_fwd, has_flanks_fwd, seq1 = validate_read(
                seq_with_flanks_fwd, L_fwd, R_fwd
            )
            in_bcs_with_flanks_rev, has_flanks_rev, seq2 = validate_read(
                seq_with_flanks_rev, L_rev, R_rev, rev=True
            )

            if seq1 != rev_comp(seq2):
                continue

            if (
                in_bcs_with_flanks_fwd
                and in_bcs_with_flanks_rev
                and has_flanks_fwd
                and has_flanks_rev
            ):
                counts[seq1] += 1

            elif has_flanks_fwd and has_flanks_rev:
                counts[seq1 + "*"] += 1

    def process_single_end(
        reads, start_pos, len_L, len_R, L_flank, R_flank, reverse=False
    ):
        for record in reads:
            if "N" in record:
                continue

            seq_with_flanks = record[start_pos : start_pos + len_L + bc_len + len_R]
            in_bcs_with_flanks, has_flanks, seq = validate_read(
                seq_with_flanks, L_flank, R_flank, rev=reverse
            )

            if reverse:
                seq = rev_comp(seq)

            if in_bcs_with_flanks and has_flanks:
                counts[seq] += 1
            elif has_flanks:
                counts[seq + "*"] += 1

    if reads1 and reads2:  # Paired-end processing
        process_paired_end(
            reads1,
            reads2,
            L_fwd_start,
            L_fwd_len,
            R_fwd_len,
            L_rev_start,
            L_rev_len,
            R_rev_len,
            L_fwd,
            R_fwd,
            L_rev,
            R_rev,
        )
    elif reads1:  # Single-end processing, fwd orientation
        process_single_end(reads1, L_fwd_start, L_fwd_len, R_fwd_len, L_fwd, R_fwd)
    elif reads2:  # Single-end processing, reverse complement
        process_single_end(
            reads2, L_rev_start, L_rev_len, R_rev_len, L_rev, R_rev, reverse=True
        )

    return (counts, len(reads1) if reads1 else len(reads2))


def main(args):
    from rich.style import Style

    console = Console(stderr=True, highlight=True, style=Style(color="white"))
    console_error = Console(
        stderr=True, highlight=True, style=Style(color="red", bold=True)
    )
    console.log("Initializing heuristic barcode counting...")

    reads1 = args.file1
    reads2 = args.file2

    num_threads = cpu_count() // 2

    # Reading FASTA File
    console.log(
        "Reading barcodes...",
    )

    try:
        validate_barcodes(args.fasta_file)
        barcodes = read_fasta(args.fasta_file)
        bcs_rev = {rev_comp(barcode) for barcode in barcodes}
        bc_len = len(next(iter(barcodes)))
        is_paired = bool(args.file2)

        # Reading reads Files
        console.log("Sampling initial reads for orientation...")

        (
            new_reads_sampled,
            bc_start1,
            bc_start2,
            sample1,
            sample2,
            unique_barcodes,
            need_swap,
        ) = sample_data(reads1, reads2, barcodes, is_paired)
        # console.log(f"Sampled {new_reads_sampled:,} unique reads. Found {processed_count1:,} valid matches, including {safe_len(sample1):,} unique forward matches and {safe_len(sample2):,} unique reverse matches...")
        console.log(
            f"Sampled {new_reads_sampled:,} reads and found {safe_len(unique_barcodes):,} barcodes in {safe_len(sample1):,} forward and {safe_len(sample2):,} reverse matches..."
        )

        if need_swap:
            console.log("Swapping orientation...")

        console.log("Identifying flanking sequences...")

        # Find flanking sequences
        if sample1 is not None:
            L_fwd, R_fwd = find_flanks(sample1, bc_start1, bc_len)
            L_fwd_start = bc_start1 - len(L_fwd) if L_fwd else 0
        else:
            L_fwd, R_fwd = None, None
            L_fwd_start = None

        if sample2 is not None:
            L_rev, R_rev = find_flanks(sample2, bc_start2, bc_len)
            L_rev_start = bc_start2 - len(L_rev) if L_rev else 0
        else:
            L_rev, R_rev = None, None
            L_rev_start = None

        # Calculate reverse complements
        L_rev_rev = rev_comp(L_rev) if L_rev else None
        R_rev_rev = rev_comp(R_rev) if R_rev else None

        # Check if the fwd and reverse flanking sequences are reverse complements
        L_min_len = R_min_len = 0

        error_messages = set()

        if L_fwd and R_rev_rev:
            L_min_len = min(len(L_fwd), len(R_rev_rev))
            if L_fwd[-L_min_len:] != R_rev_rev[:L_min_len]:
                error_messages.add("Flank complementarity violation")

        if R_fwd and L_rev_rev:
            R_min_len = min(len(R_fwd), len(L_rev_rev))
            if R_fwd[:R_min_len] != L_rev_rev[:R_min_len]:
                error_messages.add("Flank complementarity violation")

        if error_messages:
            console.log(
                "\n[bold yellow]Warning[/ bold yellow]: Expected complemtarity not found for flanking sequences..."
            )

            if L_fwd and R_rev_rev:
                console.log(
                    f"[bold]Left flank[/bold]:  '{L_fwd}' (./{reads1}) != '{R_rev_rev}' (./{reads2} revcomp)"
                )
            if R_fwd and L_rev_rev:
                console.log(
                    f"[bold]Right flank[/bold]: '{R_fwd}' (./{reads1}) != '{L_rev_rev}' (./{reads2} revcomp)"
                )

            console.log(

"""
In the analysis, guides are presumed to be in the same position across all reads. 
If this assumption disrupts your pipeline, we welcome your feedback on our GitHub page: 

[bold]https://github.com/ryandward/barcoder/issues[/bold]

You might want to consider running the analysis with a single set of reads. 

Please note, if flanks could not be identified, the localization of barcode positions may be unsuccessful with just one read. 

Without the second read, this error cannot be triggered, so please proceed with caution.

However, if one read is of superior quality, it might be advantageous to proceed with that read alone.

Start by rerunning the analysis on the read that contains the identifiable flank.
"""
            )
            raise ValueError("A critical error occurred: " + ", ".join(list(error_messages)))

        def add_flank(barcodes, L_flank=None, R_flank=None):
            L_flank, R_flank = (L_flank or ""), (R_flank or "")
            return {L_flank + barcode + R_flank for barcode in barcodes}

        # Add flanks to the ends of the barcodes, if present

        bcs_with_flanks_fwd = add_flank(barcodes, L_fwd, R_fwd)
        bcs_with_flanks_rev = add_flank(bcs_rev, L_rev, R_rev)

        console.log("Executing high-throughput read analysis...")
        chunk_generator = read_in_chunks(args.file1, args.file2 if is_paired else None)

        # Create argument generator
        args_generator = (
            (
                chunk,
                bcs_with_flanks_fwd,
                bcs_with_flanks_rev,
                L_fwd_start,
                L_rev_start,
                bc_len,
                L_fwd,
                R_fwd,
                L_rev,
                R_rev,
                need_swap,
            )
            for chunk in chunk_generator
        )

        # Execute multiprocessing
        with Pool(num_threads) as pool:
            results = pool.starmap(process_chunk, args_generator)
        # Collating Results
        console.log("Finishing up and collating results!")

        doc_bcs = Counter()
        undoc_bcs = Counter()

        total_reads = 0

        for result, chunk_size in results:
            total_reads += chunk_size
            for barcode, count in result.items():
                if barcode.endswith("*"):
                    undoc_bcs[barcode] += count
                else:
                    doc_bcs[barcode] += count

        if is_paired:
            file1_filename = (
                os.path.basename(args.file1)
                if not need_swap
                else os.path.basename(args.file2)
            )
            file2_filename = (
                os.path.basename(args.file2)
                if not need_swap
                else os.path.basename(args.file1)
            )
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
            caption_style="white",
            header_style="bold bright_white",
            border_style="bold bright_white",
            highlight=True,
            show_header=True,
        )
        # Define columns with justifications
        combined_table.add_column(
            os.path.basename(sys.argv[0]), justify="right", style="white", min_width=30
        )
        combined_table.add_column("Summary", justify="right", min_width=20)

        # Input & Configuration Sub-heading
        combined_table.add_section()
        combined_table.add_row(
            "[bold bright_magenta]Input & Config[/bold bright_magenta]", ""
        )

        if args.fasta_file:
            combined_table.add_row(
                "Barcodes", f"[bold]{os.path.basename(args.fasta_file)}[/bold]"
            )
        if file1_filename:
            combined_table.add_row("Forward Reads", f"[bold]{file1_filename}[/bold]")
        if file2_filename:
            combined_table.add_row("Reverse Reads", f"[bold]{file2_filename}[/bold]")
        if num_threads:
            combined_table.add_row("Threads", f"[bold]{num_threads}[/bold]")
        if platform.system:
            combined_table.add_row(
                "Operating System", f"[bold]{platform.system()}[/bold]"
            )

        # Heuristic Statistics Sub-heading
        combined_table.add_section()
        combined_table.add_row("[bold][bright_blue]Heuristics[/bright_blue][/bold]", "")

        if bc_len:
            combined_table.add_row("Barcode Length", f"[bold]{bc_len}[/bold]")
        if bc_start1:
            combined_table.add_row("Forward Offset", f"[bold]{bc_start1}[/bold]")
        if bc_start2:
            combined_table.add_row("Reverse Offset", f"[bold]{bc_start2}[/bold]")
        if L_fwd or R_fwd:
            combined_table.add_row("Forward Flanks", f"[bold]{L_fwd}...{R_fwd}[/bold]")
        if L_rev or R_rev:
            combined_table.add_row("Reverse Flanks", f"[bold]{L_rev}...{R_rev}[/bold]")

        # Numeric Statistics Sub-heading
        combined_table.add_section()
        # combined_table.add_row("[bold][bright_green]Barcode Alignment Stats[/bright_green][/bold]", "")
        combined_table.add_row(
            "[bold]Total Reads[/bold]", f"[bold]{total_reads:,}[/bold]"
        )
        combined_table.add_row(
            "Documented Barcode Reads", f"[bold]{sum(doc_bcs.values()):,}[/bold]"
        )
        combined_table.add_row(
            "Undocumented Barcode Reads", f"[bold]{sum(undoc_bcs.values()):,}[/bold]"
        )
        combined_table.add_section()
        combined_table.add_row("[bold]Documented Barcodes[/bold]", f"{len(barcodes):,}")
        combined_table.add_row(
            "Seen Documented Barcodes", f"[bold]{len(doc_bcs):,}[/bold]"
        )
        combined_table.add_row(
            "Unseen Documented Barcodes",
            f"[bold]{(len(barcodes) - len(doc_bcs)):,}[/bold]",
        )
        combined_table.add_section()
        combined_table.add_row(
            "[bold]Uncocumented Barcodes[/bold]", f"{len(undoc_bcs):,}"
        )
        combined_table.add_section()
        combined_table.add_row(
            "[bold]Barcoded Reads Fraction[/bold]",
            f"[bold]{((sum(doc_bcs.values()) + sum(undoc_bcs.values())) / total_reads if total_reads != 0 else 0):.3f}[/bold]",
        )
        combined_table.add_row(
            "Documented Fraction",
            f"[bold]{(sum(doc_bcs.values()) / total_reads if total_reads != 0 else 0):.3f}[/bold]",
        )
        combined_table.add_row(
            "Undocumented Fraction",
            f"[bold]{(sum(undoc_bcs.values()) / total_reads if total_reads != 0 else 0):.3f}[/bold]",
            end_section=True,
        )

        # Sequence Information Sub-heading
        combined_table.add_section()
        top_doc_bcs = min(5, len(doc_bcs))
        combined_table.add_row(
            f"[bold bright_green]Top {top_doc_bcs} Documented Barcodes[/bold bright_green]",
            "",
        )
        for idx, (barcode, count) in enumerate(doc_bcs.most_common(top_doc_bcs)):
            end_section = idx == (top_doc_bcs - 1)
            combined_table.add_row(barcode, f"{count:,}", end_section=end_section)

        combined_table.add_section()
        top_undoc_bcs = min(5, len(undoc_bcs))
        combined_table.add_row(
            f"[bold bright_red]Top {top_undoc_bcs} Undocumented Barcodes[/bold bright_red]",
            f"",
        )
        # combined_table.add_row("[bold]Found Undocumented Barcodes[bold]", f"[bold]{len(undoc_bcs):,}[/bold]")
        for idx, (barcode, count) in enumerate(undoc_bcs.most_common(top_undoc_bcs)):
            end_section = idx == (top_undoc_bcs - 1)
            combined_table.add_row(barcode, f"{count:,}", end_section=end_section)

        # Print the combined table
        console.log(combined_table)

        for barcode, count in doc_bcs.items():
            print("\t".join([barcode, str(count)]))

        for barcode, count in undoc_bcs.items():
            print("\t".join([barcode, str(count)]))

    except ValueError as ve:
        console_error.log(str(ve))
        console.log("\n" + parser.format_help())

    except Exception as e:
        console_error.log(f"An unexpected error occurred: {str(e)}")
        console.log("\n" + parser.format_help())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Barcodes.")
    parser.add_argument(
        "fasta_file", type=str, help="List or FASTA file containing barcodes."
    )
    parser.add_argument("file1", type=str, help="First reads file: FASTQ or raw reads.")
    parser.add_argument(
        "file2",
        type=str,
        nargs="?",
        default=None,
        help="Second reads file: FASTQ or raw reads (optional).",
    )
    args = parser.parse_args()
    main(args)
