import argparse
import copy
import gzip
import hashlib
import os
import platform
import re
import subprocess
import sys
from collections import Counter
from datetime import datetime
from multiprocessing import cpu_count

import pandas as pd
import pysam
import rich
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord
from rich.console import Console
from rich.table import Table


# Utility function to open files
def open_file(file, mode):
    return gzip.open(file, mode) if file.endswith(".gz") else open(file, mode)

# Create a working directory
def create_working_directory(dir_name="working_directory"):
  
    os.makedirs(dir_name, exist_ok=True)
    return dir_name

def convert_genbank_to_fasta(genbank_file, fasta_file, overhang_length=0):
    modified_records = []
    open_func = gzip.open if genbank_file.endswith(".gz") else open

    with open_func(genbank_file, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            if record.annotations.get('topology', None) == 'circular':
                # Create new sequence with overhang
                overhang_length = 10000
                new_seq = record.seq + record.seq[:overhang_length]
                modified_record = SeqRecord(Seq(str(new_seq)), id=record.id, description=record.description)
                modified_records.append(modified_record)

    # Write the modified records to a new FASTA file
    with open(fasta_file, "w") as output_handle:
        SeqIO.write(modified_records, output_handle, "fasta")

# Convert fasta to fake fastq
def convert_fasta_to_fake_fastq(fasta_file, fastq_file):
    # if it is gzipped, open with gzip.open otherwise open with open
    if fasta_file.endswith(".gz"):
        with gzip.open(fasta_file, "rt") as input_handle, open(fastq_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")
    else:
        with open(fasta_file, "rt") as input_handle, open(fastq_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")

# This function runs Bowtie, which is a tool for aligning short DNA sequences
def run_bowtie(sgrna_fastq_file, fasta_file, sam_file, num_mismatches, num_threads):

    index_prefix = "genome_index"

    with open(os.devnull, "w") as devnull:
        subprocess.run(
            ["bowtie-build", fasta_file, index_prefix],
            stdout=devnull,
            stderr=devnull,
        )
        subprocess.run(
            ["bowtie", "--all", "--nomaqround", "-p", str(num_threads), "--tryhard", "-v", str(num_mismatches), "-S", index_prefix, sgrna_fastq_file, sam_file],
            stdout=devnull,
            stderr=devnull,
        )
        
        # Remove the index files, which names start with the index_prefix and end with ebwt
        for file in os.listdir("."):
            if file.startswith(index_prefix) and file.endswith(".ebwt"):
                os.remove(file)

def create_locus_map(genbank_file):
    locus_map, overhangs, organisms, seq_lens, topologies = {}, {}, {}, {}, {}

    with open(genbank_file, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            overhangs[record.id] = 0
            organisms[record.id] = record.annotations.get('organism', None)
            seq_lens[record.id] = len(record.seq)
            topologies[record.id] = record.annotations.get('topology', None)

            overhang_length = 10000 if topologies[record.id] == 'circular' else 0

            for feature in record.features:
                if feature.type == "gene":
                    for part_location in feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]:
                        
                        # Skip the gene segment that wraps around the genome
                        if part_location.end == seq_lens[record.id]:
                            continue

                        for position in range(int(part_location.start), int(part_location.end)):
                            key = (record.id, position)
                            locus_map.setdefault(key, []).append(
                                (feature.qualifiers.get("locus_tag", [None])[0],
                                int(part_location.start),
                                int(part_location.end),
                                feature.strand)
                            )

                        # Handle special case: gene that wraps around the genome
                        if part_location.start == 0 and feature.location.parts[-1].end == seq_lens[record.id]:
                            # Store the end of the first segment for later use as overhang
                            overhangs[record.id] = feature.location.parts[0].end
                            # Add the gene that wraps around to the locus_map
                            for position in range(feature.location.parts[-1].start, part_location.end + seq_lens[record.id]):
                                key = (record.id, position)
                                locus_map.setdefault(key, []).append(
                                    (feature.qualifiers.get("locus_tag", [None])[0],
                                    int(feature.location.parts[-1].start),
                                    int(feature.location.parts[-1].end + part_location.end),
                                    feature.strand)
                                )

    # Add overhang entries to the locus_map
    with open(genbank_file, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "gene":
                    for part_location in feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]:
                        for position in range(int(part_location.start), int(part_location.end)):
                            if position > overhangs[record.id] and position <= overhang_length:
                                key = (record.id, position + seq_lens[record.id])
                                locus_map.setdefault(key, []).append(
                                    (feature.qualifiers.get("locus_tag", [None])[0],
                                    int(part_location.start) + seq_lens[record.id],
                                    int(part_location.end) + seq_lens[record.id],
                                    feature.strand)
                                )

    return locus_map, organisms, seq_lens, topologies

# Reconstruct the target sequence using the cigar string
def reconstruct_target(read):
    reference_seq = read.get_reference_sequence()
    return reference_seq

def get_chrom_lengths(gb_file):
    chrom_lengths = {}
    with open(gb_file, "r") as f:
        for rec in SeqIO.parse(f, "genbank"):
            chrom_lengths[rec.id] = len(rec.seq)
    return chrom_lengths

def get_fasta_lengths(fasta_path):
    fasta_lengths = {}
    with pysam.FastaFile(fasta_path) as fasta:
        for chrom in fasta.references:
            fasta_lengths[chrom] = len(fasta.fetch(reference=chrom))
    return fasta_lengths

def extract_pam(tar_start, tar_end, chrom, fasta_path, dir, chrom_lengths):
    # Get the length of the chromosome, default to None if not found
    chrom_length = chrom_lengths.get(chrom, None)
    
    # Open the FASTA file for reading
    with pysam.FastaFile(fasta_path) as fasta:
        # If the direction is forward ('F')
        if dir == 'F':
            # Fetch the PAM sequence, considering the chromosome's circularity
            pam = fasta.fetch(reference=chrom, 
                              start=tar_end % chrom_length, 
                              end=(tar_end + 3) % chrom_length).upper()
        # If the direction is reverse ('R')
        elif dir == 'R':
            # Fetch the PAM sequence and reverse complement it
            pam = fasta.fetch(reference=chrom, 
                              start=(tar_start - 3) % chrom_length, 
                              end=tar_start % chrom_length).upper()
            pam = str(Seq(pam).reverse_complement())
        # Fallback case: Should not occur under normal circumstances
        else:
            pam = None

    return pam

def get_diff(spacer, target):
    # Initialize an empty list to store the differences
    differences = []
    
    # Iterate through the nucleotides in the target and spacer sequences
    for i, (target_nt, spacer_nt) in enumerate(zip(target, spacer)):
        # Check if the nucleotides differ at the current position
        if target_nt != spacer_nt:
            # Record the difference in the format "target_nucleotide position spacer_nucleotide"
            diff_string = f"{target_nt}{i + 1}{spacer_nt}"
            differences.append(diff_string)
    
    # Join the list of differences into a single string, separated by commas
    diff_result = ",".join(differences)
    
    # If there are no differences, return a hyphen "-"
    if not diff_result:
        return "-"
    
    return diff_result

def get_coords(tar_start, tar_end, chrom_length):
    # Take mod to account for circularity of chromosome
    start_circular = tar_start % chrom_length
    end_circular = tar_end % chrom_length if tar_end % chrom_length != 0 else chrom_length
    
    # Condition when start and end positions are in the same circular span
    if start_circular <= end_circular:
        return f"{start_circular}..{end_circular}"
    # Condition when the end position loops back to the start of the chromosome
    elif start_circular > end_circular:
        return f"({start_circular}..{chrom_length}, 0..{end_circular})"
    # Fallback case: Should not occur under normal circumstances
    else:
        return None

def get_offset(target_dir, tar_start, tar_end, feature_start, feature_end):
    # If the direction of the target gene is forward ("F")
    if target_dir == "F":
        # The offset is calculated as the start position of the target minus the start position of the feature
        return tar_start - feature_start
    # If the direction of the target gene is reverse ("R")
    elif target_dir == "R":
        # The offset is calculated as the end position of the feature minus the end position of the target
        return feature_end - tar_end
    # Fallback case: Direction is neither forward nor reverse, return None
    else:
        return None

def hash_row(row):
    return hashlib.sha256(str(row).encode()).hexdigest()

def parse_sam_output(sam_file_path, gene_locus_map, reference_fasta_path, reference_gb_path):
    unique_rows = {}
    chromosomal_lengths = get_chrom_lengths(reference_gb_path)
    samfile = pysam.AlignmentFile(sam_file_path, "r")
    unmapped = 0

    for read in samfile.fetch():
        if read.is_unmapped: 
            unmapped+=1
            continue

        row_data = {
            'name': read.query_name,
            'spacer': read.query_sequence if not read.is_reverse else str(Seq(read.query_sequence).reverse_complement()),
            'len': read.query_length,
            'chr': read.reference_name,
            'tar_start': read.reference_start,
            'tar_end': read.reference_end,
            'target': reconstruct_target(read) if not read.is_reverse else str(Seq(reconstruct_target(read)).reverse_complement()),
            'mismatches': sum(a != b for a, b in zip(read.query_sequence, reconstruct_target(read))),
            'sp_dir': "F" if not read.is_reverse else "R",
            'pam': extract_pam(read.reference_start, read.reference_end, read.reference_name, reference_fasta_path, "F" if not read.is_reverse else "R", chromosomal_lengths),
            'diff': get_diff(read.query_sequence, reconstruct_target(read)),
            'coords': get_coords(read.reference_start, read.reference_end, chromosomal_lengths.get(read.reference_name, None)),
        }

        aligned_genes = {gene_info for pos in range(row_data['tar_start'], row_data['tar_end'])
                        for gene_info in gene_locus_map.get((row_data['chr'], pos), [])}

        # if not re.match('.GG$', row_data['pam']): #NGG pam
        #     unmapped+=1
        #     continue

        if not aligned_genes:
            row_data_copy = copy.deepcopy(row_data)
            row_data_copy.update({'locus_tag': None, 'offset': None, 'tar_dir': None})

            unique_key = hash_row(row_data_copy.values)
            unique_rows[unique_key] = row_data_copy
            
            continue

        for locus_tag, feature_start, feature_end, feature_strand in aligned_genes:
            row_data_copy = copy.deepcopy(row_data)
            target_orientation = "F" if feature_strand == 1 else "R" if feature_strand == -1 else None
            offset_value = get_offset(target_orientation, row_data['tar_start'], row_data['tar_end'], feature_start, feature_end)

            row_data_copy.update({'locus_tag': locus_tag, 'offset': offset_value, 'tar_dir': target_orientation})
            unique_key = hash_row(row_data_copy.values)
            unique_rows[unique_key] = row_data_copy
                  
    samfile.close()
    return unique_rows, unmapped

# def main(sgrna_file, gb_file, num_mismatches):
def main(args):
    console=Console(file=sys.stderr)
    console.log("[bold red]Initializing barcode target seeker[/bold red]")

    num_threads = cpu_count()

    working_dir = create_working_directory()
    
    fasta_file = os.path.join(
        working_dir, os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta")
    
    sgrna_fastq_file = os.path.join(
        working_dir, os.path.splitext(os.path.basename(args.sgrna_file ))[0] + ".fastq")
    
    output_folder = "results"

    sam_file = os.path.join(output_folder,
        f"{os.path.splitext(os.path.basename(args.sgrna_file))[0]}_{os.path.splitext(os.path.basename(args.genome_file))[0]}.sam")

    os.makedirs(output_folder, exist_ok=True)

    console.log("Generating topological coordinate maps...")
    locus_map, organisms, seq_lens, topologies = create_locus_map(args.genome_file)

    console.log("Generating topological gene maps...")
    convert_genbank_to_fasta(args.genome_file, fasta_file)
    # Delete existing .fai file if it exists
    fai_file = fasta_file + ".fai"
    if os.path.exists(fai_file):
        os.remove(fai_file)

    console.log("Annotating regions to identify...")
    convert_fasta_to_fake_fastq(args.sgrna_file, sgrna_fastq_file)
    if not os.path.exists(fasta_file):
        console.log(f"[bold red]File \"{fasta_file}\" does not exist.[/bold red]")
        sys.exit(1)

    console.log("Aligning annotations to genome...")
    run_bowtie(sgrna_fastq_file, fasta_file, sam_file, args.mismatches, num_threads)

    with pysam.FastaFile(fasta_file) as _:
        pass

    try:
        console.log("Finding matches...")    
        results, unmapped = parse_sam_output(sam_file, locus_map, fasta_file, args.genome_file)
        results = pd.DataFrame.from_dict(results, orient='index')
        column_order = ['name', 'spacer', 'pam', 'len', 'locus_tag', 'chr', 'target', 'mismatches', 'diff', 'coords', 'offset', 'sp_dir', 'tar_dir']
        results = results[column_order]
        integer_cols = ['len', 'mismatches', 'offset']
        results[integer_cols] = results[integer_cols].astype('Int64')

    except FileNotFoundError:
        console.log(f"[bold red]Trouble with Bowtie aligner. User a lower number of mismatches.[/bold red]")
        sys.exit(1)

    except KeyError:
        console.log(f"[bold red]None of the proposed barcodes map to the genomes.[/bold red]")
        sys.exit(1)



    

    console.log(f"Cleaning up...")

    # Delete fastq file, fasta file, and sam file, since they are no longer needed
    os.remove(sam_file)
    os.remove(fasta_file)
    os.remove(sgrna_fastq_file)
    os.remove(fasta_file + ".fai")

    # Delete working directory, if it is empty
    if not os.listdir(working_dir):
        os.rmdir(working_dir)

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
    combined_table.add_row("Barcodes", f"[bold]{os.path.basename(args.sgrna_file)}[/bold]")
    combined_table.add_row("Genbank Genome File", f"[bold]{os.path.basename(args.genome_file)}[/bold]")
    combined_table.add_row("Number of Mismatches", f"[bold]{args.mismatches}[/bold]")
    combined_table.add_row("Threads", f"[bold]{num_threads}[/bold]")
    combined_table.add_row("Operating System", f"[bold]{platform.system()}[/bold]")

    # Heuristic Statistics Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_blue]Heuristics[/bold bright_blue]", "")

    unique_organisms = set(organisms.values())

    if len(unique_organisms) == 1:
        organism_entry = f"[bold]{next(iter(unique_organisms))}[/bold]"
    elif len(unique_organisms) > 1:
        organism_entry = f"[bold]{', '.join(unique_organisms)}[/bold]"
    else:
        organism_entry = "[bold]Unknown[/bold]"

    combined_table.add_row("Organism", organism_entry)

    ambiguous_coordinates = {pos for pos, count in Counter(key[1] for key in locus_map.keys()).items() if count > 1}
    overlapping_genes = {tag for _, values in locus_map.items() if len(values) > 1 for value in values for tag in [value[0]]}
    combined_table.add_row("Chromosomes", f"[bold]{len(set(key[0] for key in locus_map.keys()))}[/bold]")
    combined_table.add_row("Total Genes", f"[bold]{len(set(value[0][0] for value in locus_map.values()))}[/bold]")
    
 
    combined_table.add_row("Overlapping Genes", f"[bold]{len(overlapping_genes)}[/bold]")
    combined_table.add_row("Ambiguous Coordinates", f"[bold]{len(ambiguous_coordinates)}[/bold]")

    spanning_genes = {}

    # for (chromosome, position), values in locus_map.items():
    #     if position % seq_lens[chromosome] == 0:
    #         for value in values:
    #             tag, _, _, _ = value
    #             if chromosome not in spanning_genes:
    #                 spanning_genes[chromosome] = []
    #             spanning_genes[chromosome].append(tag)

    # for chromosome, genes in spanning_genes.items():
    #     gene_list = ", ".join(genes)
    #     combined_table.add_row(f"Origin-Spanning: {chromosome}", f"[bold]{gene_list}[/bold]")

    # Barcode Mapping Stats Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_green]Barcode Mapping Stats[/bold bright_green]", "")

    unique_tags = len(set(x for x in results['locus_tag'] if x is not None))
    combined_table.add_row("Genes Targeted", f"[bold]{unique_tags}[/bold]")

    combined_table.add_row("Targeting Barcodes", f"[bold]{len(results)}[/bold]")

    unique_spacers_per_mismatch = results.groupby('mismatches')['spacer'].apply(set)

    for mismatch, unique_spacers in unique_spacers_per_mismatch.items():
        count = len(unique_spacers)
        combined_table.add_row(f"{mismatch} Mismatch Barcodes", f"[bold]{count}[/bold]")

    none_targeting_unique_spacers = len(set(spacer for spacer, tag in zip(results['spacer'], results['locus_tag']) if tag is None))
    combined_table.add_row("Intergenic Barcodes", f"[bold]{none_targeting_unique_spacers}[/bold]")

    spacer_counts = Counter(results['spacer'])
    duplicates = {spacer: count for spacer, count in spacer_counts.items() if count > 1}
    num_duplicates = len(duplicates)
    combined_table.add_row("Ambiguous Barcodes", f"[bold]{num_duplicates}[/bold]")

    combined_table.add_row("Non-targeting Barcodes", f"[bold]{unmapped}[/bold]")



    
    # Print the combined table
    console.log(combined_table)

    results.to_csv(sys.stdout, sep='\t', index=False, na_rep='None')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map barcodes to a circular genome")
    parser.add_argument("sgrna_file", help="Path to sgrna_fasta_file", type=str)
    parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
    parser.add_argument("mismatches", help="Number of allowed mismatches", type=int)

    args = parser.parse_args()

    main(args)