import argparse
import copy
import gzip
import hashlib
import os
import platform
import re
import subprocess
import sys
from collections import Counter, defaultdict
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


def open_file(file, mode):
    return gzip.open(file, mode) if file.endswith(".gz") else open(file, mode)

def create_working_directory(dir_name="working_directory"):
  
    os.makedirs(dir_name, exist_ok=True)
    return dir_name

def create_topological_fasta(genbank_file_name, topological_fasta_file_name, overhang_length=0):
    topological_records = []
    open_func = gzip.open if genbank_file_name.endswith(".gz") else open

    with open_func(genbank_file_name, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            if record.annotations.get('topology', None) == 'circular':
                overhang_length = 100_000

            new_seq = record.seq + record.seq[:overhang_length]
            topological_record = SeqRecord(Seq(str(new_seq)), id=record.id, description=record.description)
            topological_records.append(topological_record)

    with open(topological_fasta_file_name, "w") as output_handle:
        SeqIO.write(topological_records, output_handle, "fasta")


def create_fake_topological_fastq(topological_fasta_file_name, fastq_file):
    if topological_fasta_file_name.endswith(".gz"):
        with gzip.open(topological_fasta_file_name, "rt") as input_handle, open(fastq_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")
    else:
        with open(topological_fasta_file_name, "rt") as input_handle, open(fastq_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")

def run_bowtie(sgrna_fastq_file_name, topological_fasta_file_name, sam_file_name, num_mismatches, num_threads):
    index_prefix = "genome_index"

    with open(os.devnull, "w") as devnull:
        subprocess.run(
            ["bowtie-build", topological_fasta_file_name, index_prefix],
            stdout=devnull,
            stderr=devnull,
        )
        subprocess.run(
            ["bowtie", "--all", "--nomaqround", "-p", str(num_threads), "--tryhard", "-v", str(num_mismatches), "-S", index_prefix, sgrna_fastq_file_name, sam_file_name],
            stdout=devnull,
            stderr=devnull,
        )
        
        for file in os.listdir("."):
            if file.startswith(index_prefix) and file.endswith(".ebwt"):
                os.remove(file)

def create_locus_map(genbank_file_name):
    locus_map, overhangs, organisms, seq_lens, topologies = {}, {}, {}, {}, {}

    with open(genbank_file_name, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            overhangs[record.id] = 0
            organisms[record.id] = record.annotations.get('organism', None)
            seq_lens[record.id] = len(record.seq)
            topologies[record.id] = record.annotations.get('topology', None)

            overhang_length = 100_000 if topologies[record.id] == 'circular' else 0

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
    with open(genbank_file_name, "rt") as input_handle:
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

def get_true_chrom_lengths(gb_file):
    chrom_lengths = {}
    with open(gb_file, "r") as f:
        for rec in SeqIO.parse(f, "genbank"):
            chrom_lengths[rec.id] = len(rec.seq)
    return chrom_lengths

def get_topological_chrom_lengths(fasta_file):
    chrom_lengths = {}
    with open(fasta_file, "r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            chrom_lengths[rec.id] = len(rec.seq)
    return chrom_lengths

def extract_pam(pam, tar_start, tar_end, chrom, topological_fasta_file_name, dir, true_chrom_lengths, topological_chrom_lengths):
    true_chrom_length = true_chrom_lengths.get(chrom, None)
    topological_chrom_length = topological_chrom_lengths.get(chrom, None)

    if None in (pam, tar_start, tar_end, chrom, topological_fasta_file_name, dir, true_chrom_length, topological_chrom_length):
        return None

    with pysam.FastaFile(topological_fasta_file_name) as fasta:
        if dir == 'F':
            if tar_end + len(pam) > topological_chrom_length:
                return None
            pam = fasta.fetch(reference=chrom, start=tar_end, end=tar_end + len(pam)).upper()

        elif dir == 'R':
            if tar_start - len(pam) < 0:
                return None
            pam = fasta.fetch(reference=chrom, start=tar_start - len(pam), end=tar_start).upper()
            pam = str(Seq(pam).reverse_complement())

        else:
            return None

    return pam


def get_diff(spacer, target):
    differences = []
    
    for i, (target_nt, spacer_nt) in enumerate(zip(target, spacer)):
        if target_nt != spacer_nt:
            diff_string = f"{target_nt}{i + 1}{spacer_nt}"
            
            differences.append(diff_string)
    
    diff_result = ",".join(differences)
    
    if not diff_result:
        return None
    
    return diff_result

def get_coords(tar_start, tar_end, chrom_length):
    start_circular = tar_start % chrom_length
    end_circular = tar_end % chrom_length if tar_end % chrom_length != 0 else chrom_length

    if start_circular > end_circular:
        return f"({start_circular}..{chrom_length}, 0..{end_circular})"    

    return f"{start_circular}..{end_circular}"


def get_offset(target_dir, tar_start, tar_end, feature_start, feature_end):
    if target_dir == "F":
        return tar_start - feature_start
    elif target_dir == "R":
        return feature_end - tar_end
    else:
        return None

def hash_row(row_data):
    canonical_str = "".join([f"{k}:{v}|" for k, v in sorted(row_data.items(), key=lambda x: x[0])])
    return hashlib.md5(canonical_str.encode()).hexdigest()

def pam_matches(pam_pattern, extracted_pam):
    # Convert N to . for regex matching
    regex_pattern = pam_pattern.replace('N', '.')
    return bool(re.match(regex_pattern, extracted_pam))

def parse_sam_output(sam_file_name, locus_map, topological_fasta_file_name, gb_file_name, pam):

    unique_rows = {}
    true_chrom_lengths = get_true_chrom_lengths(gb_file_name)
    topological_chrom_lengths = get_topological_chrom_lengths(topological_fasta_file_name)

    with pysam.AlignmentFile(sam_file_name, "r") as samfile:
        for read in samfile.fetch():
            
            extracted_pam = extract_pam(
                pam, read.reference_start, read.reference_end, read.reference_name, topological_fasta_file_name, 
                "F" if not read.is_reverse else "R", true_chrom_lengths, topological_chrom_lengths)
            
            row_data = {
                'name': read.query_name,
                'spacer': read.query_sequence if not read.is_reverse else str(Seq(read.query_sequence).reverse_complement()),
                'len': len(read.query_sequence),
           }

            # exclude reads that don't match the PAM
            if not read.is_unmapped and not pam_matches(pam, extracted_pam):
                read.is_unmapped = True
            
            # skip reads where topological map ran off the end of the chromosome
            if not read.is_unmapped and extracted_pam is None:
                continue
        
            if read.is_unmapped:
                row_data.update({
                    'type': 'non-targeting',
                    'chr': None,
                    'tar_start': None,
                    'tar_end': None,
                    'target': None,
                    'mismatches': None,
                    'sp_dir': None,
                    'pam': None,
                    'diff': None,
                    'coords': None,
                    'locus_tag': None, 
                    'offset': None, 
                    'tar_dir': None
                })

                unique_key = hash_row(row_data)
                unique_rows[unique_key] = row_data
                continue
            
            if not read.is_unmapped:
                row_data.update({
                    'target': read.get_reference_sequence() if not read.is_reverse else str(Seq(read.get_reference_sequence()).reverse_complement()),
                    'mismatches': read.get_tag('NM'),
                    'chr': read.reference_name,
                    'tar_start': read.reference_start % true_chrom_lengths.get(read.reference_name, None),
                    'tar_end': read.reference_end % true_chrom_lengths.get(read.reference_name, None),
                    'sp_dir': "F" if not read.is_reverse else "R",
                    'pam': extracted_pam,
                    'coords': get_coords(read.reference_start, read.reference_end, true_chrom_lengths.get(read.reference_name, None))
                })    

                row_data.update({
                    'type': 'mismatch' if row_data['mismatches'] > 0 else 'perfect',
                    'diff': get_diff(row_data['spacer'], row_data['target'])
                })

                aligned_genes = {gene_info for pos in range(row_data['tar_start'], row_data['tar_end'])
                                for gene_info in locus_map.get((row_data['chr'], pos), [])}

                if not aligned_genes:
                    row_data.update({
                        'locus_tag': None, 
                        'offset': None, 
                        'tar_dir': None
                    })
                    unique_rows[hash_row(row_data)] = row_data                
                    continue
                
                for locus_tag, feature_start, feature_end, feature_strand in aligned_genes:
                    row_data_copy = row_data.copy()
                    target_orientation = "F" if feature_strand == 1 else "R" if feature_strand == -1 else None
                    offset_value = get_offset(
                        target_orientation, 
                        row_data_copy['tar_start'], 
                        row_data_copy['tar_end'], 
                        feature_start, 
                        feature_end)

                    row_data_copy.update({
                        'locus_tag': locus_tag, 
                        'offset': offset_value, 
                        'tar_dir': target_orientation
                    })
                    
                    unique_rows[hash_row(row_data_copy)] = row_data_copy

    
    def filter_offtargets_by_pam(unique_rows):
        targeting_guides = {row_data['spacer'] for row_data in unique_rows.values() if row_data['target']}
        # Only keep non-targeting guides if they don't appear in the targeting guides set.
        return {key: val for key, val in unique_rows.items() if not (val['target'] is None and val['spacer'] in targeting_guides)}

    unique_rows = filter_offtargets_by_pam(unique_rows)

    note = {}

    # Create a dictionary to store unique coordinates for each name
    coords_by_spacer = {}
    for row in unique_rows.values():
        spacer = row['spacer']
        coords = row['coords']
        if spacer not in coords_by_spacer:
            coords_by_spacer[spacer] = []
        coords_by_spacer[spacer].append(coords)

    # Update notes based on the number of unique coordinates and ambiguous annotations
    for spacer, coords in coords_by_spacer.items():
        unique_coords = set(coords)
        if len(unique_coords) > 1:
            if spacer not in note:
                note[spacer] = []
            note[spacer].append(f"{len(unique_coords)} targets")
        
        # Add note for ambiguous annotations
        for coord in unique_coords:
            count = coords.count(coord)
            if count > 1:
                if spacer not in note:
                    note[spacer] = []
                note[spacer].append(f"{count} annotations")

    # Join the notes and add them to the corresponding rows in unique_rows
    for row_data in unique_rows.values():
        spacer = row_data['spacer']
        if spacer in note:
            row_data['note'] = ', '.join(note[spacer])

    # return unique_rows
    return unique_rows



# def main(sgrna_file, gb_file, num_mismatches):
def main(args):
    console=Console(file=sys.stderr)
    console.log("[bold red]Initializing barcode target seeker[/bold red]")

    num_threads = cpu_count() // 2

    working_dir = create_working_directory()
    
    topological_fasta_file_name = os.path.join(
        working_dir, os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta")
    
    sgrna_fastq_file_name = os.path.join(
        working_dir, os.path.splitext(os.path.basename(args.sgrna_file ))[0] + ".fastq")
    
    output_folder = "results"

    sam_file_name = os.path.join(output_folder,
        f"{os.path.splitext(os.path.basename(args.sgrna_file))[0]}_{os.path.splitext(os.path.basename(args.genome_file))[0]}.sam")

    os.makedirs(output_folder, exist_ok=True)

    console.log("Generating topological coordinate maps...")
    locus_map, organisms, seq_lens, topologies = create_locus_map(args.genome_file)

    console.log("Generating topological gene maps...")
    create_topological_fasta(args.genome_file, topological_fasta_file_name)
    # Delete existing .fai file if it exists
    fai_file = topological_fasta_file_name + ".fai"
    if os.path.exists(fai_file):
        os.remove(fai_file)

    console.log("Annotating regions to identify...")
    create_fake_topological_fastq(args.sgrna_file, sgrna_fastq_file_name)
    if not os.path.exists(topological_fasta_file_name):
        console.log(f"[bold red]Unable to create topological map of the genome. [/bold red]")
        sys.exit(1)

    console.log("Aligning annotations to genome...")
    run_bowtie(sgrna_fastq_file_name, topological_fasta_file_name, sam_file_name, args.mismatches, num_threads)

    with pysam.FastaFile(topological_fasta_file_name) as _:
        pass

    try:
        console.log("Finding matches...")    
        results = parse_sam_output(sam_file_name, locus_map, topological_fasta_file_name, args.genome_file, args.pam)
        results = pd.DataFrame.from_dict(results, orient='index')
        column_order = ['name', 'spacer', 'pam', 'chr', 'locus_tag', 'target', 'type', 'mismatches', 'diff', 'coords', 'offset', 'sp_dir', 'tar_dir', 'note']
        results = results[column_order]
        integer_cols = ['mismatches', 'offset']
        
        results[integer_cols] = results[integer_cols].astype('Int64')

    except FileNotFoundError:
        console.log(f"[bold red]Trouble with Bowtie aligner. User a lower number of mismatches.[/bold red]")
        sys.exit(1)

    except KeyError:
        console.log(f"[bold red]None of the proposed barcodes map to the genomes.[/bold red]")
        sys.exit(1)

    console.log(f"Cleaning up...")

    # Delete fastq file, fasta file, and sam file, since they are no longer needed
    os.remove(sam_file_name)
    os.remove(topological_fasta_file_name)
    os.remove(sgrna_fastq_file_name)
    os.remove(topological_fasta_file_name + ".fai")

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
    combined_table.add_row("PAM", f"[bold]{args.pam}[/bold]")
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

    unique_topologies = set(topologies.values())
    combined_table.add_row("Topology", f"[bold]{', '.join(unique_topologies)}[/bold]")

    unique_seq_lens = set(seq_lens.values())
    combined_table.add_row("Sequence Length", f"[bold]{'; '.join(format(seq_len, ',') for seq_len in unique_seq_lens)}[/bold]")

    ambiguous_coordinates = {(chrom, pos % seq_lens[chrom]) for chrom, pos in locus_map if len(locus_map[(chrom, pos)]) > 1}
    ambiguous_locus_tags = {entry[0] for chrom, pos in ambiguous_coordinates for entry in locus_map[(chrom, pos)]}
  
    combined_table.add_row("Chromosomes", f"[bold]{len(set(key[0] for key in locus_map.keys()))}[/bold]")
    combined_table.add_row("Total Genes", f"[bold]{len(set(value[0][0] for value in locus_map.values())):,}[/bold]")
    
    combined_table.add_row("Ambiguous Coordinates", f"[bold]{len(ambiguous_coordinates):,}[/bold]")
    combined_table.add_row("Overlapping Genes", f"[bold]{len(ambiguous_locus_tags):,}[/bold]")

    # Barcode Mapping Stats Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_green]Barcode Mapping Stats[/bold bright_green]", "")

    unique_tags = len(set(x for x in results['locus_tag'] if x is not None))
    combined_table.add_row("Genes Targeted", f"[bold]{unique_tags:,}[/bold]")

    # spacers targeting overlapping genes
    overlapping_spacers = results[results['locus_tag'].isin(ambiguous_locus_tags)]['locus_tag'].nunique()
    combined_table.add_row("Overlapping Genes Targeted", f"[bold]{overlapping_spacers:,}[/bold]")

    unique_barcodes = results['spacer'].nunique()

    combined_table.add_row("Unique Barcodes", f"[bold]{unique_barcodes:,}[/bold]")

    unique_spacers_per_mismatch = results.groupby('mismatches')['spacer'].apply(set)

    for mismatch, unique_spacers in unique_spacers_per_mismatch.items():
        count = len(unique_spacers)
        combined_table.add_row(f"{mismatch} Mismatch Barcodes", f"[bold]{count:,}[/bold]")

    intergenic_spacers = results[(results['locus_tag'].isnull()) & (results['chr'].notnull())]['spacer'].nunique()
    combined_table.add_row("Intergenic Barcodes", f"[bold]{intergenic_spacers:,}[/bold]")

    # spacers targeting multiple coordinates
    off_target_spacers = results[results['target'].notnull()].groupby('spacer')['coords'].apply(set).apply(len).gt(1).sum()
    combined_table.add_row("Off-targeting Barcodes", f"[bold]{off_target_spacers:,}[/bold]")

    nontargeting_spacers = results[results['target'].isnull()]['spacer'].nunique()
    combined_table.add_row("Non-targeting Barcodes", f"[bold]{nontargeting_spacers:,}[/bold]")

    # Print the combined table
    console.log(combined_table)

    results.to_csv(sys.stdout, sep='\t', index=False, na_rep='None')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map barcodes to a circular genome")
    parser.add_argument("sgrna_file", help="Path to sgrna_fasta_file", type=str)
    parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
    parser.add_argument("pam", help="PAM sequence", type=str)
    parser.add_argument("mismatches", help="Number of allowed mismatches", type=int)

    args = parser.parse_args()

    main(args)