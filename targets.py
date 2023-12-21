import argparse
import copy
from collections import Counter, defaultdict
from datetime import datetime
from multiprocessing import cpu_count

import gzip
import hashlib
import os
import pandas as pd
import platform
import pysam
import re
import rich
import rich.table
import subprocess
import sys
import tempfile

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


def create_locus_map(genbank_file_name):
    locus_map, overhang_continue, organisms, seq_lens, topologies = {}, {}, {}, {}, {}

    with open(genbank_file_name, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            overhang_continue[record.id] = 0
            organisms[record.id] = record.annotations.get('organism', None)
            seq_lens[record.id] = len(record.seq)
            topologies[record.id] = record.annotations.get('topology', None)

            overhang_length = 100_000 if topologies[record.id] == 'circular' else 0

            for feature in record.features:
                if feature.type == "gene":
                    
                    if isinstance(feature.location, CompoundLocation) and any(part.start == 0 or part.end == len(record.seq) for part in feature.location.parts):
                        
                        # Find segments that wrap around the genome
                        genome_end_segment = next(part for part in feature.location.parts if part.end == len(record.seq))
                        genome_start_segment = next(part for part in feature.location.parts if part.start == 0)
                        
                        # Adjust positions
                        adj_start = int(genome_end_segment.start)
                        adj_end = int(genome_start_segment.end + len(record.seq))
                        overhang_continue[record.id] = int(genome_start_segment.end)

                        # Iterate over the adjusted range
                        for position in range(adj_start, adj_end):
                            key = (record.id, position)
                            locus_map.setdefault(key, []).append(
                                (feature.qualifiers.get("locus_tag", [None])[0],
                                adj_start,
                                adj_end,
                                feature.strand)
                            )

                    else:
                        # normal genes
                        for part_location in feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]:
                            for position in range(int(part_location.start), int(part_location.end)):
                                key = (record.id, position)
                                locus_map.setdefault(key, []).append(
                                    (feature.qualifiers.get("locus_tag", [None])[0],
                                    int(part_location.start),
                                    int(part_location.end),
                                    feature.strand)
                                )

                                # add the rest of the overhang
                                if overhang_continue[record.id] <= position < overhang_length:
                                    key = (record.id, position + len(record.seq))
                                    locus_map.setdefault(key, []).append(
                                        (feature.qualifiers.get("locus_tag", [None])[0],
                                        int(part_location.start) + len(record.seq),
                                        int(part_location.end) + len(record.seq),
                                        feature.strand)
                                    )

    return locus_map, organisms, seq_lens, topologies

# Reconstruct the target sequence using the cigar string
def reconstruct_target(read):
    reference_seq = read.get_reference_sequence()
    return reference_seq

# Get the real length of chromosomes in the genbank file
def get_true_chrom_lengths(gb_file):
    chrom_lengths = {}
    with open(gb_file, "r") as f:
        for rec in SeqIO.parse(f, "genbank"):
            chrom_lengths[rec.id] = len(rec.seq)
    return chrom_lengths

# Get length of chromosomes in topological map
def get_topological_chrom_lengths(fasta_file):
    chrom_lengths = {}
    with open(fasta_file, "r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            chrom_lengths[rec.id] = len(rec.seq)
    return chrom_lengths

# Get the differences between the spacer and target
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

# Get the coordinates of the target, accounting for circularity
def get_coords(tar_start, tar_end, chrom_length):
    start_circular = tar_start % chrom_length
    end_circular = tar_end % chrom_length if tar_end % chrom_length != 0 else chrom_length

    if start_circular > end_circular:
        return f"({start_circular}..{chrom_length}, 0..{end_circular})"    
    return f"{start_circular}..{end_circular}"

# Get the offset of the target from the feature
def get_offset(target_dir, tar_start, tar_end, feature_start, feature_end):
    if target_dir == "F":
        return tar_start - feature_start
    elif target_dir == "R":
        return feature_end - tar_end
    else:
        return None

# Get the overlap of the target and feature
def get_overlap(tar_start, tar_end, feature_start, feature_end):
    overlap_start = max(tar_start, feature_start)
    overlap_end = min(tar_end, feature_end)

    # Check if there's any overlap
    if overlap_start < overlap_end:
        return overlap_end - overlap_start
    else:
        return 0

# Check if the extracted PAM matches the PAM pattern
def pam_matches(pam_pattern, extracted_pam):
    # Convert N to . for regex matching
    regex_pattern = pam_pattern.replace('N', '.')
    if extracted_pam is None:
        return False
    return bool(re.match(regex_pattern, extracted_pam))

# Extract the PAM of the target sequence from the topological map
def extract_pam(pam, tar_start, tar_end, chrom, fasta, dir, true_chrom_lengths, topological_chrom_lengths):
    true_chrom_length = true_chrom_lengths.get(chrom, None)
    topological_chrom_length = topological_chrom_lengths.get(chrom, None)

    if None in (pam, tar_start, tar_end, chrom, fasta, dir, true_chrom_length, topological_chrom_length):
        return None

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

# Parse the SAM file and extract the relevant information
def parse_sam_output(samfile, locus_map, topological_fasta_file_name, gb_file_name, pam):
    rows_list = []  # Initialize an empty list
    true_chrom_lengths = get_true_chrom_lengths(gb_file_name)
    topological_chrom_lengths = get_topological_chrom_lengths(topological_fasta_file_name)

    with pysam.FastaFile(topological_fasta_file_name) as fasta:
        for read in samfile.fetch():
            extracted_pam = extract_pam(
                pam, read.reference_start, read.reference_end, read.reference_name, fasta, 
                "F" if not read.is_reverse else "R", true_chrom_lengths, topological_chrom_lengths)
            
            if read.query_sequence is None:
                continue

            rows = {
                'name': read.query_name,
                'spacer': read.query_sequence if not read.is_reverse else str(Seq(read.query_sequence).reverse_complement()),
                'len': len(read.query_sequence),
           }    

            if not read.is_unmapped:
                # exclude reads that don't match the PAM
                if not pam_matches(pam, extracted_pam):
                    read.is_unmapped = True

                # skip reads where topological map ran out
                if extracted_pam is None:
                    continue
        
            if read.is_unmapped:
                rows.update({
                    'type': 'control',
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
                    'overlap': None,
                    'tar_dir': None
                })

                rows_list.append(rows)
                continue
            
            if not read.is_unmapped and read.reference_start is not None and read.reference_end is not None:
                
                try:
                    target = read.get_reference_sequence() if not read.is_reverse else str(Seq(read.get_reference_sequence()).reverse_complement())
                    mismatches = int(read.get_tag('NM'))
                    chr = read.reference_name
                    tar_start = read.reference_start % true_chrom_lengths.get(chr, None)
                    tar_end = read.reference_end % true_chrom_lengths.get(chr, None)
                    sp_dir = "F" if not read.is_reverse else "R"
                    coords = get_coords(tar_start, tar_end, true_chrom_lengths.get(chr, None))
                    type = 'mismatch' if mismatches > 0 else 'perfect'
                    diff = get_diff(rows['spacer'], target)
                
                except ValueError as e:
                    print(e, file=sys.stderr)
                    sys.exit(1)


                rows.update({
                    'target': target,
                    'mismatches': mismatches,
                    'chr': chr,
                    'tar_start': tar_start,
                    'tar_end': tar_end,
                    'sp_dir': sp_dir,
                    'pam': extracted_pam,
                    'coords': coords,
                    'type': type,
                    'diff': diff
                })

                aligned_genes = {gene_info for pos in range(tar_start, tar_end)
                                for gene_info in locus_map.get((chr, pos), [])}

                if not aligned_genes:
                    rows.update({
                        'locus_tag': None, 
                        'offset': None, 
                        'overlap': None,
                        'tar_dir': None
                    })
                    rows_list.append(rows)

                else:
                    for locus_tag, feature_start, feature_end, feature_strand in aligned_genes:
                        target_orientation = "F" if feature_strand == 1 else "R" if feature_strand == -1 else None
                        offset = get_offset(target_orientation, tar_start, tar_end, feature_start, feature_end)
                        overlap = get_overlap(tar_start, tar_end, feature_start, feature_end)

                        rows_copy = rows.copy()
                        rows_copy.update({
                            'locus_tag': locus_tag, 
                            'offset': offset,
                            'overlap': overlap,
                            'tar_dir': target_orientation
                        })
                        rows_list.append(rows_copy)

   
    return rows_list


# Run bowtie and parse the output using parse_sam_output and temp files
def run_bowtie_and_parse(sgrna_fastq_file_name, topological_fasta_file_name, locus_map, num_mismatches, num_threads):
    # Create a temporary directory for the bowtie index files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a temporary name for the genome index
        genome_index_temp_name = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False).name
        index_prefix = os.path.join(temp_dir, genome_index_temp_name)

        with open(os.devnull, "w") as devnull:
            bowtie_build_command = ["bowtie-build", topological_fasta_file_name, index_prefix]
            bowtie_build_process = subprocess.Popen(bowtie_build_command, stdout=devnull, stderr=devnull)
            bowtie_build_process.wait()

            # Create a temporary file for the bowtie output
            with tempfile.NamedTemporaryFile(delete=False) as bowtie_output_temp_file:
                bowtie_command = ["bowtie", "-k 100", "--nomaqround", "-p", str(num_threads), "--tryhard", "-v", str(num_mismatches), "-S", "-x", index_prefix, sgrna_fastq_file_name]
                bowtie_process = subprocess.Popen(bowtie_command, stdout=bowtie_output_temp_file, stderr=devnull)
                bowtie_process.wait()

                with pysam.AlignmentFile(bowtie_output_temp_file.name, "r") as samfile:
                    results = parse_sam_output(samfile, locus_map, topological_fasta_file_name, args.genome_file, args.pam)

                # Delete the temporary file after we're done with it
                os.remove(bowtie_output_temp_file.name)

    # Return the results of parse_sam_output
    return results

# Filter out spacers that don't match the PAM
def filter_offtargets_by_pam(df):
    targeting_spacers = df[df['target'].notna()]['spacer'].unique()
    return df[~((df['target'].isna()) & (df['spacer'].isin(targeting_spacers)))]

# Create a note for each spacer
def create_note(row):
    parts = []
    if row['sites'] > 0:
        parts.append(f"{row['sites']} {'site' if row['sites'] == 1 else 'sites'}")
        if row['genes'] > 0:
            parts.append(f"{row['genes']} {'gene' if row['genes'] == 1 else 'genes'}")
        if row['intergenic'] > 0:
            parts.append(f"{row['intergenic']} intergenic")
    else:
        parts.append("non-targeting")
    return ', '.join(parts)

# Main function
def main(args):
    console=Console(file=sys.stderr)
    console.log("[bold red]Initializing barcode target seeker[/bold red]")

    num_threads = cpu_count() // 2

    with tempfile.TemporaryDirectory() as working_dir:
        topological_fasta_file_name = os.path.join(
            working_dir, os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta")

        console.log("Annotating regions to identify...")
        
        base_name = os.path.basename(args.sgrna_file)
        if '.fastq' in base_name:
            sgrna_fastq_file_name = args.sgrna_file
        elif base_name.endswith('.fasta') or base_name.endswith('.fasta.gz'):
            ext = '.fasta.gz' if base_name.endswith('.fasta.gz') else '.fasta'
            base_name = base_name[:-len(ext)]
            sgrna_fastq_file_name = os.path.join(working_dir, base_name + ".fastq")
            create_fake_topological_fastq(args.sgrna_file, sgrna_fastq_file_name)
        else:
            console.log(f"[bold red]File extension not recognized. [/bold red]")
            sys.exit(1)

        console.log("Generating topological coordinate maps...")
        locus_map, organisms, seq_lens, topologies = create_locus_map(args.genome_file)

        create_topological_fasta(args.genome_file, topological_fasta_file_name)

        console.log("Aligning annotations to genome...")
        results = run_bowtie_and_parse(sgrna_fastq_file_name, topological_fasta_file_name, locus_map, args.mismatches, num_threads)

        with pysam.FastaFile(topological_fasta_file_name) as _:
            pass

    try:
        console.log("Finding matches...")

        # Convert your data to a DataFrame, then drop duplicates, which are caused by overhangs
        results = pd.DataFrame(results).drop_duplicates()

        # Use the filter_offtargets_by_pam function
        results = filter_offtargets_by_pam(results)
        
        # Create a 'min_tar' column that is the minimum of 'tar_start' and 'tar_end'
        results['min_tar'] = results[['tar_start', 'tar_end']].min(axis=1)

        # Sort the DataFrame by 'chr', 'min_tar', and 'spacer'
        results = results.sort_values(by=['chr', 'min_tar', 'spacer'])

        print(results, file=sys.stderr)

        # Create a 'site' column only for rows that have a 'target'
        results.loc[results['target'].notnull(), 'site'] = results['chr'].astype(str) + '_' + results['coords'].astype(str)

        # Count the number of unique sites, genes, and intergenic regions for each spacer
        site_counts = results.groupby('spacer')['site'].nunique()
        gene_counts = results.loc[results['locus_tag'].notnull(), 'spacer'].value_counts()
        intergenic_counts = results.loc[results['locus_tag'].isnull() & results['target'].notnull(), 'spacer'].value_counts()

        # Combine the counts into a DataFrame
        note = pd.DataFrame({
            'sites': site_counts,
            'genes': gene_counts,
            'intergenic': intergenic_counts
        })
        
        console.log("Annotating results...")

        # Replace NaN values with 0 and convert the counts to integers
        note = note.fillna(0).astype(int)

        # Create a 'note' column by applying the create_note function to each row
        note['note'] = note.apply(create_note, axis=1)

        # Merge the 'note' DataFrame with the 'results' DataFrame based on the 'spacer' column
        results = results.merge(note, left_on='spacer', right_index=True, how='left')
        
        column_order = ['name', 'spacer', 'pam', 'chr', 'locus_tag', 'target', 'mismatches', 'coords', 'offset', 'overlap', 'sp_dir', 'tar_dir', 'note']

        # Ensure all columns in column_order exist in the DataFrame, fill with None if absent
        for column in column_order:
            if column not in results.columns:
                results[column] = None

        # Reorder the DataFrame columns according to column_order
        final_results = results[column_order]

        integer_cols = ['mismatches', 'offset', 'overlap']
        
        final_results.loc[:, integer_cols] = final_results[integer_cols].astype('Int64')


    except FileNotFoundError:
        console.log(f"[bold red]Trouble with Bowtie aligner. User a lower number of mismatches.[/bold red]")
        sys.exit(1)

    except KeyError as e:
        console.log(f"{e}, [bold red]None of the proposed barcodes map to the genomes.[/bold red]")
        sys.exit(1)

    console.log(f"Cleaning up...")


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

    unique_chrs = results['chr'].nunique()
    combined_table.add_row("Chromosomes Targeted", f"[bold]{unique_chrs:,}[/bold]")

    unique_tags = results['locus_tag'].nunique()
    combined_table.add_row("Genes Targeted", f"[bold]{unique_tags:,}[/bold]")

    # spacers targeting overlapping genes
    overlapping_genes_targeted = results.loc[results['genes'] > 1, 'locus_tag'].nunique()
    combined_table.add_row("Overlapping Genes Targeted", f"[bold]{overlapping_genes_targeted:,}[/bold]")

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

    final_results.to_csv(sys.stdout, sep='\t', index=False, na_rep='None')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map barcodes to a circular genome")
    parser.add_argument("sgrna_file", help="Path to sgrna_fasta_file", type=str)
    parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
    parser.add_argument("pam", help="PAM sequence", type=str)
    parser.add_argument("mismatches", help="Number of allowed mismatches", type=int)

    args = parser.parse_args()

    main(args)