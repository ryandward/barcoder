import argparse
from datetime import datetime
from multiprocessing import cpu_count

import gzip
import os
import json
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


def create_topological_fasta(
    genbank_file_name, topological_fasta_file_name, overhang_length=0
):
    topological_records = []
    open_func = gzip.open if genbank_file_name.endswith(".gz") else open

    with open_func(genbank_file_name, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            if record.annotations.get("topology", None) == "circular":
                overhang_length = 100_000

            new_seq = record.seq + record.seq[:overhang_length]
            topological_record = SeqRecord(
                Seq(str(new_seq)), id=record.id, description=record.description
            )
            topological_records.append(topological_record)

    with open(topological_fasta_file_name, "w") as output_handle:
        SeqIO.write(topological_records, output_handle, "fasta")


def create_fake_topological_fastq(topological_fasta_file_name, fastq_file):
    if topological_fasta_file_name.endswith(".gz"):
        with gzip.open(topological_fasta_file_name, "rt") as input_handle, open(
            fastq_file, "w"
        ) as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")
    else:
        with open(topological_fasta_file_name, "rt") as input_handle, open(
            fastq_file, "w"
        ) as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")


def create_locus_map(genbank_file_name):
    locus_map, overhang_continue, organisms, seq_lens, topologies, all_genes = (
        {},
        {},
        {},
        {},
        {},
        {},
    )

    with open(genbank_file_name, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            overhang_continue[record.id] = 0
            organisms[record.id] = record.annotations.get("organism", None)
            seq_lens[record.id] = len(record.seq)
            topologies[record.id] = record.annotations.get("topology", None)

            gene_count = 0

            overhang_length = 100_000 if topologies[record.id] == "circular" else 0

            for feature in record.features:
                if feature.type == "gene":
                    gene_count += 1
                    locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                    gene_name = feature.qualifiers.get("gene", [None])[0]

                    if isinstance(feature.location, CompoundLocation) and any(
                        part.start == 0 or part.end == len(record.seq)
                        for part in feature.location.parts
                    ):

                        # Find segments that wrap around the genome
                        genome_end_segment = next(
                            part
                            for part in feature.location.parts
                            if part.end == len(record.seq)
                        )
                        genome_start_segment = next(
                            part for part in feature.location.parts if part.start == 0
                        )

                        # Adjust positions
                        adj_start = int(genome_end_segment.start)
                        adj_end = int(genome_start_segment.end + len(record.seq))
                        overhang_continue[record.id] = int(genome_start_segment.end)

                        # Iterate over the adjusted range
                        for position in range(adj_start, adj_end):
                            key = (record.id, position)
                            locus_map.setdefault(key, []).append(
                                (
                                    locus_tag,
                                    gene_name,
                                    int(adj_start),
                                    int(adj_end),
                                    feature.location.strand,
                                )
                            )

                    else:
                        # normal genes
                        # Check the type of feature.location before the loop
                        if isinstance(feature.location, CompoundLocation):
                            locations = feature.location.parts
                        else:
                            locations = [feature.location]

                        # Now you can iterate over locations without creating a new list in each iteration
                        for part_location in locations:
                            # rest of your code
                            for position in range(
                                int(part_location.start), int(part_location.end)
                            ):
                                key = (record.id, position)
                                locus_map.setdefault(key, []).append(
                                    (
                                        locus_tag,
                                        gene_name,
                                        int(part_location.start),
                                        int(part_location.end),
                                        feature.location.strand,
                                    )
                                )

                                # add the rest of the overhang
                                if (
                                    overhang_continue[record.id]
                                    <= position
                                    < overhang_length
                                ):
                                    key = (record.id, position + len(record.seq))
                                    locus_map.setdefault(key, []).append(
                                        (
                                            locus_tag,
                                            gene_name,
                                            int(part_location.start) + len(record.seq),
                                            int(part_location.end) + len(record.seq),
                                            feature.location.strand,
                                        )
                                    )

                all_genes[record.id] = gene_count

    return locus_map, organisms, seq_lens, topologies, all_genes


# Reconstruct the CRISPRtTarget sequence using the cigar string
def reconstruct_CRISPRtTarget(read):
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


# Get the differences between the spacer and CRISPRtTarget
def get_diff(spacer, CRISPRtTarget):
    differences = []

    for i, (CRISPRtTarget_nt, spacer_nt) in enumerate(zip(CRISPRtTarget, spacer)):
        if CRISPRtTarget_nt != spacer_nt:
            diff_string = f"{CRISPRtTarget_nt}{i + 1}{spacer_nt}"

            differences.append(diff_string)

    diff_result = ",".join(differences)

    if not diff_result:
        return None

    return diff_result


# Get the coordinates of the CRISPRtTarget, accounting for circularity
def get_coords(targStart, targEnd, chrom_length):
    start_circular = targStart % chrom_length
    end_circular = (
        targEnd % chrom_length if targEnd % chrom_length != 0 else chrom_length
    )

    if start_circular > end_circular:
        return f"({start_circular}..{chrom_length}, 0..{end_circular})"
    return f"{start_circular}..{end_circular}"


# Get the offset of the CRISPRtTarget from the feature
# def get_offset(CRISPRtTargetDir, targStart, targEnd, featStart, featEnd):
#     if CRISPRtTargetDir == "F":
#         return targStart - featStart
#     elif CRISPRtTargetDir == "R":
#         return featEnd - targEnd
#     else:
#         return None


# Get the overlap of the CRISPRtTarget and feature
def get_overlap(targStart, targEnd, featStart, featEnd):
    overlapStart = max(targStart, featStart)
    overlapEnd = min(targEnd, featEnd)

    # Check if there's any overlap
    if overlapStart < overlapEnd:
        return overlapEnd - overlapStart
    else:
        return 0


# Check if the extracted PAM matches the PAM pattern
def pam_matches(PAMPattern, extractedPam):
    # Convert N to . for regex matching
    if extractedPam is None:
        return False

    if PAMPattern == "N" * len(PAMPattern) or not PAMPattern:
        return True

    regexPattern = PAMPattern.replace("N", ".")
    return bool(re.match(regexPattern, extractedPam))


# Extract the downstream PAM of the CRISPRtTarget sequence from the topological map
def extract_downstream_pam(
    pam,
    targStart,
    targEnd,
    chrom,
    fasta,
    dir,
    true_chrom_lengths,
    topological_chrom_lengths,
):
    true_chrom_length = true_chrom_lengths.get(chrom, None)
    topological_chrom_length = topological_chrom_lengths.get(chrom, None)

    if pam == "":
        return None

    if None in (
        pam,
        targStart,
        targEnd,
        chrom,
        fasta,
        dir,
        true_chrom_length,
        topological_chrom_length,
    ):
        return None

    if dir == "F":
        if targEnd + len(pam) > topological_chrom_length:
            return None
        extractedPam = fasta.fetch(
            reference=chrom, start=targEnd, end=targEnd + len(pam)
        ).upper()

    elif dir == "R":
        if targStart - len(pam) < 0:
            return None
        extractedPam = fasta.fetch(
            reference=chrom, start=targStart - len(pam), end=targStart
        ).upper()
        extractedPam = str(Seq(extractedPam).reverse_complement())

    else:
        return None

    return extractedPam


# Extract the upstream PAM of the CRISPRtTarget sequence from the topological map
def extract_upstream_pam(
    pam,
    targStart,
    targEnd,
    chrom,
    fasta,
    dir,
    true_chrom_lengths,
    topological_chrom_lengths,
):
    true_chrom_length = true_chrom_lengths.get(chrom, None)
    topological_chrom_length = topological_chrom_lengths.get(chrom, None)

    if pam == "":
        return None

    if None in (
        pam,
        targStart,
        targEnd,
        chrom,
        fasta,
        dir,
        true_chrom_length,
        topological_chrom_length,
    ):
        return None

    if dir == "F":
        if targStart - len(pam) < 0:
            return None
        extractedPam = fasta.fetch(
            reference=chrom, start=targStart - len(pam), end=targStart
        ).upper()

    elif dir == "R":
        if targEnd + len(pam) > topological_chrom_length:
            return None
        extractedPam = fasta.fetch(
            reference=chrom, start=targEnd, end=targEnd + len(pam)
        ).upper()
        extractedPam = str(Seq(extractedPam).reverse_complement())

    else:
        return None

    return extractedPam


# Parse the SAM file and extract the relevant information
def parse_sam_output(
    samfile, locus_map, topological_fasta_file_name, gb_file_name, pam, pam_direction
):
    """
    Parses the SAM output file and extracts relevant information for each read.

    Args:
        samfile (pysam.AlignmentFile): The SAM file object.
        locus_map (dict): A dictionary mapping chromosome positions to gene information.
        topological_fasta_file_name (str): The filename of the topological FASTA file.
        gb_file_name (str): The filename of the GenBank file.
        pam (str): The PAM sequence.
        pam_direction (str): The direction of the PAM sequence relative to the read.

    Returns:
        list: A list of dictionaries, where each dictionary represents a read and contains
        information such as read name, spacer sequence, length, CRISPRtTarget sequence, mismatches,
        chromosome, CRISPRtTarget start and end positions, spacer direction, PAM sequence, coordinates,
        type of alignment, and difference between spacer and CRISPRtTarget.

    Raises:
        ValueError: If there is an error in retrieving reference sequence or calculating
        CRISPRtTarget start and end positions.
    """
    rows_list = []  # Initialize an empty list
    true_chrom_lengths = get_true_chrom_lengths(gb_file_name)
    topological_chrom_lengths = get_topological_chrom_lengths(
        topological_fasta_file_name
    )

    with pysam.FastaFile(topological_fasta_file_name) as fasta:
        for read in samfile.fetch():
            if read.query_sequence is None:
                continue

            extractedPam = None

            if pam:
                if pam_direction == "downstream":
                    extractedPam = extract_downstream_pam(
                        pam,
                        read.reference_start,
                        read.reference_end,
                        read.reference_name,
                        fasta,
                        "F" if not read.is_reverse else "R",
                        true_chrom_lengths,
                        topological_chrom_lengths,
                    )
                elif pam_direction == "upstream":
                    extractedPam = extract_upstream_pam(
                        pam,
                        read.reference_start,
                        read.reference_end,
                        read.reference_name,
                        fasta,
                        "F" if not read.is_reverse else "R",
                        true_chrom_lengths,
                        topological_chrom_lengths,
                    )

                if read.is_mapped:
                    if extractedPam is None:
                        continue
                    # exclude reads that don't match the PAM
                    if not pam_matches(pam, extractedPam):
                        read.is_unmapped = True
                        extractedPam = None

            rows = {
                "name": read.query_name,
                "spacer": (
                    read.query_sequence
                    if not read.is_reverse
                    else str(Seq(read.query_sequence).reverse_complement())
                ),
                "len": len(read.query_sequence),
            }

            if read.is_unmapped:
                rows_list.append(rows)
                continue

            if (
                read.is_mapped
                and read.reference_start is not None
                and read.reference_end is not None
            ):

                try:
                    CRISPRtTarget = (
                        read.get_reference_sequence()
                        if not read.is_reverse
                        else str(
                            Seq(read.get_reference_sequence()).reverse_complement()
                        )
                    )
                    chrom = read.reference_name

                    targStart = read.reference_start % true_chrom_lengths.get(chrom, None)
                    targEnd = read.reference_end % true_chrom_lengths.get(chrom, None)

                    # Adjust targStart if it spans the origin
                    if targEnd < targStart:
                        targStart -= true_chrom_lengths.get(chrom, None)

                    insDirection = "F" if not read.is_reverse else "R"
                    coords = get_coords(
                        targStart, targEnd, true_chrom_lengths.get(chrom, None)
                    )
                    
                    if insDirection == "F":
                    # the actual insert is  49 nucleotides downstream of  the end
                        insSite = (targEnd + 49) % true_chrom_lengths.get(chrom, None)
                    elif insDirection == "R":
                        insSite = (targStart - 49) % true_chrom_lengths.get(chrom, None)

                except ValueError as e:
                    print(e, file=sys.stderr)
                    sys.exit(1)

                rows.update(
                    {
                        "CRISPRtTarget": CRISPRtTarget,
                        "targStart": targStart,
                        "targEnd": targEnd,
                        "coords": coords,
                        "chrom": chrom,
                        "insSite": insSite,
                        "insDirection": insDirection,
                        "pam": extractedPam,
                    }
                )

                aligned_genes = {
                    gene_info
                    for pos in range(targStart, targEnd)
                    for gene_info in locus_map.get((chrom, pos), [])
                }

                if not aligned_genes:
                    rows.update(
                        {
                            "locus_tag": None,
                            "offset": None,
                            "overlap": None,
                            "targDir": None,
                        }
                    )

                    rows_list.append(rows)

                else:
                    for (
                        locus_tag,
                        gene_name,
                        featStart,
                        featEnd,
                        feature_strand,
                    ) in aligned_genes:
                        CRISPRtTarget_orientation = (
                            "F"
                            if feature_strand == 1
                            else "R" if feature_strand == -1 else None
                        )
                        # offset = get_offset(
                        #     CRISPRtTarget_orientation,
                        #     targStart,
                        #     targEnd,
                        #     featStart,
                        #     featEnd,
                        # )
                        overlap = get_overlap(
                            targStart, targEnd, featStart, featEnd
                        )

                        rows_copy = rows.copy()
                        rows_copy.update(
                            {
                                "locus_tag": locus_tag,
                                "gene": gene_name if gene_name else locus_tag,
                                # "offset": offset,
                                "overlap": overlap,
                                "targDir": CRISPRtTarget_orientation,
                            }
                        )

                        rows_list.append(rows_copy)

    return rows_list


# Run bowtie and parse the output using parse_sam_output and temp files
def run_bowtie_and_parse(
    sgrna_fastq_file_name,
    topological_fasta_file_name,
    locus_map,
    num_mismatches,
    num_threads,
):
    results = []
    # Create a temporary directory for the bowtie index files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a temporary name for the genome index
        genome_index_temp_name = tempfile.NamedTemporaryFile(
            dir=temp_dir, delete=False
        ).name
        index_prefix = os.path.join(temp_dir, genome_index_temp_name)

        with open(os.devnull, "w") as devnull:

            bowtie_build_command = [
                "bowtie-build",
                topological_fasta_file_name,
                index_prefix,
            ]

            bowtie_build_process = subprocess.Popen(
                bowtie_build_command, stdout=devnull, stderr=devnull
            )
            bowtie_build_process.wait()

            # Create a temporary file for the bowtie output
            with tempfile.NamedTemporaryFile(delete=False) as bowtie_output_temp_file:

                bowtie_command = [
                    "bowtie",
                    "-S",
                    "-k 100",
                    "--nomaqround",
                    "-p",
                    str(num_threads),
                    "--tryhard",
                    "-v",
                    str(num_mismatches),
                    "-x",
                    index_prefix,
                    sgrna_fastq_file_name,
                ]

                bowtie_process = subprocess.Popen(
                    bowtie_command,
                    stdout=subprocess.PIPE,
                    stderr=devnull,
                    universal_newlines=True,
                )

        if bowtie_process.stdout is None:
            raise RuntimeError("Bowtie was unable to start. Check your installation.")

        if bowtie_process.stdout is not None:
            with pysam.AlignmentFile(bowtie_process.stdout, "r") as samfile:
                results = parse_sam_output(
                    samfile,
                    locus_map,
                    topological_fasta_file_name,
                    args.genome_file,
                    args.pam,
                    args.pam_direction,
                )

            bowtie_process.wait()

            # Delete the temporary file after we're done with it
            os.remove(bowtie_output_temp_file.name)

    # Return the results of parse_sam_output
    if results is []:
        raise RuntimeError(
            "No results were returned from the Bowtie process. Check your input files and parameters."
        )
    else:
        return results


# Filter out spacers that don't match the PAM
def filter_offCRISPRtTargets_by_pam(df):
    CRISPRtTargeting_spacers = df[df["CRISPRtTarget"].notna()]["spacer"].unique()
    return df[~((df["CRISPRtTarget"].isna()) & (df["spacer"].isin(CRISPRtTargeting_spacers)))]


# Create a note for each spacer
def create_note(row):
    parts = []
    if row["sites"] > 0:
        parts.append(f"{row['sites']} {'site' if row['sites'] == 1 else 'sites'}")
        if row["genes"] > 0:
            parts.append(f"{row['genes']} {'gene' if row['genes'] == 1 else 'genes'}")
        if row["intergenic"] > 0:
            parts.append(f"{row['intergenic']} intergenic")
    else:
        parts.append("non-CRISPRtTargeting")
    return ", ".join(parts)


# Main function
def main(args):
    console = Console(file=sys.stderr)
    console.log("[bold red]Initializing barcode CRISPRtTarget seeker[/bold red]")

    num_threads = cpu_count() // 2

    with tempfile.TemporaryDirectory() as working_dir:
        topological_fasta_file_name = os.path.join(
            working_dir,
            os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta",
        )

        console.log("Annotating regions to identify...")

        base_name = os.path.basename(args.sgrna_file)
        if ".fastq" in base_name:
            sgrna_fastq_file_name = args.sgrna_file
        elif base_name.endswith(".fasta") or base_name.endswith(".fasta.gz"):
            ext = ".fasta.gz" if base_name.endswith(".fasta.gz") else ".fasta"
            base_name = base_name[: -len(ext)]
            sgrna_fastq_file_name = os.path.join(working_dir, base_name + ".fastq")
            create_fake_topological_fastq(args.sgrna_file, sgrna_fastq_file_name)
        else:
            console.log(f"[bold red]File extension not recognized. [/bold red]")
            sys.exit(1)

        console.log("Generating topological coordinate maps...")
        locus_map, organisms, seq_lens, topologies, all_genes = create_locus_map(
            args.genome_file
        )

        create_topological_fasta(args.genome_file, topological_fasta_file_name)

        console.log("Aligning annotations to genome...")
        results = run_bowtie_and_parse(
            sgrna_fastq_file_name,
            topological_fasta_file_name,
            locus_map,
            args.mismatches,
            num_threads,
        )

        with pysam.FastaFile(topological_fasta_file_name) as _:
            pass

    try:
        console.log("Finding matches...")

        # Convert your data to a DataFrame, then drop duplicates, which are caused by overhangs
        results = pd.DataFrame(results).drop_duplicates()

        # Use the filter_offCRISPRtTargets_by_pam function
        try:
            results = filter_offCRISPRtTargets_by_pam(results)
        except KeyError as e:
            console.log(
                f"[bold red]The following critical attribute is missing for every barcode:[/bold red] {e}"
            )
            console.log(
                f"[bold yellow]Are you sure your PAM should be[/bold yellow] '{args.pam}' '{args.pam_direction}' [bold yellow]of the spacer?[/bold yellow]"
            )
            console.log(
                f"[bold yellow]The genome file[/bold yellow] '{args.genome_file}' [bold yellow]has definitions for these chromosomes:[/bold yellow]"
            )
            console.log(json.dumps(organisms, indent=4))
            sys.exit(1)

        # Create a 'min_tar' column that is the minimum of 'targStart' and 'targEnd'
        def adjust_min_tar(row):
            if row["targStart"] > row["targEnd"]:
                return row["targStart"] - seq_lens[row["chrom"]]
            else:
                return row["targStart"]

        results["min_tar"] = results.apply(adjust_min_tar, axis=1)

        # Sort the DataFrame by 'chr', 'min_tar', and 'spacer'
        results = results.sort_values(by=["chrom", "min_tar", "spacer"])

        # count the number of times each spacer was seen after grouping by 'name' and 'spacer'
        spacers_seen = (
            results[["name", "spacer"]].drop_duplicates().groupby("spacer").size()
        )

        # now drop the name column
        results = results.drop("name", axis=1).drop_duplicates()

        # Create a 'site' column only for rows that have a 'CRISPRtTarget'
        results.loc[results["CRISPRtTarget"].notnull(), "site"] = (
            results["chrom"].astype(str) + "_" + results["coords"].astype(str)
        )

        # Count the number of unique sites, genes, and intergenic regions for each spacer
        site_counts = results.groupby("spacer")["site"].nunique()
        gene_counts = results.loc[
            results["locus_tag"].notnull(), "spacer"
        ].value_counts()
        intergenic_counts = results.loc[
            results["locus_tag"].isnull() & results["CRISPRtTarget"].notnull(), "spacer"
        ].value_counts()

        # get the lengths of the spacers into a set
        spacer_lengths = set(results["len"])

        # convert into a string, separated by commas of all the lengths
        if len(spacer_lengths) == 1:
            spacer_len_range = str(next(iter(spacer_lengths)))
        else:
            spacer_len_range = ",".join(
                str(spacer_len) for spacer_len in sorted(spacer_lengths)
            )

        # Combine the counts into a DataFrame
        note = pd.DataFrame(
            {
                "count": spacers_seen,
                "sites": site_counts,
                "genes": gene_counts,
                "intergenic": intergenic_counts,
            }
        )

        console.log("Annotating results...")

        # Replace NaN values with 0 and convert the counts to integers
        note = note.fillna(0).astype(int)

        # Create a 'note' column by applying the create_note function to each row
        note["note"] = note.apply(create_note, axis=1)

        # Merge the 'note' DataFrame with the 'results' DataFrame based on the 'spacer' column
        results = results.merge(note, left_on="spacer", right_index=True, how="left")

        column_order = ["spacer", "locus_tag", "gene", "chrom"]

        if not (results["count"] == 1).all():
            column_order.append("count")

        if not (results["pam"].isnull().all() or results["pam"].nunique() == 1):
            column_order.append("pam")

        if not (results["mismatches"] == 0).all():
            column_order.append("mismatches")

        column_order.extend(
            [
                "CRISPRtTarget",
                "targStart",
                "targEnd",
                "offset",
                "overlap",
                "insDirection",
                "insSite",
                "targDir",
                "note",
            ]
        )

        # Reorder the DataFrame columns according to column_order
        final_results = results.reindex(columns=column_order)

        integer_cols = ["mismatches", "offset", "overlap", "targStart", "targEnd"]

        # Convert the columns in integer_cols to 'Int64', which supports NaN values
        for col in integer_cols:
            if col in final_results.columns:
                final_results[col] = final_results[col].astype("Int64")

    except FileNotFoundError:
        console.log(
            f"[bold red]Trouble with Bowtie aligner. Try using a lower number of mismatches.[/bold red]"
        )
        sys.exit(1)

    except KeyError as e:
        console.log(
            f"[bold red]All of the proposed barcodes are missing some key attributes[/bold red]: {e}"
        )
        pass
        # sys.exit(1)

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
        show_header=True,
    )
    # Define columns with justifications
    combined_table.add_column(
        os.path.basename(sys.argv[0]), justify="right", style="white", min_width=30
    )
    combined_table.add_column(
        "Summary", justify="right", style="bold bright_white", min_width=20
    )

    # Input & Configuration Sub-heading
    combined_table.add_section()
    combined_table.add_row(
        "[bold bright_magenta]Input & Config[/bold bright_magenta]", ""
    )

    # Rows for Input & Configuration
    combined_table.add_row(
        "Barcodes", f"[bold]{os.path.basename(args.sgrna_file)}[/bold]"
    )
    combined_table.add_row(
        "Genbank Genome File", f"[bold]{os.path.basename(args.genome_file)}[/bold]"
    )

    combined_table.add_row("pam", f"[bold]{args.pam}[/bold]")

    if args.pam_direction == "downstream":
        combined_table.add_row("PAM Direction", f"[bold]Downstream[/bold]")
    elif args.pam_direction == "upstream":
        combined_table.add_row("PAM Direction", f"[bold]Upstream[/bold]")

    combined_table.add_row("Number of Mismatches", f"[bold]{args.mismatches}[/bold]")
    combined_table.add_row("Threads", f"[bold]{num_threads}[/bold]")
    combined_table.add_row("Operating System", f"[bold]{platform.system()}[/bold]")

    # Heuristic Statistics Sub-heading
    combined_table.add_section()
    combined_table.add_row("[bold bright_blue]Heuristics[/bold bright_blue]", "")

    # barcode lengths
    combined_table.add_row("Spacer Lengths", f"[bold]{spacer_len_range}[/bold]")

    systematic_name = None

    if args.pam_direction == "downstream":
        # call it spacer_range + PAM
        systematic_name = f"{spacer_len_range}-{args.pam}"
    elif args.pam_direction == "upstream":
        # call it PAM + spacer_range
        systematic_name = f"{args.pam}-{spacer_len_range}"

    if systematic_name:
        combined_table.add_row("Systematic Name", f"[bold]{systematic_name}[/bold]")

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
    combined_table.add_row(
        "Sequence Length",
        f"[bold]{'; '.join(format(seq_len, ',') for seq_len in unique_seq_lens)}[/bold]",
    )

    ambiguous_coordinates = {
        (chrom, pos % seq_lens[chrom])
        for chrom, pos in locus_map
        if len(locus_map[(chrom, pos)]) > 1
    }
    ambiguous_locus_tags = {
        entry[0]
        for chrom, pos in ambiguous_coordinates
        for entry in locus_map[(chrom, pos)]
    }

    # number of chromosomes from seq_lens

    combined_table.add_row("Chromosomes", f"[bold]{len(seq_lens)}[/bold]")
    combined_table.add_row("Total Genes", f"[bold]{sum(all_genes.values()):,}[/bold]")

    combined_table.add_row(
        "Overlapping Genes", f"[bold]{len(ambiguous_locus_tags):,}[/bold]"
    )
    combined_table.add_row(
        "Ambiguous Coordinates", f"[bold]{len(ambiguous_coordinates):,}[/bold]"
    )

    # Barcode Mapping Stats Sub-heading
    combined_table.add_section()
    combined_table.add_row(
        "[bold bright_green]Barcode Mapping Stats[/bold bright_green]", ""
    )

    unique_chrs = results["chrom"].nunique()
    combined_table.add_row("Chromosomes CRISPRtTargeted", f"[bold]{unique_chrs:,}[/bold]")

    unique_tags = results["locus_tag"].nunique()
    combined_table.add_row("Genes CRISPRtTargeted", f"[bold]{unique_tags:,}[/bold]")

    # spacers CRISPRtTargeting overlapping genes
    overlapping_genes_CRISPRtTargeted = results.loc[
        results["genes"] > 1, "locus_tag"
    ].nunique()
    combined_table.add_row(
        "Overlapping Genes CRISPRtTargeted", f"[bold]{overlapping_genes_CRISPRtTargeted:,}[/bold]"
    )

    unique_barcodes = results["spacer"].nunique()

    combined_table.add_row("Unique Barcodes", f"[bold]{unique_barcodes:,}[/bold]")

    if "mismatches" in final_results.columns:
        unique_spacers_per_mismatch = final_results.groupby(["mismatches"])[
            "spacer"
        ].nunique()

        for mismatch, count in unique_spacers_per_mismatch.items():
            combined_table.add_row(
                f"{mismatch} Mismatch Barcodes", f"[bold]{count:,}[/bold]"
            )

    intergenic_spacers = results[
        (results["locus_tag"].isnull()) & (results["chrom"].notnull())
    ]["spacer"].nunique()
    combined_table.add_row(
        "Intergenic Barcodes", f"[bold]{intergenic_spacers:,}[/bold]"
    )

    # spacers CRISPRtTargeting multiple coordinates
    off_CRISPRtTarget_spacers = (
        results[results["CRISPRtTarget"].notnull()]
        .groupby("spacer")["coords"]
        .apply(set)
        .apply(len)
        .gt(1)
        .sum()
    )
    combined_table.add_row(
        "Off-CRISPRtTargeting Barcodes", f"[bold]{off_CRISPRtTarget_spacers:,}[/bold]"
    )

    nonCRISPRtTargeting_spacers = results[results["CRISPRtTarget"].isnull()]["spacer"].nunique()
    combined_table.add_row(
        "Non-CRISPRtTargeting Barcodes", f"[bold]{nonCRISPRtTargeting_spacers:,}[/bold]"
    )

    # Print the combined table
    console.log(combined_table)

    final_results.to_csv(sys.stdout, sep="\t", index=False, na_rep="None")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map barcodes to a circular genome")
    parser.add_argument("sgrna_file", help="Path to sgrna_fasta_file", type=str)
    parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
    parser.add_argument("pam", help="PAM sequence", type=str)
    parser.add_argument("mismatches", help="Number of allowed mismatches", type=int)
    parser.add_argument(
        "--pam_direction",
        choices=["upstream", "downstream"],
        default="downstream",
        help="Direction of the PAM sequence",
    )

    args = parser.parse_args()

    main(args)
