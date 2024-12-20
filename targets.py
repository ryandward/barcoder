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
                Seq(str(new_seq)),
                id=record.id,
                description=record.description,
                name=args.genome_file,
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
                        genome_end_segment = next(
                            part
                            for part in feature.location.parts
                            if part.end == len(record.seq)
                        )
                        genome_start_segment = next(
                            part for part in feature.location.parts if part.start == 0
                        )
                        adj_start = int(genome_end_segment.start)
                        adj_end = int(genome_start_segment.end + len(record.seq))
                        overhang_continue[record.id] = int(genome_start_segment.end)

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
                        locations = (
                            feature.location.parts
                            if isinstance(feature.location, CompoundLocation)
                            else [feature.location]
                        )
                        for part_location in locations:
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


def get_diff(spacer, target):
    differences = [
        f"{target_nt}{i + 1}{spacer_nt}"
        for i, (target_nt, spacer_nt) in enumerate(zip(target, spacer))
        if target_nt != spacer_nt
    ]
    return ",".join(differences) if differences else None


def get_coords(tar_start, tar_end, chrom_length):
    start_circular = tar_start % chrom_length
    end_circular = (
        tar_end % chrom_length if tar_end % chrom_length != 0 else chrom_length
    )
    return (
        f"({start_circular}..{chrom_length}, 0..{end_circular})"
        if start_circular > end_circular
        else f"{start_circular}..{end_circular}"
    )


def get_offset(target_dir, tar_start, tar_end, feature_start, feature_end):
    if target_dir == "F":
        return tar_start - feature_start
    if target_dir == "R":
        return feature_end - tar_end
    return None


def get_overlap(tar_start, tar_end, feature_start, feature_end):
    overlap_start = max(tar_start, feature_start)
    overlap_end = min(tar_end, feature_end)
    return overlap_end - overlap_start if overlap_start < overlap_end else 0


def pam_matches(pam_pattern, extracted_pam):
    if not extracted_pam:
        return False
    if pam_pattern == "N" * len(pam_pattern) or not pam_pattern:
        return True
    return bool(re.match(pam_pattern.replace("N", "."), extracted_pam))


def extract_downstream_pam(
    pam,
    tar_start,
    tar_end,
    chrom,
    fasta,
    dir,
    true_chrom_lengths,
    topological_chrom_lengths,
):
    true_chrom_length = true_chrom_lengths.get(chrom, None)
    topological_chrom_length = topological_chrom_lengths.get(chrom, None)

    if not pam or None in (
        tar_start,
        tar_end,
        chrom,
        fasta,
        dir,
        true_chrom_length,
        topological_chrom_length,
    ):
        return None

    if dir == "F":
        if tar_end + len(pam) > topological_chrom_length:
            return None
        extracted_pam = fasta.fetch(
            reference=chrom, start=tar_end, end=tar_end + len(pam)
        ).upper()
    elif dir == "R":
        if tar_start - len(pam) < 0:
            return None
        extracted_pam = fasta.fetch(
            reference=chrom, start=tar_start - len(pam), end=tar_start
        ).upper()
        extracted_pam = str(Seq(extracted_pam).reverse_complement())
    else:
        return None
    return extracted_pam


def extract_upstream_pam(
    pam,
    tar_start,
    tar_end,
    chrom,
    fasta,
    dir,
    true_chrom_lengths,
    topological_chrom_lengths,
):
    true_chrom_length = true_chrom_lengths.get(chrom, None)
    topological_chrom_length = topological_chrom_lengths.get(chrom, None)
    if not pam or None in (
        tar_start,
        tar_end,
        chrom,
        fasta,
        dir,
        true_chrom_length,
        topological_chrom_length,
    ):
        return None

    if dir == "F":
        if tar_start - len(pam) < 0:
            return None
        extracted_pam = fasta.fetch(
            reference=chrom, start=tar_start - len(pam), end=tar_start
        ).upper()
    elif dir == "R":
        if tar_end + len(pam) > topological_chrom_length:
            return None
        extracted_pam = fasta.fetch(
            reference=chrom, start=tar_end, end=tar_end + len(pam)
        ).upper()
        extracted_pam = str(Seq(extracted_pam).reverse_complement())
    else:
        return None
    return extracted_pam


def parse_sam_output(
    samfile, locus_map, topological_fasta_file_name, gb_file_name, pam, pam_direction
):
    rows_list = []
    true_chrom_lengths = get_true_chrom_lengths(gb_file_name)
    topological_chrom_lengths = get_topological_chrom_lengths(
        topological_fasta_file_name
    )

    with pysam.FastaFile(topological_fasta_file_name) as fasta:
        for read in samfile.fetch():
            if read.query_sequence is None:
                continue

            extracted_pam = None

            if pam:
                if pam_direction == "downstream":
                    extracted_pam = extract_downstream_pam(
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
                    extracted_pam = extract_upstream_pam(
                        pam,
                        read.reference_start,
                        read.reference_end,
                        read.reference_name,
                        fasta,
                        "F" if not read.is_reverse else "R",
                        true_chrom_lengths,
                        topological_chrom_lengths,
                    )

                if read.is_mapped and not pam_matches(pam, extracted_pam):
                    read.is_unmapped = True
                    extracted_pam = None

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

            if read.is_mapped:
                try:
                    target = (
                        read.get_reference_sequence()
                        if not read.is_reverse
                        else str(
                            Seq(read.get_reference_sequence()).reverse_complement()
                        )
                    )
                    mismatches = int(read.get_tag("NM"))
                    chr = read.reference_name

                    tar_start = read.reference_start % true_chrom_lengths.get(chr, None)
                    tar_end = read.reference_end % true_chrom_lengths.get(chr, None)

                    if tar_end < tar_start:
                        tar_start -= true_chrom_lengths.get(chr, None)

                    sp_dir = "F" if not read.is_reverse else "R"
                    coords = get_coords(
                        tar_start, tar_end, true_chrom_lengths.get(chr, None)
                    )
                    type = "mismatch" if mismatches > 0 else "perfect"
                    diff = get_diff(rows["spacer"], target)

                except ValueError as e:
                    print(e, file=sys.stderr)
                    sys.exit(1)

                rows.update(
                    {
                        "target": target,
                        "mismatches": mismatches,
                        "chr": chr,
                        "tar_start": tar_start,
                        "tar_end": tar_end,
                        "sp_dir": sp_dir,
                        "pam": extracted_pam,
                        "coords": coords,
                        "type": type,
                        "diff": diff,
                    }
                )

                aligned_genes = {
                    gene_info
                    for pos in range(tar_start, tar_end)
                    for gene_info in locus_map.get((chr, pos), [])
                }

                if not aligned_genes:
                    rows.update(
                        {
                            "locus_tag": None,
                            "offset": None,
                            "overlap": None,
                            "tar_dir": None,
                        }
                    )
                    rows_list.append(rows)
                else:
                    for (
                        locus_tag,
                        gene_name,
                        feature_start,
                        feature_end,
                        feature_strand,
                    ) in aligned_genes:
                        target_orientation = (
                            "F"
                            if feature_strand == 1
                            else "R" if feature_strand == -1 else None
                        )
                        offset = get_offset(
                            target_orientation,
                            tar_start,
                            tar_end,
                            feature_start,
                            feature_end,
                        )
                        overlap = get_overlap(
                            tar_start, tar_end, feature_start, feature_end
                        )

                        rows_copy = rows.copy()
                        rows_copy.update(
                            {
                                "locus_tag": locus_tag,
                                "gene": gene_name if gene_name else locus_tag,
                                "offset": offset,
                                "overlap": overlap,
                                "tar_dir": target_orientation,
                            }
                        )
                        rows_list.append(rows_copy)

    return rows_list


def run_bowtie_and_parse(
    sgrna_fastq_file_name,
    topological_fasta_file_name,
    locus_map,
    num_mismatches,
    num_threads,
    pam,
    pam_direction,
):
    results = []
    with tempfile.TemporaryDirectory() as temp_dir:
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
                    pam,
                    pam_direction,
                )

            bowtie_process.wait()
            os.remove(bowtie_output_temp_file.name)

    if not results:
        raise RuntimeError(
            "No results were returned from the Bowtie process. Check your input files and parameters."
        )
    return results


def filter_offtargets_by_pam(df):
    targeting_spacers = df[df["target"].notna()]["spacer"].unique()
    return df[~((df["target"].isna()) & (df["spacer"].isin(targeting_spacers)))]


def create_note(row):
    parts = []
    if row["sites"] > 0:
        parts.append(f"{row['sites']} {'site' if row['sites'] == 1 else 'sites'}")
        if row["genes"] > 0:
            parts.append(f"{row['genes']} {'gene' if row['genes'] == 1 else 'genes'}")
        if row["intergenic"] > 0:
            parts.append(f"{row['intergenic']} intergenic")
    else:
        parts.append("non-targeting")
    return ", ".join(parts)


def main(args):
    console = Console(file=sys.stderr)
    console.log("[bold red]Initializing barcode target seeker[/bold red]")
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
            args.pam,
            args.pam_direction,
        )

        with pysam.FastaFile(topological_fasta_file_name) as _:
            pass

    try:
        console.log("Finding matches...")
        results = pd.DataFrame(results).drop_duplicates()

        try:
            results = filter_offtargets_by_pam(results)
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

        def adjust_min_tar(row):
            if row["tar_start"] > row["tar_end"]:
                return row["tar_start"] - seq_lens[row["chr"]]
            return row["tar_start"]

        results["min_tar"] = results.apply(adjust_min_tar, axis=1)
        results = results.sort_values(by=["chr", "min_tar", "spacer"])

        spacers_seen = (
            results[["name", "spacer"]].drop_duplicates().groupby("spacer").size()
        )
        results = results.drop("name", axis=1).drop_duplicates()
        results.loc[results["target"].notnull(), "site"] = (
            results["chr"].astype(str) + "_" + results["coords"].astype(str)
        )

        site_counts = results.groupby("spacer")["site"].nunique()
        gene_counts = results.loc[
            results["locus_tag"].notnull(), "spacer"
        ].value_counts()
        intergenic_counts = results.loc[
            results["locus_tag"].isnull() & results["target"].notnull(), "spacer"
        ].value_counts()

        spacer_lengths = set(results["len"])
        spacer_len_range = (
            str(next(iter(spacer_lengths)))
            if len(spacer_lengths) == 1
            else ",".join(str(spacer_len) for spacer_len in sorted(spacer_lengths))
        )

        note = pd.DataFrame(
            {
                "count": spacers_seen,
                "sites": site_counts,
                "genes": gene_counts,
                "intergenic": intergenic_counts,
            }
        )

        console.log("Annotating results...")

        note = note.fillna(0).astype(int)
        note["note"] = note.apply(create_note, axis=1)
        results = results.merge(note, left_on="spacer", right_index=True, how="left")

        column_order = ["spacer", "locus_tag", "gene", "chr"]
        if not (results["count"] == 1).all():
            column_order.append("count")
        if not (results["pam"].isnull().all() or results["pam"].nunique() == 1):
            column_order.append("pam")
        if not (results["mismatches"] == 0).all():
            column_order.append("mismatches")
        column_order.extend(
            [
                "target",
                "tar_start",
                "tar_end",
                "offset",
                "overlap",
                "sp_dir",
                "tar_dir",
                "note",
            ]
        )

        final_results = results.reindex(columns=column_order)
        integer_cols = ["mismatches", "offset", "overlap", "tar_start", "tar_end"]
        for col in integer_cols:
            if col in final_results.columns:
                final_results[col] = final_results[col].astype("Int64")

        if args.json:  # ADDED JSON OUTPUT
            console.log("Writing to JSON...")
            print(final_results.to_json(orient="records", indent=4))
        else:  # DEFAULT TSV OUTPUT
            console.log("Writing to TSV...")
            final_results.to_csv(sys.stdout, sep="\t", index=False, na_rep="None")

    except FileNotFoundError:
        console.log(
            f"[bold red]Trouble with Bowtie aligner. Try using a lower number of mismatches.[/bold red]"
        )
        sys.exit(1)
    except KeyError as e:
        console.log(
            f"[bold red]All of the proposed barcodes are missing some key attributes[/bold red]: {e}"
        )
        sys.exit(1)

    console.log(f"Cleaning up...")

    combined_table = Table(
        box=rich.table.box.SIMPLE_HEAVY,
        caption="Finished at [u]{}[/u]".format(datetime.now()),
        title_style="bold bright_white",
        caption_style="bold white",
        header_style="bold bright_white",
        border_style="bold bright_white",
        show_header=True,
    )
    combined_table.add_column(
        os.path.basename(sys.argv[0]), justify="right", style="white", min_width=30
    )
    combined_table.add_column(
        "Summary", justify="right", style="bold bright_white", min_width=20
    )

    combined_table.add_section()
    combined_table.add_row(
        "[bold bright_magenta]Input & Config[/bold bright_magenta]", ""
    )

    combined_table.add_row(
        "Barcodes", f"[bold]{os.path.basename(args.sgrna_file)}[/bold]"
    )
    combined_table.add_row(
        "Genbank Genome File", f"[bold]{os.path.basename(args.genome_file)}[/bold]"
    )

    combined_table.add_row("PAM", f"[bold]{args.pam}[/bold]")

    if args.pam_direction == "downstream":
        combined_table.add_row("PAM Direction", f"[bold]Downstream[/bold]")
    elif args.pam_direction == "upstream":
        combined_table.add_row("PAM Direction", f"[bold]Upstream[/bold]")

    combined_table.add_row("Number of Mismatches", f"[bold]{args.mismatches}[/bold]")
    combined_table.add_row("Threads", f"[bold]{num_threads}[/bold]")
    combined_table.add_row("Operating System", f"[bold]{platform.system()}[/bold]")
    combined_table.add_section()
    combined_table.add_row("[bold bright_blue]Heuristics[/bold bright_blue]", "")
    combined_table.add_row("Spacer Lengths", f"[bold]{spacer_len_range}[/bold]")

    systematic_name = None

    if args.pam_direction == "downstream":
        systematic_name = f"{spacer_len_range}-{args.pam}"
    elif args.pam_direction == "upstream":
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

    combined_table.add_row("Chromosomes", f"[bold]{len(seq_lens)}[/bold]")
    combined_table.add_row("Total Genes", f"[bold]{sum(all_genes.values()):,}[/bold]")
    combined_table.add_row(
        "Overlapping Genes", f"[bold]{len(ambiguous_locus_tags):,}[/bold]"
    )
    combined_table.add_row(
        "Ambiguous Coordinates", f"[bold]{len(ambiguous_coordinates):,}[/bold]"
    )
    combined_table.add_section()
    combined_table.add_row(
        "[bold bright_green]Barcode Mapping Stats[/bold bright_green]", ""
    )

    unique_chrs = results["chr"].nunique()
    combined_table.add_row("Chromosomes Targeted", f"[bold]{unique_chrs:,}[/bold]")

    unique_tags = results["locus_tag"].nunique()
    combined_table.add_row("Genes Targeted", f"[bold]{unique_tags:,}[/bold]")

    overlapping_genes_targeted = results.loc[
        results["genes"] > 1, "locus_tag"
    ].nunique()
    combined_table.add_row(
        "Overlapping Genes Targeted", f"[bold]{overlapping_genes_targeted:,}[/bold]"
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
        (results["locus_tag"].isnull()) & (results["chr"].notnull())
    ]["spacer"].nunique()
    combined_table.add_row(
        "Intergenic Barcodes", f"[bold]{intergenic_spacers:,}[/bold]"
    )

    off_target_spacers = (
        results[results["target"].notnull()]
        .groupby("spacer")["coords"]
        .apply(set)
        .apply(len)
        .gt(1)
        .sum()
    )
    combined_table.add_row(
        "Off-targeting Barcodes", f"[bold]{off_target_spacers:,}[/bold]"
    )

    nontargeting_spacers = results[results["target"].isnull()]["spacer"].nunique()
    combined_table.add_row(
        "Non-targeting Barcodes", f"[bold]{nontargeting_spacers:,}[/bold]"
    )

    console.log(combined_table)


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
    parser.add_argument(
        "--json",
        action="store_true",
        default=False,
        help="Output results in JSON format",
    )

    args = parser.parse_args()

    main(args)
