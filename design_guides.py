from numpy import full
from targets import create_topological_fasta
import argparse
import sys
from rich.console import Console
from rich.highlighter import JSONHighlighter
import os
import tempfile
from multiprocessing import cpu_count
from Bio import SeqIO
import re
import subprocess
import pandas as pd
from io import StringIO
import json


def is_dna(sequence):
    return all(base in "GATC" for base in sequence)


def find_sequences_with_barcode_and_pam(
    topological_fasta_file_name, barcode_length, pam
):
    matching_sequences = set()
    pam_regex = re.compile(pam.replace("N", "[ATGC]"))

    with open(topological_fasta_file_name, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Consider both the original sequence and its reverse complement
            for sequence in [str(record.seq), str(record.seq.reverse_complement())]:
                for i in range(len(sequence) - barcode_length - len(pam) + 1):
                    # If PAM is downstream
                    if args.pam_direction == "downstream" and pam_regex.match(
                        sequence[i + barcode_length : i + barcode_length + len(pam)]
                    ):
                        spacer = sequence[i : i + barcode_length]
                        if is_dna(spacer):
                            matching_sequences.add(spacer)
                    # If PAM is upstream
                    elif args.pam_direction == "upstream" and pam_regex.match(
                        sequence[i - len(pam) : i]
                    ):
                        spacer = sequence[i : i + barcode_length]
                        if is_dna(spacer):
                            matching_sequences.add(spacer)

    return matching_sequences


# create a sgRNA fasta file such as >sequence\nsequence\n
def create_sgRNA_fasta(matching_sequences, sgRNA_fasta_file_name):
    with open(sgRNA_fasta_file_name, "wt") as handle:
        for sequence in matching_sequences:
            handle.write(f">{sequence}\n{sequence}\n")


def main(args):

    console = Console(file=sys.stderr)
    json_console = Console(file=sys.stderr, highlighter=JSONHighlighter())

    with tempfile.TemporaryDirectory() as working_dir:
        console.log("[bold red]Initializing barcode target builder[/bold red]")
        console.log("Parameters:")
        json_console.log(json.dumps(vars(args), indent=4))

        topological_fasta_file_name = os.path.join(
            working_dir,
            os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta",
        )

        sgRNA_fasta_file_name = os.path.join(working_dir, "sgRNA.fasta")

        create_topological_fasta(args.genome_file, topological_fasta_file_name)

        matching_sequences = find_sequences_with_barcode_and_pam(
            topological_fasta_file_name, args.barcode_length, args.pam
        )

        create_sgRNA_fasta(matching_sequences, sgRNA_fasta_file_name)

        console.log(f"Found {len(matching_sequences):,} potential guides in the genome")
        console.log(
            f"Stay tuned... running 'targets.py' to find guides for {args.genome_file} with {args.barcode_length}bp barcodes and {args.pam} PAM sequence"
        )

        result = subprocess.run(
            [
                "python",
                "targets.py",
                sgRNA_fasta_file_name,
                args.genome_file,
                args.pam,
                "--pam_direction",
                args.pam_direction,
                str(args.mismatches),
            ],
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )

        targets = pd.read_csv(StringIO(result.stdout), sep="\t")

        console.log(f"Found {len(targets):,} guides")
        console.log(f"Removing guides based on settings")

        targets["target"] = targets["target"].str.upper()

        if args.orientation == "forward":
            # omit everything where sp_dir and tar_dir are not the same
            targets = targets.loc[targets["sp_dir"] == targets["tar_dir"]]

        elif args.orientation == "reverse":
            # omit everything where sp_dir and tar_dir are the same
            targets = targets.loc[targets["sp_dir"] != targets["tar_dir"]]

        if args.omit_offtargets:

            console.log("[bold red]Removing guides with off-targets[/bold red]")
            len_before = len(targets)
            # Extract the number of sites from the 'note' column
            targets.loc[:, "sites"] = (
                targets["note"].str.extract(r"(\d+) site", expand=False).astype(int)
            )
            # Create a mask that is True for rows where 'sites' is 1
            mask = targets["sites"] == 1

            # Apply the mask to the DataFrame
            targets = targets[mask]
            console.log(f"Removed {(len_before - len(targets)):,} guides")
        # remove everything where mismatches > 0
        if args.mismatches > 0:
            console.log(
                "[bold red]Removing guides with mismatches.\nThere shouldn't be any if offtargets are omitted![/bold red]"
            )
            len_before = len(targets)

            targets = targets.loc[targets["mismatches"] == 0]

            console.log(f"Removed {(len_before - len(targets)):,} guides")

        if args.omit_ambiguous:
            console.log(
                "[bold red]Removing ambiguous guides, this will be a lot![/bold red]"
            )

            # Check if 'note' column exists
            if "note" in targets.columns:
                # Extract the number of sites, genes, and intergenic regions from the 'note' column
                targets["sites"] = (
                    targets["note"]
                    .str.extract(r"(\d+) site", expand=False)
                    .fillna(0)
                    .astype(int)
                )

                targets["genes"] = (
                    targets["note"]
                    .str.extract(r"(\d+) gene", expand=False)
                    .fillna(0)
                    .astype(int)
                )

                targets["intergenic"] = (
                    targets["note"]
                    .str.extract(r"(\d+) intergenic", expand=False)
                    .fillna(0)
                    .astype(int)
                )
            else:
                targets["sites"] = 0
                targets["genes"] = 0
                targets["intergenic"] = 0

            # Create a mask that is True for rows where 'sites', 'genes', are 1 and 'intergenic' is 0
            mask = (
                (targets["sites"] == 1)
                & (targets["genes"] == 1)
                & (targets["intergenic"] == 0)
            )

            len_before = len(targets)

            # Apply the mask to the DataFrame
            targets = targets[mask]

            console.log(f"Removed {(len_before - len(targets)):,} guides")

        if args.omit_intergenic:
            console.log("[bold red]Removing intergenic regions[/bold red]")

            # Create a mask that is True for rows where 'note' does not contain "intergenic"
            mask = ~targets["note"].str.contains("intergenic")

            len_before = len(targets)

            # Apply the mask to the DataFrame
            targets = targets[mask]

            console.log(f"Removed {(len_before - len(targets)):,} guides")

        if args.full_overlap:
            console.log(
                "[bold red]Removing guides that don't fully overlap with the gene[/bold red]"
            )
            # Sort the DataFrame by 'target', 'spacer' and 'locus_tag'

            len_before = len(targets)

            # Create a new DataFrame where 'overlap' is equal to args.barcode_length
            overlap_df = targets.loc[targets["overlap"] == args.barcode_length]

            # Extract the 'spacer' column from overlap_df and convert it to a set
            overlap_spacers = set(overlap_df["spacer"])

            # Create a mask that is True for rows where 'spacer' is in overlap_spacers
            mask = targets["spacer"].isin(overlap_spacers)

            # Apply the mask to the DataFrame
            targets = targets[mask]

            console.log(f"Removed {(len_before - len(targets)):,} guides")

        # 0 or -1 will result in all barcodes being selected
        if args.tile_size > 0:
            # Sort the DataFrame by 'locus_tag' and 'offset'
            targets = targets.sort_values(["locus_tag", "offset"])

            # Group the DataFrame by 'locus_tag'
            grouped = targets.groupby("locus_tag")

            # Initialize a set to store the selected spacers
            selected_spacers = set()

            # Iterate over each group
            for name, group in grouped:

                # Here, last_offset is the offset of the last spacer that was added to selected_spacers

                if args.full_overlap:
                    filtered_df = group["offset"].loc[
                        group["overlap"] == args.barcode_length
                    ]
                    if not filtered_df.empty:
                        last_offset = filtered_df.iloc[0]
                    else:
                        # Handle the case when the DataFrame is empty
                        # For example, set last_offset to a default value
                        last_offset = None
                else:
                    last_offset = group["offset"].iloc[0]

                if last_offset is not None:
                    selected_spacers.add(
                        group["spacer"].loc[group["offset"] == last_offset].iloc[0]
                    )

                # iterate through the rest of the group
                for index, row in group.iterrows():

                    # If the current offset is at least tile_size away from the last offset, add the spacer to the set
                    if (
                        last_offset is not None
                        and row["offset"] >= last_offset + args.tile_size
                    ):
                        selected_spacers.add(row["spacer"])
                        last_offset = row["offset"]

            # Create a mask that is True for rows where 'spacer' is in selected_spacers
            mask = targets["spacer"].isin(selected_spacers)

            # Apply the mask to the DataFrame
            targets = targets[mask]

        if args.keep_top > 0:
            console.log(
                f"[bold red]Keeping only the top {args.keep_top} guides for each gene[/bold red]"
            )

            len_before = len(targets)

            # Sort the DataFrame by 'locus_tag', 'offset', and 'overlap'
            targets = targets.sort_values(["locus_tag", "offset", "overlap"])

            # Group the DataFrame by 'locus_tag' and keep the top n rows
            top_targets = targets.groupby("locus_tag").head(args.keep_top)

            # Get the spacers that are in the top n rows for each gene
            top_spacers = top_targets["spacer"].unique()

            # Filter the original DataFrame to keep only the rows that are in the top n rows for each gene or that target multiple genes
            targets = targets[targets["spacer"].isin(top_spacers)]

            console.log(f"Removed {(len_before - len(targets)):,} guides")

    # Convert all numeric columns to integers
    targets = targets.apply(
        lambda col: (
            pd.to_numeric(col, errors="coerce").fillna(0).astype(int)
            if col.dtypes != object
            else col
        )
    )

    targets = targets.sort_values(
        ["chr", "tar_start", "tar_end", "locus_tag", "offset", "overlap"]
    )

    # print csv to output
    targets.to_csv(sys.stdout, sep="\t", index=False, na_rep="None")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map barcodes to a circular genome",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
    parser.add_argument("pam", help="PAM sequence", type=str)
    parser.add_argument("barcode_length", help="Length of the barcode", type=int)
    parser.add_argument(
        "--orientation",
        choices=["forward", "reverse", "both"],
        default="forward",
        help="Orientation of the barcode with respect to the gene. Default is both.",
    )
    parser.add_argument(
        "--mismatches",
        type=int,
        default=1,
        metavar="(0-2)",
        help="Number of mismatches to constitute an offtarget.",
    )
    parser.add_argument(
        "--pam_direction",
        choices=["upstream", "downstream"],
        default="downstream",
        help="Direction of the PAM sequence",
    )
    parser.add_argument(
        "--omit_intergenic",
        action="store_true",
        default=True,
        help="Omit intergenic regions",
    )
    parser.add_argument(
        "--omit_offtargets",
        action="store_true",
        default=False,
        help="Omit all guides that have off-targeting",
    )
    parser.add_argument(
        "--omit_ambiguous",
        action="store_true",
        default=False,
        help="Target only sites that have a single gene annotated.",
    )
    parser.add_argument(
        "--keep-top",
        type=int,
        default=10,
        metavar="(1-n)",
        help="Keep the top n guides for each gene",
    )
    parser.add_argument(
        "--tile_size",
        type=int,
        default=None,
        metavar="(1-n)",
        help="Tile size for the genome (defaults to barcode length)",
    )
    parser.add_argument(
        "--full-overlap",
        action="store_true",
        default=False,
        help="Require full overlap of the guide with the gene. This can get messy if genes overlap!",
    )

    args = parser.parse_args()

    if not args.tile_size:
        args.tile_size = args.barcode_length

    if args.omit_ambiguous:
        args.omit_offtargets = True

    main(args)
