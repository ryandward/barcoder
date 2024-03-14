from Bio.Seq import Seq


class CRISPRiLibrary:
    def __init__(self, pyranges_df, pam_finder):
        self.targets_df = pyranges_df
        self.pam_finder = pam_finder
        self._annotate_targets()
        self.source_unique_targets = self._get_source_unique_targets()
        self.mapped_targets = self._get_mapped_targets()
        # ⬇️ This is probably what you want
        self.unique_targets = self._get_unique_targets()
        self.unambiguous_targets = self._get_unambiguous_targets()

    def _annotate_targets(self):
        self.targets_df["PAM"] = self.targets_df.apply(
            lambda row: self.pam_finder.get_pam_seq(row), axis=1
        )
        self.targets_df["Targeting"] = self.targets_df["PAM"].apply(
            lambda x: self.pam_finder.pam_matches(x)
        )

    def _get_source_unique_targets(self):
        """
        Returns a DataFrame containing the unique source targets.

        The source targets are objects that contain information that is completely agnostic to features such as genes.
        Instead, they only contain a locus map for each chromosome.
        This method filters the targets based on their type, targeting status, and mapping status.
        The resulting DataFrame contains only the unique source targets.

        Note: The term "source" comes from NCBI GenBank annotations, where "source" refers to the chromosomal coordinate information.

        Returns:
            pandas.DataFrame: A DataFrame containing the unique source targets.
        """
        return (
            self.targets_df[
                (self.targets_df["Type"] == "source")
                & (self.targets_df["Targeting"] == True)
                & (self.targets_df["Mapped"] == True)
            ]
            .loc[lambda df: ~df.duplicated(subset=["Barcode"])]
            .reset_index(drop=True)
        )

    def _get_mapped_targets(self):
        """
        Returns a DataFrame containing the mapped targets.

        This method filters the targets DataFrame to include only the targets that meet the following criteria:
        - The target type is not "source", i.e. it is a gene or other feature
        - The target is marked as targeting
        - The target is marked as mapped

        The method then calculates the offset and overlap for each target and adds them as new columns to the DataFrame.

        Returns:
            pandas.DataFrame: DataFrame containing the mapped targets with additional columns for offset and overlap.
        """
        return (
            self.targets_df[
                (self.targets_df["Type"] != "source")
                & (self.targets_df["Targeting"] == True)
                & (self.targets_df["Mapped"] == True)
            ]
            .assign(
                Offset=lambda df: df.apply(
                    lambda row: {
                        "+": row.Start - row.Start_b,
                        "-": row.End_b - row.End,
                    }.get(row.Strand_b, None),
                    axis=1,
                ),
                Overlap=lambda df: df.apply(
                    lambda row: max(
                        min(row.End, row.End_b) - max(row.Start, row.Start_b), 0
                    ),
                    axis=1,
                ),
            )
            .reset_index(drop=True)
        )

    def _get_unique_targets(self):
        """
        Returns a list of unique guide targets based on their location in the genome.

        This method filters the mapped targets based on their barcode and sorts them by chromosome, start, and end positions.
        The resulting list contains only the targets that are unique with respect to their location in the genome.

        Returns:
            pandas.DataFrame: A DataFrame containing the unique guide targets.
        """
        mapped_targets = self._get_mapped_targets()
        return (
            mapped_targets[
                mapped_targets["Barcode"].isin(self.source_unique_targets.Barcode)
            ]
            .sort_values(["Chromosome", "Start", "End"])
            .reset_index(drop=True)
        )

    def _get_unambiguous_targets(self):
        """
        Returns a list of unambiguous guide targets based on their barcode.

        This method filters the unique targets based on their barcode and removes any duplicates.
        The resulting list contains only the targets that are unambiguous with respect to their barcode.

        Note: This method is probably excessively stringent, as it does not allow for any overlapping genes,
        which are often biologically relevant. Performing the unambiguous option will only make post-experimental
        deconvolution easier but at the potential expense of biologically meaningful information.

        Returns:
            pandas.DataFrame: A DataFrame containing the unambiguous guide targets.
        """
        return self.unique_targets[
            ~self.unique_targets.duplicated(subset=["Barcode"]).reset_index(drop=True)
        ]
