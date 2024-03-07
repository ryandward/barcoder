class CRISPRiLibrary:
    def __init__(self, targets_df, pam_finder):
        self.targets_df = targets_df
        self.pam_finder = pam_finder
        self._annotate_targets()
        self.source_unique_targets = self.get_source_unique_targets()
        self.feature_targets = self.get_feature_targets()
        self.feature_unique_targets = (
            self.get_feature_unique_targets()
        )  # This is probably what you want
        self.feature_unambiguous_targets = self.get_feature_unambiguous_targets()

    def _annotate_targets(self):
        self.targets_df["PAM"] = self.targets_df.apply(
            lambda row: self.pam_finder.get_pam_seq(row), axis=1
        )
        self.targets_df["Targeting"] = self.targets_df["PAM"].apply(
            lambda x: self.pam_finder.pam_matches(x)
        )

    def get_source_unique_targets(self):
        return (
            self.targets_df[
                (self.targets_df["Type"] == "source")
                & (self.targets_df["Targeting"] == True)
                & (self.targets_df["Mapped"] == True)
            ]
            .loc[lambda df: ~df.duplicated(subset=["Barcode"])]
            .reset_index(drop=True)
        )

    def get_feature_targets(self):
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

    def get_feature_unique_targets(self):
        feature_targets = self.get_feature_targets()
        return (
            feature_targets[
                feature_targets["Barcode"].isin(self.source_unique_targets.Barcode)
            ]
            .sort_values(["Chromosome", "Start", "End"])
            .reset_index(drop=True)
        )

    def get_feature_unambiguous_targets(self):
        return self.feature_unique_targets[
            ~self.feature_unique_targets.duplicated(subset=["Barcode"]).reset_index(
                drop=True
            )
        ]
