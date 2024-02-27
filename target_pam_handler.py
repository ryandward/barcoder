# Extract the downstream PAM of the target sequence from the topological map
import re
from Bio.Seq import Seq

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

    if pam == "":
        return None

    if None in (
        pam,
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


# Extract the upstream PAM of the target sequence from the topological map
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

    if pam == "":
        return None

    if None in (
        pam,
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

# Filter out spacers that don't match the PAM
def filter_offtargets_by_pam(df):
    targeting_spacers = df[df["target"].notna()]["spacer"].unique()
    return df[~((df["target"].isna()) & (df["spacer"].isin(targeting_spacers)))]

# Check if the extracted PAM matches the PAM pattern
def pam_matches(pam_pattern, extracted_pam):
    # Convert N to . for regex matching
    if extracted_pam is None:
        return False

    if pam_pattern == "N" * len(pam_pattern) or not pam_pattern:
        return True

    regex_pattern = pam_pattern.replace("N", ".")
    return bool(re.match(regex_pattern, extracted_pam))