# BarCoder Toolkit: The Pinnacle of Genomic Experimentation

Unlock the full potential of genome-scale experiments with the BarCoder Toolkit, a state-of-the-art suite designed for peak performance, precision, and adaptability. This toolkit is the epitome of cutting-edge genomic data processing, offering an array of utilities from sequence alignment to data compression and read analysis, all while adhering to the SOLID principles for extreme interoperability and dependency injection.

## Installation

### Conda/Mamba
Set up the BarCoder environment with ease using Conda or Mamba:

1. **Streamlined Dependencies**: The `environment.yml` file includes all necessary dependencies, such as Bowtie, for a hassle-free setup.

2. **Environment Setup**:
    - With Mamba:
      ```bash
      mamba env create -f environment.yml
      ```
    - With Conda:
      ```bash
      conda env create -f environment.yml
      ```

3. **Ready to Go**: Post-installation, dive into the toolkit's functionalities. Our classes are crafted to output JSON for valid genomic matches, streamlining downstream analysis and integration.

### Pipenv (For Development)
For developers, Pipenv provides a seamless dependency management experience:

- **Note**: Bowtie is not included in the Pipenv environment and should be managed separately.
- **Setup**:
  ```bash
  pipenv install --dev
  ```

## Counting barcodes: `heuristicount2.py`

`heuristicount2.py` counts how many reads support each barcode/sgRNA in a pooled
library. It is the successor to `heuristicount.py`; the original is kept unchanged for
reproducibility, but new analyses should use v2.

```bash
# single-end
python heuristicount2.py library.tsv reads.fastq.gz --column spacer > counts.tsv

# paired-end
python heuristicount2.py library.tsv R1.fastq.gz R2.fastq.gz --column spacer > counts.tsv
```

The library may be a `.tsv`/`.csv` (pick the barcode column with `--column`, default
`spacer`; add a guide-ID column with `--id-column`), a FASTA, or one sequence per line.

What it does, and how it differs from v1:

- **Flank-anchored localization.** It detects the constant primer flanks around the
  barcode and finds them per read (within a bounded window), so **staggered / phased
  primers** (variable spacer length) are counted correctly instead of undercounted.
- **Quality-gated matching** (`-q`, default Phred 20). A base called at or above the
  cutoff must match a guide **exactly**; only sub-threshold (low-confidence) bases may be
  wildcarded. So a low-quality sequencing error is recovered to the right guide, but a
  *high-confidence* difference is never forgiven - which is what keeps a designed
  single-mismatch guide from absorbing reads of its perfect-match sibling. A read whose
  low-quality base can't disambiguate between two guides is skipped as *ambiguous*, never
  miscounted. `-m` caps how many low-confidence bases may be forgiven (default 2);
  `-q 0` gives strict exact matching.
- **Complete output.** One tab-separated row per library barcode, sorted, with a header,
  **including zero-count guides** (a dropped-out guide is real screen signal). There is
  no "undocumented" bin - unmatched reads are reported only as a summary fraction.
- Reads `.fastq`/`.reads`, optionally `.gz`/`.zst`; single- or paired-end; auto-detects
  read orientation; streams via multiprocessing; pipeline-safe exit codes.

### Soft (EM) mode: `--soft`

The default (hard) mode skips a read as ambiguous when a low-quality base could equally
have come from two near-neighbor guides. With `--soft`, those contested reads are instead
apportioned between the candidate guides by posterior probability (per-read mixture
maximum likelihood, solved by EM), so the abundant true source gets the read rather than
losing it. Counts are then expected values and can be fractional. This is the statistically
ideal estimator; the default hard mode is its threshold approximation, exact when neighbors
are similar in abundance. The full derivation, the bias of naive counting, and the EM are
in [`heuristicount2_model.md`](heuristicount2_model.md).

```bash
python heuristicount2.py library.tsv reads.fastq.gz --column spacer --soft > counts.tsv
```

A per-run summary (mode, detected flanks/offsets, assigned/contested/unmatched counts,
guides seen vs. dropped out) is printed to stderr; only the count table goes to stdout.

Regression tests: `python test_heuristicount2.py` (or `pytest test_heuristicount2.py`).

## Evolving with Precision

The BarCoder Toolkit is in constant evolution, embracing more sophisticated data structures and seamless integration with a plethora of bioinformatics tools. Stay tuned for our latest updates and enhancements that push the boundaries of genomic research.

## Contributing

We invite you to contribute to the BarCoder Toolkit's journey towards excellence. Your expertise can help shape the future of genomic experimentation. Engage with us through the project's issue tracker for bug reports, feature proposals, or pull requests.

