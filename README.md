# BarCoder Toolkit
A toolkit for genome-scale experiments
## Development note:
* Reworking internals to use interval trees.
* Bowtie is now part of the Conda/Mamba package.
* Install using `mamba env create -f environment.yml` or `conda env create -f environment.yml`
* Implemented `classes.py` to output `json` for valid matches for genomic targets.

#  First pass at creating a simple approach for guide design in a CRISPR-flavored experiment.
```python
barcodes = BarCodeLibrary(E_coli_library.tsv", column="spacer")
genbank = GenBankParser("GCA_000005845.2.gb")

with BowtieRunner() as bowtie:
    bowtie.write_fasta(genbank.records)
    bowtie.write_fastq(barcodes.barcodes)
    bowtie.create_index()
    bowtie.align(num_mismatches=1, num_threads=12)
    sam = PySamParser(bowtie.sam_path)
    overlaps = find_overlaps(sam, genbank)
    print(json.dumps(overlaps, indent=4))


```
  ```json
 [
    {
    "read": {
      "chromosome": "U00096.3",
      "start": 1952787,
      "end": 1952807,
      "sequence": "GGATATGATCCAGCGTTCCG",
      "is_reverse": true,
      "mis_matches": 0
    },
    "feature": {
      "chromosome": "U00096.3",
      "start": 1952701,
      "end": 1953445,
      "locus_tag": "b1870",
      "gene_name": "cmoA",
      "strand": 1
    }
  },
  {
    "read": {
      "chromosome": "U00096.3",
      "start": 1952759,
      "end": 1952779,
      "sequence": "TTGATGAACGGGTAGCTGAA",
      "is_reverse": true,
      "mis_matches": 0
    },
    "feature": {
      "chromosome": "U00096.3",
      "start": 1952701,
      "end": 1953445,
      "locus_tag": "b1870",
      "gene_name": "cmoA",
      "strand": 1
    }
  }
]
  ```
