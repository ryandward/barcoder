# BarCoder Toolkit
A toolkit for genome-scale experiments
## Development note:
* Reworking internals to use interval trees.
* Bowtie is now part of the Conda/Mamba package.
* Install using `mamba env create -f environment.yml` or `conda env create -f environment.yml`
* Implemented `classes.py` to output `json` for valid matches for genomic targets.

#  First pass at creating a simple approach for guide design in a CRISPR-flavored experiment.
```python
barcodes = BarCodeLibrary("E_coli_library.tsv", column="spacer")
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
      "matches":"thousands of matches"
   },
    {
        "read": {
            "chromosome": "CP023719.1",
            "start": 21912,
            "end": 21944,
            "sequence": "ATTGTTGAAGCTTTGGCTATTCATACTGATCA",
            "is_reverse": false,
            "mismatches": 0
        },
        "feature": {
            "chromosome": "CP023719.1",
            "start": 21894,
            "end": 22275,
            "type": "gene",
            "locus_tag": "ZMO1_ZMOp39x021",
            "gene_name": null,
            "strand": 1
        },
        "match": {
            "orientation": "sense",
            "offset": 18,
            "overlap": 32
        }
    },
    {
        "read": {
            "chromosome": "CP023719.1",
            "start": 22024,
            "end": 22056,
            "sequence": "ATTATAACGATCTTATACAAGAATCTGCTGCA",
            "is_reverse": false,
            "mismatches": 0
        },
        "feature": {
            "chromosome": "CP023719.1",
            "start": 21894,
            "end": 22275,
            "type": "gene",
            "locus_tag": "ZMO1_ZMOp39x021",
            "gene_name": null,
            "strand": 1
        },
        "match": {
            "orientation": "sense",
            "offset": 130,
            "overlap": 32
        }
    },
   {
      "matches":"thousands more matches"
   }
]
  ```
