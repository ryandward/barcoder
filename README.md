# BarCoder Toolkit
A toolkit for genome-scale experiments

## Installation notes
* Bowtie is packaged as part of `environment.yml`, so it should be managed for you.
* Install using `mamba env create -f environment.yml` or `conda env create -f environment.yml`
* Implemented `classes.py` to output `json` for valid matches for genomic targets.

## Development notes:
* Reworking internals to use ~interval trees~. [Feb 09]
* Reworking internals to use **pyranges**. [Feb 12]
* Working to implement elasticsearch integration [Feb 12]
