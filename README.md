# BarCoder Toolkit

The BarCoder Toolkit is designed to facilitate genome-scale experiments with a focus on efficiency and accuracy. This toolkit offers a range of utilities for genomic data processing, including sequence alignment, data compression, and read analysis.

## Important Disclaimer

* The libraries hosted in the [Example Libraries](Example_Libraries/) directory are _specifically designed as edge cases for stress testing_ the BarCoder Toolkit. They aim to evaluate the toolkit's performance under various conditions and are **not** intended for direct application.



## Installation

### Conda/Mamba
Follow these steps to set up the BarCoder environment using Conda or Mamba:

1. **Dependencies**: Bowtie is included as part of the `environment.yml` file, ensuring that dependency management is streamlined.

2. **Environment Setup**:
    - With Mamba:
      ```bash
      mamba env create -f environment.yml
      ```
    - With Conda:
      ```bash
      conda env create -f environment.yml
      ```

3. **Post-Installation**: Once the environment is set up, you are ready to use the toolkit. `classes.py` has been implemented to output JSON for valid matches for genomic targets, facilitating downstream analysis and integration.

### Pipenv (For Development)
Pipenv is supported for development environments, offering an easy and efficient way to manage dependencies:

- **Note**: The Pipenv environment does not include Bowtie. You will need to manage this dependency separately.
- **Setup**:
  To set up your development environment using Pipenv, run:
  ```bash
  pipenv install --dev

## Development Updates

Our toolkit is continually evolving to incorporate more efficient data structures and integration with other bioinformatics tools:

- **[Feb 09]** Initiated reworking of internals to incorporate interval trees for enhanced data handling and analysis.
- **[Feb 12]** Transitioned to using **pyranges**, a Python library for genomic ranges manipulation, to improve performance and compatibility.
- **[Feb 12]** Currently implementing Elasticsearch integration for advanced data querying and management capabilities.

## Contributing

The BarCoder Toolkit is under active development, and contributions are welcome. Whether it's adding new features, enhancing existing ones, or reporting bugs, your input is valuable. Please refer to the project's issue tracker to submit issues or propose pull requests.

