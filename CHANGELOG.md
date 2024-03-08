# nf-core/metapep: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [2022-01-20]

First release of [nf-core/metapep](https://nf-co.re/metapep/dev), created based on [nf-core](https://nf-co.re) standards and [nf-core/tools](https://nf-co.re/tools) template version 1.13.1.

### `Added`

- Relational datamodel to handle large amounts of data
- Various Input options:
  - Assembly
  - Metagenomic bins
  - Taxids
- Download of Proteins via Entrez
- Protein prediction via Prodigal
- Dynamic chunking for efficient prediction and data handling
- Various options for epitope prediction:
  - SYFPEITHI
  - MHCFlurry
  - MHC-nuggets-class 1
  - MHC-nuggets-class 2
- Downstream visualizations of binding affinities and entity binding ratios
- Summarize workflow using MultiQc

- Additional subworkflow to fetch possible model and peptide lengths for the prediction tools
