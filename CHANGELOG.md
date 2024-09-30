# nf-core/metapep: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - [2022-01-20]

First release of [nf-core/metapep](https://nf-co.re/metapep), created based on [nf-core](https://nf-co.re) standards and [nf-core/tools](https://nf-co.re/tools) template version 1.14.1.

### `Added`

- Download proteins for input type taxa from [Entrez](https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html).
- Predict proteins for input type assembly or bins using [Prodigal](https://github.com/hyattpd/Prodigal).
- Generate peptides from proteins.
- Split peptide files into chunks for parallel prediction and report stats.
- Predict epitopes for given alleles and peptides using [SYFPEITHI](http://www.syfpeithi.de), [MHCflurry](https://github.com/openvax/mhcflurry) or [MHCnuggets](https://github.com/KarchinLab/mhcnuggets).
- Downstream visualizations between conditions (different microbiomes assemblies, bins, taxids or same input class with different weights) given within samplesheet
  - binding affinities
  - entity binding ratios
- Summarize workflow using MultiQC

- Relational datamodel to handle large amounts of data
  Tables defined within model:

  - `alleles.tsv`
  - `condition_alleles.tsv`
  - `conditions.tsv`
  - `entities_proteins.tsv`
  - `entities.tsv`
  - `microbiomes_entities.no_weights.tsv`
  - `microbiome_entities.nucl.tsv`
  - `microbiomes_entities.tsv`
  - `microbiomes.tsv`
  - `preptides.tsv.gz`
  - `predictions.tsv.gz`
  - `proteins_peptides.tsv`
  - `proteins.tsv.gz`
  - `stats.txt`

- Additional subworkflow to fetch possible model and peptide lengths for the prediction tools
