# nf-core/metapep: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [2022-01-20]

Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)

### `Added`

- [#12](https://github.com/nf-core/metapep/pull/12) - Updated documentation

### `Fixed`

- [#11](https://github.com/nf-core/metapep/pull/11) - Template update for nf-core/tools version 2.3
- [#13](https://github.com/nf-core/metapep/pull/13) - Update modules custom/dumpsoftwareversions and prodigal

### `Dependencies`

### `Deprecated`

## v0.9dev - [date]

Initial release of nf-core/metapep, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#1](https://github.com/skrakau/metapep/pull/1) - Download of proteins from NCBI's Entrez databases
- [#2](https://github.com/skrakau/metapep/pull/2) - Add peptide generationfrom proteins
- [#9](https://github.com/skrakau/metapep/pull/9),[#10](https://github.com/skrakau/metapep/pull/10) - Add protein prediction using `Prodigal v2.6.3`
- [#11](https://github.com/skrakau/metapep/pull/11) - Compute protein weights based on taxonomic abundances or contig depths [#4](https://github.com/skrakau/metapep/issues/4)
- [#20](https://github.com/skrakau/metapep/pull/20), [#26](https://github.com/skrakau/metapep/pull/26) - Add creation of main db files for new data model (published in `results/db_tables`)
- [#19](https://github.com/skrakau/metapep/pull/19) - Integrate epitope prediction into the pipeline
- [#20](https://github.com/skrakau/metapep/pull/20) - New input format allowing handling of multiple conditions
- [#36](https://github.com/skrakau/metapep/pull/36) - Add generation of first plot: score distributions for different conditions (and alleles)
- [#44](https://github.com/skrakau/metapep/pull/44) - Adjusted process-specific resource requirements
- [#47](https://github.com/skrakau/metapep/pull/47) - Add peptides subsampling parameter `--sample_n`
- [#60](https://github.com/skrakau/metapep/pull/60) - Add entities to data model
- [#60](https://github.com/skrakau/metapep/pull/60) - Add option for `bins` input type
- [#60](https://github.com/skrakau/metapep/pull/60) - Add plotting of entity-wise binder ratios
- [#64](https://github.com/skrakau/metapep/pull/64) - Add process to collect stats (`results/db_tables/stats.txt`)

### `Fixed`

- [#41](https://github.com/skrakau/metapep/pull/41) - Allow `assembly` input without weights
- [#53](https://github.com/skrakau/metapep/pull/53) - Add buffering of predictions and chunk-wise merging to avoid sbatch error due to too many input files [#52](https://github.com/skrakau/metapep/issues/52)
- [#3](https://github.com/nf-core/metapep/pull/3) - Fix generation of `entities_proteins.entrez.tsv` for Entrez download

### `Dependencies`

### `Deprecated`
