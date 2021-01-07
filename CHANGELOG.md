# nf-core/metapep: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [date]

Initial release of nf-core/metapep, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#1](https://github.com/skrakau/metapep/pull/1) - Download of proteins from NCBI's Entrez databases
- [#2](https://github.com/skrakau/metapep/pull/2) - Add peptide generationfrom proteins
- [#9](https://github.com/skrakau/metapep/pull/9),[#10](https://github.com/skrakau/metapep/pull/10) - Add protein prediction using `Prodigal v2.6.3`
- [#11](https://github.com/skrakau/metapep/pull/11) - Compute protein weights based on taxonomic abundances or contig depths [#4](https://github.com/skrakau/metapep/issues/4)
- [#20](https://github.com/skrakau/metapep/pull/20), [#26](https://github.com/skrakau/metapep/pull/26) - Add creation of main db files for new data model (published in `results/db_tables`)
- [#20](https://github.com/skrakau/metapep/pull/20) - New input format allowing handling of multiple conditions
- [#36](https://github.com/skrakau/metapep/pull/36) - Add generation of first plot: score distributions for different conditions (and alleles)

### `Fixed`

- [#41](https://github.com/skrakau/metapep/pull/41) - Allow `assembly` input without weights

### `Dependencies`

### `Deprecated`
