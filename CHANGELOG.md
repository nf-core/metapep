# nf-core/metapep: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [2022-01-20]

Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)

### `Added`

- [#12](https://github.com/nf-core/metapep/pull/12) - Updated documentation
- [#29](https://github.com/nf-core/metapep/pull/29) - Added data model figure to `output.md`
- [#45](https://github.com/nf-core/metapep/pull/45) - Added support for multiple weights tables for one bin (i.e. co-assembly input)
- [#56](https://github.com/nf-core/metapep/pull/56),[#61](https://github.com/nf-core/metapep/pull/61) - Updated documentation
- [#67](https://github.com/nf-core/metapep/pull/67) - Added parameters to adjust binder/non-binder calling. Additional documentation on scoring by `SYPEITHI`, `MHCflurry`and `MHCnuggets`
- [#70](https://github.com/nf-core/metapep/pull/70),[#75](https://github.com/nf-core/metapep/pull/75) - Added check for supported models and model lengths and functionality to output all supported models.

### `Changed`

- [#26](https://github.com/nf-core/metapep/pull/26) - Optimized memory usage of `PREPARE_ENTITY_BINDING_RATIOS`: peptide_ids are processed chunk-wise now. Added `ds_prep_chunk_size` parameter.
- [#32](https://github.com/nf-core/metapep/pull/32) - Optimized memory usage of `PREPARE_SCORE_DISTRIBUTION`: peptide_ids are processed chunk-wise now. Peptide `count`s are used for `weight_sum` computation.
- [#33](https://github.com/nf-core/metapep/pull/33) - Output already binned scores in `PREPARE_SCORE_DISTRIBUTION` to reduce resources needed for `PLOT_SCORE_DISTRIBUTION`.
- [#36](https://github.com/nf-core/metapep/pull/36) - Replaced the tool `csvtk` by custom script for optimized TSV file concatenation in process `MERGE_PREDICTIONS_BUFFER` and `MERGE_PREDICTIONS`.
- [#44](https://github.com/nf-core/metapep/pull/44) - Optimized memory usage of `COLLECT_STATS` process
- [#46](https://github.com/nf-core/metapep/pull/46) - Input type `taxa` changed from TXT or TSV format to TSV only and increase format checks to ensure passing of optional `abundance` column
- [#55](https://github.com/nf-core/metapep/pull/55) - Updated `PREDICT_EPITOPES` to use latest epytope version (3.3.0) and adjusted the script `predict_epitopes.py` accordingly. Exchanged the container to official epytope container.
- [#64](https://github.com/nf-core/metapep/pull/64), [#65](https://github.com/nf-core/metapep/pull/65) - Updated the data model figure and using a white background to ensure readability in both dark and light themes.

### `Fixed`

- [#11](https://github.com/nf-core/metapep/pull/11) - Template update for nf-core/tools version 2.3
- [#13](https://github.com/nf-core/metapep/pull/13) - Update modules custom/dumpsoftwareversions and prodigal
- [#21](https://github.com/nf-core/metapep/pull/21) - Optimized peptide processing and Pandas joining in process `SPLIT_PRED_TASK` to reduce memory usage.
- [#24](https://github.com/nf-core/metapep/pull/24) - Optimized peptide generation in process `GENERATE_PEPTIDES` to reduce memory usage.
- [#40](https://github.com/nf-core/metapep/pull/40) - Fix the bins processing workflow, after the co-assembly feature excluded parts of the bins in `GENERATE_PROTEIN_AND_ENTITY_IDS`
- [#58](https://github.com/nf-core/metapep/pull/58) - Ensured deterministic microbiome_id and entity_id assignments.
- [#62](https://github.com/nf-core/metapep/pull/62) - nf-core module prodigal is updated

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
- [#30](https://github.com/nf-core/metapep/pull/30) - Fix redundant protein prediction for co-assembly inputs (one assembly, multiple weights)

### `Dependencies`

### `Deprecated`
