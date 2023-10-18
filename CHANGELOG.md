# nf-core/metapep: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0dev - [2022-01-20]

Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)

### `Added`

- [#12](https://github.com/nf-core/metapep/pull/12), [#56](https://github.com/nf-core/metapep/pull/56), [#61](https://github.com/nf-core/metapep/pull/61), [#88](https://github.com/nf-core/metapep/pull/88), [#99](https://github.com/nf-core/metapep/pull/99) - Updated documentation, Citations and Changelog
- [#29](https://github.com/nf-core/metapep/pull/29), [#64](https://github.com/nf-core/metapep/pull/64), [#65](https://github.com/nf-core/metapep/pull/65) - Added data model figure to `output.md`
- [#40](https://github.com/nf-core/metapep/pull/40), [#45](https://github.com/nf-core/metapep/pull/45) - Added support for multiple weights tables for one bin (i.e. co-assembly input)
- [#67](https://github.com/nf-core/metapep/pull/67) - Added parameters to adjust binder/non-binder calling. Additional documentation on scoring by `SYPEITHI`, `MHCflurry`and `MHCnuggets`.
- [#78](https://github.com/nf-core/metapep/pull/78) - Added parameter `memory_usage_log_deep` for pandas memory usage logging.
- [#78](https://github.com/nf-core/metapep/pull/78), [#96](https://github.com/nf-core/metapep/pull/96) - Added chunk size parameters (`--prediction_chunk_size`, `--pred_chunk_size_scaling`, `--downstream_chunk_size`, `--pred_buffer_files`) to optimise the number and size of jobs for large datasets
- [#70](https://github.com/nf-core/metapep/pull/70) - Added check for supported models and functionality to output all supported models (`--show_supported_models`).
- [#75](https://github.com/nf-core/metapep/pull/75), [#89](https://github.com/nf-core/metapep/pull/85), [#92](https://github.com/nf-core/metapep/pull/92), [#95](https://github.com/nf-core/metapep/pull/95) - Added module `UNIFY_MODEL_LENGTHS`: Peptide lengths are checked if supported generally and unified if models are not available
- [#80](https://github.com/nf-core/metapep/pull/80) - Add mean comparison to entity binding ratios plots.
- [#84](https://github.com/nf-core/metapep/pull/84) - Added module `MultiQC` again with updated method description and updated references

### `Changed`

#### `Resource usage optimization`
- [#26](https://github.com/nf-core/metapep/pull/26) - Optimized memory usage of `PREPARE_ENTITY_BINDING_RATIOS`: peptide_ids are processed chunk-wise now.
- [#32](https://github.com/nf-core/metapep/pull/32) - Optimized memory usage of `PREPARE_SCORE_DISTRIBUTION`: peptide_ids are processed chunk-wise now. Peptide `count`s are used for `weight_sum` computation.
- [#21](https://github.com/nf-core/metapep/pull/21) - Optimized peptide processing and Pandas joining in process `SPLIT_PRED_TASK` to reduce memory usage.
- [#24](https://github.com/nf-core/metapep/pull/24) - Optimized peptide generation in process `GENERATE_PEPTIDES` to reduce memory usage.
- [#44](https://github.com/nf-core/metapep/pull/44) - Optimized memory usage of `COLLECT_STATS` process
- [#33](https://github.com/nf-core/metapep/pull/33) - Output already binned scores in `PREPARE_SCORE_DISTRIBUTION` to reduce resources needed for `PLOT_SCORE_DISTRIBUTION`.
- [#36](https://github.com/nf-core/metapep/pull/36) - Replaced the tool `csvtk` by custom script for optimized TSV file concatenation in process `MERGE_PREDICTIONS_BUFFER` and `MERGE_PREDICTIONS`.

#### `Miscellaneous`

- [#46](https://github.com/nf-core/metapep/pull/46) - Input type `taxa` changed to TSV only and increase format checks to ensure passing of optional `abundance` column
- [#76](https://github.com/nf-core/metapep/pull/76) - Restructure and rename `CHECK_INPUT` subworkflow to `PREPROCESS_INPUT` and move some channel logic into the subworkflow. Also `CHECK_SAMPLESHEET` module is renamed to `INPUT_TO_DATAMODEL` to have a more descriptive name and circumvent clashing with the nf-core template.
- [#82](https://github.com/nf-core/metapep/pull/82) - Add label `error_retry` to process `DOWNLOAD_PROTEINS`
- [#91](https://github.com/nf-core/metapep/pull/91) - Added new parameter (`--allow_inconsistent_pep_lengths`) to allow all input peptide lengths to be predicted, and also necessary checks to circumvent errors during predictions
- [#94](https://github.com/nf-core/metapep/pull/94) - Updated workflow overview figure
- [#93](https://github.com/nf-core/metapep/pull/93), [#98](https://github.com/nf-core/metapep/pull/98) - Added new tests and replaced tests to have a broader test landscape

### `Fixed`

- [#58](https://github.com/nf-core/metapep/pull/58) - Ensured deterministic microbiome_id and entity_id assignments.

### `Dependencies`

- [#11](https://github.com/nf-core/metapep/pull/11), [#15](https://github.com/nf-core/metapep/pull/15), [#16](https://github.com/nf-core/metapep/pull/16), [#19](https://github.com/nf-core/metapep/pull/19), [#28](https://github.com/nf-core/metapep/pull/28), [#31](https://github.com/nf-core/metapep/pull/31), [#69](https://github.com/nf-core/metapep/pull/69), [#87](https://github.com/nf-core/metapep/pull/87), [#97](https://github.com/nf-core/metapep/pull/97)   - Template update for nf-core/tools up to version 2.10
- [#13](https://github.com/nf-core/metapep/pull/13), [#62](https://github.com/nf-core/metapep/pull/62) - Update modules custom/dumpsoftwareversions and prodigal
- [#55](https://github.com/nf-core/metapep/pull/55), [#73](https://github.com/nf-core/metapep/pull/73) - Updated `PREDICT_EPITOPES` to use latest epytope version (3.3.1) and adjusted the script `predict_epitopes.py` accordingly. Exchanged the container to official epytope container.
- [#79](https://github.com/nf-core/metapep/pull/79) - Replaced multi-package BioContainers by single-package pandas BioContainer and updated containers with old pandas versions

### `Deprecated`

- [#77](https://github.com/nf-core/metapep/pull/77) - Remove unused github workflow files to push dockerimage to dockerhub
- [#81](https://github.com/nf-core/metapep/pull/81) - Removed subsampling option
- [#90](https://github.com/nf-core/metapep/pull/90) - Removed script headers and unused imports

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
