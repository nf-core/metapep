# nf-core/metapep: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Create data model](#data-model) - Create tables according to the relational data model.
* [Download proteins](#download-proteins) - Download proteins for input type taxa from Entrez.
* [Prodigal](#prodigal) - Predict proteins for input type assembly or bins.
* [Generate peptides](#generate-peptides) - Generate peptides from proteins.
* [Report stats](#report-stats) - Report some statistics on proteins and peptides.
* [Epitope prediction](#epitope-prediction) - Predict epitopes for given alleles and peptides.
* [Plot results](#plot-results) - Produce plots that summarize results.
* [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
* [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Data model

<details markdown="1">
<summary>Output files</summary>

* `db_tables/`
    * `alleles.tsv`: contains allele_id and allele_name for all unique alleles used for epitope prediction.
    * `conditions.tsv`: contains condition_id, condition_name and microbiome_id for all unique conditions.
    * `entities.tsv`: contains entity_id and entity_name for all unique entities. An entity can be a contig (for input type assembly and bins) or a taxon (for input type taxa).
    * `microbiomes_entities.nucl.tsv`: matches entities to microbiomes. Contains entity_name, microbiome_id and entity_weight for all entities of input types assembly and bins.
    * `microbiomes.tsv`: contains microbiome_id, microbiome_path, microbiome_type and weights_path for all unique microbiomes (combination of path, type and weights).
    * `proteins.tsv.gz`: contains protein_id (new unique id), protein_orig_id and protein_sequence for all unique proteins.
    * `conditions_alleles.tsv`: matches alleles to conditions. Contains condition_id and allele_id for all unique condition - allele combinations.
    * `entities_proteins.tsv`: matches proteins to entities. Contains entity_id	and protein_id for all unique entity - protein combinations.
    * `microbiomes_entities.no_weights.tsv`: matches entities to microbiomes. Contains microbiome_id and entity_id for all unique microbiome - entity combinations.
    * `microbiomes_entities.tsv`: matches entities and their weights to microbiomes. Contains microbiome_id, entity_id and entity_weight for all unique microbiome - entity combinations.
    * `proteins_peptides.tsv`: matches peptides to proteins. Contains protein_id, peptide_id and count (number of occurences of peptide in respective protein) for all unique protein - peptide combinations.

</details>

Metapep uses a relational data model that consists of tables that can describe relationships between different objects.

### Download proteins

<details markdown="1">
<summary>Output files</summary>

* `entrez_data/`
    * `entities_proteins.entrez.tsv`: matches temporary protein id given by Entrez to entities. Contains protein_tmp_id and entity_name.
    * `microbiomes_entities.entrez.tsv`: matches entities (taxa) and their weights to microbiomes. Contains microbiome_id, entity_id and entity_weight for unique microbiome - entity combinations downloaded from Entrez.
    * `proteins.entrez.tsv.gz`: contains protein_tmp_id (protein id given by Entrez) and protein_sequence for all proteins downloaded from Entrez.
    * `taxa_assemblies.tsv`: matches taxon id to assembly id.

</details>

Proteins are downloaded for input type taxa from Entrez.

### Prodigal

<details markdown="1">
<summary>Output files</summary>

* `prodigal/`
    * `*.gff`: contains proteins predicted by Prodigal in gff format.
    * `proteins.pred_*.tsv.gz`: contains proteins predicted by Prodigal in tsv format. The columns are protein_tmp_id (<contig-id_suffix>) and protein_sequence.

</details>

Proteins are predicted for input type assembly and bins.
### Generate peptides

<details markdown="1">
<summary>Output files</summary>

* `db_tables/`
    * `peptides.tsv.gz`: contains peptide_id and peptide_sequence for all unique peptides.

</details>

Peptides are generated for downloaded or predicted proteins.

### Report stats

<details markdown="1">
<summary>Output files</summary>

* `db_tables/`
    * `stats.txt`: contains statistics: unique protein counts, total peptide counts, unique peptide counts, unique peptides across all conditions.

</details>

Some statistics on protein and peptide number are calculated.

### Epitope prediction

<details markdown="1">
<summary>Output files</summary>

* `db_tables/`
    * `predictions.tsv.gz`: contains peptide_id, prediction_score (epitope prediction score) and allele_id for all unique peptide - allele combinations.
* `logs/`
    * `prediction_warnings.log`: contains warnings that occured during epitope prediction.

</details>

Epitopes are predicted for unique peptide - allele combinations.

### Plot results

<details markdown="1">
<summary>Output files</summary>

* `figures/`
    * `entity_binding_ratios/`
        * `entity_binding_ratios.allele_*.tsv`: Contains condition_name, binding_rate and entity_weight. The binding rate is calculated per entity as number of binders divided by total number of peptides. Multiple occurences of peptides within one protein are not counted.
    * `entity_binding_ratios.*.pdf`: box plot showing the binding ratios per condition and entity.
    * `entity_binding_ratios.with_points.*.pdf`: box plot showing the binding ratios per condition and entity. Each point corresponds to one entity (contig or taxon, depending on input type).
    * `prediction_scores/`
        * `prediction_scores.allele_*.tsv`: Contains prediction_score, condition_name and weight_sum. The weight_sum is calculated as the sum of all weights that belong to the entites the peptide is contained in.
    * `prediction_score_distribution.*.pdf`: weighted violin plot showing the distribution of prediction scores per condition.

</details>

Results are summarized by plots.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

* `multiqc/`
    * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
    * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
    * `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

* `pipeline_info/`
    * Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
    * Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
    * Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
