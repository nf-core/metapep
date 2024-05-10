# nf-core/metapep: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/metapep/usage](https://nf-co.re/metapep/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 5 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

A final samplesheet file may look something like the one below.

```csv title="samplesheet.csv"
condition,type,microbiome_path,alleles,weights_path
cond_1,taxa,testdata/taxids.txt,A*01:01,
cond_2,taxa,testdata/taxids.tiny.txt,A*01:01 B*07:02,
cond_3,taxa,testdata/taxids.tiny.txt,A*01:01,
cond_4,assembly,testdata/test_minigut.contigs.fa.gz,A*01:01,testdata/test_minigut.contig_weights.tsv
```

| Column            | Description                                                                                                                                                                                                |
| ----------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `condition`       | The condition name for this entry. Conditions have to be unique and describe a combination of a microbiome, alleles and weights.                                                                           |
| `type`            | Input type, can be one of "assembly", "bins" or "taxa".                                                                                                                                                    |
| `microbiome_path` | Full path to microbiome file, the format of which can vary with type: fasta, folder, compressed folder or tsv file for taxon ids (taxon_id ["\\t" assembly_ids "\\t" abundance]).                          |
| `alleles`         | List of alleles to predict epitopes for.                                                                                                                                                                   |
| `weights_path`    | Full path to a tab-separated file contataining weights. Currently allowed are contig weights for input types assembly and bins. Please use "contig_name" or "bin_basename" and "weight" as column headers. |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Input type taxa

The input type taxa allows the user to specify a taxid on strain level, for which all proteins are downloaded. The microbiome path corresponds to a tsv file containing one column for taxon_id and optionally specific assembly_ids and/or the abundance of a specific strain (see below). If only the taxid(s) (and optionally abundance) is provided, the pipeline will automatically download the largest assembly of the given taxid(s). If a specific assembly_id is provided it will download proteins of the given assembly_id. The input allows for a mixed assignment of specific assembly_ids and unspecific taxon_ids, but this is only recommended for specific use cases.

As the pipeline is only generating reproducible results for this type of input if a specific assembly_id is chosen, it is highly recommended to use this option.
If taxids without assembly ids were chosen as input, the pipeline results can be reproduced in following runs using the `./outdir/entrez_data/taxa_assemblies.tsv` reference file. If additional abundances were given for each taxon_id, the `input.csv` and `taxa_assemblies.tsv` files can be merged by matching the taxon_id column.

| Column        | Description                                                           |
| ------------- | --------------------------------------------------------------------- |
| `taxon_id`    | Chosen Taxids for the microbiome condition (Must be on strain level). |
| `assembly_id` | Specific assembly id for a strain level taxid.                        |
| `abundance`   | Abundance of strain level taxid and/or assembly id.                   |

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/metapep --input ./samplesheet.csv --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/metapep -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
pred_method: 'syfpeithi'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Memory considerations

The pipeline needs to handle large amounts of data, depending on the size and number of microbiomes the user has defined in the input. To handle these data the pipeline mainly uses python scripts and the python module pandas. As the data needs to be compared memory consumption is currently peaking at around 150 GB for the full-size test, but can easily be higher depending on the input.

If the memory is still an issue one can try to reduce the chunk sizes for high memory consuming processes. The parameters are: `--chunk_size <INTEGER>` and the scaling factor `--chunk_size_scaling <INTEGER>` which are used for the preprocessing of the peptides prior to the epitope prediction in `SPLIT_PRED_TASKS` and the downstream processes `MERGE_PREDICTIONS`, `PREPARE_ENTITY_BINDING_RATIOS` and `PREPARE_SCORE_DISTRIBUTION`. For for the epitope prediction process `PREDICT_EPITOPES` the chunk size equals the unscaled parameter `--chunk_size <INTEGER>`.

### Supported allele models

The pipeline predicts epitopes for specific peptide lengths and for specific alleles of MHC class I or class II. As the prediction is performed by external tools, the user is restricted to the corresponding combinations the external tools are offering. Therefore, the metapep pipeline comes with a functionality to output all supported alleles and supported lengths of the supported external tools, which is invoked by:

```bash
nextflow run nf-core/metapep -profile <YOURPROFILE> --outdir <OUTDIR> --show_supported_models
```

More on the output can be found at https://nf-co.re/metapep/dev/output#supported-allele-models

Moreover, the pipeline checks if a supported prediction model (combination of allele and peptide length) is available if a PSSMs method like SYFPEITHI is chosen and reduces the peptide lengths to a common denominator for further analysis if models are not available.

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/metapep
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/metapep releases page](https://github.com/nf-core/metapep/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

If the pipeline input contains microbiomes of type `taxa` it may not generate the same results, as for each taxid the largest assembly is chosen, which might change. Therefore, using a specific assembly_id as explained in the section `Input type taxa` is highly recommended. If only a taxid is used for input the pipeline generates a linking file for taxid and assembly id (`./outdir/entrez_data/taxa_assembly.tsv`) which can be used for following pipeline runs.

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data but needs your own NCBI credentials
- `test_assembly_only`
  - A small test profile, whith complete configuration for automated testing
  - Includes links to test data and needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
