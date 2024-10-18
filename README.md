<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-metapep_logo_dark.png">
    <img alt="nf-core/metapep" src="docs/images/nf-core-metapep_logo_light.png">
  </picture>
</h1>

**From metagenomes to peptides**.

[![GitHub Actions CI Status](https://github.com/nf-core/metapep/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/metapep/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/metapep/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/metapep/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/metapep/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/metapep)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23metapep-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/metapep)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/metapep** is a bioinformatics best-practice analysis pipeline for epitope prediction specifically designed for metagenomes. It integrates multiple types of input (proteins, taxa, assemblies and bins), generates peptides and predicts their MHC-/HLA-affinity.

<p align="center">
    <img src="docs/images/metapep_workflow_figure.png" alt="nf-core/metapep workflow overview" width="90%">
</p>

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/metapep/results).

## Pipeline summary

1. Download proteins for input type taxa from [Entrez](https://www.ncbi.nlm.nih.gov/Web/Search/entrezfs.html).
2. Predict proteins for input type assembly or bins using [Prodigal](https://github.com/hyattpd/Prodigal).
3. Generate peptides from proteins.
4. Split peptide files into chunks for parallel prediction and report stats.
5. Predict epitopes for given alleles and peptides using [SYFPEITHI](http://www.syfpeithi.de), [MHCflurry](https://github.com/openvax/mhcflurry) or [MHCnuggets](https://github.com/KarchinLab/mhcnuggets).
6. Produce plots.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

```bash
nextflow run nf-core/metapep \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/metapep/usage) and the [parameter documentation](https://nf-co.re/metapep/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/metapep/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/metapep/output).

## Credits

nf-core/metapep was written by [Léon Kuchenbecker](https://github.com/lkuchenb) at the [Kohlbacher lab](https://kohlbacherlab.org), [Sabrina Krakau](https://github.com/skrakau) and [Till Englert](https://github.com/tillenglert) at the [QBiC](http://qbic.life/).

We thank the following people for their extensive assistance in the development of this pipeline:

- [Antonia Schuster](https://github.com/AntoniaSchuster)
- [Daniel Straub](https://github.com/d4straub)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#metapep` channel](https://nfcore.slack.com/channels/metapep) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/metapep for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
