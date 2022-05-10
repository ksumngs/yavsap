# ![yavsap](docs/images/nf-core/yavsap_logo_light.png#gh-light-mode-only) ![yavsap](docs/images/nf-core/yavsap_logo_dark.png#gh-dark-mode-only)

<!--markdownlint-disable line-length -->

[![GitHub Actions CI Status](https://github.com/ksumngs/yavsap/actions/workflows/ci.yml/badge.svg)](https://github.com/ksumngs/yavsap/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/ksumngs/yavsap/actions/workflows/linting.yml/badge.svg)](https://github.com/ksumngs/yavsap/actions/workflows/linting.yml)
<!--[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/yavsap/results)-->
<!--[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)-->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://apptainer.org/docs/)

<!--
[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23yavsap-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/yavsap)
[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)
[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)
-->

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://ksumngs.github.io/yavsap)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/ksumngs/yavsap?label=version)](https://github.com/ksumngs/yavsap/blob/master/CHANGELOG.md)
[![GitHub license](https://img.shields.io/github/license/ksumngs/yavsap)](https://github.com/ksumngs/yavsap/blob/master/LICENSE)

> This project follows the [semver](https://semver.org) _pro forma_ and uses the [git-flow](https://nvie.com/posts/a-successful-git-branching-model) branching model.

## Introduction

**yavsap** is a bioinformatics best-practice analysis pipeline for identifying and analyzing viral haplotypes in metagenomic NGS reads.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
<!--
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/yavsap/results).
-->

## Pipeline summary

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)/[`NanoStat`](https://github.com/wdecoster/nanostat))
2. Read trimming ([`Trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic)/[`NanoFilt`](https://github.com/wdecoster/nanofilt/))
3. Host read filtering ([`Kraken2`](https://github.com/DerrickWood/kraken2/wiki)+[`krakentools`](https://github.com/jenniferlu717/KrakenTools))
4. Consensus sequence generation
    1. Reference genome download ([`entrez-direct`](https://www.ncbi.nlm.nih.gov/books/NBK179288/))
    2. Read alignment ([`minimap2`](https://lh3.github.io/minimap2/))
    3. Variant calling ([`CliqueSNV`](https://github.com/vtsyvina/CliqueSNV)/[`HapLink.jl`](https://ksumngs.github.io/HapLink.jl))
    4. Consensus sequence generation ([`CliqueSNV`](https://github.com/vtsyvina/CliqueSNV)/[`HapLink.jl`](https://ksumngs.github.io/HapLink.jl))
5. Strain identification ([`BLAST+`](https://www.ncbi.nlm.nih.gov/books/NBK569839/))
6. Variant calling
    1. Read alignment ([`minimap2`](https://lh3.github.io/minimap2/))
    2. Variant calling ([`CliqueSNV`](https://github.com/vtsyvina/CliqueSNV)/[`HapLink.jl`](https://ksumngs.github.io/HapLink.jl))
7. Haplotype calling ([`CliqueSNV`](https://github.com/vtsyvina/CliqueSNV)/[`HapLink.jl`](https://ksumngs.github.io/HapLink.jl))
8. Phylogenetic tree generation
    1. Multiple sequence alignment ([`MAFFT`](https://mafft.cbrc.jp/alignment/software/))
    2. Maximum-likelihood phylogenetic trees ([`RAxML-ng`](https://github.com/amkozlov/raxml-ng))
9. Output visualization
    - Haplotypes table ([`BioJulia`](https://biojulia.net))
    - Read QC results ([`MultiQC`](http://multiqc.info/))
    - Metagenomic classifications ([`Krona`](https://github.com/marbl/Krona/wiki/KronaTools))
    - Alignments ([`IGV`](https://igv.org/))
    - Phylogenetic tree ([`phylotree.js`](https://github.com/veg/phylotree.js))

```mermaid
flowchart TD
    A([--input]) --> B[Quality analysis]
    A --> C[Quality trimming]
    C --> D[Read classification]
    D --> E[Host read removal]
    E --> F[Alignment]
    F --> G[Consensus Sequence]

    REF.A([--genome]) --> REF.B[(NCBI Download)]
    REF.B --> F

    CL.A([--genome_list]) --> CL.B[(NCBI Download)]
    CL.B --> CL.C[Make BLAST database]

    G --> H[BLAST]
    CL.C --> H

    E --> I[Realignment]
    CL.B --Closest reference--> I
    H --> I

    I --> J[Variant Calling]
    J --> K[Haplotype Calling]
    K --> L[Multiple sequence alignment]
    L --> M[Phylogenetic tree]
```

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```console
   nextflow run ksumngs/yavsap -profile test,YOUR_PROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOUR_PROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   ```console
   nextflow run ksumngs/yavsap -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input . --outdir <OUTDIR> --platform illumina --kraken2_db https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20210517.tar.gz --keep_taxid classified
   ```

## Documentation

The nf-core/yavsap pipeline comes with documentation about the pipeline [usage](https://ksumngs.github.io/yavsap/usage), [parameters](https://ksumngs.github.io/yavsap/parameters) and [output](https://ksumngs.github.io/yavsap/output).

## Credits

nf-core/yavsap was originally written by [Thomas A. Christensen II](https://millironx.com), under the supervision of [Rachel Palinski](https://www.vet.k-state.edu/academics/dmp/faculty-staff/faculty/palinski/) at the [Kansas State University Veterinary Diagnostic Laboratory](http://www.ksvdl.org/).

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

<!-- cspell:disable -->
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: 10.1038/s41587-020-0439-x.
> In addition, references of tools and data used in this pipeline are as follows:

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/yavsap for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->
<!-- cspell:enable -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.
