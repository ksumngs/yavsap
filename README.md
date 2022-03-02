<!-- markdownlint-disable -->
<p align="center">
  <img src="assets/yavsap-logo.svg" alt="logo" width="270">
</p>

# YAVSAP (Yet Another Viral Subspecies Analysis Pipeline)

[![Testing](https://github.com/ksumngs/yavsap/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/ksumngs/yavsap/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://ksumngs.github.io/yavsap)
[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/ksumngs/yavsap?label=version)](https://github.com/ksumngs/yavsap/blob/master/CHANGELOG.md)
[![GitHub license](https://img.shields.io/github/license/ksumngs/yavsap)](https://github.com/ksumngs/yavsap/blob/master/LICENSE)

<!-- markdownlint-enable -->

A [Nextflow] pipeline for studying viral populations within a single sample,
tuned for [Japanese Encephalitis Virus]. :dna::computer::chart_with_upwards_trend:
Yeah, we're still looking for a better name. :shrug:

> This project follows the [semver] _pro forma_ and uses the [git-flow]
> branching model.

## Installation

1. Install [Nextflow] (>= 21.10.6)
2. Install [Conda]
3. Install one or more of
   - [Singularity] (**Recommended**)
   - [Podman]
   - [Docker]
4. Download a [Kraken2 database]

Check out the [Installation] docs for a more nuanced take on the requirements.

## Usage

### Syntax

```bash
nextflow run ksumngs/yavsap              \
  -profile <singularity,podman,docker>   \
  --platform <illumina,nanopore>         \
  --kraken2_db /path/to/kraken2/database \
  [--input /path/to/reads/folder]        \
  [--genome accession_number]            \
  [--keep_taxid list]                    \
  [--outdir /path/to/output]
```

### Example: Illumina reads with a Kraken2 database containing the host

```bash
nextflow run ksumngs/yavsap \
  -profile singularity      \
  --platform illumina       \
  --kraken2_db /databases/kraken2/nt
```

### Example: Nanopore reads with a viral-only Kraken2 database

```bash
nextflow run ksumngs/yavsap             \
  -profile podman                       \
  --platform nanopore                   \
  --kraken2_db /databases/kraken2/viral \
  --keep_taxid '10239'
```

### Example: Illumina reads aligned against a different reference genome

```bash
nextflow run ksumngs/yavsap                                \
  -profile docker                                          \
  --platform illumina                                      \
  --kraken2_db /databases/kraken2/refseq-complete_unmasked \
  --genome 'KT957423.1'
```

There are _way_ more parameters than listed here. For a more complete
description, please read the docs on [Usage] and [Parameters].

## Process Summary

Here's what happens to your reads in the pipeline.

 1. Quality analysis ([FastQC])
 2. Quality trimming ([Trimmomatic]/[NanoFilt])
 3. Read classification ([Kraken2])
 4. Host read removal ([KrakenTools])
 5. _de novo_ assembly of viral reads ([SPAdes]/[Canu])
 6. Alignment of assemblies against the reference genome ([minimap2])
 7. Alignment of reads against the reference genome ([minimap2])
 8. Variant calling ([iVar])
 9. Haplotype calling ([CliqueSNV]/custom)
10. Haplotype alignment ([MAFFT])
11. Phylogenetic tree generation ([raxml-ng])
12. Alignments and phylogenetics output to browser ([IGV]+[phylotree.js])

[Canu]: https://canu.readthedocs.io
[CliqueSNV]: https://github.com/vtsyvina/CliqueSNV
[Conda]: https://conda.io/miniconda.html
[Docker]: https://docs.docker.com/engine/installation
[FastQC]: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
[git-flow]: https://nvie.com/posts/a-successful-git-branching-model
[IGV]: https://igv.org/
[Installation]: https://ksumngs.github.io/yavsap/install
[iVar]: https://andersen-lab.github.io/ivar/html/manualpage.html
[Japanese Encephalitis Virus]: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=11072
[Kraken2 database]: https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases
[Kraken2]: https://github.com/DerrickWood/kraken2/wiki
[KrakenTools]: https://github.com/jenniferlu717/KrakenTools
[MAFFT]: https://mafft.cbrc.jp/alignment/software/
[minimap2]: https://lh3.github.io/minimap2/
[NanoFilt]: https://github.com/wdecoster/nanofilt/
[Nextflow]: https://nextflow.io
[Parameters]: https://ksumngs.github.io/yavsap/parameters
[phylotree.js]: https://github.com/veg/phylotree.js
[Podman]: https://podman.io
[raxml-ng]: https://github.com/amkozlov/raxml-ng
[semver]: https://semver.org
[Singularity]: https://www.sylabs.io/guides/3.8/user-guide
[SPAdes]: cab.spbu.ru/spades
[Trimmomatic]: www.usadellab.org/cms/?page=trimmomatic
[Usage]: https://ksumngs.github.io/yavsap/usage
