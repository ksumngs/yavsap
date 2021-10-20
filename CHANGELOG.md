# Changelog

All notable changes to the jev-analysis-pipeline will be documented in this
file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic
Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## Added

- Made a changelog (365b3f1964f24c686aedc8a8019a14287a99a0f7)

## Changed

- Documentation more fully-fleshed out
  (e9604d29edacb828154ba50d9f17cd7e0072d2ad)
- Visualizer made more modular and more complete
  (c408cc00531cde39c9df8daa8eed8fc95d20a3f5)
- Visualizer switched to https://github.com/veg/phylotree.js for phylogenetic
  tree viewing (3d50edc91d4c7479e15552ff251d0a10834dc55e)
- Debug profile output increased (ee2b3d0552b5487302ae475816a9b38ee4e97d2a)

## [0.1.0-alpha] - 2021-10-12

## Added

- Standard out logs are now printed via built-in cowsay (0defca3)
- Sphinx-based documentation (b1398a1)
- Dynamically calculated reference genome size (96524ed)
- Ability to alter the assembly mode of SPAdes (3e2c5d7)
- [MultiQC](https://multiqc.info/) support (c6e7787)
- CI pipeline testing via simulated reads and GitHub actions
- _N_-dimensional haplotype finding for ONT reads (e56899a)

## Changed

- Resource allocations have been increased (45e4ea5)
- The visualization generator code is now a part of the pipeline (d5a788f)
- Filtering and haplotyping have been moved to subworkflows
- Visualizer is now written in Pug

## Removed

- The Julia-based ONT haplotype finder no longer outputs a CSV of the linkage
  statistics

## Fixed

- Canu can now use its full memory allocation (75780a6)
- Parameter names are now consistent internally (517d0e3)
- The Kraken2 database no longer requires custom binding point configs (6b0b03f)
- Failure to build a phylogenetic tree due to few haplotypes no longer kills the
  pipeline (3d6e38e)
- Failure to _de novo_ assemble no longer kills the pipeline (097ac95)
- Remove `-k` flag from gzip commands (963125d)
- Contig files from SPAdes now contain the sample name (ad810bc)

## [0.0.1] - 2021-09-15

> **Notice:** jev-analysis-pipeline is forked from another project,
> viral-metagenomics-pipeline. While this was the first release, there are still
> 'changes' and 'removals' from that pipeline that are addressed from the point
> of the fork, and not from the complete beginning of the project.

### Added

- Pipeline published in usable format to GitHub via nextflow.config manifest
- _de novo_ assembly via [SPAdes](http://cab.spbu.ru/spades) (Illumina) and
  [Canu](https://canu.readthedocs.io) (Oxford Nanopore)
- Read alignment via [minimap2](https://github.com/lh3/minimap2)
- Haplotype calling via [CliqueSNV](https://github.com/vtsyvina/CliqueSNV)
  (Illumina) and a combination of
  [iVar](https://andersen-lab.github.io/ivar/html/manualpage.html) and custom
  Julia scripts (Oxford Nanopore)

### Changed

- viral-metagenomics-pipeline converted to jev-analysis-pipeline
- Converted to [Nextflow DSL2](https://nextflow.io/docs/latest/dsl2.html)
- Container and conda directives are based on labels and not process names

### Removed

- Krona graphs of Kraken2 output
- BLAST of assemblies and unclassified reads

[Unreleased]:
https://github.com/MillironX/jev-analysis-pipeline/compare/v0.1.0-alpha...HEAD
[0.1.0-alpha]: https://github.com/MillironX/jev-analysis-pipeline/compare/v0.0.1...v0.1.0-alpha
[0.0.1]: https://github.com/MillironX/jev-analysis-pipeline/releases/tag/v0.0.1
