# Changelog

All notable changes to YAVSAP will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.1] - 2022-08-01

## Fixed

- Errors from incorrect CliqueSNV VCF output

## [0.9.0] - 2022-07-29

## Added

- Multi-strand reference genome support ([#46](https://github.com/ksumngs/yavsap/pull/46)
  fixes [#25](https://github.com/ksumngs/yavsap/issues/25))
- Genome annotations downloading ([#50](https://github.com/ksumngs/yavsap/pull/50))
- Table of variant calls to report ([#51](https://github.com/ksumngs/yavsap/pull/51)
  fixes [#42](https://github.com/ksumngs/yavsap/issues/42))

## Changed

- nf-core template updated to v1.4.1
- Update e-direct containers ([#48](https://github.com/ksumngs/yavsap/pull/48))

## [0.8.0] - 2022-06-06

### Changed

- Visualizer created entirely in MultiQC ([#46](https://github.com/ksumngs/yavsap/pull/46))
- nf-core template updated to v1.4.0 ([#46](https://github.com/ksumngs/yavsap/pull/46))
- Consensus sequences now called using the variant callers for Illumina/Nanopore
  reads

## [0.7.1] - 2022-04-25

### Fixed

- No visualizer output for Nanopore reads ([#41](https://github.com/ksumngs/yavsap/issues/41))

## [0.7.0] - 2022-04-19

### Changed

- Visualizer changed to single, static page ([#40](https://github.com/ksumngs/yavsap/pull/40))

## [0.6.5-alpha] - 2022-03-24

### Fixed

- Bug that prevented external genome tables again

## [0.6.4-alpha] - 2022=03-23

### Fixed

- Kraken2 memory allocates based on the database size again

## [0.6.3-alpha] - 2022-03-23

### Fixed

- Bug introduced for internal genome table files in v0.6.2

## [0.6.2-alpha] - 2022-03-07

### Fixed

- Workflow no longer fails when passed an external genome table file

## [0.6.1-alpha] - 2022-03-07

### Fixed

- `test` profile parameters are no longer applied by default

## [0.6.0-alpha] - 2022-03-03

This is a major overhaul of YAVSAP to make it use nf-core's DSL2 modules.

### Added

- Interleaved and Single-end Illumina read support (#34)
- nf-core compliance
- SARS-CoV2 genome preset
- Support for tarballed and/or remote Kraken2 database

### Changed

- Multiple sequence alignments and phylogenetic trees now include all samples
  collectively

### Removed

- `--pe` and `--ont` parameters (#34)
- _de novo_ assembly steps and parameters

## [0.5.0-alpha] - 2022-02-09

### Added

- Samplesheet compatability (#21/#31)
- Ability to use user-created genome lists (#23)
- Nanopore read QC via NanoStat and MultiQC (#33)

### Fixed

- Oversubscription error in RAxML-NG process
- Reroot functionality in phylogenetic tree viewer (#12/#29)

## [0.4.0-alpha] - 2022-01-11

### Added

- VCF output for Illumina reads via CliqueSNV (#18)
- Parameter schema in
  [Nextflow schema](https://help.tower.nf/pipeline-schema/overview/) format
  (#19)
- Phylogenetic tree quality cutoff parameter (#20)

### Changes

- HapLink.jl updated to 0.4 (#16)
- ONT reads now simulated using https://github.com/yukiteruono/pbsim2 (#17)
- Docs now stay in sync with docstrings/`nextflow_schema.json` (#10)

### Fixed

- It is now nearly impossible to receive cooldown timeouts from NCBI (#15)

## [0.3.3-alpha] - 2021-12-02

### Fixed

- Container for realigning reads fixed (#9)

## [0.3.2-alpha] - 2021-12-02

### Fixed

- File copy no longer fails under Podman profile/container engine (#8)

## [0.3.1-alpha] - 2021-12-01

### Changed

- HapLink.jl scripts outsourced to container (#7)

## [0.3.0-alpha] - 2021-12-01

### Added

- Kraken2 database download instructions to docs (547f67b)
- Parameter option for number of phylogenetic bootstrap trees (474652f)
- Options to skip read trimming (`--skip_trimming`) and read QC (`--skip_qc`)
  (#5)
- Realignment to closest reference genome based on BLAST (#6)

### Changed

- Dynamically calculate Kraken2's memory requirements (5d9a873)
- Processes using gzip sped up by removing `-9` flag (e3c2a13)
- ONT Haplotyping is now performed by the external library
  https://github.com/ksumngs/HapLink.jl (#2)
- Test suite updated to more predicable haplotypes (3a3f30b)

### Removed

- All `conda`-dependent processes (#3)

## [0.2.1-alpha] - 2021-11-23

### Fixed

- Kraken2 memory is now allocated correctly even when running under a `test`
  profile (12aea9c)

## [0.2.0-alpha] - 2021-11-04

### Added

- Made a changelog (365b3f1)
- Made a markdownlint config file (817f7e2)
- FastQC analysis added to pipeline and visualizer (6fb498f)

### Changed

- Documentation more fully-fleshed out (e9604d2)
- Visualizer made more modular and more complete (c408cc0)
- Visualizer switched to https://github.com/veg/phylotree.js for phylogenetic
  tree viewing (3d50edc)
- Debug profile output increased (ee2b3d0)
- Config files made more modular (f82aedc)
- Pipeline rename to YAVSAP

### Fixed

- Trimming ONT reads no longer removes variant features (152eafe)

## [0.1.0-alpha] - 2021-10-12

### Added

- Standard out logs are now printed via built-in cowsay (0defca3)
- Sphinx-based documentation (b1398a1)
- Dynamically calculated reference genome size (96524ed)
- Ability to alter the assembly mode of SPAdes (3e2c5d7)
- [MultiQC](https://multiqc.info/) support (c6e7787)
- CI pipeline testing via simulated reads and GitHub actions
- _N_-dimensional haplotype finding for ONT reads (e56899a)

### Changed

- Resource allocations have been increased (45e4ea5)
- The visualization generator code is now a part of the pipeline (d5a788f)
- Filtering and haplotyping have been moved to subworkflows
- Visualizer is now written in Pug

### Removed

- The Julia-based ONT haplotype finder no longer outputs a CSV of the linkage
  statistics

### Fixed

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
> viral-metagenomics-pipeline (now v-met). While this was the first release,
> there are still 'changes' and 'removals' from that pipeline that are addressed
> from the point of the fork, and not from the complete beginning of the
> project.

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

[unreleased]: https://github.com/ksumngs/yavsap/compare/v0.9.1...HEAD
[0.9.1]: https://github.com/ksumngs/yavsap/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/ksumngs/yavsap/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/ksumngs/yavsap/compare/v0.7.1...v0.8.0
[0.7.1]: https://github.com/ksumngs/yavsap/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/ksumngs/yavsap/compare/v0.6.5-alpha...v0.7.0
[0.6.5-alpha]: https://github.com/ksumngs/yavsap/compare/v0.6.4-alpha...v0.6.5-alpha
[0.6.4-alpha]: https://github.com/ksumngs/yavsap/compare/v0.6.3-alpha...v0.6.4-alpha
[0.6.3-alpha]: https://github.com/ksumngs/yavsap/compare/v0.6.2-alpha...v0.6.3-alpha
[0.6.2-alpha]: https://github.com/ksumngs/yavsap/compare/v0.6.1-alpha...v0.6.2-alpha
[0.6.1-alpha]: https://github.com/ksumngs/yavsap/compare/v0.6.0-alpha...v0.6.1-alpha
[0.6.0-alpha]: https://github.com/ksumngs/yavsap/compare/v0.5.0-alpha...v0.6.0-alpha
[0.5.0-alpha]: https://github.com/ksumngs/yavsap/compare/v0.4.0-alpha...v0.5.0-alpha
[0.4.0-alpha]: https://github.com/ksumngs/yavsap/compare/v0.3.3-alpha...v0.4.0-alpha
[0.3.3-alpha]: https://github.com/ksumngs/yavsap/compare/v0.3.2-alpha...v0.3.3-alpha
[0.3.2-alpha]: https://github.com/ksumngs/yavsap/compare/v0.3.1-alpha...v0.3.2-alpha
[0.3.1-alpha]: https://github.com/ksumngs/yavsap/compare/v0.3.0-alpha...v0.3.1-alpha
[0.3.0-alpha]: https://github.com/ksumngs/yavsap/compare/v0.2.1-alpha...v0.3.0-alpha
[0.2.1-alpha]: https://github.com/ksumngs/yavsap/compare/v0.2.0-alpha...v0.2.1-alpha
[0.2.0-alpha]: https://github.com/ksumngs/yavsap/compare/v0.1.0-alpha...v0.2.0-alpha
[0.1.0-alpha]: https://github.com/ksumngs/yavsap/compare/v0.0.1...v0.1.0-alpha
[0.0.1]: https://github.com/ksumngs/yavsap/releases/tag/v0.0.1
