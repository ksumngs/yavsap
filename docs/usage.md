<!--markdownlint-disable line-length -->
# Usage

Basic jist:

```bash
nextflow run ksumngs/yavsap \
  -profile <docker,podman,singularity,INSTITUTE> \
  --platform <illumina,nanopore> --kraken2_db <kraken2_db>
```

## Input preparation

### Using a Directory as Input

If the {ref}`--input <Input/Output Options>` parameter is not passed, then
YAVSAP assumes reads files are located in the current directory. You can also
pass a directory that contains reads files to `--input`. YAVSAP will not look
in subdirectories of `--input`. When running YAVSAP on a directory, Reads
files must have one of one of the following extensions:

- .fastq
- .fq
- .fastq.gz
- .fq.gz

There must be only one file per sample for single-end and inteleaved paired-end
reads, and one pair of files (two files) for regular paired-end reads. Regular
paired-end files must have either the characters `_1` and `_2` or `_R1`
and `_R2` in the filename to be identified as a pair. The sample name is taken
as the all the characters _before_ the first underscore in the file name.

For example, a folder with the following files

```text
.
├── pig-serum_S4_L001_R1_001.fastq.gz
├── pig-serum_S4_L001_R2_001.fastq.gz
├── pig-feces_S5_L001_R1_001.fastq.gz
├── pig-feces_S5_L001_R2_001.fastq.gz
├── mosquito1_S6_L001_R1_001.fastq.gz
├── mosquito1_S6_L001_R2_001.fastq.gz
├── mosquito2_S7_L001_R1_001.fastq.gz
└── mosquito2_S7_L001_R2_001.fastq.gz
```

will produce a sample list of

- pig-serum
- pig-feces
- mosquito1
- mosquito2

Paired-endedness is determined based on the {ref}`--paired <Input/Output Options>` flag. By default `--paired` is `true` when {ref}`--platform <Input/Output Options>` is `illumina` and `false` for `nanopore`, but this
can be overridden. To use interleaved paired-end reads, use the
{ref}`--interleaved <Input/Output Options>` flag.

This file naming system is very common for Illumina output, however for Oxford
Nanopore reads, the default file structure breaks this in several ways, and it is
often better to use a samplesheet.

### Using a Samplesheet as Input

A samplesheet allows you to pass multiple reads files in as
a single sample, which can be useful for e.g. Oxford Nanopore reads. When using
a samplesheet, simply pass the path to the sheet to {ref}`--input <Input/Output Options>`.

Samplesheets are tab-delimited files. A header is optional, and must be
delineated by beginning with the pound (`#`) symbol. The first column contains
sample names, and each subsequent column contains a path to reads files. All
sample names will be cleaned of any shell metacharacters, and will be truncated
at the first underscore. Paths are resolved relative to the current directory.
Path names can contain wildcard characters, and will resolve as glob patterns
identical to how they would in Nextflow's {ref}`Channel.fromPath <channel-path>`
constructor. For paired-end reads, forward reads should be specified in the
even-number columns (`2,4,6`), and reverse reads should be specified in the
odd-number columns (`3,5,7`).

#### Single-end example

| #samplename |                                                                                        |                                                                                         |
| ----------- | -------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| pig-serum   | `/data/run/fastq{pass,fail}/barcode01/FAP01234_{pass,fail}_barcode01_abcde01_*.fastq*` | `/data/run2/fastq{pass,fail}/barcode07/FAP01234_{pass,fail}_barcode07_abcde01_*.fastq*` |
| pig-feces   | `/data/run/fastq{pass,fail}/barcode02/FAP01234_{pass,fail}_barcode02_abcde01_*.fastq*` |                                                                                         |
| mosquito1   | `/data/run/fastq{pass,fail}/barcode03/FAP01234_{pass,fail}_barcode03_abcde01_*.fastq*` |                                                                                         |
| mosquito2   | `` ./seq-results/mosquito2/*.fastq*` ``                                                |                                                                                         |

#### Paired-end example

| #Sample   | Forward1                                 | Reverse1                                 | Forward2                              | Reverse2                              |
| --------- | ---------------------------------------- | ---------------------------------------- | ------------------------------------- | ------------------------------------- |
| pig-serum | `/basespace/run/PIG-SERUM*_R1*.fastq.gz` | `/basespace/run/PIG-SERUM*_R2*.fastq.gz` | `/dragen/run/PIG-SERUM*_R1*.fastq.gz` | `/dragen/run/PIG-SERUM*_R2*.fastq.gz` |
| pig-feces | `/basespace/run/PIG-FECES*_R1*.fastq.gz` | `/basespace/run/PIG-FECES*_R2*.fastq.gz` |                                       |                                       |
| mosquito1 | `/basespace/run/MOSQUITO1*_R1*.fastq.gz` | `/basespace/run/MOSQUITO1*_R1*.fastq.gz` |                                       |                                       |
| mosquito2 | `./seq-results/mosquito2/*_R1*.fastq.gz` | `./seq-results/mosquito2/*_R2*.fastq.gz` |                                       |                                       |

Once the samplesheet is constructed, pass it on the command line as:

```bash
nextflow run ksumngs/yavsap \
  -profile <docker,podman,singularity,INSTITUTE> \
  --input /path/to/sheet.tsv \
  --platform <illumina,nanopore> \
  --kraken2_db <kraken2_db>
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run ksumngs/yavsap -profile docker --input . --outdir <OUTDIR> --platform illumina --kraken2_db https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20210517.tar.gz --keep_taxid classified
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work                # Directory containing the nextflow working files
<OUTIDR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```console
nextflow pull ksumngs/yavsap
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [yavsap releases page](https://github.com/ksumngs/yavsap/releases) and find the latest version number (eg. `v1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r v1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/)<!-- and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/)-->.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility<!--, however when this is not possible, Conda is also supported-->.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

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
- ~~`conda`~~
  - ~~A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.~~ **YAVSAP doesn't currently have conda support, but we're working on it!. Don't use this profile!**
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Mandatory pipeline parameters

See {doc}`the page on parameters <parameters>` for the complete lowdown on
parameters.

The pipeline is pretty much set up to run itself. So long as you have your input
reads formatted correctly, it doesn't need much input from you. (Provided, the
defaults are for analyzing Japanese Encephalitis Virus, so your results might be
way off if you don't have any JEV in your samples, see {ref}`Genome Preparation`.) These are the bare minimum parameters that you must provide on
the command-line for the pipeline to complete. Note that there has been mixed
success with placing these parameters in a `nextflow.config` file, so keeping
them on the command-line is best.

--kraken2_db

The path to a Kraken2 database. See {ref}`--kraken2_db <Kraken2 Options>`.

--platform

Must be set to `illumina` or `nanopore`, depending on the type of reads
you are analyzing. See {ref}`--platform <Input/Output Options>`.

## Genome Preparation

YAVSAP needs some idea of what kind of viruses it's looking for. You will need
to provide it with

1. A reference genome of the virus of interest
2. Genome examples of different strains or relatives of the virus of interest

YAVSAP needs these in the form of NCBI GenBank Accession Numbers. To pick a new
reference genome, simply pass the accession number to the {ref}`--genome <Reference Genome Options>` parameter, e.g. `nextflow run ksumngs/yavsap --genome 'NC_001437.1'`.

Close relative examples must be provided as a tab-delimited document with the
name in the first column, and the NCBI accession number in the second column. No
header is allowed in the comparison genomes file.

Nodes will be grouped in the output phylogenetic tree based on the name before
any underscore, to allow for easy visualization of known related groups. In
other words, `Nakayama\tEF571853.1` and `SA14\tMH258848.1` are perfectly
valid entries in the relatives list, but `GIII_Nakayama\tEF571853.1` and
`GIII_SA14\tMH258848.1` will give these two entries the same color so you can
easily see the Genotype III node.

There **MUST** be one entry whose name begins with `ROOT`. The genome of this
entry will be used as the outgroup for the resulting phylogenetic tree.

Once a comparison genomes file is prepared, the path to it can be passed to the
{ref}`--genome_list <Reference Genome Options>` parameter.

YAVSAP comes with some example comparison genomes files. These can be referred
to by name, rather than by path(e.g. `nextflow run ksumngs/yavsap --genome_list jev`). They are

| Virus                                                                                                                        | Name  | Recommended Reference                                           |
| ---------------------------------------------------------------------------------------------------------------------------- | ----- | --------------------------------------------------------------- |
| [Japanese Encephalitis Virus (JEV)](https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/11072/)                                   | `jev` | [NC_001437.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001437.1) |
| [Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)](https://www.ncbi.nlm.nih.gov/labs/data-hub/taxonomy/2697049/) | `sc2` | [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) |

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```

## Setting up for HPC job schedulers

YAVSAP comes preconfigured for local use only. Yes, that's about as ridiculous
as it sounds. What's even more ridiculous is trying to make a configuration that
can easily be adapted to multiple HPCs and job-scheduler frameworks. There is a
compromise, however.

### Process Labels

Rather than provide hard-coded configurations that will certainly
break, there are several Nextflow 'labels,' that can be used for assigning
processes to specific node queues or partitions. The labels are

process_low

: For processes with low resource usage that take a short time

process_medium

: For processes with moderate resource usage that take a moderate amount of time

process_high

: For processes with high resource usage that take a moderately high amount of
time

process_long

: For processes that take a long time

process_high_memory

: For processes that use a lot of memory

run_local

: For processes that have to be run on the login node. This label was created
specially for processes that download resources from the internet on HPCs
where compute nodes do not have internet access

Using a custom `nextflow.config` and these process labels, you can construct a
setup for your own HPC needs.

### Example

As an example, here is a guide on how to set up a configuration for the
[USDA's SCINet Ceres cluster](https://scinet.usda.gov/guide/ceres/), using the
publically available info on their website.

First, we see that Ceres uses SLURM and Singularity. Excellent.
Let's set Nextflow up to use SLURM:

```groovy
process {
    executor = 'slurm'
}
```

Some SLURM systems require an account to be passed with every job submission.
Let's add ours just to be safe.

```groovy
process {
    executor       = 'slurm'
    clusterOptions = '--account=millironx'
}
```

For this example, I don't think we'll need to do anything special with the low,
medium, and high processes, but let's make sure that the long and high memory
processes get submitted to partitions that can handle them.

```groovy
process {
    executor       = 'slurm'
    clusterOptions = '--account=millironx'
    module         = 'singularity'
    withLabel: process_long {
        queue      = 'long'
    }
    withLabel: process_high_memory {
        queue      = 'mem'
            }
}
```

Now, we can place this file in `$HOME/.nextflow/nextflow.config`, and these
settings will be applied every time we run the pipeline.
