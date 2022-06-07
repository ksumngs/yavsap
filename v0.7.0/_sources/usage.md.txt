# Usage

Basic jist:

```bash
nextflow run ksumngs/yavsap \
  -profile <docker,podman,singularity,INSTITUTE> \
  --platform <illumina,nanopore> --kraken2_db <kraken2_db>
```

## Input Preparation

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

<!-- markdownlint-disable line-length -->

| #samplename |                                                                                        |                                                                                         |
| ----------- | -------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| pig-serum   | `/data/run/fastq{pass,fail}/barcode01/FAP01234_{pass,fail}_barcode01_abcde01_*.fastq*` | `/data/run2/fastq{pass,fail}/barcode07/FAP01234_{pass,fail}_barcode07_abcde01_*.fastq*` |
| pig-feces   | `/data/run/fastq{pass,fail}/barcode02/FAP01234_{pass,fail}_barcode02_abcde01_*.fastq*` |                                                                                         |
| mosquito1   | `/data/run/fastq{pass,fail}/barcode03/FAP01234_{pass,fail}_barcode03_abcde01_*.fastq*` |                                                                                         |
| mosquito2   | `` ./seq-results/mosquito2/*.fastq*` ``                                                |                                                                                         |

<!-- markdownlint-enable line-length -->

#### Paired-end example

<!-- markdownlint-disable line-length -->

| #Sample   | Forward1                                 | Reverse1                                 | Forward2                              | Reverse2                              |
| --------- | ---------------------------------------- | ---------------------------------------- | ------------------------------------- | ------------------------------------- |
| pig-serum | `/basespace/run/PIG-SERUM*_R1*.fastq.gz` | `/basespace/run/PIG-SERUM*_R2*.fastq.gz` | `/dragen/run/PIG-SERUM*_R1*.fastq.gz` | `/dragen/run/PIG-SERUM*_R2*.fastq.gz` |
| pig-feces | `/basespace/run/PIG-FECES*_R1*.fastq.gz` | `/basespace/run/PIG-FECES*_R2*.fastq.gz` |                                       |                                       |
| mosquito1 | `/basespace/run/MOSQUITO1*_R1*.fastq.gz` | `/basespace/run/MOSQUITO1*_R1*.fastq.gz` |                                       |                                       |
| mosquito2 | `./seq-results/mosquito2/*_R1*.fastq.gz` | `./seq-results/mosquito2/*_R2*.fastq.gz` |                                       |                                       |

<!-- markdownlint-enable line-length -->

Once the samplesheet is constructed, pass it on the command line as:

```bash
nextflow run ksumngs/yavsap \
  -profile <docker,podman,singularity,INSTITUTE> \
  --input /path/to/sheet.tsv \
  --platform <illumina,nanopore> \
  --kraken2_db <kraken2_db>
```

## Profile Selection

Profiles allow for unique combinations of settings within a Nextflow pipeline.
For the purposes of YAVSAP, they reconfigure the pipeline to run on a particular
container engine. Whichever engine you choose must be installed and available
(e.g. `module load`) on each node that pipeline processes are running on. The
available options are

(none)

: Don't use a container engine. Requires that every tool and the right version
of every tool be installed on each machine the pipeline is running on. Don't
use this one.

docker

: Use [Docker](https://docker.com) as the container engine. Note that Docker
often requires root or nearly-root permissions that usually aren't available
on HPCs, and has a weird license that might forbid its use in commercial
settings. Works well on local machines, though.

podman

: Use [Podman](https://podman.io) instead of Docker. Podman is similar enough
to Docker that they can often be used interchangably, but doesn't require root
permissions and has a free license. Some files might not be accessible via
container on RHEL-based distros thanks to their particular SELinux
implementation.

singularity

: **Recommended**

Use the {doc}`Singularity <singularity:index>` container engine. This engine
was built with HPC use in mind, and doesn't require any special permissions to
run. Singularity's downfall is how it will expose your home directory to the
container, resulting in rare, but difficult to explain bugs when files
conflict. Every effort has been made to minimize the effects of Singularity's
file mounting in this pipeline.

INSTITUTE

: If your computing center is listed in [nf-core configs](https://github.com/nf-core/configs/)
then you can pass that name to `-profile` to have the container engine and
resource limits configured automatically.

test

: Download a test dataset from nf-core and run a sample to ensure your machine
configuration is correct. Must be used with one of the container engine
profiles.

test_nanopore

: Download a MinION test dataset and run just like `test`.

gh

: Used to limit resource usage during continuous integration on GitHub Actions.
You should never have to use these for real-life analysis.

To select a profile, you must pass the desired profile name to Nextflow's
`-profile` flag. Note that this is a Nextflow flag, and not a pipeline flag,
so it is a single dash (`-profile`), not a double dash (`--profile`).

## Mandatory Parameters

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

-p

rofile

So, this isn't really a parameter, but the container engine needs to be set
using Nextflow's `-profile` flag. See {ref}`Profile Selection`.

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

<!-- markdownlint-disable line-length -->

| Virus                                                                                                                        | Name  | Recommended Reference                                           |
| ---------------------------------------------------------------------------------------------------------------------------- | ----- | --------------------------------------------------------------- |
| [Japanese Encephalitis Virus (JEV)](https://www.ncbi.nlm.nih.gov/data-hub/taxonomy/11072/)                                   | `jev` | [NC_001437.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001437.1) |
| [Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)](https://www.ncbi.nlm.nih.gov/labs/data-hub/taxonomy/2697049/) | `sc2` | [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) |

<!-- markdownlint-enable line-length -->

## Setting up for HPC Job Schedulers

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
