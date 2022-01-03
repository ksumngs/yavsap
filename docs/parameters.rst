Complete Parameter Reference
============================

Every single parameter in the entire pipeline is described here. If you need to
tweak every little detail of how your JEV samples are analyzed, you've come to
the right place! Note that *every* parameter is described here, even those that
shouldn't be used, so proceed with caution!

Input Options
-------------

--input
^^^^^^^

======== ======
Type     String
======== ======
Required No
Default  .
======== ======

The path to the folder containing input reads. For Illumina (paired-end) reads,
the file names must be identical until the ending underscore with a read number,
e.g. 'sample1_S10_L001_R1_001.fastq.gz' and 'sample1_S10_L001_R2_001.fastq.gz'.
The read number must be designated using either '_1' and '_2' or '_R1' and
'_R2'. For Nanopore reads, each fastq file is assumed to be a unique sample, so,
e.g. 'FAP01234_pass_barcode05_abcd01234_0.fastq.gz' and
'FAP01234_pass_barcode05_abcd01234_1.fastq.gz' are assumed to be different
samples even though they are from the same barcode. All read files must be
gzipped, and have the extension '.fastq.gz' or '.fq.gz'.

Defaults to the current working directory.

--sra
^^^^^

======== ======
Type     Flag
======== ======
Required No
Default  false
======== ======

This flag switches the meaning of :ref:`--input` to mean an NCBI Short Read Archive
(SRA) accession number, then pulls the associated files directly from NCBI and
runs the pipeline on them. To use this flag, you **must** have an `NCBI API key
<https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/>`_
and it **must** be exported in the shell environment as ``NCBI_API_KEY``, e.g.
``NCBI_API_KEY=0123456789abcdef``. This feature is currenly broken due to an
upstream issue in Nextflow, but should be fixed by Nextflow version 21.10.0.

--platform
^^^^^^^^^^

======== ======
Type     String
======== ======
Required Yes
Default  n/a
======== ======

Defines which type of reads to analyze. Must be either 'illumina' or 'nanopore'.


--pe
^^^^

======== ======
Type     Boolean
======== ======
Required No
Default  false
======== ======

.. note:: This parameter is hidden and should **never** be used on the command
    line! Bad things will happen if you do!

Internally-used parameter to check if the pipeline is processing Illumina reads.
Automatically set by :ref:`--platform`

--ont
^^^^^

======== ======
Type     Boolean
======== ======
Required No
Default  false
======== ======

.. note:: This parameter is hidden and should **never** be used on the command
    line! Bad things will happen if you do!

Internally-used parameter to check if the pipeline is processing Nanopore reads.
Automatically set by :ref:`--platform`

--help
^^^^^^

======== ======
Type     Boolean
======== ======
Required No
Default  false
======== ======

Print a friendly help message that is both less daunting and less complete than
the document you're reading right now, and exit.

Output Options
--------------

--outdir
^^^^^^^^

======== ======
Type     String
======== ======
Required No
Default  ./results
======== ======

Where to store the results files. For more info, see
:doc:`the page on output <output>`.

--tracedir
^^^^^^^^^^

======== ======
Type     String
======== ======
Required No
Default  :ref:`--outdir`/.trace
======== ======

Where to store pipeline performance reports and diagnostic info. If you keep
this as the default, then you can access these documents through the visualizer.

--publish_dir_mode
^^^^^^^^^^^^^^^^^^

======== ======
Type     String
======== ======
Required No
Default  copy
======== ======

How to take files out of the ``work`` dirctories they were generated in and put
them into :ref:`--outdir`. Supports every mode that the
:ref:`Nextflow publishDir directive <process-publishdir>`
does, which as of Nextflow 21.04, includes

* symlink
* relink
* link
* copy
* copyNoFollow
* move

Reference Genome Options
------------------------

--genome
^^^^^^^^

======== ======
Type     String
======== ======
Required No
Default  NC_001437.1
======== ======

The NCBI genbank accession number of the reference genome to align reads
against.

Defaults to the accession number of the Japanese Encephalitis Virus RefSeq
record.

--genome_list
^^^^^^^^^^^^^

======== ======
Type     String
======== ======
Required No
Default  jev
======== ======

The name of the file that contains NCBI accession numbers of related strains for
checking closest reference and constructing phylogenetic trees.

Kraken2 Options
---------------

--kraken2_db
^^^^^^^^^^^^

======== ======
Type     String
======== ======
Required Yes
Default  n/a
======== ======

The path to a `Kraken2 database
<https://github.com/DerrickWood/kraken2/wiki/Manual#custom-databases>`_ that
will be used to filter out host reads in the pipeline. This path will be
automatically mounted into the container environments if a containerized profile
is used.

Corresponds to the
`--db <https://github.com/DerrickWood/kraken2/wiki/Manual#classification>`_
option of Kraken2.

--keep_taxid
^^^^^^^^^^^^
======== ======
Type     String
======== ======
Required No
Default  0 10239
======== ======

A space-separated list (use quotes on the command line), of the taxonomic ids to
keep based on Kraken2's classification.

Defaults to keeping all unclassified reads and all viral reads. Note that this
requires the host to be present in the Kraken2 database. When dealing with
animals and the databases available from ``kraken2-build``, this is not the
case, and this parameter should be modified.

Read Trimming Options
---------------------

Common Options
^^^^^^^^^^^^^^

--trim_minlen
"""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  100/300 (Illumina/Nanopore)
======== ======

Remove reads that are shorter than this length in bases.

Corresponds to the `MINLEN <http://www.usadellab.org/cms/?page=trimmomatic>`_
option of Trimmomatic for Illumina reads.

Corresponds to the `--length <https://github.com/wdecoster/nanofilt/#usage>`_
option of NanoFilt for Nanopore reads.

--trim_headcrop
"""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  0
======== ======

The number of bases to remove from the start of the read.

Corresponds to the `HEADCROP <http://www.usadellab.org/cms/?page=trimmomatic>`_
option of Trimmomatic for Illumina reads.

Corresponds to the `--headcrop <https://github.com/wdecoster/nanofilt/#usage>`_
option of NanoFilt for Nanopore reads.


Illumina-Specific (Trimmomatic) Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

--trim_adapters
"""""""""""""""

======== ======
Type     String
======== ======
Required No
Default  NexteraPE-PE.fa
======== ======

Illumina adapters to be removed during trimming.

Due to the way the container is built, custom adapters cannot be used, and this
option **must** be set to one of the following

* NexteraPE-PE.fa
* TruSeq2-PE.fa
* TruSeq3-PE-2.fa
* TruSeq3-PE.fa

If left blank (i.e. ``--trim_adapters ''``), then adapter trimming is disabled.

Corresponds to the first
`ILLUMINACLIP <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of
Trimmomatic.

--trim_mismatches
"""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  2
======== ======

The maximum mismatch count which will still allow a full adapter match to be
performed. If set to ``0``, then adapter trimming is disabled.

Corresponds to the second
`ILLUMINACLIP <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of
Trimmomatic.

--trim_pclip
""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  30
======== ======

``pclip``: palindrome clip. How accurate the match between the two adapter
ligated reads must be for paired-end palindrome read alignment. If set to ``0``,
then adapter trimming is disabled.

Corresponds to the third
`ILLUMINACLIP <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of
Trimmomatic.

--trim_clip
"""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  10
======== ======

How accurate the match between any adapter sequence must be against a read. If
set to ``0``, then adapter trimming is disabled.

Corresponds to the final
`ILLUMINACLIP <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of
Trimmomatic.

--trim_winsize
""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  50
======== ======

The number of bases to average quality accross during sliding window trimming.
If set to ``0``, then sliding window trimming is disabled.

Corresponds to the first
`SLIDINGWINDOW <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of
Trimmomatic.

--trim_winqual
""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  15
======== ======

The minimum average quality within the sliding window to keep a read. If set to
``0``, then sliding window trimming is disabled.

Corresponds to the second
`SLIDINGWINDOW <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of
Trimmomatic.

--trim_leading
""""""""""""""

======== ======
Type     Integer
======== ======
Required NoFloat
Default  15
======== ======

The minimum quality to keep a base in the leading end of a read. If set to
``0``, LEADING trimming is disabled.

Corresponds to the
`LEADING <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of
Trimmomatic.

--trim_trailing
"""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  15
======== ======

The minimum quality to keep a base in the trailing end of a read. If set to
``0``, TRAILING trimming is disabled.

Corresponds to the
`TRAILING <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of
Trimmomatic.

--trim_crop
"""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  0
======== ======

The number of bases to keep from the start of the read. If set to ``0``, CROP
trimming is disabled.

Corresponds to the
`CROP <http://www.usadellab.org/cms/?page=trimmomatic>`_ option of Trimmomatic.

Nanopore-Specific (NanoFilt) Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

--trim_maxlen
"""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  0
======== ======

Remove reads that are longer than this number of bases.

Corresponds to the
`--maxlength <https://github.com/wdecoster/nanofilt/#usage>`_ option of
NanoFilt.

--trim_meanqual
"""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  7
======== ======

Remove reads where the average basecall quality is lower than this.

Corresponds to the
`--quality <https://github.com/wdecoster/nanofilt/#usage>`_ option of NanoFilt.

--trim_mingc
"""""""""""""""

======== ======
Type     Float
======== ======
Required No
Default  0
======== ======

Remove reads that don't have at least this fraction of GC content.

Corresponds to the
`--minGC <https://github.com/wdecoster/nanofilt/#usage>`_ option of NanoFilt.

--trim_maxgc
"""""""""""""""

======== ======
Type     Float
======== ======
Required No
Default  0
======== ======

Remove reads that have more than this fraction of GC content. If set to ``0``,
then there is no upper limit of GC content.

Corresponds to the
`--maxGC <https://github.com/wdecoster/nanofilt/#usage>`_ option of NanoFilt.

--trim_tailcrop
"""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  0
======== ======

Remove this many bases from the end of each read.

Corresponds to the
`--tailcrop <https://github.com/wdecoster/nanofilt/#usage>`_ option of NanoFilt.

*de novo* Assembly Options
--------------------------

Illumina-Specific (SPAdes) Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

--spades_mode
"""""""""""""

======== ======
Type     String
======== ======
Required No
Default  rnaviral
======== ======

If provided, this parameter is turned into a flag and passed as the 'mode' to
the SPAdes assembly, e.g.
``nextflow run ksumngs/yavsap --spades_mode 'metaviral'`` will run SPAdes
as ``spades.py --metaviral``. The available modes are

* meta
* plasmid
* metaplasmid
* metaviral
* rna
* rnaviral

Due to parameter mismatches, the ``isolate`` and ``bio`` modes normally provided
by SPAdes are unavailable in the pipeline.

See `SPAdes command line options <https://cab.spbu.ru/files/release3.15.3/manual.html#sec3.2>`_
for more info on what each of these mean.

Nanopore-Specific (Canu) Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

--canu_corrected_error_rate
"""""""""""""""""""""""""""

======== ======
Type     Float
======== ======
Required No
Default  0.144
======== ======

How dissimilar overlap between two reads can be and still be assembled together.

Corresponds to the :ref:`correctedErrorRate` option of Canu.

--canu_min_read_length
""""""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  1000
======== ======

The shortest length allowed for a read to be considered in the assembly.

Corresponds to the :ref:`minReadLength` option of Canu.

--canu_min_overlap_length
"""""""""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  :ref:`--canu_min_read_length` ÷ 2
======== ======

The shortest length allowed for reads to overlap to be considered in the assembly.

Corresponds to the :ref:`minOverlapLength` option of Canu.

--canu_stop_on_low_coverage
"""""""""""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  10
======== ======

Lowest depth allowable for assembly to proceed.

Corresponds to the :ref:`stopOnLowCoverage` option of Canu.

Variant Calling Options
-----------------------

--variant_quality
^^^^^^^^^^^^^^^^^

======== ======
Type     Integer
======== ======
Required No
Default  30/21 (Illumina/Nanopore)
======== ======

Minimum average quality (PHRED score) of a basecall to consider it a variant.

--variant_depth
^^^^^^^^^^^^^^^

======== ======
Type     Integer
======== ======
Required No
Default  10/15 (Illumina/Nanopore)
======== ======

Minimum depth to consider a variant

--variant_position
^^^^^^^^^^^^^^^^^^

======== ======
Type     Float
======== ======
Required No
Default  0.1
======== ======

Remove variants that occur only in positions within this percentage of the end.

--variant_frequency
^^^^^^^^^^^^^^^^^^^

======== ======
Type     Float
======== ======
Required No
Default  0.05
======== ======

Minimum frequency at which a variant must appear.

--variant_significance
^^^^^^^^^^^^^^^^^^^^^^

======== ======
Type     Float
======== ======
Required No
Default  1e-3
======== ======

The highest p-value that will be considered a significant variant based on
Fisher's Exact test.

Haplotyping Options
-------------------

--haplotype_significance
^^^^^^^^^^^^^^^^^^^^^^^^

======== ======
Type     Float
======== ======
Required No
Default  0.05
======== ======

The highest p-value that will be considered a significant haplotype based on
a Χ-squared test of linkage disequilibrium.

--haplotype_depth
^^^^^^^^^^^^^^^^^

======== ======
Type     Integer
======== ======
Required No
Default  10
======== ======

The minimum number of times a particular haplotype has to occur for it to be
considered real and processed downstream and output.

--haplotype_frequency
^^^^^^^^^^^^^^^^^^^^^

======== ======
Type     Integer
======== ======
Required No
Default  10
======== ======

**The minimum haplotype frequency**

Only applies to Illumina/PE reads. Corresponds to the
`-tf <https://github.com/vtsyvina/CliqueSNV#parameters>`_ option of CliqueSNV.

Workflow Options
----------------

These options allow you to skip entire sections of the pipeline. They can
significantly speed up pipeline execution if you know they are not needed, e.g.
skipping read trimming from reads already trimmed in CLC Genomic Workbench. All
of these options are boolean flags that are disabled by default.

--skip_filtering
^^^^^^^^^^^^^^^^

Treats the input reads as trimmed reads and skips using Trimmomatic or Filtlong
on the reads.

--skip_assembly
^^^^^^^^^^^^^^^

Bypasses *de novo* assembly entirely. Note that phylogenetic trees may not be
generated if no assemblies are input.

--skip_haplotype
^^^^^^^^^^^^^^^^

Basically negate the entire purpose of the pipeline by performing no variant
calling, haplotype calling, or phylogenetic analysis on the samples. Can be
useful for debugging.

Resource Allocation Options
---------------------------

These are the maximum resources allowed for a single process within the
pipeline. Place these in a ``nextflow.config`` in a central location on your HPC
to ensure that your pipeline is not cancelled for requesting too many resources.

================== =======
Parameter          Default
================== =======
``--max_memory``   750.GB
``--max_cpus``     72
``--max_time``     240.h
================== =======
