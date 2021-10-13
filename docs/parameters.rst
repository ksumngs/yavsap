Complete Parameter Reference
============================

Every single parameter in the entire pipeline is described here. If you need to
tweak every little detail of how your JEV samples are analyzed, you've come to
the right place! Note that *every* parameter is described here, even those that
shouldn't be used, so proceed with caution!

Input Options
-------------

``--input``
^^^^^^^^^^^

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

``--sra``
^^^^^^^^^

======== ======
Type     Flag
======== ======
Required No
Default  false
======== ======

This flag switches the meaning of ``--input`` to mean an NCBI Short Read Archive
(SRA) accession number, then pulls the associated files directly from NCBI and
runs the pipeline on them. To use this flag, you **must** have an `NCBI API key
<https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/>`_
and it **must** be exported in the shell environment as ``NCBI_API_KEY``, e.g.
``NCBI_API_KEY=0123456789abcdef``. This feature is currenly broken due to an
upstream issue in Nextflow, but should be fixed by Nextflow version 21.10.0.

``--platform``
^^^^^^^^^^^^^^

======== ======
Type     String
======== ======
Required Yes
Default  n/a
======== ======

Defines which type of reads to analyze. Must be either 'illumina' or 'nanopore'.


``--pe``
^^^^^^^^

======== ======
Type     Boolean
======== ======
Required No
Default  false
======== ======

.. note:: This parameter is hidden and should **never** be used on the command
    line! Bad things will happen if you do!

Internally-used parameter to check if the pipeline is processing Illumina reads.
Automatically set by ``--platform``

``--ont``
^^^^^^^^^

======== ======
Type     Boolean
======== ======
Required No
Default  false
======== ======

.. note:: This parameter is hidden and should **never** be used on the command
    line! Bad things will happen if you do!

Internally-used parameter to check if the pipeline is processing Nanopore reads.
Automatically set by ``--platform``

Reference Genome Options
------------------------

``--genome``
^^^^^^^^^^^^

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

Kraken2 Options
---------------

``--kraken2_db``
^^^^^^^^^^^^^^^^

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

Corresponds to the |--db option of Kraken2|_.

.. |--db option of Kraken2| replace:: ``--db`` option of Kraken2
.. _--db option of Kraken2: https://github.com/DerrickWood/kraken2/wiki/Manual#classification

``--keep_taxid``
^^^^^^^^^^^^^^^^
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

``--trim_minlen``
"""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  100/300 (Illumina/Nanopore)
======== ======

Remove reads that are shorter than this length in bases.

Corresponds to the |MINLEN option of Trimmomatic|_ for Illumina reads.

Corresponds to the |--min_length option of Filtlong|_ for Nanopore reads.

.. |MINLEN option of Trimmomatic| replace:: ``MINLEN:`` option of Trimmomatic
.. _MINLEN option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. |--min_length option of Filtlong| replace:: ``--min_length`` option of Filtlong
.. _--min_length option of Filtlong: https://github.com/rrwick/Filtlong#full-usage

``--trim_winsize``
""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  50/250 (Illumina/Nanopore)
======== ======

The number of bases to average quality accross during sliding window trimming.

Corresponds to the |first SLIDINGWINDOW option of Trimmomatic|_ for Illumina reads.

Corresponds to the |--window_size option of Filtlong|_ for Nanopore reads.

.. |first SLIDINGWINDOW option of Trimmomatic| replace:: ``first SLIDINGWINDOW`` option of Trimmomatic
.. _first SLIDINGWINDOW option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. |--window_size option of Filtlong| replace:: ``--window_size`` option of Filtlong
.. _--window_size option of Filtlong: https://github.com/rrwick/Filtlong#full-usage

``--trim_winqual``
""""""""""""""""""

======== ======
Type     Integer/Float (Illumina/Nanopore)
======== ======
Required No
Default  15/0.9 (Illumina/Nanopore)
======== ======

The minimum average quality within the sliding window to keep a read. Note that
this value is the minmum PHRED score when trimming Illumina reads, but it is a
percentage score when trimming Nanopore reads.

Corresponds to the |second SLIDINGWINDOW option of Trimmomatic|_ for Illumina reads.

Corresponds to the |--min_mean_q option of Filtlong|_ for Nanopore reads.

.. |second SLIDINGWINDOW option of Trimmomatic| replace:: ``second SLIDINGWINDOW`` option of Trimmomatic
.. _second SLIDINGWINDOW option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
.. |--min_mean_q option of Filtlong| replace:: ``--min_mean_q`` option of Filtlong
.. _--min_mean_q option of Filtlong: https://github.com/rrwick/Filtlong#full-usage


Illumina-Specific (Trimmomatic) Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--trim_adapers``
""""""""""""""""""

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

Corresponds to the |first ILLUMINACLIP option of Trimmomatic|_.

.. |first ILLUMINACLIP option of Trimmomatic| replace:: first ``ILLUMINACLIP`` option of Trimmomatic
.. _first ILLUMINACLIP option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

``--trim_mismatches``
"""""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  2
======== ======

The maximum mismatch count which will still allow a full adapter match to be
performed.

Corresponds to the |second ILLUMINACLIP option of Trimmomatic|_.

.. |second ILLUMINACLIP option of Trimmomatic| replace:: second ``ILLUMINACLIP`` option of Trimmomatic
.. _second ILLUMINACLIP option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

``--trim_pclip``
""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  30
======== ======

``pclip``: palindrome clip. How accurate the match between the two adapter
ligated reads must be for paired-end palindrome read alignment.

Corresponds to the |third ILLUMINACLIP option of Trimmomatic|_.

.. |third ILLUMINACLIP option of Trimmomatic| replace:: third ``ILLUMINACLIP`` option of Trimmomatic
.. _third ILLUMINACLIP option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

``--trim_clip``
"""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  10
======== ======

How accurate the match between any adapter sequence must be against a read.

Corresponds to the |final ILLUMINACLIP option of Trimmomatic|_.

.. |final ILLUMINACLIP option of Trimmomatic| replace:: final ``ILLUMINACLIP`` option of Trimmomatic
.. _final ILLUMINACLIP option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

``--trim_leading``
""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  15
======== ======

The minimum quality to keep a base in the leading end of a read. If set to
``0``, LEADING trimming is disabled.

Corresponds to the |LEADING option of Trimmomatic|_.

.. |LEADING option of Trimmomatic| replace:: ``LEADING:`` option of Trimmomatic
.. _LEADING option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

``--trim_trailing``
"""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  15
======== ======

The minimum quality to keep a base in the trailing end of a read. If set to
``0``, TRAILING trimming is disabled.

Corresponds to the |TRAILING option of Trimmomatic|_.

.. |TRAILING option of Trimmomatic| replace:: ``TRAILING:`` option of Trimmomatic
.. _TRAILING option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

``--trim_crop``
"""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  0
======== ======

The number of bases to keep from the start of the read. If set to ``0``, CROP
trimming is disabled.

Corresponds to the |CROP option of Trimmomatic|_.

.. |CROP option of Trimmomatic| replace:: ``CROP:`` option of Trimmomatic
.. _CROP option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

``--trim_headcrop``
"""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  0
======== ======

The number of bases to remove from the start of the read. If set to ``0``,
HEADCROP trimming is disabled.

Corresponds to the |HEADCROP option of Trimmomatic|_.

.. |HEADCROP option of Trimmomatic| replace:: ``HEADCROP:`` option of Trimmomatic
.. _HEADCROP option of Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic

Nanopore-Specific (Filtlong) Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--trim_keep_percent``
"""""""""""""""""""""""

======== ======
Type     Float
======== ======
Required No
Default  0.9
======== ======

The number of reads to keep as a percentage of the total reads, discarding the
lowest-quality reads.

Corresponds to the |--keep_percent option of Filtlong|_.

.. |--keep_percent option of Filtlong| replace:: ``--keep_percent`` option of Filtlong
.. _--keep_percent option of Filtlong: https://github.com/rrwick/Filtlong#full-usage

``--trim_target_bases``
"""""""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  500,000,000
======== ======

The maximum number of bases to keep, discarding the lowest quality reads until
this threshold is met. Can be effectively disabled by using a very high number
of bases.

Corresponds to the |--target_bases option of Filtlong|_.

.. |--target_bases option of Filtlong| replace:: ``--target_bases`` option of Filtlong
.. _--target_bases option of Filtlong: https://github.com/rrwick/Filtlong#full-usage

*de novo* Assembly Options
--------------------------

Illumina-Specific (SPAdes) Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``--spades_mode``
"""""""""""""""""

======== ======
Type     String
======== ======
Required No
Default  rnaviral
======== ======

If provided, this parameter is turned into a flag and passed as the 'mode' to
the SPAdes assembly, e.g.
``nextflow run jev-analysis-pipeline --spades_mode 'metaviral'`` will run SPAdes
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

``--canu_corrected_error_rate``
"""""""""""""""""""""""""""""""

======== ======
Type     Float
======== ======
Required No
Default  0.144
======== ======

How dissimilar overlap between two reads can be and still be assembled together.

Corresponds to the |correctedErrorRate option of Canu|_.

.. |correctedErrorRate option of Canu| replace:: ``correctedErrorRate`` option of Canu
.. _correctedErrorRate option of Canu: https://canu.readthedocs.io/en/latest/parameter-reference.html#correctederrorrate

``--canu_min_read_length``
"""""""""""""""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  1000
======== ======

The shortest length allowed for a read to be considered in the assembly.

Corresponds to the |minReadLength option of Canu|_.

.. |minReadLength option of Canu| replace:: ``minReadLength`` option of Canu
.. _minReadLength option of Canu: https://canu.readthedocs.io/en/latest/parameter-reference.html#minreadlength

``--canu_min_overlap_length``
"""""""""""""""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  ``--canu_min_read_length`` ÷ 2
======== ======

The shortest length allowed for reads to overlap to be considered in the assembly.

Corresponds to the |minOverlapLength option of Canu|_.

.. |minOverlapLength option of Canu| replace:: ``minOverlapLength`` option of Canu
.. _minOverlapLength option of Canu: https://canu.readthedocs.io/en/latest/parameter-reference.html#minoverlaplength

``--canu_stop_on_low_coverage``
"""""""""""""""""""""""""""""""

======== ======
Type     Integer
======== ======
Required No
Default  10
======== ======

Lowest depth allowable for assembly to proceed.

Corresponds to the |stopOnLowCoverage option of Canu|_.

.. |stopOnLowCoverage option of Canu| replace:: ``stopOnLowCoverage`` option of Canu
.. _stopOnLowCoverage option of Canu: https://canu.readthedocs.io/en/latest/parameter-reference.html#stoponlowcoverage

Haplotyping Options
-------------------

``--haplotype_significance``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

======== ======
Type     Float
======== ======
Required No
Default  0.05
======== ======

The highest p-value that will be considered a significant haplotype based on
linkage disequilibrium and proportional equivalence.

``--haplotype_minimum``
^^^^^^^^^^^^^^^^^^^^^^^

======== ======
Type     Integer
======== ======
Required No
Default  10
======== ======

The minimum number of times a particular haplotype has to occur for it to be
considered real and processed downstream and output.
