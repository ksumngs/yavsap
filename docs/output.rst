Output
======

You ran the pipeline, but now what? Here is the description of the pipeline
outputs and where to find them, but first...

The Visualizer
--------------

If you like pictures, hate juggling files into various Java applications, and
are OCD about having everything in one spot, then we have you covered. The
pipeline outputs its own web application with every run that allows you to view,
interact with, and download nearly every output file.

To use, you'll need Node.js v14 (LTS) or later. I highly recommend installing
Node.js via `n <https://github.com/tj/n>`_, which keeps your Node.js packages
clean and organized. You can use n, Node.js, and the Visualizer on any computer,
including Mac and Windows computers. Once you have Node.js installed, point a
terminal to the output directory, and run

::

    npm install && npm start

A link should appear in the terminal. Follow the link in your browser to access
the Visualizer.

The Visualizer contains links to the following:

* Latest MultiQC report
* Nextflow reports
* Nextflow timelines
* Nextflow traces
* Nextflow process graph
* For each sample
   * Interactive alignment viewer
      * Reads alignment
      * Assembly alignment
   * Downloadable alignment in BAM format
   * Interactive phylogenetic tree

Once you are done, press ``CTRL+C`` in the terminal, and the Visualizer will
shut down. If you need it again, just run ``npm start`` in the results folder
and follow the link again.

File Output
-----------

You don't like all the fluff, huh? Here's the lowdown on the file structure you
can expect from each pipeline run.

Visualizer Data
^^^^^^^^^^^^^^^

All the following files and folders are strictly for powering the Visualizer,
and have no scientific merit. You may safely ignore them.

::

    results
    ├── 📁 css
    ├── 📁 views
    ├── 📝 index.js
    ├── 📝 package.json
    ├── 📝 package-lock.json
    └── 📝 favicon.ico

MultiQC Report
^^^^^^^^^^^^^^

MultiQC is a tool that provides overall diagnostic insights into metagenomic
pipelines. Although its use is limited in this pipeline, the report does
contain information on trimmed reads from Trimmomatic and filtered reads from
Kraken2. The reports can be viewed via :ref:`The Visualizer`, or they can be
viewed alone. The report overwrites itself with every pipeline run, and is named
``multiqc_report.html``. Its sidecar data files are located in the
``multiqc_data`` folder.

::

    results
    ├── 📁 multiqc_data
    └── 📝 multiqc_report.html

Reference Genome
^^^^^^^^^^^^^^^^

The pulled reference genome in fasta format can be found in the ``data`` folder,
as well as its samtools index. The name of the files are given by the NCBI
description of the reference genome, with some replacements to make it a safe
name for all the processes that use it.

::

    results
    └── 📁 data
        ├── 📝 Japanese_Encephalitis_Virus.fasta
        └── 📝 Japanese_Encephalitis_Virus.fasta.fai

Alignments
^^^^^^^^^^

Alignments for each sample can be found in the ``data`` folder in the BAM
format, as well as their sidecar index files. Files are always in the
``<samplename>.bam`` and ``<samplename>.bam.bai`` naming scheme.

::

    results
    └── 📁 data
        ├── 📝 pig-serum.bam
        ├── 📝 pig-serum.bam.bai
        ├── 📝 pig-feces.bam
        └── 📝 pig-feces.bam.bai

Variant Calls
^^^^^^^^^^^^^

Variant calls for each sample are output directly to the results folder when
analyzing Nanopore reads. The filename is always ``<samplename>.filtered.tsv``.

::

    results
    ├── 📝 pig-serum.filtered.tsv
    └── 📝 pig-feces.filtered.tsv

Haplotypes
^^^^^^^^^^

Haplotypes for each sample are output directly to the results folder. There is a
a fasta file containing the mutated sequences, and also a data file describing
the haplotypes. The data file is in JSON format for Illumina reads and YAML
format for Nanopore reads. The filename is always
``<samplename>.haplotypes.<ext>``.

::

    results
    ├── 📝 pig-serum.haplotypes.fasta
    ├── 📝 pig-serum.haplotypes.yaml
    ├── 📝 pig-feces.haplotypes.fasta
    └── 📝 pig-feces.haplotypes.yaml


Multiple Alignments
^^^^^^^^^^^^^^^^^^^

The alignments of the reference genome, assembly, and haplotypes in fasta format
are contained in the ``<samplename>.haplotypes.fas`` file in the results folder.

::

    results
    ├── 📝 pig-serum.haplotypes.fas
    └── 📝 pig-feces.haplotypes.fas

Phylogenetic Trees
^^^^^^^^^^^^^^^^^^

Phylogenetic trees of the haplotypes in Newick format are contained in the
``data`` folder in files with the name ``<samplename>.tree``.

::

    results
    └── 📁 data
        ├── 📝 pig-serum.tree
        └── 📝 pig-feces.tree
