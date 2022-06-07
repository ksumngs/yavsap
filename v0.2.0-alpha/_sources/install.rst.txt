Installation
============

Prerequisites
-------------

Here's what you'll need to get started. Unless otherwise noted, every program
must be available on your ``PATH`` in bash (other shells don't count).

* Git v2.7.0 or higher (🌐/🍺/🐍)
* Curl v2.41.0 or higher (🌐/🍺/🐍)
* Java Runtime v8-v15 (🌐/🍺/🐍)
* Nextflow v20.10.0 or higher (🌐/🐍)
* One or more of the following container engines
   * Docker v20.10.0 or higher (🌐)
   * Podman v3.0 or higher (🌐)
   * Singularity v3.7 or higher (🌐)
* Conda v3.7 or higher (🌐/🍺)
   * Full Anaconda and Miniconda both work

Most of these programs can be installed without root privileges in your
homedirectory using either `Homebrew <https://brew.sh>`_ or
`Conda <https://docs.conda.io/en/latest/miniconda.html>`_. Sources legend:

* 🌐: Distro/Vendor Repos
* 🍺: Linuxbrew
* 🐍: Bioconda/conda-forge

Getting the Pipeline
--------------------

The recommended way to run it is to pull and run at the same time

.. code-block::

    nextflow run ksumngs/yavsap -latest ...

If you need a specific version, use Nextflow's ``-v`` option with the version
tag.

.. code-block::

    nextflow run ksumngs/yavsap -v v0.1.0-alpha ...

If you *really* want to download the pipeline before use, you can run

.. code-block::

    nextflow pull ksumngs/yavsap

If you want to tweak the code and run it locally, you can clone the repo and run
from the code itself.

.. code-block::

    git clone https://github.com/ksumngs/yavsap.git
    ./yavsap/main.nf ...
