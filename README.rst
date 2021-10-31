Longbow
"""""""

.. |pypi_v| image:: https://img.shields.io/pypi/v/maslongbow
.. _pypi_v: https://pypi.org/project/maslongbow/
.. |pypi_dm| image:: https://img.shields.io/pypi/dm/maslongbow
.. _pypi_dm: https://pypi.org/project/maslongbow/

Longbow is a command line tool to process MAS-ISO-seq data. Longbow employs a generative modelling approach to accurately annotate and segment MAS-ISO-seq's concatenated full-length transcript isoforms from single-cell or bulk long read RNA sequencing libraries.

Documentation for all ``longbow`` commands can be found on the `Longbow documentation page <https://broadinstitute.github.io/longbow/>`_.

Installation
------------

``pip`` is recommended for Longbow installation.

::

   pip install maslongbow

To install from github source for development, the following commands can be run.

::

   git clone https://github.com/broadinstitute/longbow.git
   pip install -e longbow/

Getting Started
---------------

The commands below illustrate the Longbow workflow on a small library of SIRVs (Spike-in RNA Variant Control Mixes). MAS-ISO-seq concatenated transcripts are annotated, segmented, and filtered using the `mas15` model.  A number of statistics and QC images are generated along the way.  Final filtered transcripts can then be aligned using standard splice-aware long read mappers (e.g. minimap2). More detail for each command can be found in the `full documentation <https://broadinstitute.github.io/longbow/commands.html>`_.

::

    # Download a tiny test dataset (less than 300K)
    wget https://github.com/broadinstitute/longbow/raw/main/tests/test_data/mas15_test_input.bam
    wget https://github.com/broadinstitute/longbow/raw/main/tests/test_data/mas15_test_input.bam.pbi
    wget https://github.com/broadinstitute/longbow/raw/main/tests/test_data/resources/SIRV_Library.fasta

    # Annotate reads according to the mas15 model
    longbow annotate -m mas15 -o ann.bam mas15_test_input.bam

    # Generate annotation stats and automatic QC figures in .png and .svg format
    longbow stats -o stats ann.bam

    # Segment reads according to the model annotations
    longbow segment -o seg.bam ann.bam

    # Filter out improperly-constructed arrays (will generate filter_passed.bam and filter_failed.bam files)
    longbow filter -o filter seg.bam

    # Align reads with long read aligner (e.g. minimap2, pbmm2)
    samtools fastq filter_passed.bam | \
        minimap2 -ayYL --MD -x splice:hq SIRV_Library.fasta - | \
        samtools sort > align.bam &&
        samtools index align.bam


Getting help
------------

The `Longbow documentation page <https://broadinstitute.github.io/longbow/>`_ provides detailed descriptions of command line options and algorithmic details. If you encounter bugs or have questions/comments/concerns, please file an issue on our `Github page <https://github.com/broadinstitute/longbow/issues>`_.

Developers' guide
-----------------

For information on contributing to Longbow development, visit our `developer documentation <DEVELOP.md>`_.