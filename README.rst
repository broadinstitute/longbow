Longbow
"""""""

|GitHub release| |Generic badge| |PyPI version maslongbow|

.. |GitHub release| image:: https://img.shields.io/github/release/broadinstitute/longbow.svg
   :target: https://github.com/broadinstitute/longbow/releases/

.. |Generic badge| image:: https://img.shields.io/badge/Docker-v0.5.19-blue.svg
   :target: https://console.cloud.google.com/gcr/images/broad-dsp-lrma/US/lr-longbow

.. |PyPI version maslongbow| image:: https://img.shields.io/pypi/v/maslongbow.svg
   :target: https://pypi.python.org/pypi/maslongbow/

Longbow is a command line tool to process MAS-ISO-seq data. Longbow employs a generative modelling approach to accurately annotate and segment MAS-ISO-seq's concatenated full-length transcript isoforms from single-cell or bulk long read RNA sequencing libraries.

Documentation for all ``longbow`` commands can be found on the `Longbow documentation page <https://broadinstitute.github.io/longbow/>`_.

Installation
------------

``pip`` is recommended for Longbow installation.

::

   pip install maslongbow

For a pre-built version including all dependencies, access our Docker image.

::

   docker pull us.gcr.io/broad-dsp-lrma/lr-longbow:0.5.19

To install from Github source for development, the following commands can be run.

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

    # Basic processing workflow
    longbow annotate -m mas15v2 mas15_test_input.bam | \  # Annotate reads according to the mas15v2 model
      tee ann.bam | \                                     # Save annotated BAM for later
      longbow filter | \                                  # Filter out improperly-constructed arrays
      longbow segment | \                                 # Segment reads according to the model
      longbow extract -o filter_passed.bam                # Extract adapter-free cDNA sequences

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
