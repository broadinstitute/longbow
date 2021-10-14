---
layout: default
title: annotate
description: "Annotate adapters and cDNA sequences in MAS-seq reads."
nav_order: 1
parent: Commands
---

# Annotate

## Description

The `annotate` command is the first in the chain of commands for processing MAS-seq data. Given an array model (e.g. 'mas15'), `annotate` determines an optimal annotation (Viterbi path) for the adapters and unknown sequences in each read. The Viterbi path is automatically computed for both the forward and reverse-complement orientations of the read, and the path with the highest overall maximum likelihood is considered to be the correct one.  Longbow then appends this information to the read's auxillary BAM tags for use with downstream Longbow tools (e.g. `filter`, `segment`).

Several models are available for use with the `--model` argument.  They are:

| name          | description                                 | version        |
|---------------|---------------------------------------------|----------------|
| mas15         | The standard MAS-seq 15-element array model | 1.0.0          |
| mas10         | The MAS-seq 10-element array model          | 1.0.0          |
| mas8prototype | The prototype MAS-seq 8-element array       | 1.0.0          |
| slide-seq     | The Slide-seq 15-element array model        | 0.0.1          |
| bulk15        | The 15-element bulk RNA model               | (experimental) |

By default, the `annotate` command streams its results to `stdout` for use with other Longbow commands or tools (e.g. `samtools`) via Unix pipes.

This command is parallelizable on a per-read level.  On a system with N cores, it will attempt to use N-1 cores by default.

Optionally, this command can make use of a PacBio index (.pbi) file, which specifies the total number of reads in the file and can therefore be useful for accurate progress reporting.  The .pbi file has no effect on speed or accuracy of results; it is solely for convenience in progress logging.

## Command help

```shell
$ longbow annotate --help
Usage: longbow annotate [OPTIONS] INPUT_BAM

  Annotate reads in a BAM file with segments from the model.

Options:
  -v, --verbosity LVL    Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -p, --pbi PATH         BAM .pbi index file
  -t, --threads INTEGER  number of threads to use (0 for all)  [default: 11]
  -o, --output-bam PATH  annotated bam output  [default: stdout]
  -m, --model TEXT       The model to use for annotation.  If the given value
                         is a pre-configured model name, then that model will
                         be used.  Otherwise, the given value will be treated
                         as a file name and Longbow will attempt to read in
                         the file and create a LibraryModel from it.  Longbow
                         will assume the contents are the configuration of a
                         LibraryModel as per LibraryModel.to_json().
                         [default: mas15]
                         
  -c, --chunk TEXT       Process a single chunk of data (e.g. specify '2/4' to
                         process the second of four equally-sized chunks
                         across the dataset)

  --help                 Show this message and exit.
```

## Example

```shell
$ longbow annotate -o annotated.bam tests/test_data/mas15_test_input.bam
[INFO 2021-08-12 09:55:07 annotate] Invoked via: longbow annotate -o annotated.bam tests/test_data/mas15_test_input.bam
[INFO 2021-08-12 09:55:07 annotate] Running with 11 worker subprocess(es)
[INFO 2021-08-12 09:55:07 annotate] Using The standard MAS-seq 15 array element model.
[INFO 2021-08-12 09:55:10 annotate] Annotating 8 reads
Progress: 100%|████████████████████████████████████████████████████████████████████████████████| 8/8 [00:02<00:00,  3.26 read/s]
[INFO 2021-08-12 09:55:12 annotate] Annotated 8 reads with 557 total sections.
[INFO 2021-08-12 09:55:12 annotate] Done. Elapsed time: 5.35s. Overall processing rate: 1.50 reads/s.
```