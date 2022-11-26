---
layout: default
title: peek
description: "Guess the best pre-built array model to use for annotation."
nav_order: 11
parent: Commands
---

# Convert

## Description

Guess the best pre-built array model to use for annotation.

## Command help

```shell
$ longbow peek --help
Usage: longbow peek [OPTIONS] INPUT_BAM

  Guess the best pre-built array model to use for annotation.

Options:
  -v, --verbosity LVL             Either CRITICAL, ERROR, WARNING, INFO or
                                  DEBUG
  -p, --pbi PATH                  BAM .pbi index file
  -t, --threads INTEGER           number of threads to use (0 for all)
                                  [default: 7]
  -o, --output-model PATH         model name output  [default: stdout]
  -c, --chunk TEXT                Process a single chunk of data (e.g. specify
                                  '2/4' to process the second of four equally-
                                  sized chunks across the dataset)
  -n, --num-reads INTEGER         Number of reads to examine when trying to
                                  guess the array model
  -l, --min-length INTEGER        Minimum length of a read to process.  Reads
                                  shorter than this length will not be
                                  annotated.  [default: 0]
  -L, --max-length INTEGER        Maximum length of a read to process.  Reads
                                  longer than this length will not be
                                  annotated.  [default: 30000]
  --min-rq FLOAT                  Minimum ccs-determined read quality for a
                                  read to be annotated.  CCS read quality
                                  range is [-1,1].  [default: -2.0]
  -f, --force                     Force overwrite of the output files if they
                                  exist.  [default: False]
  -d, --include-deprecated-models
                                  Examine the deprecated built-in models as
                                  well  [default: False]
  --help                          Show this message and exit.
```