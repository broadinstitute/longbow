---
layout: default
title: annotate
description: "Annotate adapters and cDNA sequences in MAS-seq reads."
nav_order: 1
parent: Commands
---

# Annotate

## Description

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

  --help                 Show this message and exit.
```

## Example
