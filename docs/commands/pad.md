---
layout: default
title: pad
description: "Pad tag by specified number of adjacent bases from the read."
nav_order: 9
parent: Commands
---

# Convert

## Description

Pad tag by specified number of adjacent bases from the read.

## Command help

```shell
$ longbow pad --help
Usage: longbow pad [OPTIONS] INPUT_BAM

  Pad tag by specified number of bases from the read.

Options:
  -v, --verbosity LVL     Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -t, --threads INTEGER   number of threads to use (0 for all)  [default: 7]
  -o, --output-bam PATH   annotated bam output  [default: stdout]
  -m, --model TEXT        The model to use for annotation.  If not specified,
                          it will be autodetected from the BAM header.  If the
                          given value is a pre-configured model name, then
                          that model will be used.  Otherwise, the given value
                          will be treated as a file name and Longbow will
                          attempt to read in the file and create a
                          LibraryModel from it.  Longbow will assume the
                          contents are the configuration of a LibraryModel as
                          per LibraryModel.to_json().
  -f, --force             Force overwrite of the output files if they exist.
                          [default: False]
  -b, --barcode-tag TEXT  The barcode tag to adjust
  -e, --expand INTEGER    Expand tag by specified number of bases
  --help                  Show this message and exit.


```