---
layout: default
title: convert
description: "Convert reads from fastq{,.gz} files for use with `annotate`."
nav_order: 2
parent: Commands
---

# Convert

## Description

Convert reads from fastq{,.gz} files for use with `annotate`.

## Command help

```shell
$ longbow convert --help
Usage: longbow convert [OPTIONS] INPUT_SPEC

  Convert reads from fastq{,.gz} files for use with `annotate`.

Options:
  -v, --verbosity LVL       Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -o, --output-bam PATH     bam output  [default: stdout]
  -q, --default-rq FLOAT    Default read quality for ingested reads (ranging
                            from [-1,1])  [default: -1.0]
  -g, --read-group-id TEXT  Read group ID to set
  -s, --sample-name TEXT    Sample name to set in the read group
  -l, --min-length INTEGER  Minimum length of reads to process
  -L, --max-length INTEGER  Maximum length of reads to process
  -f, --force               Force overwrite of the output files if they exist.
                            [default: False]
  --help                    Show this message and exit.
```