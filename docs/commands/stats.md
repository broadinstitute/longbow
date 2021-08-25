---
layout: default
title: version
description: "Calculate stats on input file."
nav_order: 9
parent: Commands
---

# Stats

## Description

Stats produces several output files that are useful for quality control.  Examples of files produced are:

- summary statistics
- ligation heatmap plots

All figures are created as both `.png` and `.svg` files.

## Command help

```shell
$ longbow stats --help
Usage: longbow stats [OPTIONS] INPUT_BAM

  Calculate and produce stats on the given input bam file.

Options:
  -v, --verbosity LVL       Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -o, --output-prefix TEXT  prefix to give to output files
  -m, --model TEXT          The model to use for annotation.  If the given
                            value is a pre-configured model name, then that
                            model will be used.  Otherwise, the given value
                            will be treated as a file name and Longbow will
                            attempt to read in the file and create a
                            LibraryModel from it.  Longbow will assume the
                            contents are the configuration of a LibraryModel
                            as per LibraryModel.to_json().  [default: mas15]
  --help                    Show this message and exit.

```

## Example

```shell
$ longbow stats input.bam
[INFO 2021-08-25 15:01:37    stats] Invoked via: longbow stats input.bam
[INFO 2021-08-25 15:01:37    stats] Using The standard MAS-seq 15 array element model.
Progress: 0 read [00:00, ? read/s]
[INFO 2021-08-25 15:01:40    stats] Processing statistics...
[INFO 2021-08-25 15:01:40    stats] Writing summary stats file...
[INFO 2021-08-25 15:01:40    stats] Writing complete ligation matrix...
[INFO 2021-08-25 15:01:50    stats] Writing reduced ligation matrix...
[INFO 2021-08-25 15:01:53    stats] Done. Elapsed time: 16.36s.

```