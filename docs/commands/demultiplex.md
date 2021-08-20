---
layout: default
title: demultiplex
description: "Demultiplex reads from different array structures."
nav_order: 6
parent: Commands
---

# Demultiplex

## Description

If multiple MAS-seq libraries with different array designs are pooled into a single sequencing run, Longbow can demultiplex the reads into separate BAM files.  This is accomplished with the `demultiplex` command, which tests multiple MAS-seq models per read (see `annotate`) and chooses the one with the highest overall likelihood.

Note that in Longbow parlance, the act of splitting a MAS-seq read array into its constituent elements is referred to as "segmentation", not "demultiplexing". If you're trying to split the reads into individual transcripts, see the `segment` command.

## Command help

```shell
$ longbow demultiplex --help
Usage: longbow demultiplex [OPTIONS] INPUT_BAM

  Separate reads into files based on which model they fit best.

  Resulting reads will be annotated with the model they best fit as well as
  the score and segments for that model.

Options:
  -v, --verbosity LVL       Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -p, --pbi PATH            BAM .pbi index file
  -t, --threads INTEGER     number of threads to use (0 for all)  [default:
                            11]

  -o, --out-base-name TEXT  base name for output files  [default:
                            longbow_demultiplexed]

  --help                    Show this message and exit.
```

## Example

```shell
$ longbow demultiplex -t4 -v INFO multiplexed_data.bam -o MAS-seq_demuxed
[INFO 2021-08-20 17:59:16 demultiplex] Invoked via: longbow demultiplex -t4 -v INFO multiplexed_data.bam -o MAS-seq_demuxed
[INFO 2021-08-20 17:59:16 demultiplex] Running with 4 worker subprocess(es)
[INFO 2021-08-20 17:59:16 demultiplex] Annotating 996 reads
Progress: 100%|████████████████████████████| 996/996 [10:24<00:00,  1.59 read/s]
[INFO 2021-08-20 18:09:45 demultiplex] Annotated 996 reads with 46004 total sections.
[INFO 2021-08-20 18:09:45 demultiplex] Model mas15 annotated 818 reads.
[INFO 2021-08-20 18:09:45 demultiplex] Model mas10 annotated 178 reads.
[INFO 2021-08-20 18:09:45 demultiplex] Done. Elapsed time: 629.44s. Overall processing rate: 1.58 reads/s.
```
