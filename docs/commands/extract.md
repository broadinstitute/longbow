---
layout: default
title: extract
description: "Extract reads."
nav_order: 4
parent: Commands
---

# Extract

## Description

Extract coding segments from the reads in the given bam. The main coding
segments are assumed to be labeled as `random` segments. Uses known
segments flanking the region to be extracted as markers to indicate the
start and end of what to extract.

## Command help

```shell
$ longbow extract --help
Usage: longbow extract [OPTIONS] INPUT_BAM

  Extract coding segments from the reads in the given bam. The main coding
  segments are assumed to be labeled as `random` segments. Uses known
  segments flanking the region to be extracted as markers to indicate the
  start and end of what to extract.

Options:
  -v, --verbosity LVL         Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -p, --pbi PATH              BAM .pbi index file
  -o, --out-file TEXT         Output file name.  [required]
  --force                     Force overwrite of the output files if they
                              exist.  [default: False]

  -b, --base-padding INTEGER  Number of bases to include on either side of the
                              extracted region(s).  [default: 2]

  --leading-adapter TEXT      Adapter preceding the region to extract.
                              [default: 10x_Adapter]

  --trailing-adapter TEXT     Adapter following the region to extract.
                              [default: Poly_A]

  --start-offset INTEGER      Number of bases to ignore from the extracted
                              region start.  These bases will not be included
                              in the extracted sequences.  [default: 26]

  --help                      Show this message and exit.
```

## Example

```shell
$ longbow extract -o extracted.bam filtered.bam
[INFO 2021-08-20 17:53:02  extract] Invoked via: longbow extract -o extracted.bam filtered.bam
[INFO 2021-08-20 17:53:02  extract] Writing extracted read segments to: extracted.bam
[INFO 2021-08-20 17:53:02  extract] Extracting `random` segments between 10x_Adapter and Poly_A.
[INFO 2021-08-20 17:53:02  extract] Ignoring the first 26 bases from extracted read segments.
[INFO 2021-08-20 17:53:02  extract] Including 2 flanking bases.
Progress: 8 read [00:00, 465.82 read/s]
[INFO 2021-08-20 17:53:02  extract] Done. Elapsed time: 0.02s.
[INFO 2021-08-20 17:53:02  extract] Total # Reads Processed: 8
[INFO 2021-08-20 17:53:02  extract] # Reads Containing Extracted Segments: 8 (100.00%)
[INFO 2021-08-20 17:53:02  extract] Total # Segments Extracted: 113
[INFO 2021-08-20 17:53:02  extract] Total # Segments Skipped: 0
[INFO 2021-08-20 17:53:02  extract] # Segments extracted per read: 14.12
```
