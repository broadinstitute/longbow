---
layout: default
title: extract
description: "Extract reads."
nav_order: 4
parent: Commands
---

# Extract

## Description

The `extract` command is a read post-processing tool.  `Extract` 
removes leading and trailing adapter sequences, as well as poly-A 
tails from the given reads so that they can be aligned to the 
genome / transcriptome without library prep artifacts affecting 
the alignment process.  This effectively extracts the cDNA 
portion of each read.

The cDNA portion of each read is expected to be in sections that
Longbow has labeled as `random`.  The leading and trailing 
sequences that `extract` uses as markers are configurable and 
default to the regions that most commonly flank the cDNA regions
in the current models (`10x_adapter` and `Poly_A`).

`extract` can ignore a number of bases from the start of the 
`random` / cDNA section to account for barcodes prepended to 
reads (by default, 26).  In addition, a configurable number
of bases that extend beyond the `random` segment are included in
the output to mitigate potential off-by-one errors in the Longbow
model (by default, 2).  

Input reads are expected to be either `annotate`d or `segment`ed.

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
