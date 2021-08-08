---
layout: default
title: extract
description: "Extract reads."
nav_order: 4
parent: Commands
---

# Extract

`
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
`