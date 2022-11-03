---
layout: default
title: extract
description: "Extract reads."
nav_order: 14
parent: Commands
---

# Tagfix

## Description

Update longbow read tags after alignment.

According to the SAM spec, aligners reverse-complement the read sequence in the event that a read maps to the reverse strand.
However, the aligners do not know about the tags that `Longbow` adds to the reads and therefore certain tags will be incorrect
for reverse-stranded reads after the alignment process (e.g. the `SG` tag).

`Tagfix` checks each read in the input bam file for its strandedness and in the event of a reverse strand read, the position-based
tags that have been added to the read by `Longbow` are corrected to reflect the sequence of bases in the BAM file itself.

## Command help

```shell
$ longbow tagfix --help
Usage: longbow tagfix [OPTIONS] INPUT_BAM

  Update longbow read tags after alignment.

Options:
  -v, --verbosity LVL    Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -t, --threads INTEGER  number of threads to use (0 for all)  [default: 7]
  -o, --output-bam PATH  annotated bam output  [default: stdout]
  -f, --force            Force overwrite of the output files if they exist.
                         [default: False]
  --help                 Show this message and exit.
```
