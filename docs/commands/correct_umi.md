---
layout: default
title: correct_umi 
description: "Get info on built-in models."
nav_order: 4
parent: Commands
---

# Correct_UMI

## Description

Correct UMIs with Set Cover algorithm.

Corrects all UMIs in the given bam file.  

Algorithm originally developed by Victoric Popic.

### Data Requirements:

- Bam file should be aligned and annotated with genes and transcript equivalence classes prior to running.
- It is critical that you give the proper input for the `--pre-extracted` flag.
  - If the file has been run through `longbow extract`, use `--pre-extracted`
  - If the file has _NOT_ been run through `longbow extract` _DO NOT USE_ `--pre-extracted`

The following tags are required in the input file:

- `CB` Cell Barcode
- `JX` (Adjusted UMI - Configurable)
- `eq` (Equivalence class assignment - Configurable)
- `XG` (Gene assignment)
- `rq` (Read Quality: [-1.0, 1.0])
- `JB` (Back / UMI trailing segment Smith-Waterman alignment score - Configurable)

## Command help

```shell
$ longbow correct_umi --help
Usage: longbow correct_umi [OPTIONS] INPUT_BAM

  Correct UMIs with Set Cover algorithm.

Options:
  -v, --verbosity LVL             Either CRITICAL, ERROR, WARNING, INFO or
                                  DEBUG
  --max-ccs-edit-dist INTEGER     Maximum edit distance between a CCS read and
                                  a target to consider them both part of a
                                  group of UMIs  [default: 2]
  --max-clr-edit-dist INTEGER     Maximum edit distance between a CLR read and
                                  a target to consider them both part of a
                                  group of UMIs  [default: 3]
  --max-ccs-length-diff INTEGER   Maximum length difference between a CCS read
                                  and a target to consider them both part of a
                                  group of UMIs  [default: 50]
  --max-clr-length-diff INTEGER   Maximum length difference between a CLR read
                                  and a target to consider them both part of a
                                  group of UMIs  [default: 150]
  --max-ccs-gc-diff FLOAT         Maximum GC content difference between a CCS
                                  read and a target to consider them both part
                                  of a group of UMIs  [default: 0.05]
  --max-clr-gc-diff FLOAT         Maximum GC content difference between a CLR
                                  read and a target to consider them both part
                                  of a group of UMIs  [default: 0.15]
  --max-ccs-umi-length-delta INTEGER
                                  Maximum length difference between the UMI of
                                  a CCS read and the given UMI length to be
                                  included in the loci for possible
                                  correction.  [default: 3]
  --max-clr-umi-length-delta INTEGER
                                  Maximum length difference between the UMI of
                                  a CLR read and the given UMI length to be
                                  included in the loci for possible
                                  correction.  [default: 4]
  --max-final-ccs-umi-length-delta INTEGER
                                  Maximum length difference between the final,
                                  corrected UMI of a CCS read and the given
                                  UMI length to be included in the final good
                                  results file.  [default: 3]
  --max-final-clr-umi-length-delta INTEGER
                                  Maximum length difference between the final,
                                  corrected UMI of a CLR read and the given
                                  UMI length to be included in the final good
                                  results file.  [default: 3]
  --min-back-seg-score INTEGER    Minimum score of the back alignment tag for
                                  a read to be included in the final good
                                  results file.  [default: 10]
  --umi-tag TEXT                  Tag from which to read in UMIs from the
                                  input bam file.  [default: JX]
  --eq-class-tag TEXT             Tag from which to read in read transcript
                                  equivalence classes (i.e. transcript IDs)
                                  from the input bam file.  [default: eq]
  --final-umi-tag TEXT            Tag into which to put final, corrected UMI
                                  values.  [default: BX]
  --umi-corrected-tag TEXT        Tag into which to put whether a given UMI
                                  was actually corrected.  [default: UX]
  -l, --umi-length INTEGER        Length of the UMI for this sample.
                                  [default: 10]
  -o, --output-bam PATH           Corrected UMI bam output [default: stdout].
  -x, --reject-bam PATH           Filtered bam output (failing reads only).
  -f, --force                     Force overwrite of the output files if they
                                  exist.
  --pre-extracted                 Whether the input file has been processed
                                  with `longbow extract`
  --help                          Show this message and exit.
```


