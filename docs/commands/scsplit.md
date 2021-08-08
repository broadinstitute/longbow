---
layout: default
title: scsplit 
description: "Split single-cell reads."
nav_order: 7
parent: Commands
---

# Scsplit

`
$ longbow scsplit --help
Usage: longbow scsplit [OPTIONS] INPUT_BAM

  Create files for use in `alevin` for single-cell analysis. This tool
  coerces a set of reads from a single source into a format that `alevin`
  can ingest.

  Segment names are assumed to be those in the default model
  (utils/model.py).

  INPUT_BAM should contain reads that have been processed by `longbow
  segment`.

  The output from this tool consists of several files:
  OUTPUT_BASE_NAME_mates_1.fastq:         A file containing partial
  sequences for all reads in the given input file.  These partial reads
  consist of the          dummy cell barcode + detected UMI for each read in
  the given input file.     OUTPUT_BASE_NAME_mates_2.fastq:         A file
  containing partial sequences for all reads in the given input file.  These
  partial reads consist of the          transcript sequences for all reads
  in the given input file.  Transcript sequences include data after the UMI
  and before the Poly-A tail.  All bases outside of this range are excluded
  from the output.     OUTPUT_BASE_NAME_whitelist.txt:         A whitelist
  file for alevin containing the given dummy cell barcode.

Options:
  -v, --verbosity LVL          Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -t, --threads INTEGER        number of threads to use (0 for all)  [default:
                               11]

  -o, --output-base-name TEXT  base name for output files [default: scsplit]
  -c, --cell-barcode TEXT      dummy cell barcode to use for the dataset
                               [default: CTGCCTAACCTGATCC, length: 16]

  -u, --umi-length INTEGER     length of the UMI from this library  [default:
                               10]

  -b, --write-bam              Write out an annotated bam file in addition to
                               the mates files.  [default: False]

  --force                      Force scsplit to run on the input bam without
                               checking for compatibility.  [default: False]

  -m, --model TEXT             The model to use for annotation.  If the given
                               value is a pre-configured model name, then that
                               model will be used.  Otherwise, the given value
                               will be treated as a file name and Longbow will
                               attempt to read in the file and create a
                               LibraryModel from it.  Longbow will assume the
                               contents are the configuration of a
                               LibraryModel as per LibraryModel.to_json().
                               [default: mas15]

  --help                       Show this message and exit.
`
