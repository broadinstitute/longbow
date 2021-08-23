---
layout: default
title: scsplit 
description: "Split single-cell reads."
nav_order: 7
parent: Commands
---

# Scsplit

## Description

_*This is an experimental tool and should not be used in most applications.*_

The Single Cell Split (`scsplit`) tool is designed to convert MAS-seq data into
a format that is readily ingested by the `alevin` single-cell analysis tool.  

`scsplit` assumes all data are from the same cell, and takes in a dummy cell
barcode as one of its inputs.  The primary output is two [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) files.
Both of these output files have information from each read and are in the 
same order (read 1 of the input bam will be represented as the first entry
in both output files).

`OUTPUT_BASE_NAME_mates_1.fastq` contains barcode information for each 
input read.  This barcode information consists of the given dummy cell
barcode and the UMI for each read.

`OUTPUT_BASE_NAME_mates_2.fastq` contains the transcript sequences for
each read in the given input bam file.  These transcript sequences are
calculated to be any bases that occur after the UMI and before the 
Poly-A tail for each read.

These two files should be considered "linked" and should not be reordered
before processing with `alevin`.

The third and final output is a text file containing the dummy cell 
barcode given to `scsplit` so that `alevin` can know what the cell barcode
is for the data.

As mentioned above, this tool is experimental and not yet mature enough for
use on real datasets. 

More information on `alevin` can be found [here](https://combine-lab.github.io/alevin-tutorial/).

## Command help

```shell
$ longbow scsplit --help
Usage: longbow scsplit [OPTIONS] INPUT_BAM

  Create files for use in `alevin` for single-cell analysis. This tool
  coerces a set of reads from a single source into a format that `alevin`
  can ingest.

  Segment names are assumed to be those in the default model (utils/model.py).

  INPUT_BAM should contain reads that have been processed by `longbow segment`.

  The output from this tool consists of several files:
    OUTPUT_BASE_NAME_mates_1.fastq:  A file containing partial sequences for all
                                     reads in the given input file. These partial
                                     reads consist of the dummy cell barcode + detected
                                     UMI for each read in the given input file.
    OUTPUT_BASE_NAME_mates_2.fastq:  A file containing partial sequences for all
                                     reads in the given input file. These partial
                                     reads consist of the transcript sequences for
                                     all reads in the given input file. Transcript
                                     sequences include data after the UMI and before
                                     the Poly-A tail. All bases outside of this range
                                     are excluded from the output.
    OUTPUT_BASE_NAME_whitelist.txt:  A whitelist file for alevin containing the given
                                     dummy cell barcode.

Options:
  -v, --verbosity LVL          Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -t, --threads INTEGER        number of threads to use (0 for all)  [default: 11]

  -o, --output-base-name TEXT  base name for output files [default: scsplit]
  -c, --cell-barcode TEXT      dummy cell barcode to use for the dataset
                               [default: CTGCCTAACCTGATCC, length: 16]

  -u, --umi-length INTEGER     length of the UMI from this library  [default: 10]

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
```

## Example

```shell
$ longbow scsplit -o single_cell_split segmented.bam
[INFO 2021-08-20 18:14:49  scsplit] Invoked via: longbow scsplit -o single_cell_split tmp.segment.bam
[INFO 2021-08-20 18:14:49  scsplit] Running with 7 worker subprocess(es)
[INFO 2021-08-20 18:14:49  scsplit] Using The standard MAS-seq 15 array element model.
[INFO 2021-08-20 18:14:52  scsplit] Processed 108 reads.
[INFO 2021-08-20 18:14:52  scsplit] CBC length: 16.
[INFO 2021-08-20 18:14:52  scsplit] UMI length: 10.
[INFO 2021-08-20 18:14:52  scsplit] Done. Elapsed time: 3.01s.
```
