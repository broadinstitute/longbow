---
layout: default
title: sift
description: "Sift reads."
nav_order: 13
parent: Commands
---

# Sift

## Description

Filter segmented reads by conformation to expected cDNA design.

## Command help

```shell
$ longbow sift --help
Usage: longbow sift [OPTIONS] INPUT_BAM

  Filter segmented reads by conformation to expected cDNA design.

Options:
  -v, --verbosity LVL       Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -p, --pbi PATH            BAM .pbi index file
  -o, --output-bam PATH     filtered bam output (passing reads only)
                            [default: stdout]
  -x, --reject-bam PATH     Filtered bam output (failing reads only)
                            [default: /dev/null]
  -m, --model TEXT          The model to use for annotation.  If not
                            specified, it will be autodetected from the BAM
                            header.  If the given value is a pre-configured
                            model name, then that model will be used.
                            Otherwise, the given value will be treated as a
                            file name and Longbow will attempt to read in the
                            file and create a LibraryModel from it.  Longbow
                            will assume the contents are the configuration of
                            a LibraryModel as per LibraryModel.to_json().
  -f, --force               Force overwrite of the output files if they exist.
  -s, --stats PATH          Table describing the ways in which the reads do
                            not conform to expectation (failing reads
                            only)[default: /dev/null]
  -u, --summary-stats PATH  Table containing summary statistics for the sifted
                            reads that are output.[default: /dev/null]
  -k, --ignore-list PATH    Txt file containing a list of read names to
                            ignore.
  --help                    Show this message and exit.
```

## Example

```shell
$ longbow sift -f -v DEBUG -m mas_15+sc_10x5p -s sift_stats.tsv -u sift_summary.tsv -o sift_out.bam -x sift_rejects.bam tests/test_data/segment/mas_15+sc_10x5p.expected.bam
[INFO 2022-11-03 10:19:14     sift] Invoked via: longbow sift -f -v DEBUG -m mas_15+sc_10x5p -s sift_stats.tsv -u sift_summary.tsv -o sift_out.bam -x sift_rejects.bam tests/test_data/segment/mas_15+sc_10x5p.expected.bam
[E::idx_find_and_load] Could not retrieve index file for 'tests/test_data/segment/mas_15+sc_10x5p.expected.bam'
[INFO 2022-11-03 10:19:15     sift] Using mas_15+sc_10x5p: 15-element MAS-ISO-seq array, single-cell 10x 5' kit
[WARNING 2022-11-03 10:19:15 bam_utils] Output file exists: sift_out.bam.  Overwriting.
[WARNING 2022-11-03 10:19:15 bam_utils] Output file exists: sift_rejects.bam.  Overwriting.
[INFO 2022-11-03 10:19:15     sift] Writing reads that conform to the model to: sift_out.bam
[INFO 2022-11-03 10:19:15     sift] Writing reads that do not conform to the model to: sift_rejects.bam
mas_15+sc_10x5p
Progress: 0 read [00:00, ? read/s][DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001001/ccs
[INFO 2022-11-03 10:19:15     sift] Sifted 0 reads.
[DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001002/ccs
[DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001003/ccs
[DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001004/ccs
[DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001005/ccs
[DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001006/ccs
[DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001007/ccs
[DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001008/ccs
[DEBUG 2022-11-03 10:19:15     sift] Read is mas_15+sc_10x5p valid: m64020_201213_022403/1000000001009/ccs
...
Progress: 112 read [00:00, 4092.25 read/s]
[INFO 2022-11-03 10:19:15     sift] Total Reads Processed:	112
[INFO 2022-11-03 10:19:15     sift] # Reads Passing Model Filter:	112	100.0000
[INFO 2022-11-03 10:19:15     sift] # Reads Failing Model Filter:	0	0.0000
[INFO 2022-11-03 10:19:15     sift] # Reads Ignored:	0	0.0000
[INFO 2022-11-03 10:19:15     sift] Done. Elapsed time: 0.43s.
```
