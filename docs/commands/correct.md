---
layout: default
title: correct
description: "Correct tag to values provided in barcode allowlist."
nav_order: 3
parent: Commands
---

# Correct

## Description

Correct tag to values provided in barcode allowlist.

This tag correction scheme requires one key input files aside from the reads:
- A barcode allow list (i.e. a list of all possible/allowed barcodes).
  
In addition, a frequency barcode list can be given.  This frequency barcode list should be a TSV file of the form `BARCODE    COUNT`.  It should ideally be obtained by an orthogonal sequencing method (i.e. short reads).  This frequency list will weight barcodes based on the rate at which they are known to appear in the data.  The more often a barcode appears, the more likely that a read's barcode will be corrected to it.  

## Command help

```shell
$ longbow correct --help
Usage: longbow correct [OPTIONS] INPUT_BAM

  Correct tag to values provided in barcode allowlist.

Options:
  -v, --verbosity LVL          Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -p, --pbi PATH               BAM .pbi index file
  -t, --threads INTEGER        number of threads to use (0 for all)  [default:
                               7]
  -o, --output-bam PATH        annotated bam output  [default: stdout]
  -m, --model TEXT             The model to use for annotation.  If not
                               specified, it will be autodetected from the BAM
                               header.  If the given value is a pre-configured
                               model name, then that model will be used.
                               Otherwise, the given value will be treated as a
                               file name and Longbow will attempt to read in
                               the file and create a LibraryModel from it.
                               Longbow will assume the contents are the
                               configuration of a LibraryModel as per
                               LibraryModel.to_json().
  -f, --force                  Force overwrite of the output files if they
                               exist.  [default: False]
  -r, --restrict-to-allowlist  Restrict barcode correction possibilities to
                               only those on the allowlist.  [default: True]
  -b, --barcode-tag TEXT       The tag from which to read the uncorrected
                               barcode.  [default: CR]
  -c, --corrected-tag TEXT     The tag in which to store the corrected
                               barcode.  [default: CB]
  -a, --allow-list PATH        List of allowed barcodes for specified tag
                               (.txt, .txt.gz).  [required]
  --barcode-freqs PATH         TSV file containing barcodes and the
                               frequencies associated with them in the data
                               (BARCODE      FREQ).  If not provided, barcode
                               freqs will be uniformly seeded by the barcode
                               whitelist.
  --max-hifi-dist INTEGER      Maximum levenshtein distance to allow for
                               hifi/CCS reads during correction.  [default: 2]
  --max-clr-dist INTEGER       Maximum levenshtein distance to allow for CLR
                               (ccs uncorrected) reads during correction.
                               [default: 3]
  --help                       Show this message and exit.

```

## Example
```shell
$ longbow correct -a 737K-august-2016.txt -m mas15v2 --barcode-freqs read_barcode_freqs.tsv tagged_array_elements.bam -o array_elements_corrected_barcodes.bam
[INFO 2022-03-01 14:08:03  correct] Invoked via: longbow correct -f -a 737K-august-2016.txt -m mas15v2 --barcode-freqs read_barcode_freqs.tsv tagged_array_elements.bam -o array_elements_corrected_barcodes.bam
[INFO 2022-03-01 14:08:03  correct] Running with 7 worker subprocess(es)
[INFO 2022-03-01 14:08:08  correct] Using mas15v2: The standard MAS-seq 15 array element model.
[INFO 2022-03-01 14:08:09  correct] Corrected tags in 570 reads of 796 total (71.61%).
[INFO 2022-03-01 14:08:09  correct] Num reads with barcodes: 796/796 (100.00%)
[INFO 2022-03-01 14:08:09  correct] Num reads without barcodes: 0/796 (0.00%)
[INFO 2022-03-01 14:08:09  correct] Done. Elapsed time: 6.09s. Overall processing rate: 130.69 reads/s.
```