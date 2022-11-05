---
layout: default
title: extract
description: "Extract reads."
nav_order: 6
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

If a given model has a `coding_region` defined, `extract` will 
use any sections labeled with the same name as thhe `coding_region`
as the section of the read to retain during the extraction process.  

For models that do not have a `coding_region` defined, the 
cDNA portion of each read is expected to be in sections that
Longbow has labeled as `random`.  The leading and trailing 
sequences that `extract` uses as markers are configurable and 
default to the regions that most commonly flank the cDNA regions
in the current models (`BOREAS` and `Poly_A`).

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
  segments are assumed to be labeled as `random` segments or labeled with the
  same name as the `coding_region` in the model (if specified).

  For `random` segment extraction: Uses known segments flanking the region to
  be extracted as markers to indicate the start and end of what to extract.

  For `coding_region` segment extraction: Looks for any section of the reads
  labeled with the same name as the `coding_region` in the model, regardless
  of position.

Options:
  -v, --verbosity LVL         Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -p, --pbi PATH              BAM .pbi index file
  -o, --output-bam PATH       extracted bam output.  [default: stdout]
  -f, --force                 Force overwrite of the output files if they
                              exist.  [default: False]
  --create-barcode-conf-file  Create a barcode confidence score file based on
                              the barcodes in the given model.  This only
                              applies for models that have annotation_segments
                              where one such segment is annotated into the raw
                              barcode field (XC)
  -b, --base-padding INTEGER  Number of bases to include on either side of the
                              extracted region(s).  [default: 2]
  --leading-adapter TEXT      Adapter preceding the region to extract.
                              Required if the given model does not name a
                              `coding_region`.
  --trailing-adapter TEXT     Adapter following the region to extract.
                              Required if the given model does not name a
                              `coding_region`.
  --start-offset INTEGER      Number of bases to ignore from the extracted
                              region start.  These bases will not be included
                              in the extracted sequences.  Required if the
                              given model does not name a `coding_region`.
  -m, --model TEXT            The model to use for annotation.  If not
                              specified, it will be autodetected from the BAM
                              header.  If the given value is a pre-configured
                              model name, then that model will be used.
                              Otherwise, the given value will be treated as a
                              file name and Longbow will attempt to read in
                              the file and create a LibraryModel from it.
                              Longbow will assume the contents are the
                              configuration of a LibraryModel as per
                              LibraryModel.to_json().
  --help                      Show this message and exit.
```

## Examples

```shell
$ longbow extract --model mas15 --leading-adapter 10x_Adapter --trailing-adapter Poly_A --start-offset 26 -o extracted.bam filtered.bam
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

```shell
[INFO 2021-11-30 10:48:47  extract] Invoked via: longbow extract --model mas15threeP --force --create-barcode-conf-file -o extracted.bam annotated.bam
[INFO 2021-11-30 10:48:47  extract] Writing extracted read segments to: extracted.bam
[INFO 2021-11-30 10:48:47  extract] Including 2 flanking bases.
[INFO 2021-11-30 10:48:48  extract] Using mas15threeP: The 3' kit MAS-seq 15 array element model.
[INFO 2021-11-30 10:48:48  extract] Extracting coding region from model mas15threeP: cDNA
[INFO 2021-11-30 10:48:48  extract] Creating barcode confidence file: barcode_confidence_scores.txt
Progress: 926 read [00:02, 317.51 read/s]
[INFO 2021-11-30 10:48:50  extract] Done. Elapsed time: 2.92s.
[INFO 2021-11-30 10:48:50  extract] Total # Reads Processed: 926
[INFO 2021-11-30 10:48:50  extract] # Reads Containing Extracted Segments: 905 (97.73%)
[INFO 2021-11-30 10:48:50  extract] Total # Segments Extracted: 7899
[INFO 2021-11-30 10:48:50  extract] Total # Segments Skipped: 68
[INFO 2021-11-30 10:48:50  extract] # Segments extracted per read: 8.53
```
