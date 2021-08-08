---
layout: default
title: segment
description: "Segment reads."
nav_order: 2
parent: Commands
---

# Segment

```
> longbow segment --help
Usage: longbow segment [OPTIONS] INPUT_BAM

  Segment pre-annotated reads from an input BAM file.

Options:
  -v, --verbosity LVL        Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -t, --threads INTEGER      number of threads to use (0 for all)  [default:
                             11]

  -o, --output-bam PATH      segment-annotated bam output  [default: stdout]
  -s, --do-simple-splitting  Do splitting of reads based on splitter
                             delimiters, rather than whole array structure.
                             This splitting will cause delimiter sequences to
                             be repeated in each read they bound.

  -m, --model TEXT           The model to use for annotation.  If the given
                             value is a pre-configured model name, then that
                             model will be used.  Otherwise, the given value
                             will be treated as a file name and Longbow will
                             attempt to read in the file and create a
                             LibraryModel from it.  Longbow will assume the
                             contents are the configuration of a LibraryModel
                             as per LibraryModel.to_json().  [default: mas15]

  --help                     Show this message and exit.
```