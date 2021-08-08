---
layout: default
title: inspect 
description: "Inspect reads."
nav_order: 5
parent: Commands
---

# Inspect

```
> longbow inspect --help
Usage: longbow inspect [OPTIONS] INPUT_BAM

  Inspect the classification results on specified reads.

Options:
  -v, --verbosity LVL          Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -r, --read-names TEXT        read names (or file(s) of read names) to
                               inspect

  -p, --pbi PATH               BAM .pbi index file
  -f, --file-format [png|pdf]  Image file format
  -o, --outdir PATH            Output directory
  -m, --model TEXT             The model to use for annotation.  If the given
                               value is a pre-configured model name, then that
                               model will be used.  Otherwise, the given value
                               will be treated as a file name and Longbow will
                               attempt to read in the file and create a
                               LibraryModel from it.  Longbow will assume the
                               contents are the configuration of a
                               LibraryModel as per LibraryModel.to_json().
                               [default: mas15]

  --seg-score                  Display alignment score for annotated segments.
                               [default: False]

  --help                       Show this message and exit.
```