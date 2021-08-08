---
layout: default
title: filter
description: "Filter reads."
nav_order: 3
parent: Commands
---

# Filter

```
> longbow filter --help
Usage: longbow filter [OPTIONS] INPUT_BAM

  Filter reads by whether they conform to expected segment order.

Options:
  -v, --verbosity LVL    Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -p, --pbi PATH         BAM .pbi index file
  -o, --out-prefix TEXT  Output file prefix  [required]
  -m, --model TEXT       The model to use for annotation.  If the given value
                         is a pre-configured model name, then that model will
                         be used.  Otherwise, the given value will be treated
                         as a file name and Longbow will attempt to read in
                         the file and create a LibraryModel from it.  Longbow
                         will assume the contents are the configuration of a
                         LibraryModel as per LibraryModel.to_json().
                         [default: mas15]

  --force                Force overwrite of the output files if they exist.
                         [default: False]

  --help                 Show this message and exit.
```
