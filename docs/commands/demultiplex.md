---
layout: default
title: demultiplex
description: "Demultiplex reads from different array structures."
nav_order: 6
parent: Commands
---

# Demultiplex

## Description

## Command help

```
$ longbow demultiplex --help
Usage: longbow demultiplex [OPTIONS] INPUT_BAM

  Separate reads into files based on which model they fit best.

  Resulting reads will be annotated with the model they best fit as well as
  the score and segments for that model.

Options:
  -v, --verbosity LVL       Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -p, --pbi PATH            BAM .pbi index file
  -t, --threads INTEGER     number of threads to use (0 for all)  [default:
                            11]

  -o, --out-base-name TEXT  base name for output files  [default:
                            longbow_demultiplexed]

  --help                    Show this message and exit.
```

## Example
