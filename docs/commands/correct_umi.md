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

### Data Requirements:

- Bam file should be aligned and annotated with genes and transcript equivalence classes prior to running.

The following tags are required in the input file:

- `JX` (Adjusted UMI)
- `eq` (Equivalence class assignment)
- `XG` (Gene assignment)
- `rq` (Read Quality: [-1.0, 1.0])
- `JB` (Back / UMI trailing segment Smith-Waterman alignment score)

## Command help

```shell
$ longbow correct_umi --help
Usage: longbow correct_umi [OPTIONS] INPUT_BAM

  Correct UMIs with Set Cover algorithm.

Options:
  -v, --verbosity LVL       Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -l, --umi-length INTEGER  Length of the UMI for this sample.  [default: 10]
  -o, --output-bam PATH     Corrected UMI bam output [default: stdout].
  -x, --reject-bam PATH     Filtered bam output (failing reads only).
  -f, --force               Force overwrite of the output files if they exist.
                            [default: False]
  --pre-extracted           Whether the input file has been processed with
                            `longbow extract`  [default: True]
  --help                    Show this message and exit.
```


