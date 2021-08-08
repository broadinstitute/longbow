---
layout: default
title: Commands
description: "Longbow commands."
nav_order: 4
has_children: true
---

# Commands

$ longbow --help
Usage: longbow [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  annotate     Annotate reads in a BAM file with segments from the model.
  demultiplex  Separate reads into files based on which model they fit best.
  extract      Extract coding segments from the reads in the given bam.
  filter       Filter reads by whether they conform to expected segment...
  inspect      Inspect the classification results on specified reads.
  scsplit      Create files for use in `alevin` for single-cell analysis.
  segment      Segment pre-annotated reads from an input BAM file.
  train        Train transition and emission probabilities on real data.
  version      Print the version of longbow.