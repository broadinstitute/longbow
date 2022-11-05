---
layout: default
title: train
description: "Train model."
nav_order: 16
parent: Commands
---

# Train

## Description

Train transition and emission probabilities using the Baum-Welch algorithm.

(This command is currently being rewritten and is not recommended for use by anyone other than those who wish to develop the Longbow software package.)

## Command help

```shell
$ longbow train --help
Usage: longbow train [OPTIONS] TRAINING_BAM

  Train transition and emission probabilities on real data.

Options:
  -v, --verbosity LVL             Either CRITICAL, ERROR, WARNING, INFO or
                                  DEBUG

  -n, --num-training-samples INTEGER
                                  number of training samples to use  [default:
                                  10]

  -i, --max-training-iterations INTEGER
                                  number of training iterations to use
                                  [default: 5]

  -t, --threads INTEGER           number of threads to use (0 for all)
                                  [default: 11]

  -o, --output-yaml PATH          trained model  [required]
  -m, --model TEXT                The model to use for annotation.  If the
                                  given value is a pre-configured model name,
                                  then that model will be used.  Otherwise,
                                  the given value will be treated as a file
                                  name and Longbow will attempt to read in the
                                  file and create a LibraryModel from it.
                                  Longbow will assume the contents are the
                                  configuration of a LibraryModel as per
                                  LibraryModel.to_json().  [default: mas15]

  --help                          Show this message and exit.
```

## Example
