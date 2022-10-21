---
layout: default
title: model 
description: "Get info on built-in models."
nav_order: 9
parent: Commands
---

# Model

## Description

Get information about built-in Longbow models.

Can list all built-in models with their version and descriptions or can dump the details of a single model to several files that contain information about that model.

## Command help

```shell
$ longbow model --help
Usage: longbow model [OPTIONS]

  Get information about built-in Longbow models.

Options:
  -v, --verbosity LVL  Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  -l, --list-models    List the names of all models supported natively by this
                       version of Longbow. NOTE: This argument is mutually
                       exclusive with  arguments: [dump].
  -d, --dump TEXT      Dump the details of a given model.  This command
                       creates a set of files corresponding to the given model
                       including: json representation, dot file
                       representations, transmission matrix, state emission
                       json file. NOTE: This argument is mutually exclusive
                       with  arguments: [list_models].
  --help               Show this message and exit.

```

## Examples

```shell
$ longbow model --list-models
[INFO 2021-11-02 15:27:28    model] Invoked via: longbow model --list-models
Longbow includes the following models:
Name                    Version  Description
mas15                   1.0.0    The standard MAS-seq 15 array element model.
mas15v2                 2.0.0    The standard MAS-seq 15 array element model.
mas15threeP             2.0.0    The 3' kit MAS-seq 15 array element model.
mas15BulkWithIndices    1.0.0    A MAS-seq 15 array element model with a 10 base index just before the 3' adapter for bulk sequencing.
mas10                   1.0.0    The MAS-seq 10 array element model.
mas10v2                 2.0.0    The MAS-seq 10 array element model.
mas10threeP             2.0.0    The 3' kit MAS-seq 10 array element model.
slide-seq               0.0.1    The Slide-seq 15 array element model.
slide-seqV2             2.0.0    The Slide-seq 15 array element model.
mas8prototype           1.0.0    The prototype MAS-seq 8 array element model.
mas8prototypeV2         2.0.0    The prototype MAS-seq 8 array element model.
```

```shell
$ longbow model --dump mas15
[INFO 2021-11-02 15:31:00    model] Invoked via: longbow model --dump mas15
[INFO 2021-11-02 15:31:00    model] Dumping mas15: The standard MAS-seq 15 array element model.
[INFO 2021-11-02 15:31:00    model] Dumping dotfile: longbow_model_mas15.v1.0.0.dot
[INFO 2021-11-02 15:31:00    model] Dumping simple dotfile: longbow_model_mas15.v1.0.0.simple.dot
[INFO 2021-11-02 15:31:00    model] Dumping json model specification: longbow_model_mas15.v1.0.0.spec.json
[INFO 2021-11-02 15:31:00    model] Dumping dense transition matrix: longbow_model_mas15.v1.0.0.dense_transition_matrix.pickle
[INFO 2021-11-02 15:31:00    model] Dumping emission distributions: longbow_model_mas15.v1.0.0.emission_distributions.txt

$ ls
longbow_model_mas15.v1.0.0.json
longbow_model_mas15.v1.0.0.dot
longbow_model_mas15.v1.0.0.simple.dot
longbow_model_mas15.v1.0.0.spec.json
longbow_model_mas15.v1.0.0.dense_transition_matrix.pickle
longbow_model_mas15.v1.0.0.emission_distributions.txt
```


