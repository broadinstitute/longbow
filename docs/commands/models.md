---
layout: default
title: models
description: "Get info on built-in models."
nav_order: 9
parent: Commands
---

# Models

## Description

Get information about built-in Longbow models.

Can list all built-in models with their version and descriptions or can dump the details of a single model to several files that contain information about that model.

## Command help

```shell
$ longbow models --help
Usage: longbow models [OPTIONS]

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
$ longbow models --list-models
Longbow includes the following models:

Array models
============
Name      Version  Description
mas_15    3.0.0    15-element MAS-ISO-seq array
mas_10    3.0.0    10-element MAS-ISO-seq array
isoseq    3.0.0    PacBio IsoSeq model

cDNA models
===========
Name                Version  Description
sc_10x3p            3.0.0    single-cell 10x 3' kit
sc_10x5p            3.0.0    single-cell 10x 5' kit
bulk_10x5p          3.0.0    bulk 10x 5' kit
bulk_teloprimeV2    3.0.0    Lexogen TeloPrime V2 kit
spatial_slideseq    3.0.0    Slide-seq protocol

Specify a fully combined model via '<array model>+<cDNA model>' syntax, e.g. 'mas_15+sc_10x5p'.

```

```shell
$ longbow models --dump mas_15+sc_10x5p
[INFO 2022-11-02 20:33:28   models] Generating model: mas_15+sc_10x5p
[INFO 2022-11-02 20:33:29   models] Dumping mas_15+sc_10x5p: 15-element MAS-ISO-seq array, single-cell 10x 5' kit
[INFO 2022-11-02 20:33:29   models] Dumping dotfile: longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.dot
[INFO 2022-11-02 20:33:29   models] Dumping simple dotfile: longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.simple.dot
[INFO 2022-11-02 20:33:29   models] Dumping json model specification: longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.spec.json
[INFO 2022-11-02 20:33:30   models] Dumping dense transition matrix: longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.dense_transition_matrix.pickle
[INFO 2022-11-02 20:33:30   models] Dumping emission distributions: longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.emission_distributions.txt
[INFO 2022-11-02 20:33:30   models] Creating model graph from 1109 states...
[INFO 2022-11-02 20:33:47   models] Rendering model graph now...
[INFO 2022-11-02 20:33:57   models] Writing model graph now to longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.graph.png ...

$ ls 
longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.dense_transition_matrix.pickle  longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.graph.png
longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.dot                             longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.simple.dot
longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.emission_distributions.txt      longbow_model-mas_15+sc_10x5p-Av3.0.0_Cv3.0.0.spec.json
```


