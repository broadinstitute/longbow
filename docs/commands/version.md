---
layout: default
title: version
description: "Print version number."
nav_order: 17
parent: Commands
---

# Version

## Description

Print the version of Longbow.

Longbow uses [Semantic Versioning](https://semver.org/) (MAJOR.MINOR.PATCH).

## Command help

```shell
$ longbow version --help
Usage: longbow version [OPTIONS]

  Print the version of longbow.

Options:
  -v, --verbosity LVL  Either CRITICAL, ERROR, WARNING, INFO or DEBUG
  --help               Show this message and exit.
```

## Example

```shell
$ longbow version
0.3.0
```