---
layout: default
title: Installation
nav_order: 3
description: "Installation instructions."
---

# Installation

## From Github release

```shell
$ wget https://github.com/broadinstitute/longbow/archive/refs/tags/v0.3.0.tar.gz \
    && tar zxvf v0.3.0.tar.gz \
    && cd longbow-0.3.0/

$ python -mvenv venv \
    && . venv/bin/activate \
    && pip install -r dev-requirements.txt \
    && pip install -e .
```

## From Github source

```shell
$ git clone git@github.com:broadinstitute/longbow.git \
    && cd longbow

$ python -mvenv venv \
    && . venv/bin/activate \
    && pip install -r dev-requirements.txt \
    && pip install -e .
```

## From Docker
If you have Docker installed, you can use the pre-built Docker image which has all dependencies already installed:

```shell
$ docker run -v $PWD:/data -it "us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0"
(venv) root@3d3b800e48ca:/# longbow version
0.3.0
```

## From conda
(TODO: conda recipes to be written before 1.0 release.)


## Testing installation

Have Longbow print its version number to verify things have been installed correctly:

```shell
$ longbow version
0.3.0
```

Note: first invocation may be slow (owing to longbow commands being discovered by introspection, rather than using hardcoded imports).  Introspective imports will be replaced for 1.0 release.
