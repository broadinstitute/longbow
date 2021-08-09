---
layout: default
title: Installation
nav_order: 3
description: "Installation instructions."
---

# Installation

## From scratch

```
wget https://github.com/broadinstitute/longbow/archive/refs/tags/v0.3.0.tar.gz
tar zxvf v0.3.0.tar.gz
cd longbow-0.3.0/

python -mvenv venv
. venv/bin/activate
pip install -r dev-requirements.txt
pip install -e .
```

Note: conda/pip recipes to be written before 1.0 release.

## From Docker

```
docker run -it us.gcr.io/broad-dsp-lrma/lr-longbow:0.3.0
```

## Testing installation

Have Longbow print its version number to verify things have been installed correctly:

```
> longbow version
0.3.0
```

Note: first invocation may be slow (owing to longbow commands being discovered by introspection, rather than using hardcoded imports).  Introspective imports will be replaced for 1.0 release.
