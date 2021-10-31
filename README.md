# Longbow 
Annotation and segmentation of MAS-seq data

Current version: Current version: Current version: 0.4.5

## Documentation

[Longbow documentation.](https://broadinstitute.github.io/longbow/)

## Installation

    pip install .

## Usage

TBW

## Development

To do development in this codebase, the python3 development package must
be installed.

After installation the development environment can be set up by
the following commands:

    python3 -mvenv venv
    . venv/bin/activate
    pip install --upgrade pip
    pip install -r dev-requirements.txt
    pip install -e .

### Linting files

    # run all linting commands
    tox -e lint

    # reformat all project files
    black src tests setup.py

    # sort imports in project files
    isort -rc src tests setup.py

    # check pep8 against all project files
    flake8 src tests setup.py

    # lint python code for common errors and codestyle issues
    pylint src

### Tests

    # run all linting and test
    tox

    # run only (fast) unit tests
    tox -e unit
    
    # run only integration tests
    tox -e integration

    # run only linting
    tox -e lint

Note: If you run into "module not found" errors when running tox for testing, verify the modules are listed in test-requirements.txt and delete the .tox folder to force tox to refresh dependencies.

## Versioning

We use `bumpversion` to maintain version numbers.  It will automatically create a new tag every time you run it.
*DO NOT MANUALLY EDIT ANY VERSION NUMBERS.*

Our versions are specified by a 3 number semantic version system (https://semver.org/):

	major.minor.patch

To update the version with bumpversion do the following:

`bumpversion PART` where PART is one of:
- major
- minor
- patch

This will increase the corresponding version number by 1.

You should always bump the version from the *main branch* _after_ merging in any PRs for that version.  Then you will have to push both the bumpversion code update *AND* the tag that was created:

```
git push && git push --tags
```

## Releasing

The release process is as follows:

1. Bump the version by the appropriate PART (see *Versioning* for details).
2. Build and push the docker container (See *Docker*).  NOTE: The version will have already been updated by using `bumpversion` in step 1.

## Docker 

This repository includes a docker folder that contains a `Dockerfile` and a `Makefile`.  
*DO NOT manually modify the Makefile or the Dockerfile* (unless you have to add dependencies).  The version is automatically updated by `bumpversion`. 

### Building

To build the docker image and push it to the release location do the following from the root of this repository:

```
cd docker
make build && make push
```

