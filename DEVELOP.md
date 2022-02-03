# Longbow 
Annotation and segmentation of MAS-seq data

Current version: 0.5.13

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

## Releasing

Longbow is released automatically upon pushes to the main branch (generally by merging branches into main).  Automatic releases of Longbow consist of the following steps:

### Versioning

We use `bumpversion` to maintain version numbers.  It will automatically create a new tag each time it is run.
*DO NOT MANUALLY EDIT ANY VERSION NUMBERS.*

Our versions are specified by a 3 number semantic version system (https://semver.org/):

	major.minor.patch

By default, pushes to main will increment the patch number.  Major and minor version numbers can be incremented through special keywords in commit messages.  From the [automated-version-bump](https://github.com/marketplace/actions/automated-version-bump) Github action documentation:

> Based on the commit messages, increment the version from the latest release.
> * If the string "BREAKING CHANGE", "major" or the Attention pattern refactor!: drop support for Node 6 is found anywhere in any of the commit messages or descriptions the major version will be incremented.
> * If a commit message begins with the string "feat" or includes "minor" then the minor version will be increased. This works for most common commit metadata for feature additions: "feat: new API" and "feature: new API".
> * If a commit message contains the word "pre-alpha" or "pre-beta" or "pre-rc" then the pre-release version will be increased (for example specifying pre-alpha: 1.6.0-alpha.1 -> 1.6.0-alpha.2 or, specifying pre-beta: 1.6.0-alpha.1 -> 1.6.0-beta.0)
All other changes will increment the patch version.

Versions will always be bumped from the *main branch* _after_ merging in any PRs for that version.

### Docker image

This repository includes a Docker folder that contains a `Dockerfile`. New images will be built automatically and pushed to [Google's Docker registry](https://console.cloud.google.com/gcr/images/broad-dsp-lrma/US/lr-longbow) . 

*DO NOT manually modify the Makefile or the Dockerfile* (unless you have to add dependencies).  The version is automatically updated by `bumpversion`. 

### PyPI

A PyPI distribution will be built automatically and uploaded to https://pypi.org/project/maslongbow .
