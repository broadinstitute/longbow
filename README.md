# Python CLI Project Template
This template was created so that a new python3 CLI program could be quickly whipped up with a nice infrastructure. 

## Acknowledgements
This template is based on several projects that originated at or as part of collaborations involving The Broad Institute of MIT and Harvard.
In particular, the following people were instrumental in those projects: 

- Kiran Garimella (@kvg)
- Winni Kretzschmar (@winni2k)
- Karl Johan Westrin (@karljohanw)

## Contents / Features
Included in this template are the following features:
- `click` for argument parsing with examples
- an over-engineered logging module (ensures that all log messages have correct leading whitespace to be column-aligned)
- `tox` setup (linting, black for code formatting, testing)
- examples of tests (already connected to tox for quick starting)
- `.gitignore` for most kinds of local files
- include sort ordering defaults 
- a `readthedocs` template
- a `bumpversion` config file
- a `pylint` config file

If you check this out and run `tox` the tests will all pass.

## Make It Your Own
### Required Changes
To make this project-specific, several placeholder values need to be replaced with your content.  

| FIELD | DESCRIPTION |
| ----- | ----------- |
|  PROJECT\_NAME | The name of your new project / tool. |
| \_SHORT\_PROJECT\_DESCRIPTION\_ | A concise description of your new project / tool. |
| \_AUTHOR\_ | The name of the primary author / point of contact. |
| \_AUTHOR\_EMAIL\_ | The email address of the primary author / point of contact. |

This find/replace will be needed until github [implements template repo variables.](https://github.com/isaacs/github/issues/1716)

However, fear not!  I have provided `initialize.sh` for this very purpose.  Run it with the following parameters after creating your repo to initialize your specific repository info:

```
./initialize.sh PROJECT_NAME "SHORT PROJECT DESCRIPTION" "AUTHOR NAME" AUTHOR_EMAIL
```
I quoted the fields that are likely to contain spaces - this will be necessary.

### Optional Changes
You may also need to do the following:
- Update the LICENSE file - it defaults to the BSD 3-Clause License.

When you add development dependencies, add them to `dev-requirements.txt`.

When you add testing dependencies, add them to `test-requirements.txt`.

When you add runtime dependencies, add them to `setup.py` in the `install_requires` section.

## General Notes

When running `tox`, you'll notice that the linter runs before black8.  This is intentional.  Rather than have it blindly reformat your code, I wanted to make sure you knew what you were doing wrong (i.e. against PEP 8 standards).  You can configure this order to suit your needs.

Notable projects already using this template are:
- [Cromshell 2.0](https://github.com/broadinstitute/cromshell/tree/cromshell_2.0)
- [CARROT CLI](https://github.com/broadinstitute/carrot_cli)

Lastly, it's worth double-checking everything to make sure you want what's included here.

Below this line begins the project-specific portion of the README that should be modified after creating a new project:

---

# PROJECT\_NAME 
\_SHORT\_PROJECT\_DESCRIPTION\_

Current version: 0.0.1

## Installation

    pip install .

## Development

To do development in this codebase, the python3 development package must
be installed.

After installation the development environment can be set up by
the following commands:

    python3 -mvenv venv
    . venv/bin/activate
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

    # run only linting
    tox -e lint

Note: If you run into "module not found" errors when running tox for testing, verify the modules are listed in test-requirements.txt and delete the .tox folder to force tox to refresh dependencies.

### Versioning

We use `bumpversion` to maintain version numbers.
*DO NOT MANUALLY EDIT ANY VERSION NUMBERS.*

Our versions are specified by a 3 number semantic version system (https://semver.org/):

	major.minor.patch

To update the version with bumpversion do the following:

`bumpversion PART` where PART is one of:
- major
- minor
- patch

This will increase the corresponding version number by 1.

