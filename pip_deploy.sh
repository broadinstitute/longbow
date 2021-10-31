#!/bin/bash

set -euxo pipefail

echo 'Creating and uploading pip package...' && \
    python setup.py check && \
    python setup.py sdist && \
    python setup.py bdist_wheel --universal && \
    twine upload --repository-url https://test.pypi.org/legacy/ dist/* && \
    twine upload dist/*

