from io import open

from setuptools import find_packages, setup

# following src dir layout according to
# https://blog.ionelmc.ro/2014/05/25/python-packaging/#the-structure
version = "version = "version = "0.4.4"""
setup(
    name="maslongbow",
    version=version,
    description="Annotation and segmentation of MAS-seq data",
    url="https://broadinstitute.github.io/longbow/",
    author="Kiran V Garimella, Jonn Smith",
    author_email="kiran@broadinstitute.org, jonn@broadinstitute.org",
    license="BSD 3-Clause",
    long_description=open("README.md").read(),
    install_requires=[
        'alabaster==0.7.12',
        'appdirs==1.4.4',
        'astroid==2.5.6',
        'Babel==2.9.1',
        'black==21.5b1',
        'bleach==3.3.0',
        'bump2version==1.0.1',
        'bumpversion==0.6.0',
        'certifi==2020.12.5',
        'chardet==4.0.0',
        'click==7.1.2',
        'click-log==0.3.2',
        'colorama==version = "version = "0.4.4""',
        'construct==2.10.67',
        'cycler==0.10.0',
        'Cython==0.29.23',
        'databind.core==1.1.2',
        'databind.json==1.1.2',
        'decorator==4.4.2',
        'Deprecated==1.2.12',
        'dill==0.3.3',
        'distlib==0.3.1',
        'docspec==1.0.1',
        'docspec-python==1.0.1',
        'docutils==0.17.1',
        'filelock==3.0.12',
        'flake8==3.9.2',
        'idna==2.10',
        'imagesize==1.2.0',
        'importlib-metadata==4.0.1',
        'isort==4.3.21',
        'Jinja2==2.11.3',
        'joblib==1.0.1',
        'keyring==23.0.1',
        'kiwisolver==1.3.1',
        'lazy-object-proxy==1.6.0',
        'MarkupSafe==1.1.1',
        'maslongbow==version = "version = "0.4.4""',
        'matplotlib==3.4.2',
        'mccabe==0.6.1',
        'multiprocess==0.70.11.1',
        'mypy==0.770',
        'mypy-extensions==0.4.3',
        'networkx==2.5.1',
        'nr.fs==1.6.3',
        'nr.optional==0.1.1',
        'nr.parsing.date==1.0.2',
        'nr.preconditions==0.0.4',
        'nr.pylang.utils==0.1.3',
        'nr.stream==0.1.2',
        'nr.utils.re==0.3.1',
        'numpy==1.20.3',
        'packaging==20.9',
        'pandas==1.2.4',
        'pathos==0.2.7',
        'pathspec==0.8.1',
        'Pillow==8.2.0',
        'pkginfo==1.7.0',
        'pluggy==0.13.1',
        'pomegranate==0.14.5',
        'pox==0.2.9',
        'ppft==1.6.6.3',
        'progress==1.5',
        'py==1.10.0',
        'pycodestyle==2.7.0',
        'pydoc-markdown==4.1.5',
        'pyflakes==2.3.1',
        'Pygments==2.9.0',
        'pylint==2.8.2',
        'pyparsing==2.4.7',
        'pysam==0.16.0.1',
        'python-dateutil==2.8.1',
        'pytz==2021.1',
        'PyYAML==5.4.1',
        'readme-renderer==29.0',
        'regex==2021.4.4',
        'requests==2.25.1',
        'requests-toolbelt==0.9.1',
        'rfc3986==1.5.0',
        'scipy==1.6.3',
        'six==1.16.0',
        'snowballstemmer==2.1.0',
        'Sphinx==4.0.1',
        'sphinx-rtd-theme==0.5.2',
        'sphinxcontrib-applehelp==1.0.2',
        'sphinxcontrib-devhelp==1.0.2',
        'sphinxcontrib-htmlhelp==1.0.3',
        'sphinxcontrib-jsmath==1.0.1',
        'sphinxcontrib-qthelp==1.0.3',
        'sphinxcontrib-serializinghtml==1.1.4',
        'ssw==0.4.1',
        'toml==0.10.2',
        'tox==3.23.1',
        'tqdm==4.60.0',
        'twine==3.4.1',
        'typed-ast==1.4.3',
        'typing-extensions==3.10.0.0',
        'urllib3==1.26.4',
        'virtualenv==20.4.6',
        'watchdog==2.1.3',
        'webencodings==0.5.1',
        'wrapt==1.12.1',
        'zipp==3.4.1',
    ],
    tests_require=["coverage", "pytest"],
    python_requires=">=3.6, <=3.7.9",
    packages=find_packages("src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
    ],
    entry_points={"console_scripts": ["longbow=longbow.__main__:main_entry"]},
    include_package_data=True,
)
