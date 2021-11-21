FROM python:3.7.9-slim-buster

MAINTAINER Kiran V Garimella

# install gsutil
RUN apt-get update
RUN apt-get install -y curl git-lfs time datamash
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# Fix python malloc bug
# Detailed here: https://github.com/tensorflow/tensorflow/issues/6968 
RUN apt-get install -y libtcmalloc-minimal4
RUN echo 'export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"' >> ~/.bashrc

# Set up pip environment (see https://pythonspeed.com/articles/activate-virtualenv-dockerfile/)
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN git clone https://github.com/broadinstitute/longbow.git && \
    cd longbow && \
    git checkout f3e3dbabfb9b3b035e9c1d9c318e8e8c5f35c4ad && \
    pip install -r dev-requirements.txt && \
    pip install -e .

# activate conda environment
RUN echo "source /opt/venv/bin/activate" >> ~/.bashrc

# We want this to be interactive
CMD ["/bin/bash"]
