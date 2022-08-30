FROM python:3.7.9-slim-buster

MAINTAINER Kiran V Garimella

ARG branch

# install gsutil
RUN apt-get update --allow-releaseinfo-change
RUN apt-get update -m
RUN apt-get install -y curl git-lfs time datamash wget zlib1g-dev libssl-dev lbzip2 make gcc libbz2-dev libncurses5-dev libncursesw5-dev liblzma-dev
RUN curl https://sdk.cloud.google.com | bash

# Setup crcmodc for gsutil:
RUN apt-get install -y gcc python3-dev python3-setuptools && \
			pip3 uninstall -y crcmod && \
			pip3 install --no-cache-dir -U crcmod

# Fix python malloc bug
# Detailed here: https://github.com/tensorflow/tensorflow/issues/6968 
RUN apt-get install -y libtcmalloc-minimal4
RUN echo 'export LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libtcmalloc_minimal.so.4"' >> ~/.bashrc

# Create known hosts for github ssh info:
RUN mkdir -p ~/.ssh && touch ~/.ssh/known_hosts
RUN ssh-keyscan github.com >> ~/.ssh/known_hosts

# Install Longbow
RUN git clone http://github.com/broadinstitute/longbow.git
RUN cd longbow \
    && git checkout $branch \
    && python3 -mvenv venv \
    && . venv/bin/activate \
    && pip install -r dev-requirements.txt \
    && pip install -e .

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2 \
    && tar xjf samtools-1.14.tar.bz2 \
    && rm samtools-1.14.tar.bz2 \
    && cd samtools-1.14 \
    && ./configure \
    && make \
    && make install

# install Starcode
RUN wget -O /usr/local/bin/starcode https://github.com/gui11aume/starcode/releases/download/1.4/starcode-1.4 \
    && chmod 755 /usr/local/bin/starcode

# Install vim so we can edit files later:
RUN apt-get install -y vim

# activate conda environment
ENV PATH="/longbow/venv/bin:/root/google-cloud-sdk/bin/:$PATH"
RUN echo "source /longbow/venv/bin/activate" >> ~/.bashrc

# We want this to be interactive:
CMD ["/bin/bash"]
