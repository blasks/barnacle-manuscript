#!/bin/bash

# Use these commands to build all the requisite  docker images for running the
# analyses in this manuscript. Standard bioinformatics tool containers are 
# preferentially pulled from the developer when possible. If no container is 
# maintained by the developer, the container is pulled from biocontainers
# (see https://biocontainers.pro/). Two custom containers are also constructed:
# one with the software necessary for pre-processing fastq files, and another
# for running all Python code.

# Pull Docker containers for tools in which one is maintained by the developer.
docker pull combinelab/salmon:1.10.2

# # make and save docker image
# docker image build -t fastq-preprocess docker/preprocessing/
# docker save fastq-preprocess | gzip > fastq-preprocess.tar.gz

