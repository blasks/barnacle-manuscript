#!/bin/bash

# Use these commands to build all the requisite  Singularity images for running 
# the analyses in this manuscript. Standard bioinformatics tool containers are 
# preferentially pulled from the developer when possible. If no container is 
# maintained by the developer, the container is pulled from biocontainers
# (see https://biocontainers.pro/). Two custom containers are also constructed:
# one with the software necessary for pre-processing fastq files, and another
# for running all Python code.

# Pull Singularity containers for tools in which one is maintained by the developer.
singularity build salmon.sif docker://combinelab/salmon:1.10.2

# # build singularity image from tarball
# singularity build fastq-preprocess.sif docker-archive://fastq-preprocess.tar.gz
