#!/bin/bash

# Use this command to build a singularity image from a 
# Docker image tarball


# build singularity image from tarball
singularity build fastq-preprocess.sif docker-archive://fastq-preprocess.tar.gz
