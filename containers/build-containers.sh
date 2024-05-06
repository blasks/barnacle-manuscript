#!/bin/bash

# Use these commands to build all the requisite  containers for running the
# analyses in this manuscript. Standard bioinformatics tool containers are 
# preferentially pulled from the developer when possible. If no container is 
# maintained by the developer, the container is pulled from biocontainers
# (see https://biocontainers.pro/). Two custom containers are also constructed:
# one with the software necessary for pre-processing fastq files, and another
# for running all Python code.

#################################################
# Select Singularity or Docker containerization #
#################################################
CONTAINER=""
while [ "${CONTAINER}" == "" ]; do
    printf "\nPlease select an option:\n\n\t1 - singularity\n\t2 - docker\n\nEnter '1' or '2': "
    read SELECTION 
    if [ "${SELECTION}" == "1" ]; then
        echo "Continuing with Singularity containerized workflow"
        CONTAINER="singularity"
    elif [ "${SELECTION}" == "2" ]; then
        echo "Continuing with Docker containerized workflow"
        CONTAINER="docker"
    else 
        printf "\nInvalid selection\n"
    fi
done

####################
# Build containers #
####################
if [ "${CONTAINER}" == "singularity" ]; then
    # Pull Singularity containers for tools in which one is maintained by the developer.
    singularity build salmon.sif docker://combinelab/salmon:1.10.2
    singularity build sra-tools.sif docker://quay.io/biocontainers/sra-tools:3.1.0--h4304569_1
    singularity build trimmomatic.sif docker://quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4 
    # # build singularity image from tarball
    # singularity build barnacle.sif docker-archive://barnacle.tar.gz
elif [ "${CONTAINER}" == "docker" ]; then
    # Pull Docker containers for tools in which one is maintained by the developer.
    docker pull combinelab/salmon:1.10.2
    docker pull quay.io/biocontainers/sra-tools:3.1.0--h4304569_1
    docker pull quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4 
    # export poetry environment requirements
    if [[ ! -e requirements.txt ]]; then
        poetry export --without-hashes >> requirements.txt
    fi
    # make Docker containers from Dockerfiles
    docker image build -t barnacle .
    # # save tarball of barnacle image to convert to Singularity image
    # docker save barnacle | gzip > barnacle.tar.gz
fi

