#!/bin/bash

# exit when your script tries to use undeclared variables
set -o nounset
# exit if any pipe commands fail
set -o pipefail
# exit when a command fails
set -o errexit
# # trace what gets executed
# set -o xtrace

# data inputs
WORKDIR=$(realpath "./" )
BASEDIR=$(realpath "../../")
CONTAINERDIR="${BASEDIR}/containers"

####################################################
# 0. Select Singularity or Docker containerization #
####################################################
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

#################
# 5. Fit models #
#################
# build barnacle container
printf "\nChecking ${CONTAINER} container\n"
if [ "${CONTAINER}" == "singularity" ]; then
    if [[ ! -e "${CONTAINERDIR}/barnacle.sif" ]]; then 
        singularity build "${CONTAINERDIR}/barnacle.sif" "docker-archive://${CONTAINERDIR}/barnacle.tar.gz"
    fi
elif [ "${CONTAINER}" == "docker" ]; then
    # export poetry environment requirements
    if [[ ! -e "${CONTAINERDIR}/requirements.txt" ]]; then
        poetry export --without-hashes >> "${CONTAINERDIR}/requirements.txt"
    fi
    # make Docker containers from Dockerfiles
    docker image build -t barnacle "${CONTAINERDIR}"
fi
printf "\nRunning parameter grid search\n"
if [ "${CONTAINER}" == "singularity" ]; then
    singularity exec --bind ${BASEDIR}:${BASEDIR} --pwd "${WORKDIR}" ${CONTAINERDIR}/barnacle.sif python grid-search.py
elif [ "${CONTAINER}" == "docker" ]; then
    docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} -w "${WORKDIR}" barnacle python grid-search.py
fi
