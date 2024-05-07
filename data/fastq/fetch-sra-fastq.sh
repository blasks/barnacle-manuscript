#!/bin/bash
# AUTHOR: Stephen Blaskowski

# Fetch raw fastq data from the NCBI Sequence Read Archive (SRA)
# This script depends on the prefetch and fasterq-dump functions of 
# the sra-tools software suite. Instructions for using sra-tools can be found
# on the NCBI website, or the project repository page:
# https://github.com/ncbi/sra-tools/wiki

# Inputs:
BASEDIR=$(realpath "../../")
FASTQDIR="${BASEDIR}/data/fastq"
SRADIR="${FASTQDIR}/sra"
RAWDIR="${FASTQDIR}/raw"
CONTAINERDIR="${BASEDIR}/containers"
METADATA="SRA-Sample-Metadata.csv"
THREADS=32


#####################
# 0. Select options #
#####################
# Singularity or Docker containerization
CONTAINER=""
while [ "${CONTAINER}" == "" ]; do
    printf "\nPlease select an option:\n\n\t1 - Singularity\n\t2 - Docker\n\nEnter '1' or '2': "
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
# File cleanup
CLEANUP=""
while [ "${CLEANUP}" == "" ]; do
    printf "\nWould you like to save storage by removing SRA files following fastq extraction?\n\n\t1 - Yes\n\t2 - No\n\nEnter '1' or '2': "
    read SELECTION 
    if [ "${SELECTION}" == "1" ]; then
        echo "Will remove SRA files following fastq extraction"
        CLEANUP="True"
    elif [ "${SELECTION}" == "2" ]; then
        echo "Will not remove any files"
        CLEANUP="False"
    else 
        printf "\nInvalid selection\n"
    fi
done

##########################
# 1. Download containers #
##########################
printf "\n\t* Checking ${CONTAINER} container\n"
if [ "${CONTAINER}" == "singularity" ]; then
    if [[ ! -e "${CONTAINERDIR}/sra-tools.sif" ]]; then 
        singularity build "${CONTAINERDIR}/sra-tools.sif" docker://quay.io/biocontainers/sra-tools:3.1.0--h4304569_1
    fi
elif [ "${CONTAINER}" == "docker" ]; then
    docker pull quay.io/biocontainers/sra-tools:3.1.0--h4304569_1
fi

#############################
# 2. Download data from SRA #
#############################
if [[ ! -d ${SRADIR} ]]; then
    mkdir -p ${SRADIR}
fi
for SRA in `tail -n +2 $METADATA | cut -d "," -f 3`; do
    # prefetch SRA
    printf "\n\t* Downloading ${SRA}\n"
    if [ "${CONTAINER}" == "singularity" ]; then
        singularity exec --bind ${BASEDIR}:${BASEDIR} ${CONTAINERDIR}/sra-tools.sif prefetch \
            ${SRA} -O ${SRADIR} --max-size 1t
    elif [ "${CONTAINER}" == "docker" ]; then
        docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} \
            quay.io/biocontainers/sra-tools:3.1.0--h4304569_1 prefetch \
            ${SRA} -O ${SRADIR} --max-size 1t
    fi
done

##################################
# 3. Extract fastq data from SRA #
##################################
if [[ ! -d ${RAWDIR} ]]; then
    mkdir -p ${RAWDIR}
fi
for SRA in `tail -n +2 $METADATA | cut -d "," -f 3`; do
    # extract fastq files
    printf "\n\t* Extracting fastq files for ${SRA}\n"
    if [ "${CONTAINER}" == "singularity" ]; then
        singularity exec --bind ${BASEDIR}:${BASEDIR} --pwd ${SRADIR} ${CONTAINERDIR}/sra-tools.sif fasterq-dump \
            ${SRA} -O ${RAWDIR} --threads ${THREADS} -v
    elif [ "${CONTAINER}" == "docker" ]; then
        docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} --workdir ${SRADIR} \
            quay.io/biocontainers/sra-tools:3.1.0--h4304569_1 fasterq-dump \
            ${SRA} -O ${RAWDIR} --threads ${THREADS} -v
    fi
    # clean up
    if [ "${CLEANUP}" == "True" ]; then
        rm -rf "${SRADIR}/${SRA}"
    fi
    # gzip files in parallel every ${THREADS} files
    FILES=$( ls ${RAWDIR}/*.fastq | wc -l )
    if [ $(( ${FILES} % ${THREADS} )) == 0 ]; then
        printf "\n\t* Compressing ${THREADS} fastq files in parallel\n"
        find ${RAWDIR}/*.fastq | parallel gzip
    fi
done
# gzip remaining fastq files
FILES=$( ls ${RAWDIR}/*.fastq | wc -l )
printf "\n\t* Compressing ${FILES} fastq files in parallel\n"
find ${RAWDIR}/*.fastq | parallel gzip
# clean up
if [ "${CLEANUP}" == "True" ]; then
    rmdir "${SRADIR}"
fi

