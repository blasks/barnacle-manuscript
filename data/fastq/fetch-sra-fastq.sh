#!/bin/bash
# AUTHOR: Stephen Blaskowski

# Fetch raw fastq data from the NCBI Sequence Read Archive (SRA)
# This script depends on the prefetch and fasterq-dump functions of the sra-tools
# software. Instructions for downloading/installing sra-tools can be found
# on the NCBI website, or the project repository page:
# https://github.com/ncbi/sra-tools/wiki

# Inputs:
BASEDIR=$(realpath "../../")
FASTQDIR="${BASEDIR}/data/fastq"
SRADIR="${FASTQDIR}/sra"
CONTAINERDIR="${BASEDIR}/containers"
METADATA="SRA-Sample-Metadata.csv"
THREADS=32

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

###############################
# 1. Download fastqs from SRA #
###############################
if [[ ! -d ${SRADIR} ]]; then
    mkdir -p ${SRADIR}
fi
for SRA in `tail -n +2 $METADATA | cut -d "," -f 3`; do
    printf "\n\t* Downloading ${SRA}\n"
    if [ "${CONTAINER}" == "singularity" ]; then
        # prefetch SRA
        singularity exec --bind ${BASEDIR}:${BASEDIR} ${CONTAINERDIR}/sra-tools.sif prefetch \
            ${SRA} -O ${SRADIR} --max-size 1t
        # extract fastq files from SRA
        singularity exec --bind ${BASEDIR}:${BASEDIR} --workdir ${SRADIR} ${CONTAINERDIR}/sra-tools.sif fasterq-dump --threads ${THREADS} ${SRA}
    elif [ "${CONTAINER}" == "docker" ]; then
        # prefetch SRA
        docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} \
            quay.io/biocontainers/sra-tools:3.1.0--h4304569_1 prefetch \
            ${SRA} -O ${SRADIR} --max-size 1t
        # extract fastq files from SRA
        docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} --workdir ${SRADIR} \
            quay.io/biocontainers/sra-tools:3.1.0--h4304569_1 fasterq-dump --threads ${THREADS} ${SRA}
    fi
done

##############################################
# 2. Combine all reads from the same samples #
##############################################
# Iterate through each unique SampleID
# Pull out all SRR IDs associated with the SampleID
# Concatenate forward reads
# Concatenate reverse reads
