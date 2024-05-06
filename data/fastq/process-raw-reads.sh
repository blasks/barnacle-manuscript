#!/bin/bash
# AUTHOR: Stephen Blaskowski

# Process raw reads downloaded from SRA by combining all reads originating
# from the same sample, and then trimming all reads using Trimmomatic

# Inputs:
BASEDIR=$(realpath "../../")
FASTQDIR="${BASEDIR}/data/fastq"
RAWDIR="${FASTQDIR}/raw"
DEFRACDIR="${FASTQDIR}/defrac"
QCDIR="${FASTQDIR}/qc"
CONTAINERDIR="${BASEDIR}/containers"
METADATA="SRA-Sample-Metadata.csv"
ADAPTERS="TruSeq2-PE.fa"
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
    printf "\nWould you like to save storage by removing raw fastq files following processing?\n\n\t1 - Yes\n\t2 - No\n\nEnter '1' or '2': "
    read SELECTION 
    if [ "${SELECTION}" == "1" ]; then
        echo "Will remove raw fastq files following processing"
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
    if [[ ! -e "${CONTAINERDIR}/trimmomatic.sif" ]]; then 
        singularity build "${CONTAINERDIR}/trimmomatic.sif" docker://quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4 
    fi
elif [ "${CONTAINER}" == "docker" ]; then
    docker pull quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4 
fi

##############################################
# 2. Combine all reads from the same samples #
##############################################
if [[ ! -d ${DEFRACDIR} ]]; then
    mkdir -p ${DEFRACDIR}
fi
for SAMPLE in `tail -n +2 $METADATA | cut -d "," -f 2 | sort | uniq`; do
    printf "\n\t* Collecting reads for sample ${SAMPLE}\n"
    # pull out all the SRRs associated with each sample ID
    grep "${SAMPLE}" "${METADATA}" >> temp-srr-list.csv
    for SRR in `cut -d "," -f 3 temp-srr-list.csv`; do
        # concatenate forward reads
        cat "${RAWDIR}/${SRR}_1.fastq.gz" >> "${DEFRACDIR}/${SAMPLE}.fw.fastq.gz"
        # concatenate reverse reads
        cat "${RAWDIR}/${SRR}_2.fastq.gz" >> "${DEFRACDIR}/${SAMPLE}.rv.fastq.gz"
        # clean up
        if [ "${CLEANUP}" == "True" ]; then
            rm "${RAWDIR}/${SRA}_?.fastq.gz"
        fi
    done
    rm temp-srr-list.csv
done
# clean up
if [ "${CLEANUP}" == "True" ]; then
    rmdir "${RAWDIR}"
fi

##################################
# 3. Trim files with Trimmomatic #
##################################
if [[ ! -d ${QCDIR}/logs ]]; then
    mkdir -p ${QCDIR}/logs
fi
for SAMPLE in `tail -n +2 $METADATA | cut -d "," -f 2 | sort | uniq`; do
    printf "\n\t* Trimming reads for sample ${SAMPLE}\n"
    if [ "${CONTAINER}" == "singularity" ]; then
        singularity exec --bind ${BASEDIR}:${BASEDIR} ${CONTAINERDIR}/trimmomatic.sif \
            trimmomatic PE -threads ${THREADS} \
            "${DEFRACDIR}/${SAMPLE}.fw.fastq.gz" \
            "${DEFRACDIR}/${SAMPLE}.rv.fastq.gz" \
            "${QCDIR}/${SAMPLE}.fw.fastq.gz" \
            "${QCDIR}/${SAMPLE}.fw.unpaired.fastq.gz" \
            "${QCDIR}/${SAMPLE}.rv.fastq.gz" \
            "${QCDIR}/${SAMPLE}.rv.unpaired.fastq.gz" \
            ILLUMINACLIP:"${ADAPTERS}":2:30:10:1:true \
            MAXINFO:135:0.5 LEADING:3 TRAILING:3 MINLEN:60 AVGQUAL:20 >> "${QCDIR}/logs/${sample}.trimmomatic.log"
    elif [ "${CONTAINER}" == "docker" ]; then
        docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} \
            quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4 \
            trimmomatic PE -threads ${THREADS} \
            "${DEFRACDIR}/${SAMPLE}.fw.fastq.gz" \
            "${DEFRACDIR}/${SAMPLE}.rv.fastq.gz" \
            "${QCDIR}/${SAMPLE}.fw.fastq.gz" \
            "${QCDIR}/${SAMPLE}.fw.unpaired.fastq.gz" \
            "${QCDIR}/${SAMPLE}.rv.fastq.gz" \
            "${QCDIR}/${SAMPLE}.rv.unpaired.fastq.gz" \
            ILLUMINACLIP:"${ADAPTERS}":2:30:10:1:true \
            MAXINFO:135:0.5 LEADING:3 TRAILING:3 MINLEN:60 AVGQUAL:20 >> "${QCDIR}/logs/${sample}.trimmomatic.log"
    fi
    if [ "${CLEANUP}" == "True" ]; then
        rm "${DEFRACDIR}/${SAMPLE}.??.fastq.gz"
    fi
done
# clean up
if [ "${CLEANUP}" == "True" ]; then
    rmdir "${DEFRACDIR}"
fi
