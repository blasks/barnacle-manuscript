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
BASEDIR=$(realpath "../../")
REFS="${BASEDIR}/data/refseqs/pro-syn-virus-refseqs.genes.fna.gz"
FASTQ="${BASEDIR}/data/fastq"
CONTAINERDIR="${BASEDIR}/containers"
THREADS=48

# data output
MAPDIR="${BASEDIR}/data/2-mapping/"
MAPPINGS="${MAPDIR}/collated_salmon_data.csv.gz"

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

##################################
# 1. Compile reference sequences #
##################################
printf "\nStep 3: Compiling reference sequences\n"
if [[ ! -e ${REFS} ]]; then
    pushd $(dirname ${REFS})
    ./compile-refseqs.sh
    popd
fi

#########################################
# 2. Map fastq reads against references #
#########################################
printf "\nStep 4: Mapping processed fastq reads against reference set\n"

# build salmon container
printf "\n\t* Checking ${CONTAINER} container\n"
if [ "${CONTAINER}" == "singularity" ]; then
    if [[ ! -e "${CONTAINERDIR}/salmon.sif" ]]; then 
        singularity build "${CONTAINERDIR}/salmon.sif" docker://combinelab/salmon:1.10.2
    fi
elif [ "${CONTAINER}" == "docker" ]; then
    docker pull combinelab/salmon:1.10.2
fi

# if necessary, build salmon index files and store in reference directory
IDXDIR="$(dirname ${REFS})/salmon_index"
# build index (keep duplicates to properly divide reads across identical refs)
if [[ ! -d ${IDXDIR} ]]; then
    printf "\n\t* Building salmon index\n"
    if [ "${CONTAINER}" == "singularity" ]; then
        singularity exec --bind ${BASEDIR}:${BASEDIR} ${CONTAINERDIR}/salmon.sif salmon index \
            -t ${REFS} -i ${IDXDIR} -k 31 -p ${THREADS} --keepDuplicates
    elif [ "${CONTAINER}" == "docker" ]; then
        docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} \
            combinelab/salmon:1.10.2 salmon index \
            -t ${REFS} -i ${IDXDIR} -k 31 -p ${THREADS} --keepDuplicates
    fi
else
    printf "\n\t* Using salmon index found at ${IDXDIR}\n"
fi

# map data to reference
i=1
TOTAL=$( ls ${FASTQ}/qc/*.fw.fastq.gz | wc -l )
for R1 in ${FASTQ}/qc/*.fw.fastq.gz; do 
    SAMPLE=$(basename ${R1})
    SAMPLE=${SAMPLE%.fw.fastq.gz}
    R2="${FASTQ}/qc/${SAMPLE}.rv.fastq.gz"
    OUTDIR="${MAPDIR}/salmon/${SAMPLE}"
    # quanitify with salmon
    if [[ ! -f ${OUTDIR}/quant.sf ]]; then
        mkdir -p ${OUTDIR}
        printf "\n\t* Mapping sample ${i}/${TOTAL}: ${SAMPLE}\n\n"
        if [ "${CONTAINER}" == "singularity" ]; then
            singularity exec --bind ${BASEDIR}:${BASEDIR} ${CONTAINERDIR}/salmon.sif salmon quant \
                -i ${IDXDIR} -l A -1 ${R1} -2 ${R2} -o ${OUTDIR} --validateMappings
        elif [ "${CONTAINER}" == "docker" ]; then
            docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} \
                combinelab/salmon:1.10.2 salmon quant \
                -i ${IDXDIR} -l A -1 ${R1} -2 ${R2} -o ${OUTDIR} --validateMappings
        fi
    fi
    let i=$i+1
done

#######################
# 3. Collate mappings #
#######################
printf "\nStep 5: Collating mapping outputs\n"
# build barnacle container
printf "\n\t* Checking ${CONTAINER} container\n"
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
printf "\n\t* Collating mapping files\n"
if [[ ! -e ${MAPPINGS} ]]; then
    if [ "${CONTAINER}" == "singularity" ]; then
        singularity exec --bind ${BASEDIR}:${BASEDIR} ${CONTAINERDIR}/barnacle.sif python \
            "${BASEDIR}/analyses/2-mapping/collate_salmon_files.py" "${MAPDIR}/salmon/" "${MAPDIR}"
    elif [ "${CONTAINER}" == "docker" ]; then
        docker run --mount type=bind,source=${BASEDIR},target=${BASEDIR} barnacle python \
            "${BASEDIR}/analyses/2-mapping/collate_salmon_files.py" "${MAPDIR}/salmon/" "${MAPDIR}"
    fi
fi
