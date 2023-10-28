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
REFS="../../data/refseqs/pro-syn-virus-refseqs.genes.fna.gz"
FASTQ="../../data/fastq"
CONTAINER_DIR="../../containers"

# data output
MAPPINGS=../../data/2-mapping/collated_salmon_files.py

####################################################
# 0. Select Singularity or Docker containerization #
####################################################
CONTAINER=""
while [ "${CONTAINER}" == "" ]; do
    printf "\nPlease select an option:\n\n\t1 - singularity\n\t2 - docker\n\nEnter '1' or '2': "
    read selection 
    if [ "${selection}" == "1" ]; then
        echo "Continuing with Singularity containerized workflow"
        CONTAINER="singularity"
    elif [ "${selection}" == "2" ]; then
        echo "Continuing with Docker containerized workflow"
        CONTAINER="docker"
    else 
        printf "\nInvalid selection\n"
    fi
done

############################
# 1. Trim & QC fastq reads #
############################
printf "\nStep 1: Preprocessing raw fastq reads\n"
if [[ ! -e ${FASTQ} ]]; then
    pushd ${FASTQ}
    ./process-raw-reads.sh
    popd
fi

#############################
# 2. Combine size fractions #
#############################
printf "\nStep 2: Combining size fractions\n"
if [[ ! -e ${FASTQ} ]]; then
    pushd ${FASTQ}
    ./combine-size-fractions.sh
    popd
fi

##################################
# 3. Compile reference sequences #
##################################
printf "\nStep 3: Compiling reference sequences\n"
if [[ ! -e ${REFS} ]]; then
    pushd $(dirname ${REFS})
    ./compile-refseqs.sh
    popd
fi

#########################################
# 4. Map fastq reads against references #
#########################################
printf "\nStep 4: Mapping processed fastq reads against reference set\n"
# # build singularity image from docker image
# if [[ ! -e salmon_latest.sif ]]; then 
#     singularity build salmon_latest.sif docker://combinelab/salmon
# fi

# # if necessary, build salmon index files and store in reference directory
# for ref in $REFS; do
#     ref=$(realpath ${ref})
#     idxdir="$(dirname ${ref})/salmon_index"
#     if [[ ! -d ${idxdir} ]]
#     then
#         # build index (keep duplicates to properly divide reads across identical txs)
#         singularity exec --bind ${MOUNT} salmon_latest.sif salmon index \
#             -t ${ref} -i ${idxdir} -k 31 -p ${THREADS} --keepDuplicates
#     fi
# done

# # map gradients data to combined reference
# for ref in $REFS; do
#     idxdir="$(dirname ${ref})/salmon_index"
#     ref_set=$(basename ${ref})
#     ref_set=${ref_set%.genes.renamed.fna.gz}
#     mkdir -p data/${ref_set}
#     i=0
#     for dataset in `realpath ${TX_READS}`; do 
#         for r1 in ${dataset}/*.fw.fastq.gz
#         do 
#             sample=$(basename ${r1})
#             sample=${sample%.fw.fastq.gz}
#             r2="${dataset}/${sample}.rv.fastq.gz"
#             outdir=data/${ref_set}/${sample}
#             # quanitify with salmon
#             if [[ ! -f ${outdir}/quant.sf ]]; then
#                 mkdir -p ${outdir}
#                 printf "\nmapping sample #${i}: ${sample} vs reference ${ref_set}\n"
#                 singularity exec --bind ${MOUNT} \
#                     salmon_latest.sif salmon quant \
#                     -i ${idxdir} -l A \
#                     -1 ${r1} -2 ${r2} -o ${outdir} \
#                     --dumpEq --validateMappings
#             fi
#             let i=$i+1
#         done
#     done
#     # tar all mappings
#     tar -czvf data/${ref_set}.tar.gz data/${ref_set}
# done

#######################
# 5. Collate mappings #
#######################
printf "\nStep 5: Collating mapping outputs\n"
# outdir=${mappings}/collated/$(basename ${mapping})
# if [[ ! -e $outdir/collated_TPM_data.csv.gz ]]; then
#     mkdir -p ${outdir}
#     ../../bin/collate_salmon_files.py ${mapping} ${outdir}
# fi
