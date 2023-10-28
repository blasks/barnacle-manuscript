#!/bin/bash
# AUTHOR: Stephen Blaskowski

# The purpose of this script is to compile a single fasta containing
# all the Prochlorococcus, Synechococcus, and viral gene reference sequences 
# against which the metatranscriptomic data will be mapped. The script also
# renames reference gene sequences to match the nomenclature of 
# CyCOGs v6 (Berube et al., 2018).

# exit when your script tries to use undeclared variables
set -o nounset
# exit if any pipe commands fail
set -o pipefail
# exit when a command fails
set -o errexit
# # trace what gets executed
# set -o xtrace

# data inputs
GENOME_METADATA="../metadata/genome-metadata.csv"
GENOMES="genomes"
CYCOGS="../metadata/cycogs.tsv"

# data output
REFS=pro-syn-virus-refseqs.genes.fna.gz

# check for inputs
for input in ${GENOME_METADATA} ${GENOMES} ${CYCOGS}; do
    if [[ ! -e ${input} ]]; then
        echo "Missing data input: ${input}"
        echo "Please run data/assemble-data-inputs.ipynb to interactively assemble data inputs"
        exit
    fi
done

# standardize reference gene names to match CyCOGs
printf "\n*** Standardizing reference gene names to match CyCOGs ***\n"
while read id name; do
    newdir=${GENOMES}/renamed/${id}
    if [[ ! -d ${newdir} ]]; then
        echo "standardizing gene names for ${name}"
        mkdir -p ${newdir}
        # regex: replace all lines starting with ">" with the genome name 
        # and the gene id found in the second capture group
        sed -E "s/(^>)([^ ]*)(.*$)/>${name}_\2/g" \
            ${GENOMES}/${id}/${id}.genes.fna \
            >> ${newdir}/${id}.genes.renamed.fna
    else
        echo "${name} is all set"
    fi
done < <(tail -n +2 ${GENOME_METADATA} | cut -d , -f 1,2 | tr , ' ')

# concatenate reference gene sequences
printf "\n*** Concatenating reference gene sequences ***\n"
if [[ ! -e ${REFS} ]]; then
    mkdir -p $(dirname ${REFS})
    for id in `tail -n +2 ${GENOME_METADATA} | cut -d , -f 1`; do
        cat ${GENOMES}/renamed/${id}/${id}.genes.renamed.fna >> ${REFS%.gz}
    done
    # clean up renamed reference directory
    rm -rf ${GENOMES}/renamed
fi

# gzip reference sequences
printf "\n*** Compressing fasta of compiled reference sequences ***\n"
gzip ${REFS%.gz}
echo ${REFS}

# check that all CyCOG genes are represented in reference genes
printf "\n*** Checking that all CyCOG genes are represented in reference gene set ***\n"
# collect reference genes
zgrep "^>" ${REFS} | tr -d ">" >> ref-genes.txt 
# collect cycog genes
for gene in `tail -n +2 ${CYCOGS} | cut -f 9 | tr , ' '`; do 
    echo ${gene} >> cycog-genes.txt
done
# find the intersection of reference and cycog gene sets
cat ref-genes.txt cycog-genes.txt | sort | uniq -d >> intersection-genes.txt
# compare the intersection back to cycog genes
cat intersection-genes.txt cycog-genes.txt | sort | uniq -u >> missing-genes.txt
# clean up scratch files
rm ref-genes.txt cycog-genes.txt intersection-genes.txt
# check missing genes
if [[ -s missing-genes.txt ]]; then
    echo "The following CyCOG genes are missing from the reference sequences:"
    cat missing-genes.txt
    exit
else
    echo "All CyCOG genes present and accounted for"
    rm missing-genes.txt
fi

