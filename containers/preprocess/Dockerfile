# FROM python
FROM ubuntu:18.04
MAINTAINER Stephen Blaskowski

# Set source and bin directories
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

# Install dependencies
RUN apt-get update && \
    apt-get install -y default-jre git python3 python3-dev python3-pip unzip wget

# Install Trimmomatic
WORKDIR $SRC
RUN wget www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    echo '#!/bin/bash\njava -jar /usr/local/src/Trimmomatic-0.39/trimmomatic-0.39.jar $@' > $BIN/trimmomatic && \
    chmod +x $BIN/trimmomatic && \
    rm Trimmomatic-0.39.zip

# Install Fastqc
WORKDIR $SRC
RUN wget www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip && \
    unzip fastqc_v0.11.9.zip && \
    chmod +x FastQC/fastqc && \
    ln -s $SRC/FastQC/fastqc $BIN/ && \
    rm fastqc_v0.11.9.zip

# Install MultiQC
WORKDIR $SRC
RUN git clone https://github.com/ewels/MultiQC.git && \
    cd MultiQC && \
    pip3 install .

# Set default working directory
WORKDIR /home/

# Fix python3 error
ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
