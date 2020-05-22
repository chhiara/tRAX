FROM continuumio/miniconda:latest

# C Compiler for DESeq2
RUN apt-get update && apt-get install -yq --no-install-recommends \ 
    build-essential

# Install Core TRAX Dependencies
RUN conda install --quiet --yes -c conda-forge -c bioconda \        
    'bowtie2=2.3.5' \
    'cutadapt=1.18' \
    'infernal=1.1.2' \
    'pysam=0.15.3' \
    'samtools=1.9' \
    'seqprep=1.3.2' \
    'sra-tools=2.8.0'

RUN conda install --quiet --yes -c conda-forge -c r -c bioconda \   
    'bioconductor-deseq2=1.26.0' \
    'r-ggplot2=3.3.0' \
    'r-getopt=1.20.3' \
    'r-ggrepel=0.8.2' \
    'r-reshape2=1.4.4'

# Set working directory and copy TRAX software into docker container
COPY . /opt/trax/
#RUN chmod -R 777 trax/
ENV PATH /opt/trax:$PATH

# Add user for the container
RUN useradd -ms /bin/bash jerry

# Add empty folder for RNA database docker volumes
RUN mkdir /rnadb &&\
     chmod -R 777 /rnadb &&\
     chmod -R 777 /home

USER jerry
WORKDIR /home/jerry