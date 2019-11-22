FROM continuumio/miniconda:latest

# C Compiler for DESeq2
RUN apt-get update && apt-get install -yq --no-install-recommends build-essential

# Install Core TRAX Dependencies
RUN conda install --quiet --yes -c r -c bioconda \
    'bowtie2' \
    'cutadapt' \
    'infernal' \
    'pysam' \
    'samtools' \
    'sra-tools'

# Install R Packages
RUN conda install --quiet --yes -c r \
    'r-ggplot2' \
    'r-getopt'

# Install Bioconductor Manager and deSeq2
ADD install_biocmanager_deseq2.R /tmp/
RUN R -f /tmp/install_biocmanager_deseq2.R

# Set working directory and copy TRAX software into docker container
ADD * /opt/trax/
#RUN chmod -R 777 trax/
ENV PATH /opt/trax:$PATH

# Add user for the container
RUN useradd -ms /bin/bash jerry

# Add empty folder for RNA database docker volumes
RUN mkdir /rnadb &&\
     chmod -R 755 /rnadb &&\
     chmod -R 777 /home

USER jerry
WORKDIR /home/jerry

