FROM continuumio/miniconda:latest

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
ADD install_biocmanager.R /tmp/
RUN R -f /tmp/install_biocmanager.R
RUN Rscript -e 'BiocManager::install(c("DESeq2"))'

# Set working directory and copy TRAX software into docker container
COPY * /trax/
RUN chmod -R 777 /trax/
WORKDIR /trax/

# Add empty folder for RNA database docker volumes
RUN mkdir /rnadb && chmod -R 777 /rnadb
