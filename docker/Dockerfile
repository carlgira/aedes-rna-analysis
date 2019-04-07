FROM ubuntu:18.04

# docker build . -t carlgira/rna-analysis:latest

#######
## OS #
#######

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV TZ=Europe/Minsk
ENV PATH /opt/conda/bin:/opt/STAR/source:$PATH
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git r-base && \
    apt-get update && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

##############
## MiniConda #
##############

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.11-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

##################
## RNA-seq tools #
##################

RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda

RUN conda install -y fastqc=0.11.8
RUN conda install -y hisat2=2.1.0
RUN conda install -y samtools=1.9
RUN conda install -y subread=1.6.4
RUN conda install -y picard=2.19.0
RUN conda install -y multiqc=1.7
RUN conda install -y star=2.7.0e
RUN conda install -y pandoc=2.7.1
RUN conda install -y bioconductor-edger=3.24.1
RUN conda install -y bioconductor-biomart=2.38.0
RUN conda install -y bioconductor-noiseq=2.26.0
RUN conda install -y bioconductor-glimma=1.10.0
RUN conda install -y r-rmarkdown=1.12
RUN conda install -y r-gplots=3.0.1.1
RUN conda install -y r-ggplot2=3.1.0
