FROM debian:stretch

# docker build . -t carlgira/rna-analysis:latest

#######
## OS #
#######

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 && \
    apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/debian stretch-cran35/' && \
    apt-get update && \
    apt-get -y install r-base && \
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
    conda config --add channels bioconda && \
    conda install -y fastqc=0.11.8 && \
    conda install -y star=2.7.0e && \
    conda install -y hisat2=2.1.0 && \
    conda install -y samtools=1.9 && \
    conda install -y subread=1.6.4 && \
    conda install -y pandas=0.23.4


###############
## R packages #
###############

RUN echo "source('http://bioconductor.org/biocLite.R') \n biocLite('edgeR')" | R --no-save
#RUN echo "install.packages('ggplot2', dep = TRUE)" | R --no-save
