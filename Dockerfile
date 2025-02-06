################## BASE IMAGE #####################
FROM ubuntu:22.04

################## METADATA #######################
LABEL base_image="ubuntu:22.04"
LABEL software="RNAseq"
LABEL software.version="1.0"
LABEL about.summary="Image contenant toutes les dépendances pour le pipeline RNAseq"
LABEL about.home="http://github.com/Lipinski-B/RNAseq"
LABEL about.documentation="http://github.com/Lipinski-B/RNAseq/README.md"
LABEL about.license="GNU-3.0"

################## INSTALLATION ######################
RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y --no-install-recommends \
    curl \
    bzip2 \
    ca-certificates \
    procps \
    cmake \
    python3-pip \
    python3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y --no-install-recommends wget

################## MINICONDA ######################
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
    bash /miniconda.sh -b -p /opt/conda && \
    rm /miniconda.sh && \
    /opt/conda/bin/conda clean -a
RUN chmod +x /opt/conda/bin/conda
ENV PATH="/opt/conda/bin:$PATH"

################## KALLISTO ######################
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.51.1/kallisto_linux-v0.51.1.tar.gz && \
    tar -xvf kallisto_linux-v0.51.1.tar.gz -C /opt/
RUN chmod +x /opt/kallisto/kallisto
ENV PATH="/opt/kallisto:$PATH"

################## ENVIRONNNEMENT ######################
COPY environnement.yml /
RUN conda env create -n RNAseq -f /environnement.yml && conda clean -a
RUN pip install cget numpy
ENV PATH="/opt/conda/envs/RNAseq/bin:$PATH"
SHELL ["conda", "run", "-n", "RNAseq", "/bin/bash", "-c"]
