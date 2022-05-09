FROM docker.io/continuumio/miniconda3
ENV GTDBTK_VERSION=2.0.0
RUN conda update conda && \
    conda install python=3 && \
    conda install -y -c bioconda gtdbtk=${GTDBTK_VERSION}