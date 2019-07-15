FROM conda/miniconda3
ENV GTDBTK_VERSION=0.3.2
RUN conda update conda && \
    conda install python=2 && \
    conda install -y -c bioconda gtdbtk=${GTDBTK_VERSION}