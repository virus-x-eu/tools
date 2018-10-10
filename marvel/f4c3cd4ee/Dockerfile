FROM conda/miniconda3
ENV MARVEL_VERSION=f4c3cd4ee5abcaf5b9e2756df4ec5cae23b54e1a
RUN apt-get update && \
    apt-get install -y wget less && \
    cd /usr/local/share/ && \
    wget -e dotbytes=1M -O- "https://github.com/LaboratorioBioinformatica/MARVEL/archive/${MARVEL_VERSION}.tar.gz" | tar xz && \
    mv "MARVEL-${MARVEL_VERSION}" marvel && \
    python3.6 -m pip install numpy scipy biopython scikit-learn && \
    conda install -y -c conda-forge -c bioconda prokka && \
    conda install -y -c biocore hmmer && \
    conda clean --all -y && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*