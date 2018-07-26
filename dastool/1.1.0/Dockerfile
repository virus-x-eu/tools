FROM ubuntu:rolling
ENV DASTOOL_VERSION=1.1.0
RUN apt-get update && \
    apt-get install -y wget r-base r-cran-data.table r-cran-domc r-cran-ggplot2 ruby prodigal ncbi-blast+ && \
    cd /usr/local/bin/ && \
    wget --progress=dot:giga "http://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz" -O- | tar xz && \
    cd /opt/ && \
    wget --progress=dot:giga "https://github.com/cmks/DAS_Tool/archive/${DASTOOL_VERSION}.tar.gz" -O- | tar xz && \
    mv "DAS_Tool-${DASTOOL_VERSION}" /opt/dastool && \
    cd /opt/dastool && \
    R CMD INSTALL "./package/DASTool_${DASTOOL_VERSION}.tar.gz" && \
    unzip db.zip -d db && \
    chmod +x DAS_Tool && \
    ln -s /opt/dastool/DAS_Tool /usr/local/bin/ && \
    apt-get purge -y wget && \
    apt-get autoremove -y --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*