FROM ubuntu:rolling
ENV DASTOOL_VERSION=1.1.0
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -yq wget r-base r-cran-data.table r-cran-domc r-cran-ggplot2 ruby prodigal ncbi-blast+ && \
    cd /usr/local/bin/ && \
    wget --progress=dot:giga "http://github.com/bbuchfink/diamond/releases/download/v0.9.22/diamond-linux64.tar.gz" -O- | tar xz && \
    wget "https://github.com/bcthomas/pullseq/releases/download/1.0.2/pullseq_v1.0.2_linux64.zip" && \
    unzip pullseq_v1.0.2_linux64.zip && \
    rm pullseq_v1.0.2_linux64.zip && \
    ln -s /lib/x86_64-linux-gnu/libpcre.so.3 /usr/lib/libpcre.so.1 && \
    cd /opt/ && \
    wget --progress=dot:giga "https://github.com/cmks/DAS_Tool/archive/e46914f31f9ba036fcd73f33742dcb61d68632bd.tar.gz" -O- | tar xz && \
    mv "DAS_Tool-e46914f31f9ba036fcd73f33742dcb61d68632bd" /opt/dastool && \
    cd /opt/dastool && \
    R CMD INSTALL "./package/DASTool_${DASTOOL_VERSION}.tar.gz" && \
    unzip db.zip -d db && \
    chmod +x DAS_Tool && \
    ln -s /opt/dastool/DAS_Tool /usr/local/bin/ && \
    apt-get purge -y wget && \
    apt-get autoremove -y --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*