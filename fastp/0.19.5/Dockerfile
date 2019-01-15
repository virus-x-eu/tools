FROM ubuntu:rolling
ENV FASTP_VERSION=0.19.5
RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev && \
    git clone https://github.com/OpenGene/fastp && \
    cd fastp && \
    git checkout v${FASTP_VERSION} && \
    make && \
    make install && \
    apt-get remove -y --purge git build-essential && \
    apt-get autoremove -y --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*