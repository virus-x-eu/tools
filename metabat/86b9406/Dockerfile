FROM ubuntu:16.04
ENV METABAT_VERSION=86b9406
RUN apt-get update && \
    apt-get install -y build-essential libboost-all-dev scons git curl libncurses5-dev zlib1g-dev && \
    cd /root/ && \
    git clone https://bitbucket.org/berkeleylab/metabat.git metabat && \
    cd metabat && \
    git checkout ${METABAT_VERSION} && \
    mkdir install && \
    scons install PREFIX=./install && \
    cp -r /root/metabat/install/* /usr/local/ && \
    cd /root && \
    rm -rf /root/metabat && \
    apt-get purge -y build-essential libboost-all-dev scons git curl libncurses5-dev zlib1g-dev && \
    apt-get autoremove -y --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
