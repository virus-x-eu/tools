FROM ubuntu:rolling
RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev && \
    git clone https://github.com/OpenGene/fastp && \
    cd fastp && \
    make && \
    make install && \
    apt-get remove --purge git build-essential && \
    apt-get autoremove --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*