FROM ubuntu:16.04 as builder
ENV METABAT_VERSION=86b9406
RUN apt-get update && \
    apt-get install -y build-essential libboost-all-dev scons git curl libncurses5-dev zlib1g-dev && \
    cd /root/ && \
    git clone https://bitbucket.org/berkeleylab/metabat.git metabat && \
    cd metabat && \
    git checkout ${METABAT_VERSION} && \
    mkdir install && \
    scons install PREFIX=./install

FROM alpine
RUN apk add --no-cache bash
COPY --from=builder /root/metabat/install/ /usr/
