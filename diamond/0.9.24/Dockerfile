FROM alpine:latest
ENV DIAMOND_VERSION="0.9.24"
RUN apk add --no-cache g++ make cmake zlib-dev && \
    wget "http://github.com/bbuchfink/diamond/archive/v${DIAMOND_VERSION}.tar.gz" && \
    tar xf v${DIAMOND_VERSION}.tar.gz && \
    cd diamond-${DIAMOND_VERSION}/ && \
    mkdir bin && \
    cd bin && \
    cmake .. && \
    make -j$(nproc) && \
    mv diamond /usr/local/bin/ && \
    apk del g++ make cmake zlib-dev && \
    rm -rf /var/cache/apk/*