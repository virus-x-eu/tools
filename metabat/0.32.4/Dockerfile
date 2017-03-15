FROM alpine
RUN apk add --update wget ca-certificates bash && \
    wget "https://bitbucket.org/berkeleylab/metabat/downloads/metabat-static-binary-linux-x64_v0.32.4.tar.gz" -O- | tar --no-same-owner -xz && \
    mv metabat/* /usr/local/bin/ && \
    rmdir metabat && \
    apk del wget ca-certificates && \
    rm -rf /var/cache/apk/*
