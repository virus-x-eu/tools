FROM alpine:latest
ENV BINARY=annotatedcontigs2json
COPY ${BINARY}.xz /usr/local/bin/
RUN apk --no-cache add xz && \
    xz -d /usr/local/bin/${BINARY}.xz && \
    apk --no-cache del xz && \
    chmod +x /usr/local/bin/${BINARY}