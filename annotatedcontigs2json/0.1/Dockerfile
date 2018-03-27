FROM ubuntu:rolling
COPY annotatedcontigs2json.py /usr/local/bin/
RUN apt-get update && \
    apt-get install -y python-biopython && \
    chmod +x /usr/local/bin/annotatedcontigs2json.py && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
