FROM ubuntu:21.10
ENV DEBIAN_FRONTEND=noninteractive
COPY ena-fastq-download.py /usr/local/bin/
RUN apt-get update && \
    apt-get install -y python3-pip wget sra-toolkit && \
    pip3 install pandas && \
    chmod +x /usr/local/bin/ena-fastq-download.py && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

