FROM java:openjdk-8-jre-alpine
ENV FASTQ_VERSION=0.11.5 PATH="/FastQC/:${PATH}"
RUN apk add --no-cache wget perl fontconfig && \
    wget "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQ_VERSION}.zip" -O fastqc.zip && \
    unzip fastqc.zip && \
    chmod +x /FastQC/fastqc && \
    rm fastqc.zip