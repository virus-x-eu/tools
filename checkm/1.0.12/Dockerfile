FROM ubuntu:18.04
ENV CHECKM_VERSION=1.0.12
COPY run_checkm.sh /usr/local/bin/
COPY config_location.patch /root/
RUN apt-get update && \
    apt-get install -y hmmer prodigal wget unzip python python-pip && \
    pip install checkm-genome==${CHECKM_VERSION} && \
    wget "https://github.com/matsen/pplacer/releases/download/v1.1.alpha19/pplacer-linux-v1.1.alpha19.zip" -O /root/pplacer.zip && \
    cd /root/ && \
    unzip pplacer.zip && \
    cp -r /root/pplacer-Linux-v1.1.alpha19/* /usr/local/bin/ && \
    rm -r /root/pplacer* && \
    patch --binary -i /root/config_location.patch /usr/local/lib/python2.7/dist-packages/checkm/checkmData.py && \
    chmod +x /usr/local/bin/run_checkm.sh && \
    apt-get purge -y wget unzip && \
    apt-get autoremove -y --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*