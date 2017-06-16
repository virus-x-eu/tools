FROM ubuntu:16.04
RUN apt-get update && \
    apt-get install -y wget build-essential make perl && \
    cd /root/ && \
    wget "http://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz" -O- \
    | tar --no-same-owner -xz && \
    cpan install Moose  && \
    cpan install IPC::Run && \
    cd  PfamScan/ && \
    mv Bio/ /usr/share/perl/5.22/ && \
    mv pfam_scan.pl /usr/local/bin/ && \
    cd /root && \
    rm -rf PfamScan* && \
    apt-get purge -y build-essential make && \
    apt-get autoremove -y --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

