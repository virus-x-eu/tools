FROM ubuntu:16.04
COPY megan_taxon2lineage.pl /usr/local/bin/
RUN apt-get update && \
    apt-get install -y wget build-essential make perl && \
    cpan -fi Moose  && \
    cpan -fi IPC::Run && \
    cpan -fi Bio::LITE::Taxonomy::NCBI && \
    cpan -fi Bio::SearchIO && \
    chmod +x /usr/local/bin/megan_taxon2lineage.pl && \
    apt-get purge -y build-essential make && \
    apt-get autoremove -y --purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

