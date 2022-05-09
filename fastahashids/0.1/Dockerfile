FROM python:alpine
COPY regenerate_ids.py regenerate_contig_ids.py update_ids.py /usr/local/bin/
RUN apk --no-cache add build-base && \
    pip install biopython base32-crockford && \
    chmod +x /usr/local/bin/regenerate_ids.py && \
    chmod +x /usr/local/bin/regenerate_contig_ids.py && \
    chmod +x /usr/local/bin/update_ids.py && \
    apk del build-base


