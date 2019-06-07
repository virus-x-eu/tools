FROM python:latest
COPY sra-run-lookup.py /usr/local/bin/
RUN  pip3 --no-cache-dir --disable-pip-version-check install lxml && \
     chmod +x /usr/local/bin/sra-run-lookup.py
