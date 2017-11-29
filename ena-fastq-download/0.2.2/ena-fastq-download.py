#!/usr/bin/env python3

# (c) Christian Henke
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0

import pandas as pd
import numpy
import subprocess
import argparse

parser = argparse.ArgumentParser(description='''
Download fastq files from the European Nucleotide Archive (ENA) to the current directory.
''')
parser.add_argument("--run-accession", dest='run_accession', required=True, help='The run accession ID, e.g. ERR0123456')
parser.add_argument("--size-limit", dest='size_limit', type=int, default=0, required=False, help='Files are only downloaded if the SRA file size (not FASTQ) in mebibytes is below this limit.')
args = parser.parse_args()

raw_url_template = 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run'
data = pd.read_csv(raw_url_template % args.run_accession, sep='\t')
fastq_urls = data['fastq_ftp'][0]
wget = ['/usr/bin/wget', '-nc', '--progress=dot:giga']


if args.size_limit == 0 or numpy.isnan(data['sra_bytes'][0]) or data['sra_bytes'][0] < args.size_limit*1024*1024:
  # check if tsv field was not empty (pandas converts empty fields to NaN which is != str)
  # we also require paired ends in separate files (_1 and _2)
  if type(fastq_urls) is str and '_1.fastq.gz' in fastq_urls and '_2.fastq.gz' in fastq_urls:
    print('FastQ file URLs found. Downloading them directly ...')
    for fastq in fastq_urls.split(';'):
      url = 'http://%s' % fastq
      print(url)
      subprocess.run(wget + [url])
  else:  # if fastq is not directly available or at least not in separate files. downloads sra format and converts it to fastq
    print('No FastQ file URLs found!')
    print('Fallback: Downloading SRA file ...')
    sra_url = data['sra_ftp'][0]
    if type(sra_url) is str:
      url = 'http://%s' % sra_url
      print(url)
      subprocess.run(wget + [url])
      print('Fallback: Converting with fastq-dump ...')
      subprocess.run(['/usr/bin/fastq-dump', '-v', '--readids', '--gzip', '--minReadLen', '15', '--split-3', './%s' % args.run_accession])
    else:
      print('Fallback: Fatal error: No SRA URL found!')
else:
  print('Data size is above the specified size limit. Nothing will be downloaded.')