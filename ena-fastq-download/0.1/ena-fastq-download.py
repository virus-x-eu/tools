#!/usr/bin/env python3

# (c) Christian Henke
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0

import pandas as pd
import subprocess
import argparse

parser = argparse.ArgumentParser(description='''
Download fastq files from the European Nucleotide Archive (ENA) to the current directory.
''')
parser.add_argument("--run-accession", dest='run_accession', required=True, help='The run accession ID, e.g. ERR0123456')
args = parser.parse_args()

raw_url_template = 'https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=%s&result=read_run'
data = pd.read_csv(raw_url_template % args.run_accession, sep='\t')
for fastq in data['fastq_ftp'][0].split(';'):
  url = 'http://%s' % fastq
  print(url)
  subprocess.run(['/usr/bin/wget', url])