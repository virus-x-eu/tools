#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
  print('Usage: \n      %s SRA_ACCESSION [PREFIX]' % sys.argv[0])
  exit(1)

acc = sys.argv[1]
prefix = sys.argv[2] if len(sys.argv) > 2 else 'http://'

urls = [
  '%sftp.sra.ebi.ac.uk/vol1/fastq/%s/%s/%s/%s%s.fastq.gz' % (prefix, acc[:6], '00' + acc[-1:], acc, acc, suf)
  for suf in ['', '_1', '_2']
]

for url in urls:
  print(url)
