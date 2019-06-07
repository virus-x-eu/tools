#!/usr/bin/env python3

import argparse
import tarfile
from lxml import etree

parser = argparse.ArgumentParser(description='''
Look up all "run accessions" (SRR/ERR/DRR...) belonging to a "submission accession" (SRA/ERA/DRA...) 
using NCBI SRA Metadata XML files inside a tar archive.
''')
parser.add_argument("metadata_tar_gz", help='NCBI SRA metadata tar.gz file '
                                            '(NCBI_SRA_Metadata_Full_20XXXXXX.tar.gz)')
parser.add_argument("submission_accession", help='Run accession (starting with SRA/ERA/DRA)')
args = parser.parse_args()

with tarfile.open(args.metadata_tar_gz, 'r') as tar:
  for item in tar:
    if item.isfile and item.name == ('%s/%s.run.xml' % (args.submission_accession, args.submission_accession)):
      root = etree.XML(tar.extractfile(item).read())
      runs = root.findall('RUN')
      for run in runs:
        print(run.get('accession'))
      break
