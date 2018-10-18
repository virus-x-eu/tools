#!/usr/bin/python3

import argparse
import tarfile
from lxml import etree
import csv

parser = argparse.ArgumentParser(description='''
Extract specific fields from the NCBI SRA Metadata XML files inside a tar archive.
''')
parser.add_argument("metadata_tar_gz", help='NCBI SRA metadata tar.gz file '
                                            '(NCBI_SRA_Metadata_Full_20XXXXXX.tar.gz)')
parser.add_argument("output_tsv", help='tsv output file the selected fields will be written to')
parser.add_argument("tag", nargs='+', help='The fields to retrieve from the XML files, '
                                           'e.g. "latitude"')
args = parser.parse_args()

total_sample_count = 0

with tarfile.open(args.metadata_tar_gz, 'r') as tar, open(args.output_tsv, 'w') as tsv:
  tsv_out = csv.writer(tsv, delimiter='\t')
  tsv_out.writerow(['sample_accession'] + args.tag)
  for item in tar:
    if item.isfile and item.name.endswith('.sample.xml'):
      root = etree.XML(tar.extractfile(item).read())
      samples = root.findall('SAMPLE')
      for sample in samples:
        total_sample_count += 1
        desired_tag_values = {key: '' for key in args.tag}
        sample_attrs = sample.findall('SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')
        for sample_attr in sample_attrs:
          for tag in args.tag:
            if sample_attr is not None and sample_attr.findtext('TAG') == tag:
              if total_sample_count % 1000 == 0:
                print('Samples processed:  %s' % total_sample_count)
              desired_tag_values[tag] = sample_attr.findtext('VALUE')
        if any([False if desired_tag_values[tag] == '' else True for tag in args.tag]):
          tsv_out.writerow([sample.get('accession')] + [desired_tag_values[tag] for tag in args.tag])
  print('Total amount of samples processed:  %s' % total_sample_count)
