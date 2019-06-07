#!/usr/bin/env python3

import argparse
import tarfile
from lxml import etree
import re
import json
import gzip

parser = argparse.ArgumentParser(description='''
Generate Elasticsearch JSON using NCBI SRA Metadata XML files inside a tar archive.
''')
parser.add_argument("metadata_tar_gz", help='NCBI SRA metadata tar.gz file '
                                            '(NCBI_SRA_Metadata_Full_20XXXXXX.tar.gz)')
parser.add_argument("json_gz", help='(gzipped) JSON output file')
args = parser.parse_args()

submissions = {}


def write_output():
  with gzip.open(args.json_gz, 'wt') as gz:
    for key in submissions:
      header = {'index': {'_id': key}}
      gz.write(json.dumps(header))
      gz.write('\n')
      gz.write(json.dumps(submissions[key], separators=(',', ':')))
      gz.write('\n')


def clean(val):
  return val.lower().replace(' ', '_')


def parse_first_sample(doc):
  # use only first sample because their attributes are usually almost identical
  sample_xml = doc.find('SAMPLE')
  return {
    'title': sample_xml.findtext('TITLE') or '',
    'id': sample_xml.findall('IDENTIFIERS')[0].findtext('PRIMARY_ID') or '',
    'taxonid': sample_xml.findall('SAMPLE_NAME')[0].findtext('TAXON_ID') or '',
    'scientificname': sample_xml.findall('SAMPLE_NAME')[0].findtext('SCIENTIFIC_NAME') or '',
    'attributes': [{'tag': clean(a.findtext('TAG')), 'value': a.findtext('VALUE')} for a in
                   sample_xml.findall('SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')]
  }


def parse_first_experiment(doc):
  # use only first experiment because their attributes are usually almost identical
  exp_xml = doc.find('EXPERIMENT')
  platform_xml = exp_xml.find('PLATFORM')
  return {
    'title': exp_xml.findtext('TITLE') or '',
    'platform': list(platform_xml)[0].tag.lower() or '',
    'instrument_model': platform_xml[0].findtext('INSTRUMENT_MODEL') or '',
    'strategy': exp_xml.xpath('//DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY/text()')[0] or '',
  }


if __name__ == '__main__':
  with tarfile.open(args.metadata_tar_gz, 'r') as tar:
    for item in tar:
      if item.isfile:
        m = re.search(r'\w+/(\w+)\.(\w+)\.xml', item.name)
        if m:
          submission_id = m.group(1)
          filetype = m.group(2)
          # print('Parsing submission "%s" (type "%s")' % (submission_id, filetype))
          if submission_id not in submissions:
            submissions[submission_id] = {"id": submission_id}
          root = etree.XML(tar.extractfile(item).read())
          if filetype == 'experiment':
            submissions[submission_id]['experiment1'] = parse_first_experiment(root)
            pass
          elif filetype == 'run':
            submissions[submission_id]['run_count'] = root.xpath('count(//RUN)')
            pass
          elif filetype == 'sample':
            submissions[submission_id]['sample1'] = parse_first_sample(root)
            submissions[submission_id]['sample_count'] = root.xpath('count(//SAMPLE)')
          elif filetype == 'study':
            pass
          elif filetype == 'submission':
            pass
          elif filetype == 'analysis':
            pass
          else:
            print('Unknown type: %s' % filetype)

  write_output()
