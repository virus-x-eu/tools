#!/usr/bin/env python3

import argparse
import tarfile
from lxml import etree
import re
import json
import gzip
from datetime import datetime
import locale

locale.setlocale(locale.LC_ALL, 'en_US.UTF-8')

parser = argparse.ArgumentParser(description='''
Generate Elasticsearch JSON using NCBI SRA Metadata XML files inside a tar archive.
''')
parser.add_argument("metadata_tar_gz", help='NCBI SRA metadata tar.gz file '
                                            '(NCBI_SRA_Metadata_Full_20XXXXXX.tar.gz)')
parser.add_argument("json_gz", help='(gzipped) JSON output file')
args = parser.parse_args()

submissions = {}


def clean(val):
  return val.lower().replace(' ', '_')


date_patterns = [
  '%Y-%m-%d',
  '%m/%d/%Y',
  '%d-%b-%Y',
  '%d-%B-%Y',
  '%d. %b %Y',
  '%d. %B %Y',
  '%d-%m-%Y',
  '%d.%m.%Y',
  '%d-%b-%y',
  '%d-%B-%y',
  '%d. %b %y',
  '%d. %B %y',
  '%d.%m.%y',
  '%Y',
  '%b %Y',
  '%B %Y',
  '%b-%Y',
  '%B-%Y',
  '%m.%Y',
  '%m %Y',
  '%m-%Y',
]


def parse_date(datestr):
  if not isinstance(datestr, str):
    return None
  for p in date_patterns:
    try:
      return datetime.strptime(datestr, p).date()
    except ValueError:
      pass
  return None


def retrieve_value_for_key(sample_xml, key):
  hits = sample_xml.xpath(
    '(//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/TAG[text()="' + key + '"]/following-sibling::VALUE/text())[1]')
  if len(hits) > 0:
    return hits[0]
  return None


def find_coordinate(sample_xml, axis):
  source1 = retrieve_value_for_key(sample_xml, 'lat_lon')
  if source1:
    m = re.search(r'([0-9.]+ [NS]) ([0-9.]+ [WE])', source1)
    if m:
      if axis == 'lat':
        coord = m.group(1).split(' ')
        try:
          if coord[1] == 'S':
            return float(coord[0]) * -1
          else:
            return float(coord[0])
        except ValueError:
          return None
      elif axis == 'lon':
        coord = m.group(2).split(' ')
        try:
          if coord[1] == 'W':
            return float(coord[0]) * -1
          else:
            return float(coord[0])
        except ValueError:
          return None


def parse_sample(sample_xml):
  return {
    'title': sample_xml.findtext('TITLE') or '',
    'id': sample_xml.findall('IDENTIFIERS')[0].findtext('PRIMARY_ID') or sample_xml.get('accession') or '',
    'taxonid': sample_xml.findall('SAMPLE_NAME')[0].findtext('TAXON_ID') or '',
    'scientificname': sample_xml.findall('SAMPLE_NAME')[0].findtext('SCIENTIFIC_NAME') or '',
    'attributes': [{'tag': clean(a.findtext('TAG')), 'value': a.findtext('VALUE')} for a in
                   sample_xml.findall('SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE')],
  }


def parse_first_experiment(doc):
  # use only first experiment because their attributes are usually almost identical
  exp_xml = doc.find('EXPERIMENT')
  platform_xml = exp_xml.find('PLATFORM')
  design_description = exp_xml.xpath('//DESIGN/DESIGN_DESCRIPTION/text()')
  return {
    'id': exp_xml.get('accession') or '',
    'title': exp_xml.findtext('TITLE') or '',
    'platform': str(list(platform_xml)[0].tag).lower() or '',
    'instrument_model': platform_xml[0].findtext('INSTRUMENT_MODEL') or '',
    'strategy': exp_xml.xpath('//DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_STRATEGY/text()')[0] or '',
    'layout': exp_xml.xpath('name(//DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/*[1])').lower() or '',
    'library_source': exp_xml.xpath('//DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_SOURCE/text()')[0].lower() or '',
    'design_description': design_description[0] if len(design_description) > 0 else '',
  }


def parse_first_study(doc):
  # use only first study because generally there is only one per file
  study_xml = doc.find('STUDY')
  return {
    'id': study_xml.get('accession') or '',
    'alias': study_xml.get('alias') or '',
    'title': first_element(study_xml.xpath('DESCRIPTOR/STUDY_TITLE/text()')) or '',
    'abstract': first_element(study_xml.xpath('DESCRIPTOR/STUDY_ABSTRACT/text()')) or '',
  }


def first_element(l):
  if len(l) > 0:
    return l[0]
  else:
    return ''


if __name__ == '__main__':
  submission_count = 1
  prev_submission_id = ""
  with tarfile.open(args.metadata_tar_gz, 'r') as tar, gzip.open(args.json_gz, 'wt') as gz:
    for item in tar:
      if item.isfile:
        m = re.search(r'\w+/(\w+)\.(\w+)\.xml', item.name)
        if m:
          submission_id = m.group(1)
          if prev_submission_id == "":
            submission = {"id": submission_id}
          if submission_id != prev_submission_id and prev_submission_id != "":
            # write previous submission data
            header = {'index': {'_id': submission['id']}}
            gz.write(json.dumps(header))
            gz.write('\n')
            gz.write(json.dumps(submission, separators=(',', ':')))
            gz.write('\n')
            # create object for current submission
            submission = {"id": submission_id}
            submission_count += 1
            if submission_count % 1000 == 0:
              print(submission_count)
          prev_submission_id = submission_id
          filetype = m.group(2)
          # print('Parsing submission "%s" (type "%s")' % (submission_id, filetype))
          try:
            root = etree.XML(tar.extractfile(item).read())
            if filetype == 'experiment':
              experiment1 = parse_first_experiment(root)
              submission['experiment1'] = experiment1
              if 'title' in experiment1:
                submission['experiment1_title'] = experiment1['title']
            elif filetype == 'run':
              submission['run_count'] = int(root.xpath('count(//RUN)'))
            elif filetype == 'sample':
              # use only first sample because their attributes are usually almost identical
              sample_xml = root.find('SAMPLE')
              sample1 = parse_sample(sample_xml)
              submission['sample1'] = sample1
              if 'title' in sample1:
                submission['sample1_title'] = sample1['title']
              submission['sample_count'] = int(root.xpath('count(//SAMPLE)'))
              lat = find_coordinate(sample_xml, 'lat')
              lon = find_coordinate(sample_xml, 'lon')
              if lat and lon:
                submission['sample1_location'] = {'lat': lat, 'lon': lon}
              for attr in sample1['attributes']:
                if attr['tag'] == clean('INSDC first public') or attr['tag'] == clean('ENA-FIRST-PUBLIC'):
                  submission['date'] = attr['value']
                if attr['tag'] in ['collection_date', 'sampling_date', 'run_date'] \
                        and 'date' not in submission:
                  date = parse_date(attr['value'])
                  if date:
                    submission['date'] = str(date)
            elif filetype == 'study':
              study1 = parse_first_study(root)
              submission['study1'] = study1
              if 'title' in study1:
                submission['study1_title'] = study1['title']
            elif filetype == 'submission':
              pass
            elif filetype == 'analysis':
              pass
            else:
              print('Unknown type: %s' % filetype)
          except etree.XMLSyntaxError as xse:
            print(repr(xse))
