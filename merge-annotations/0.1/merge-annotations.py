#!/usr/bin/env python3

import argparse
import os
import json
import gzip
import csv
import tempfile
import sys
import time

MIN_ALIGNMENT_LENGTH = 15
NUMBER_OF_TOP_HITS_TO_KEEP = 5

PFAM_IDMAP_URL = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/database_files/pfamA.txt.gz'
PDB_IDMAP_URL = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/compound.idx'

TDM_SINGLE_FIELD_MAP = [('nice_name', 'description'), ('seq_count', 'seq_count'),
                        ('structure_count', 'structure_count')]
TDM_MULTI_FIELD_MAP = [('ec_numbers', 'ec'), ('taxonomyids', 'tax'), ('gene_names', 'gene'), ('go_terms', 'go')]

parser = argparse.ArgumentParser(description='Merge annotations.')
parser.add_argument('--dataset-files', help='.genes.json.gz input file(s)', nargs='+')
parser.add_argument('--3dm-dir', dest='tdm_dir', required=True)
parser.add_argument('--3dm-annotation', dest='tdm_annotation', required=True,
                    help='tab-separated file that has to be sorted by evalue (best hits first)')
parser.add_argument('--pfam-idmap', required=True, help='Download from: %s' % PFAM_IDMAP_URL)
# parser.add_argument('--pfam-annotation', required=True,
#                     help='tab-separated file that has to be sorted by evalue (best hits first)')
parser.add_argument('--pdb-idmap', required=True, help='Download from: %s' % PDB_IDMAP_URL)
parser.add_argument('--pdb-annotation', required=True,
                    help='tab-separated file that has to be sorted by evalue (best hits first)')
parser.add_argument('--out-dir', help='Output folder for .genes.json.gz output file(s)')
args = parser.parse_args()

tdm_map = {}
pfam_map = {}
pdb_map = {}

tdm_annotation_map = {}
pfam_annotation_map = {}
pdb_annotation_map = {}


def read_tdm_files():
  filelist = list(os.walk(args.tdm_dir))[0][2]
  print('Loading data from 3DM metadata files...')
  for filename in filelist:
    with open(os.path.join(args.tdm_dir, filename), 'rt') as file:
      # cat names | sed 's/_virusx_2016//' | sed 's/fam/f/'  | sed 's/sub/s/'
      tdm_id = filename.replace('.json', '').replace('_virusx_2016', '').replace('fam', 'f').replace('sub', 's')
      tdm_map[tdm_id] = json.load(file)[0]  # id => object
  print('Got %s entries after parsing %s files.' % (len(tdm_map), len(filelist)))
  print('Done.')


def read_pfam_idmap():
  print('Loading data from PFAM ID map file...')
  row_count = 0
  with gzip.open(args.pfam_idmap, 'rt', newline='') as file:
    pfam = csv.reader(file, delimiter='\t')
    for row in pfam:
      pfam_map[row[0].strip()] = (row[1].strip(), row[3].strip())  # id => (hmmname, description)
      row_count += 1
  print('Got %s entries after parsing %s rows.' % (len(pfam_map), row_count))
  print('Done.')


def read_pdb_idmap():
  print('Loading data from PDB ID map file...')
  row_count = 0
  with open(args.pdb_idmap, 'rt', newline='') as file:
    for _ in [1, 2, 3, 4]:  # skip header rows
      file.readline()
    pdb = csv.reader(file, delimiter='\t')
    for row in pdb:
      pdb_map[row[0].strip()] = row[1].strip()  # id => description
      row_count += 1
  print('Got %s entries after parsing %s rows.' % (len(pdb_map), row_count))
  print('Done.')


def gz_approximate_number_of_records(path):
  sample_amount = 1000
  with tempfile.SpooledTemporaryFile() as tmp:
    with gzip.open(tmp, 'wt') as tmpgz:
      gz_size = os.path.getsize(path)
      with gzip.open(path, 'rt') as file:
        for _ in range(0, sample_amount):
          tmpgz.write(file.readline())
        sample_size = file.tell()
    sample_gz_size = tmp.tell()
  return int(gz_size * (sample_size / sample_gz_size) * (sample_amount / sample_size))


def read_annotation(path, label, target_map):
  status_nth = (0b1 << 20) - 1
  stamp = time.time()
  print('Loading data from %s annotation file...' % label)
  approx_count = gz_approximate_number_of_records(path)
  print('Approx. number of records: %s' % approx_count)
  one_percent = int(approx_count / 100)
  with gzip.open(path, 'rt', newline='') as file:
    ann = csv.reader(file, delimiter='\t')
    row_count = 0
    for row in ann:
      row_count += 1
      if (row_count & status_nth) == status_nth:
        if time.time() > stamp + 5:
          print('%s%%' % int(row_count / one_percent), end=' ')
          stamp = time.time()
        sys.stdout.flush()
      geneid = row[0][10:]
      if not geneid in target_map:
        target_map[geneid] = []
      if int(row[3]) >= MIN_ALIGNMENT_LENGTH and len(target_map[geneid]) < NUMBER_OF_TOP_HITS_TO_KEEP:
        target_map[geneid].append({
          'id': row[1],
          'ident': float(row[2]),
          'alength': int(row[3]),
          'mismatch': int(row[4]),
          'gapopen': int(row[5]),
          'qstart': int(row[6]),
          'qend': int(row[7]),
          'sstart': int(row[8]),
          'send': int(row[9]),
          'evalue': float(row[10]),
          'bitscore': float(row[11])
        })
  print()
  print('Got %s entries after parsing %s rows.' % (len(target_map), row_count))


def extend_gene_annotation_pdb(gene):
  if gene['geneid'] in pdb_annotation_map:
    pdbs = pdb_annotation_map[gene['geneid']]
    for pdb in pdbs:
      try:
        pdb['description'] = pdb_map[pdb['id'][:4]]
      except KeyError:
        pass
    gene['x_pdb'] = pdbs


def extend_gene_annotation_tdm(gene):
  if gene['geneid'] in tdm_annotation_map:
    tdms = tdm_annotation_map[gene['geneid']]
    for tdm in tdms:
      ref = tdm_map[tdm['id']]
      # single value fields
      for source_field, target_field in TDM_SINGLE_FIELD_MAP:
        tdm[target_field] = ref[source_field]
      # multi value maps
      for source_field, target_field in TDM_MULTI_FIELD_MAP:
        tdm[target_field] = []
        if ref[source_field]:
          for k, v in ref[source_field].items():
            tdm[target_field].append({'id': k, 'perc': float(v)})
        elif 'superfamily' in ref:
          for k, v in ref['superfamily'][source_field].items():
            tdm[target_field].append({'id': k, 'perc': float(v)})
    gene['x_3dm'] = tdms


def extend_gene_annotation_pfam(gene):
  if gene['geneid'] in pfam_annotation_map:
    pfams = pdb_annotation_map[gene['geneid']]
    for pfam in pfams:
      pfam['name'], pfam['description'] = pdb_map[pfam['id']]
    gene['x_pdb'] = pfams


def deduplicate_ecs(gene):
  if 'ecs' in gene:
    gene['ecs'] = sorted(list(set(gene['ecs'])))


def extend_gene_annotation(gene):
  extend_gene_annotation_pdb(gene)
  extend_gene_annotation_tdm(gene)
  # extend_gene_annotation_pfam(gene) #TODO
  deduplicate_ecs(gene)
  return gene


def extend_dataset(dataset_file_path):
  status_nth = (0b1 << 14) - 1
  stamp = time.time()
  print("Extending dataset '%s' ..." % os.path.basename(dataset_file_path))
  approx_count = gz_approximate_number_of_records(dataset_file_path)
  print('Approx. number of records: %s' % approx_count)
  one_percent = int(approx_count / 100)
  out_file_path = os.path.join(args.out_dir, os.path.basename(dataset_file_path))
  with gzip.open(dataset_file_path, 'rt') as file, gzip.open(out_file_path, 'wt') as out:
    is_header = True
    row_count = 0
    for line in file:
      row_count += 1
      if is_header:
        out.write(line)
      else:
        json.dump(extend_gene_annotation(json.loads(line)), out, separators=(',', ':'))
        out.write('\n')
      if (row_count & status_nth) == status_nth:
        if time.time() > stamp + 5:
          print('%s%%' % int(row_count / one_percent), end=' ')
          stamp = time.time()
        sys.stdout.flush()
      is_header = not is_header
  print()
  print('Done')
  print('Output: %s' % out_file_path)


if __name__ == '__main__':
  read_tdm_files()
  read_pfam_idmap()
  read_pdb_idmap()
  read_annotation(args.tdm_annotation, '3DM', tdm_annotation_map)
  # read_annotation(args.pfam_annotation, 'Pfam', pfam_annotation_map) #TODO
  read_annotation(args.pdb_annotation, 'PDB', pdb_annotation_map)
  print('\n')
  for d in args.dataset_files:
    extend_dataset(d)
    print('\n')
