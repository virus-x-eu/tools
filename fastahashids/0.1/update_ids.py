#!/usr/bin/env python

# (c) Christian Henke
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0

from Bio import SeqIO
from Bio.SeqIO import FastaIO
import argparse

parser = argparse.ArgumentParser(description='''
Replace the IDs of the input fasta with those from the mapping file.
''')
parser.add_argument("--input", dest='input_fasta', required=True, help='fasta input file')
parser.add_argument("--output", dest='output_fasta', required=True, help='fasta output file')
parser.add_argument("--use-mapping", dest='input_mapping', required=True, help='mapping input file')
args = parser.parse_args()

def hash(data):
  h = hashlib.sha256(data).hexdigest()
  i = int(h[0:16], 16)
  return base32_crockford.encode(i).strip()

map = {}
with open(args.input_mapping, 'r') as mapping:
  for line in mapping:
    mapping_args = line.strip().split('\t')
    map[mapping_args[0]] = mapping_args[1]

with open(args.input_fasta, 'rU') as sourceFile:
  with open(args.output_fasta, 'w') as destFile:
    dest = FastaIO.FastaWriter(destFile, wrap=None)
    dest.write_header()
    for record in SeqIO.parse(sourceFile, "fasta"):
      new_id = map[record.id]
      record.id = new_id
      rest_of_description = record.description.split(' ')[1:]
      record.description = new_id + ' ' + ' '.join(rest_of_description)
      dest.write_record(record)
