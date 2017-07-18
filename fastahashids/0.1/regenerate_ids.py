#!/usr/bin/env python
from Bio import SeqIO
from Bio.SeqIO import FastaIO
import hashlib
import base32_crockford
import argparse

parser = argparse.ArgumentParser(description='''
Change the IDs of all fasta entries to hash-based IDs by using SHA256 and Crockford base32.
Output format of ID is: PREFIX _ NUCL-HASH _ OLD-ID-HASH
''')
parser.add_argument("--prefix", dest='dataset_name', required=True, help='ID prefix')
parser.add_argument("--input", dest='input_fasta', required=True, help='fasta input file')
parser.add_argument("--output", dest='output_fasta', required=True, help='fasta output file')
args = parser.parse_args()

def hash(data):
  h = hashlib.sha256(data).hexdigest()
  i = int(h[0:16], 16)
  return base32_crockford.encode(i).strip()


with open(args.input_fasta, 'rU') as sourceFile:
  with open(args.output_fasta, 'w') as destFile:
    dest = FastaIO.FastaWriter(destFile, wrap=None)
    dest.write_header()
    for record in SeqIO.parse(sourceFile, "fasta"):
      nucl_hash = hash(str(record.seq))
      old_id_hash = hash(record.id)
      new_id = args.dataset_name + '_' + nucl_hash[-6:] + '_' + old_id_hash[-4:]
      record.id = new_id
      rest_of_description = record.description.split(' ')[1:]
      record.description = new_id + ' ' + ' '.join(rest_of_description)
      dest.write_record(record)
