#!/usr/bin/env python

# (c) Christian Henke
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0

import re


def extract_unique_element_from_contigid(contigid):
  if re.search("k141_[0-9]*", contigid):
    return contigid[len("k141_"):]
  elif re.search("k[0-9]{1,3}_[0-9]*", contigid):
    return contigid
  elif re.search("NODE_[0-9]*_.*", contigid):
    parts = contigid.split("_", 3)
    return parts[1]
  elif re.search("tig([0-9]{8})_.*", contigid):
    return contigid[3:11]
  elif re.search("contig_[0-9]{1,9}_pilon", contigid):
    return contigid.replace('contig_', '').replace('_pilon', '')
  else:
    raise NotImplementedError('It looks like parsing for this type of contig ID has not been implemented yet!')


def regenerate():
  from Bio import SeqIO
  from Bio.SeqIO import FastaIO
  import argparse

  parser = argparse.ArgumentParser(description='''
  Change the IDs of all fasta entries (contigs) to prefixed, sanitized IDs.
  The output format of each new ID is: PREFIX _ 'C' _ OLD-ID
  ''')
  parser.add_argument("--prefix", dest='dataset_id', required=True, help='ID prefix')
  parser.add_argument("--input", dest='input_fasta', required=True, help='fasta input file')
  parser.add_argument("--output", dest='output_fasta', required=True, help='fasta output file')
  args = parser.parse_args()

  with open(args.input_fasta, 'r') as sourceFile:
    with open(args.output_fasta, 'w') as destFile:
      dest = FastaIO.FastaWriter(destFile, wrap=None)
      dest.write_header()
      for record in SeqIO.parse(sourceFile, "fasta"):
        new_id = args.dataset_id + '_C_' + extract_unique_element_from_contigid(record.id)
        record.id = new_id
        record.description = ""  # any comments after contig id are removed
        dest.write_record(record)


if __name__ == '__main__':
  regenerate()
