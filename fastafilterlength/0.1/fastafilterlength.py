#!/usr/bin/env python

# (c) Christian Henke
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0

from Bio import SeqIO
from Bio.SeqIO import FastaIO
import argparse

parser = argparse.ArgumentParser(description='''
Remove sequences from fasta that are shorter than the desired minimum length.
''')
parser.add_argument("--input", dest='input_fasta', required=True, help='fasta input file')
parser.add_argument("--output", dest='output_fasta', required=True, help='fasta output file')
parser.add_argument("--min-seq-length", dest='min_seq_len', type=int, required=True, help='minimum sequence length (integer)')
args = parser.parse_args()

with open(args.input_fasta, 'r') as sourceFile:
  with open(args.output_fasta, 'w') as destFile:
    dest = FastaIO.FastaWriter(destFile, wrap=None)
    dest.write_header()
    for record in SeqIO.parse(sourceFile, "fasta"):
      if len(record.seq) >= args.min_seq_len:
        dest.write_record(record)
