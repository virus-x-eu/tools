#!/usr/bin/env python

# (c) Christian Henke
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0

from Bio import SeqIO
from Bio.SeqIO import FastaIO
import argparse
import json

parser = argparse.ArgumentParser(description='''
Convert fasta to json
''')
parser.add_argument("--fasta", dest='input_fasta', required=True, help='fasta input file')
parser.add_argument("--json", dest='output_json', required=True, help='json output file')
args = parser.parse_args()

with open(args.input_fasta, 'rU') as sourceFile:
    with open(args.output_json, 'w') as destFile:
        for record in SeqIO.parse(sourceFile, "fasta"):
            obj = {"id": str(record.id), "nucleotide": str(record.seq), "length": len(str(record.seq))}
            destFile.write(json.dumps(obj))
            destFile.write('\n')
