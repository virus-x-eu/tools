#!/usr/bin/env python

# (c) Christian Henke
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0

from Bio import SeqIO
import argparse
import json

parser = argparse.ArgumentParser(description='''
Convert fasta to Elasticsearch json while including contig based annotation information
''')
parser.add_argument("--fasta", dest='input_fasta', required=True, help='fasta input file')
parser.add_argument("--virsorter-csv", dest='input_virsorter', required=True, help='VirSorter input file')
parser.add_argument("--json", dest='output_json', required=True, help='json output file')
args = parser.parse_args()

contigid_to_virsorter_phage = {}
contigid_to_virsorter_prophage = {}

with open(args.input_fasta, 'rU') as sourceFile,  open(args.input_virsorter, 'r') as virsorterFile:
    phage_confidence = ''
    prophage_confidence = ''
    state = ''
    for line in virsorterFile:
        if line.strip() == '## 1 - Complete phage contigs - category 1 (sure)':
            phage_confidence = 'sure'
            state = 'phage'
        elif line.strip() == '## 2 - Complete phage contigs - category 2 (somewhat sure)':
            phage_confidence = 'somewhat sure'
            state = 'phage'
        elif line.strip() == '## 3 - Complete phage contigs - category 3 (not so sure)':
            phage_confidence = 'not so sure'
            state = 'phage'
        elif line.strip() == '## 4 - Prophages - category 1 (sure)':
            prophage_confidence = 'sure'
            state = 'prophage'
        elif line.strip() == '## 5 - Prophages - category 2 (somewhat sure)':
            prophage_confidence = 'somewhat sure'
            state = 'prophage'
        elif line.strip() == '## 6 - Prophages - category 3 (not so sure)':
            prophage_confidence = 'not so sure'
            state = 'prophage'
        else:
            if not line.startswith('#'):
                contigid = line.split(',')[0][len('VIRSorter_'):line.index('_flag=')]
                if state == 'phage':
                    contigid_to_virsorter_phage[contigid] = phage_confidence
                elif state == 'prophage':
                    contigid_to_virsorter_prophage[contigid] = prophage_confidence

    with open(args.output_json, 'w') as destFile:
        for record in SeqIO.parse(sourceFile, "fasta"):
            id = str(record.id)
            header = {"index": {"_id": str(record.id)}}
            obj = {
                "id": str(record.id),
                "nucleotide": str(record.seq),
                "length": len(str(record.seq)),
            }
            if id in contigid_to_virsorter_phage or id in contigid_to_virsorter_prophage:
                obj["virsorter"] = {}
                if id in contigid_to_virsorter_phage:
                    obj["virsorter"]["phage"] = contigid_to_virsorter_phage[id]
                if id in contigid_to_virsorter_prophage:
                    obj["virsorter"]["prophage"] = contigid_to_virsorter_prophage[id]
            destFile.write(json.dumps(header))
            destFile.write('\n')
            destFile.write(json.dumps(obj))
            destFile.write('\n')
