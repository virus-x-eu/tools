#!/usr/bin/env python3

# (c) Christian Henke
# Licensed under Apache License 2.0
# https://www.apache.org/licenses/LICENSE-2.0

from Bio import SeqIO
import argparse
import json
import pysam
import numpy
from math import ceil
from multiprocessing import Pool, cpu_count

parser = argparse.ArgumentParser(description='''
Convert fasta to Elasticsearch json while including contig based annotation information and coverage
''')
parser.add_argument("--fasta", dest='input_fasta', required=True, help='fasta input file')
parser.add_argument("--virsorter-csv", dest='input_virsorter', required=True, help='VirSorter input file')
parser.add_argument('--sample-names', nargs='*', dest='names')
parser.add_argument('--sample-bam-files', nargs='*', dest='bams')
parser.add_argument("--json", dest='output_json', required=True, help='json output file')
parser.add_argument("--threads", type=int, dest='threads', help='number parallel threads for BAM file reader',
                    default=cpu_count())
args = parser.parse_args()

threads = args.threads

coverage_avg_span = 10
ones_val = int(coverage_avg_span / 2)

contigid_to_virsorter_phage = {}
contigid_to_virsorter_prophage = {}

coverage_map = {}


def parse_bam_for_ids(chunk, sample_name, bam_path):
  print('Parsing chunk of size %s' % len(chunk))
  coverage_map_chunk = {}
  sam_tmp = pysam.AlignmentFile(bam_path, "rb")
  for ref in chunk:
    ref_id = ref.split(' ', 1)[0]
    if ref_id not in coverage_map_chunk:
      coverage_map_chunk[ref_id] = {}
    try:
      cov_vals = [p.n for p in sam_tmp.pileup(ref)]
      if len(cov_vals) < coverage_avg_span:
        coverage_map_chunk[ref_id][sample_name] = []
      else:
        coverage_map_chunk[ref_id][sample_name] = numpy.convolve(cov_vals, numpy.ones(ones_val, ) / ones_val,
                                                                 mode='same')[
                                                  ones_val - 1:-ones_val:coverage_avg_span].astype(int).tolist()
    except ValueError:
      print('ERROR: Cannot find index file %s.bai' % bam_path)
      return False
    except:
      print('ERROR: Unknown error in thread.')
      return False
  print('Chunk done')
  return coverage_map_chunk


def merge_with_global_coverage_map(map_chunk):
  for ref_id, ref_obj in map_chunk.items():
    if ref_id in coverage_map:
      coverage_map[ref_id].update(ref_obj)
    else:
      coverage_map[ref_id] = ref_obj


if args.names and args.bams:
  named_bams = zip(args.names, args.bams)

  for name, bam in named_bams:
    print("For sample name '%s' parsing BAM file: %s" % (name, bam))
    sam = pysam.AlignmentFile(bam, "rb")
    total_ref_count = len(sam.references)
    print('File has %s refs' % total_ref_count)
    chunk_size = ceil(total_ref_count / threads)
    chunks = [sam.references[i:i + chunk_size] for i in range(0, len(sam.references), chunk_size)]
    print('Using %s threads' % threads)
    # real processes are used instead of threads to get around the CPython Global Interpreter Lock
    with Pool(processes=threads) as pool:
      results = [pool.apply_async(parse_bam_for_ids, args=(chunk, name, bam)) for chunk in chunks]
      for result in results:
        proc_result = result.get()
        if proc_result:
          merge_with_global_coverage_map(proc_result)
        else:
          exit(1)

with open(args.input_fasta, 'rU') as sourceFile, open(args.input_virsorter, 'r') as virsorterFile:
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
    print('Writing output to %s' % args.output_json)
    for record in SeqIO.parse(sourceFile, "fasta"):
      id = str(record.id)
      header = {"index": {"_id": str(record.id)}}
      obj = {
        "id": str(record.id),
        "nucleotide": str(record.seq),
        "length": len(str(record.seq))
      }
      if args.names and args.bams and str(record.id) in coverage_map:
        obj['coverage'] = [{'name': k, 'values': v} for k, v in coverage_map[str(record.id)].items()]
      if id in contigid_to_virsorter_phage or id in contigid_to_virsorter_prophage:
        obj["virsorter"] = {}
        if id in contigid_to_virsorter_phage:
          obj["virsorter"]["phage"] = contigid_to_virsorter_phage[id]
        if id in contigid_to_virsorter_prophage:
          obj["virsorter"]["prophage"] = contigid_to_virsorter_prophage[id]
      destFile.write(json.dumps(header))
      destFile.write('\n')
      destFile.write(json.dumps(obj, separators=(',', ':')))
      destFile.write('\n')
