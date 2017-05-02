#!/bin/bash
if [ -z "${1}" ]; then
  echo "Convert fastq to fasta."
  echo "Usage: $0 <reads.fastq.gz>"
else
  zcat "${1}" | paste - - - - | awk 'BEGIN{FS="\t"}{print ">"substr($1,2)"\n"$2}'
fi
