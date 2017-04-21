#!/bin/bash
if [ -z "${2}" ]; then
  echo "Merge kallisto tsv counts and TPMs."
  echo "Usage: $0 <abundance.paired.tsv> <abundance.unpaired.tsv>"
else
  paste \
    <(tail -n+2 "${1}") \
    <(paste <(tail -n+2 "${1}" | awk '{print $1 "\t" $4}') <(tail -n+2 "${2}" | awk '{print $4}') | awk '{print $2+$3}') \
    | awk -v pm=$(paste \
      <(tail -n+2 "${1}") <(paste <(tail -n+2 "${1}" | awk '{print $1 "\t" $4}') \
      <(tail -n+2 "${2}" | awk '{print $4}') | awk '{print $2+$3}') | awk '{sum += $6/($3/1000)} END  {print sum/1000000}') \
    '{print $1 "\t" $6 "\t" ($6/($3/1000))/pm}'
fi

