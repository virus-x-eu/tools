#!/usr/bin/env bash

echo -n "Deleting index: "; curl -s -XDELETE http://elasticsearch:9200/submission/ >/dev/null && echo "OK" || echo "Failed"
echo -n "Writing mapping: "; curl -sS -XPUT http://elasticsearch:9200/submission/ --data-binary @/import/submission_mapping.json -H "Content-Type: application/json" >/dev/null && echo "OK" || echo "Failed"
zcat /data/submissions.gz | parallel --gnu -j 4 --blocksize 30M -L 500 --no-notice --pipe "curl -sS -H 'Content-Type: application/x-ndjson' -XPOST http://elasticsearch:9200/submission/_bulk --data-binary @- | jq '.items[].index.error.reason | select(. != null)'"