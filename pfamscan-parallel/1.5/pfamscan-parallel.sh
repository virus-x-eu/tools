#!/usr/bin/env bash
cat ${1} | parallel --gnu -j ${2} --block 100K --recstart '>' --pipe --cat pfam_scan.pl -fasta {} -dir ${3}