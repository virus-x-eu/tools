#!/usr/bin/env bash
cat ${1} | parallel --gnu -j ${2} --block ${3} --no-notice --recstart '>' --pipe --cat pfam_scan.pl -as -fasta {} -dir ${4}