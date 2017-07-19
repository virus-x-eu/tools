#!/usr/bin/env bash
cat ${1} | parallel --gnu -j ${2} --recstart '>' --pipe --cat pfam_scan.pl -fasta {} -dir ${3}