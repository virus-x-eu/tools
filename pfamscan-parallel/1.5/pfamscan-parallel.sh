#!/usr/bin/env bash
cat ${1} | parallel --gnu -j ${2} --block ${3} --recstart '>' --pipe --cat pfam_scan.pl -fasta {} -dir ${4}