#!/bin/sh
#
# ensembl2symbol.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

mkdir -p test-output
rm -f test-output/ensembl2symbol_*.tsv

time  ensembl2symbol.py \
    -s=rat \
    -ga=test-data/gene_annotation_file_01fhl.tsv \
    -id=UniqueID \
    -e=GeneName \
    -o=test-output/ensembl2symbol_annotation_file_01fhl.tsv

if [[ -s test-output/ensembl2symbol_annotation_file_01fhl.tsv ]]; then
    echo "ensembl2symbol.py Test SUCCESS"
else
    echo "ensembl2symbol.py Test FAIL"
fi

