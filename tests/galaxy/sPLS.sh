#!/bin/sh
#
# sPLS.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

mkdir -p test-output
rm -f test-output/spls_*

time sPLS.py \
    -m=test-data/metabolite_wide_dataset_01fhl.tsv \
    -mid=UniqueID \
    -mo=generic \
    -mka=test-data/metabolite_to_keggId_link_01fhl.tsv \
    -mn=MetName \
    -g=test-data/gene_wide_dataset_01fhl.tsv \
    -gid=UniqueID \
    -go=all \
    -k=10 \
    -t=0.8 \
    -o1=test-output/spls_correlation_file_01fhl.tsv \
    -f1=test-output/spls_figure_01fhl.pdf


if [[ -s test-output/spls_correlation_file_01fhl.tsv ]] || [[ -s test-output/spls_figure_01fhl.pdf ]]; then
    echo "sPLS.py Test #1 SUCCESS"
else
    echo "sPLS.py Test #1 FAIL"
fi

