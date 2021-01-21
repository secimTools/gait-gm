#! /bin/sh
#
# all_by_all_correlation.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

mkdir -p test-output
rm -f test-output/correlation_*

time all_by_all_correlation.py \
    -g=test-data/gene_wide_dataset_01fhl.tsv \
    -gid=UniqueID \
    -ga=test-data/gene_annotation_file_01fhl.tsv \
    -gn=GeneName \
    -m=test-data/metabolite_wide_dataset_01fhl.tsv \
    -mid=UniqueID \
    -ma=test-data/metabolite_annotation_file_01fhl.tsv \
    -mn=MetName \
    -me=pearson \
    -t=0.05 \
    -o=test-output/correlation_file_01fhl.tsv \
    -c=test-output/correlation_matrix_01fhl.tsv \
    -f=test-output/correlation_figure_01fhl.pdf

if [[ -s test-output/correlation_file_01fhl.tsv ]] && [[ -s test-output/correlation_matrix_01fhl.tsv ]] && [[ -s test-output/correlation_figure_01fhl.pdf ]]; then
    echo "all_by_all_correlation.py Test SUCCESS"
else
    echo "all_by_all_correlation.py Test FAIL"
fi
