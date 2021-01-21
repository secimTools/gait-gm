#!/bin/sh
#
# split_wide_dataset.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

mkdir -p test-output
rm -f test-output/gene_* test-output/metabolite_wide_dataset_01fhl.tsv test-output/met_*tsv

# Test 1
time split_wide_dataset.py \
    -i=test-data/gene_input_dataset_01fhl.tsv \
    -p=Gene \
    -s=2,3,4,5,6,7,8,9,10,11 \
    -w=test-output/gene_wide_dataset_01fhl.tsv \
    -d=test-output/gene_design_file_01fhl.tsv \
    -a=test-output/gene_annot_file_01fhl.tsv

if [[ -s test-output/gene_wide_dataset_01fhl.tsv ]] || [[ -s test-output/gene_design_file_01fhl.tsv ]] || [[ -s test-output/gene_annot_file_01fhl.tsv ]]; then
    echo "split_wide_dataset.py Test #1 SUCCESS"
else
    echo "split_wide_dataset.py Test #1 FAIL"
fi

# Test 2
time split_wide_dataset.py \
    -i=test-data/metabolite_input_dataset_01fhl.tsv \
    -p=Met \
    -s=2,3,4,5,6,7,8,9,10,11 \
    -w=test-output/metabolite_wide_dataset_01fhl.tsv \
    -d=test-output/met_design_file_01fhl.tsv \
    -a=test-output/met_annot_file_01fhl.tsv

if [[ -s test-output/metabolite_wide_dataset_01fhl.tsv ]] || [[ -s test-output/met_design_file_01fhl.tsv ]] || [[ -s test-output/met_annot_file_01fhl.tsv ]]; then
    echo "split_wide_dataset.py Test #2 SUCCESS"
else
    echo "split_wide_dataset.py Test #2 FAIL"
fi


