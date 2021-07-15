#!/bin/bash
# split_wide_dataset.py test
# Copyright (C) 2018-2021 Oleksandr Moskalenko <om@rc.ufl.edu>
# Distributed under terms of the MIT license.

SCRIPT=$(basename "${BASH_SOURCE[0]}");
echo "script $SCRIPT"
TEST="${SCRIPT%.*}"
TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"


mkdir -p test-output
rm -f test-output/gene_* test-output/metabolite_wide_dataset.tsv test-output/met_*tsv

# Test 1 - gene 
split_wide_dataset.py \
    -i=$INPUT_DIR/gene_input_dataset.tsv \
    -p=Gene \
    -s=2,3,4,5,6,7,8,9,10,11 \
    -w=$OUTPUT_DIR/gene_wide_dataset.tsv \
    -d=$OUTPUT_DIR/gene_design_file.tsv \
    -a=$OUTPUT_DIR/gene_annotation_file.tsv

echo "### Finished gene test: ${TEST} on $(date)"

# Test 2 - metabolite
split_wide_dataset.py \
    -i=$INPUT_DIR/metabolite_input_dataset.tsv \
    -p=met \
    -s=2,3,4,5,6,7,8,9,10,11 \
    -w=$OUTPUT_DIR/metabolite_wide_dataset.tsv \
    -d=$OUTPUT_DIR/metabolite_design_file.tsv \
    -a=$OUTPUT_DIR/metabolite_annotation_file.tsv


echo "### Finished metabolite test: ${TEST} on $(date)"

