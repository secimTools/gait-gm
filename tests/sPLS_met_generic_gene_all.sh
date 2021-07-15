#!/bin/bash
# sPLS.py test
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

sPLS.py \
    -m=$INPUT_DIR/metabolite_wide_dataset.tsv \
    -mid=UniqueID \
    -mo=generic \
    -mka=$INPUT_DIR/metabolite_to_keggId_link.tsv \
    -mn=MetName \
    -g=$INPUT_DIR/gene_wide_dataset.tsv \
    -gid=UniqueID \
    -go=all \
    -k=10 \
    -t=0.8 \
    -o1=$OUTPUT_DIR/spls_correlation_generic_all.tsv \
    -f1=$OUTPUT_DIR/spls_figure_generic_all.pdf

echo "### Finished test: ${TEST} on $(date)"
