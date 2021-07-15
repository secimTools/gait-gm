#! /bin/bash
# ensembl2symbol.py test
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

ensembl2symbol.py \
    -s=rat \
    -ga=$INPUT_DIR/gene_annotation.tsv \
    -id=UniqueID \
    -e=GeneName \
    -o=$OUTPUT_DIR/ensembl2symbol_annotation.tsv

echo "### Finished test: ${TEST} on $(date)"

