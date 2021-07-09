#!/bin/sh
#
# ensembl2symbol.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#
echo "bash source  [${BASH_SOURCE[0]}]"

SCRIPT=$(basename "${BASH_SOURCE[0]}");
echo "script $SCRIPT"

TEST="${SCRIPT%.*}"

TESTDIR="testout/${TEST}"
INPUT_DIR="../galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"

ensembl2symbol.py \
    -s=rat \
    -ga=$INPUT_DIR/gene_annotation_file.tsv \
    -id=UniqueID \
    -e=GeneName \
    -o=$OUTPUT_DIR/ensembl2symbol_annotation_file.tsv

echo "### Finished test: ${TEST} on $(date)"


