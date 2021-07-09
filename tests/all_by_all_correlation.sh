#! /bin/sh
#
# all_by_all_correlation.py test
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

all_by_all_correlation.py \
    -g=$INPUT_DIR/gene_wide_dataset.tsv \
    -gid=UniqueID \
    -ga=$INPUT_DIR//gene_annotation_file.tsv \
    -gn=GeneName \
    -m=$INPUT_DIR//metabolite_wide_dataset.tsv \
    -mid=UniqueID \
    -ma=$INPUT_DIR//metabolite_annotation_file.tsv \
    -mn=MetName \
    -me=pearson \
    -t=0.05 \
    -o=$OUTPUT_DIR/correlation_file.tsv \
    -c=$OUTPUT_DIR/correlation_matrix.tsv \
    -f=$OUTPUT_DIR/correlation_figure.pdf

echo "### Finished test: ${TEST} on $(date)"

