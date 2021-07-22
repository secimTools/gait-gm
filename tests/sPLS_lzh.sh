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
    -m=${INPUT_DIR}/metabolite_wide_dataset.tsv \
    -mid=UniqueID \
    -mo=generic \
    -mka=${INPUT_DIR}/metabolite_to_keggId_link.tsv \
    -mn=MetName \
    -g=${INPUT_DIR}/gene_wide_dataset.tsv \
    -gid=UniqueID \
    -go=pana \
    -gka=${INPUT_DIR}/gene_to_keggId_link.tsv \
    -gkn=GeneSymbol \
    -p2n=${INPUT_DIR}/pathwayId2pathwayNames.tsv \
    -p2g=${INPUT_DIR}/geneKeggId2pathwayId.tsv \
    -cu 0.23 \
    -k=10 \
    -t=0.8 \
    -f=single \
    -o1=${OUTPUT_DIR}/spls_correlation.tsv \
    -f1=${OUTPUT_DIR}/spls_figure.pdf \
    -o3=${OUTPUT_DIR}/spls_pana_outputTable.tsv
ec=?!
if [[ $ec ]]; then
    echo "### Finished metabolite test: ${TEST} on $(date)"
else
    echo "Test ${TEST} failed. Exit code $ec"
fi
