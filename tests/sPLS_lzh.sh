#!/bin/sh
#
# sPLS.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

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
    -m=${INPUT_DIR}/metabolite_wide_dataset_01fhl.tsv \
    -mid=UniqueID \
    -mo=generic \
    -mka=${INPUT_DIR}/metabolite_to_keggId_link_01fhl.tsv \
    -mn=MetName \
    -g=${INPUT_DIR}/gene_wide_dataset_01fhl.tsv \
    -gid=UniqueID \
    -go=pana \
    -gka=${INPUT_DIR}/gene_to_keggId_link_01fhl.tsv \
    -gkn=GeneSymbol \
    -p2n=${INPUT_DIR}/pathwayId2pathwayNames_01fhl.tsv \
    -p2g=${INPUT_DIR}/geneKeggId2pathwayId_01fhl.tsv \
    -cu 0.23 \
    -k=10 \
    -t=0.8 \
    -f=single \
    -o1=${OUTPUT_DIR}/spls_correlation_file_01fhl.tsv \
    -f1=${OUTPUT_DIR}/spls_figure_01fhl.pdf \
    -o3=${OUTPUT_DIR}/spls_pana_outputTable.tsv

echo "### Finished metabolite test: ${TEST} on $(date)"
