#! /bin/bash
#
# add_kegg_anno_info.py galaxy test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

echo "bash source  [${BASH_SOURCE[0]}]"

SCRIPT=$(basename "${BASH_SOURCE[0]}");
echo "script $SCRIPT"

TEST="${SCRIPT%.*}"

TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"


add_kegg_anno_info.py \
    -s=rno \
    -ga=$INPUT_DIR/ensembl2symbol_annotation_file_01fhl.tsv \
    -gid=UniqueID \
    -gn=GeneSymbol \
    -ma=$INPUT_DIR/metabolite_annotation_file_01fhl.tsv \
    -mid=UniqueID \
    -mn=MetName \
    -go=$OUTPUT_DIR/gene_link_kegg_annotation_file_01fhl.tsv \
    -mo=$OUTPUT_DIR/metabolite_link_kegg_annotation_file_01fhl.tsv

echo "### Finished test: ${TEST} on $(date)"


