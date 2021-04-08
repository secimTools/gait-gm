#! /bin/sh
#
# add_pval_flags.py test
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



# Test 1
add_pval_flags.py \
    -de=$INPUT_DIR/limma_voom_gene_file_01fhl.tsv \
    -id=UniqueID \
    -p=P.Value \
    -t=0.1,0.05,0.01 \
    -o=$OUTPUT_DIR/add_flags_gene_output_file_01fhl.tsv \
    -fl=$OUTPUT_DIR/add_flags_gene_flags_file_01fhl.tsv

echo "### Finished test1: ${TEST} on $(date)"


# Test 2
add_pval_flags.py \
    -de=$INPUT_DIR/limma_voom_metabolite_file_01fhl.tsv \
    -id=UniqueID \
    -p=P.Value \
    -t="0.1,0.05,0.01" \
    -o=$OUTPUT_DIR/add_flags_metabolite_output_file_01fhl.tsv \
    -fl=$OUTPUT_DIR/add_flags_metabolite_flags_file_01fhl.tsv

echo "### Finished test2: ${TEST} on $(date)"

