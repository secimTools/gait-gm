#! /bin/sh
#
# add_kegg_pathway_info.py galaxy test
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


add_kegg_pathway_info.py \
    -sp=rno \
    -gka=$INPUT_DIR/gene_to_keggId_link.tsv \
    -gid=UniqueID \
    -gn=GeneSymbol \
    -gkid=KEGG_ID \
    -kg2p=KGEN2PATHWAY \
    -mka=$INPUT_DIR/metabolite_to_keggId_link.tsv \
    -mid=UniqueID \
    -mn=MetName \
    -mkid=KEGG_ID \
    -km2p=KMET2PATHWAY \
    -p=PATHWAYS \
    -go=$OUTPUT_DIR/kegg_pathway_gene_output.tsv \
    -mo=$OUTPUT_DIR/kegg_pathway_metabolite_output.tsv


echo "### Finished test: ${TEST} on $(date)"

