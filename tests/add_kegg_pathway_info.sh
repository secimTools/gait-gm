#! /bin/bash
# add_kegg_pathway_info.py galaxy test
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

add_kegg_pathway_info.py \
    -sp=rno \
    -gka=$INPUT_DIR/gene_to_keggId_link.tsv \
    -gid=UniqueID \
    -gn=GeneSymbol \
    -gkid=KEGG_ID \
    -mka=$INPUT_DIR/metabolite_to_keggId_link.tsv \
    -mid=UniqueID \
    -mn=MetName \
    -mkid=KEGG_ID \
    -kg2p=$OUTPUT_DIR/geneKeggId2pathwayId.tsv \
    -km2p=$OUTPUT_DIR/metaboliteKeggId2pathwayId.tsv \
    -go=$OUTPUT_DIR/gene_kegg_pathway.tsv \
    -mo=$OUTPUT_DIR/metabolite_kegg_pathway.tsv \
    -p=pathwayId2pathwayNames.tsv


echo "### Finished test: ${TEST} on $(date)"

