#!/bin/bash
#
# sPLS.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

SCRIPT=$(basename "${BASH_SOURCE[0]}");
TEST="${SCRIPT%.*}"

TESTDIR="testout/${TEST}"
INPUT_DIR="galaxy/test-data"
OUTPUT_DIR=$TESTDIR
rm -rf "${TESTDIR}"
mkdir -p "${TESTDIR}"
echo "### Starting test: ${TEST}"

sPLS.py \
    --metDataset=$INPUT_DIR/metabolite_wide_dataset_01fhl.tsv \
    --metId=UniqueID \
    --metOption=mmc \
    --metAnno=$INPUT_DIR/metabolite_annotation_file_01fhl.tsv \
    --metName=MetName \
    --metKeggAnno=$INPUT_DIR/metabolite_to_keggId_link_01fhl.tsv \
    --design=$INPUT_DIR/gene_design_file_01fhl.tsv \
    --correlation=pearson \
    --sigmaLow=0.05 \
    --sigmaHigh=0.50 \
    --sigmaNum=451 \
    --geneDataset=$INPUT_DIR/gene_wide_dataset_01fhl.tsv \
    --geneId=UniqueID \
    --geneOption=pana \
    --geneKeggAnno=$INPUT_DIR/gene_to_keggId_link_01fhl.tsv \
    --geneKeggName=Gene \
    --path2genes=$INPUT_DIR/geneKeggId2pathwayId_01fhl.tsv \
    --path2names=$INPUT_DIR/gene_kegg_pathway_file_01fhl.tsv \
    --facSel=single \
    --geneAnno=$INPUT_DIR/gene_annotation_file_01fhl.tsv \
    --geneName=GeneName \
    -keepX=5 \
    -thres=0.8 \
    --figure1=$OUTPUT_DIR/integration_mmc_pana_genes_sPLS_fig.pdf \
    --splsOut=$OUTPUT_DIR/integration_mmc_pana_genes_sPLS_output.tsv \
    --figure2=$OUTPUT_DIR/integration_mmc_pana_genes_MMC_heatmap.pdf \
    --mmcOut=$OUTPUT_DIR/integration_mmc_pana_genes_MMC_output.tsv \
    --panaOut=$OUTPUT_DIR/integration_mmc_pana_genes_PANA_output.tsv


echo "### Finished test: ${TEST} on $(date)"
