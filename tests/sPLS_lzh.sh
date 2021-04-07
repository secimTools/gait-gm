#!/bin/sh
#
# sPLS.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

PROJ=/ufgi-serrata/data/zihaoliu/galaxy/tools/gait-gm/galaxy


mkdir -p test-output
rm -f test-output/spls_*

time sPLS.py \
    -m=$PROJ/test-data/metabolite_wide_dataset_01fhl.tsv \
    -mid=UniqueID \
    -mo=generic \
    -mka=$PROJ/test-data/metabolite_to_keggId_link_01fhl.tsv \
    -mn=MetName \
    -g=$PROJ/test-data/gene_wide_dataset_01fhl.tsv \
    -gid=UniqueID \
    -go=pana \
    -gka=$PROJ/test-data/gene_to_keggId_link_01fhl.tsv \
    -gkn=GeneSymbol \
    -p2n=$PROJ/test-data/pathwayId2pathwayNames_01fhl.tsv \
    -p2g=$PROJ/test-data/geneKeggId2pathwayId_01fhl.tsv \
    -cu 0.23 \
    -k=10 \
    -t=0.8 \
    -f=single \
    -o1=test-output/spls_correlation_file_01fhl.tsv \
    -f1=test-output/spls_figure_01fhl.pdf \
    -o3=test-output/spls_pana_outputTable.tsv


if [[ -s test-output/spls_correlation_file_01fhl.tsv ]] || [[ -s test-output/spls_figure_01fhl.pdf ]]; then
    echo "sPLS.py Test #1 SUCCESS"
else
    echo "sPLS.py Test #1 FAIL"
fi

