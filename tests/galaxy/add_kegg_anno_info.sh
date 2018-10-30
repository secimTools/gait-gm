#! /bin/sh
#
# add_kegg_anno_info.py galaxy test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

mkdir -p test-output
rm -f test-output/metabolite_link_kegg_annotation_file_01fhl.tsv

time add_kegg_anno_info.py -s=rno \
    -ga=test-data/ensembl2symbol_annotation_file_01fhl.tsv \
    -gid=UniqueID \
    -gn=GeneSymbol \
    -ma=test-data/metabolite_annotation_file_01fhl.tsv \
    -mid=UniqueID \
    -mn=MetName \
    -go=gene_link_kegg_annotation_file_01fhl.tsv \
    -mo=test-output/metabolite_link_kegg_annotation_file_01fhl.tsv

if [[ -s test-output/metabolite_link_kegg_annotation_file_01fhl.tsv ]]; then
    echo "add_kegg_anno_info.py Test SUCCESS"
else
    echo "add_kegg_anno_info.py Test FAIL"
fi
