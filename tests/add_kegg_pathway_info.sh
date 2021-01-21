#! /bin/sh
#
# add_kegg_pathway_info.py galaxy test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

time add_kegg_pathway_info.py \
    -sp=rno \
    -gka=test-data/gene_to_keggId_link_01fhl.tsv \
    -gid=UniqueID \
    -gn=GeneSymbol \
    -gkid=KEGG_ID \
    -kg2p=KGEN2PATHWAY \
    -mka=test-data/metabolite_to_keggId_link_01fhl.tsv \
    -mid=UniqueID \
    -mn=MetName \
    -mkid=KEGG_ID \
    -km2p=KMET2PATHWAY \
    -p=PATHWAYS \
    -go=test-output/kegg_pathway_gene_output.tsv \
    -mo=test-output/kegg_pathway_metabolite_output.tsv


