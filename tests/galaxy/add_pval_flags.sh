#! /bin/sh
#
# add_pval_flags.py test
# Copyright (C) 2018 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#
mkdir -p test-output
rm -f test-output/add_flags*.tsv

# Test 1
add_pval_flags.py \
    -de=test-data/limma_voom_gene_file_01fhl.tsv \
    -id=UniqueID \
    -p=P.Value \
    -t=0.1,0.05,0.01 \
    -o=test-output/add_flags_gene_output_file_01fhl.tsv \
    -fl=test-output/add_flags_gene_flags_file_01fhl.tsv

if [[ -s test-output/add_flags_gene_output_file_01fhl.tsv ]] && [[ -s test-output/add_flags_gene_flags_file_01fhl.tsv ]]; then
    echo "add_pval_flags.py Test SUCCESS"
else
    echo "add_pval_flags.py Test FAIL"
fi

# Test 2
add_pval_flags.py \
    -de=test-data/limma_voom_metabolite_file_01fhl.tsv \
    -id=UniqueID \
    -p=P.Value \
    -t="0.1,0.05,0.01" \
    -o=test-output/add_flags_metabolite_output_file_01fhl.tsv \
    -fl=test-output/add_flags_metabolite_flags_file_01fhl.tsv

if [[ -s test-output/add_flags_metabolite_output_file_01fhl.tsv ]] && [[ -s test-output/add_flags_metabolite_flags_file_01fhl.tsv ]]; then
    echo "add_pval_flags.py Test #2 SUCCESS"
else
    echo "add_pval_flags.py Test #2 FAIL"
fi
