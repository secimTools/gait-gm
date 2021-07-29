#!/bin/bash

PROJ=~/projects/bioinformatics/mcintyre-cid/GAIT-GM

# TOOLS
declare -a tests=("split_wide_dataset" "ensembl2symbol" "add_pval_flags"  "sPLS" "all_by_all_correlation" "add_kegg_anno_info" "add_kegg_pathway_info")
#rm -f "testout/galaxy_test_timings"

if [[ ! $# -eq 1 ]]; then
    echo "Pass a single argument - tool name. E.g. add_pval_flags"
    exit 0
else
    test=$1
    if [[ ! "${tests[@]}"  =~ ${test} ]]; then
        echo "Unknown test"
        exit 0
    else
        OUTDIR="testout/galaxy/${test}"
        rm -rf "${OUTDIR}"
        mkdir -p "${OUTDIR}"
        echo "Running test ${test}"
        echo "START ${test}: $(date)" >> "testout/galaxy_test_timings"
        planemo test --galaxy_root "${HOME}/projects/bioinformatics/mcintyre-cid/galaxy/" \
            --skip_venv --no_cleanup --tool_dependency_dir "${PROJ}/conda" \
            --test_data "${PROJ}/gait_gm/galaxy/test-data/" --no_conda_auto_init \
            --no_conda_auto_install --conda_use_local \
            --test_output "${OUTDIR}/test_output.html" \
            --test_output_text "${OUTDIR}/test_output.txt" \
            --job_output_files "${OUTDIR}" \
            "${PROJ}/gait_gm/galaxy/${test}.xml"
        echo "END ${test}: $(date)" >> "testout/galaxy_test_timings"
    fi
fi

