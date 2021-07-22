#! /bin/sh
#
# galaxy_sPLS.sh
# Copyright (C) 2021 Oleksandr Moskalenko <om@rc.ufl.edu>
#
# Distributed under terms of the MIT license.
#

PROJ=~/projects/bioinformatics/mcintyre-cid/GAIT-GM

planemo test --galaxy_root ~/projects/bioinformatics/mcintyre-cid/galaxy/ \
    --skip_venv --no_cleanup \
    --tool_dependency_dir ${PROJ}/conda \
    --test_data ${PROJ}/gait_gm/galaxy/test-data/ \
    --no_conda_auto_init --no_conda_auto_install \
    --conda_use_local ${PROJ}/gait_gm/galaxy/sPLS.xml
