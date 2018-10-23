=======
GAIT-GM
=======


.. image:: https://img.shields.io/pypi/v/gait_gm.svg
        :target: https://pypi.python.org/pypi/gait_gm

.. image:: https://img.shields.io/travis/moskalenko/gait_gm.svg
        :target: https://travis-ci.org/moskalenko/gait_gm


Modeling Metabolites as a function of gene expression


* Free software: MIT license


Features
--------

* add_kegg_anno_info.py: This tool takes a column with gene and/or metabolites names and returns an
annotation file with linked KEGG information (Name in KEGG, KEGG ID, ...)
* add_kegg_pathway_info.py: Download the Kegg information of the selected species and associate
gene expression and/or metabolomics data with pathways.
* add_pval_flags.py: Add flags to a Differential Expression Analysis Dataset. Output shows whether
the feature belongs (1) or not (0) to a certain flag.
* all_by_all_correlation.py: Perform a correlation analysis of two datasets.
* ensembl2symbol.py: Take a Gene Expression Annotation File with ENSEMBL IDs and return an
Annotation File with Gene Symbols.
keggPeaModules.py: Extra functions for KEGG PEA Tools
pathway_enrichment_analysis.py: Perform Pathway Enrichment Analysis
split_wide_dataset.py: Take a wide dataset and append a unique identifier to it.
* sPLS.py: Take a column with the gene and/or metabolites names and return an annotation file with
KEGG information (Name in KEGG, KEGG ID, ...)
sPLS.R: Perform a sparse PLS over subsets of gene Expression and metabolite data.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

