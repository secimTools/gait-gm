#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script for gait_gm."""

from setuptools import setup, find_packages

# with open('README.rst') as readme_file:
#    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['secimtools', 'requests']

setup_requirements = []

test_requirements = []

setup(
    author="Oleksandr Moskalenko",
    author_email='om@rc.ufl.edu',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Modeling Metabolites as a function of gene expression",
    entry_points={
        'console_scripts': [
            'gait_gm=gait_gm.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    # long_description=readme + '\n\n' + history,
    long_description="Modeling Metabolites as a function of gene expression",
    include_package_data=True,
    keywords='gait_gm',
    name='gait_gm',
    packages=find_packages(include=['gait_gm']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/secimTools/gait-gm',
    version='1.0.0',
    zip_safe=False,
    scripts=[
        'src/add_kegg_anno_info.py',
        'src/add_kegg_pathway_info.py',
        'src/add_pval_flags.py',
        'src/all_by_all_correlation.py',
        'src/ensembl2symbol.py',
        'src/keggPeaModules.py',
        'src/split_wide_dataset.py',
        'src/sPLS.py',
        'src/all_by_all_correlation.R',
        'src/sPLS.R'
        ]
)
