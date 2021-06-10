#!/usr/bin/env python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
# DESCRIPTION: Take a Gene Expression Annotation File with ENSEMBL IDs and return an Annotation
# File with Gene Symbols
#######################################################################################

import os
import logging
import argparse
from argparse import RawDescriptionHelpFormatter
import mygene
import pandas as pd
import gaitGM.keggPeaModules as modules
from secimtools.dataManager import logger as sl


def getOptions():
    parser = argparse.ArgumentParser(
        description="kegg_anno", formatter_class=RawDescriptionHelpFormatter
    )
    tool = parser.add_argument_group(title="Tool Specific Inputs")
    tool.add_argument(
        "-s",
        "--species",
        dest="species",
        action="store",
        required=True,
        help="Species to download. One of rat, human, mouse or fruitfly",
    )
    tool.add_argument(
        "-ga",
        "--geneAnnot",
        dest="geneAnnot",
        action="store",
        required=True,
        help="Gene Expression Annotation File.",
    )
    tool.add_argument(
        "-id",
        "--uniqId",
        dest="uniqId",
        action="store",
        required=True,
        help="Name of the column with gene Unique Ids.",
    )
    tool.add_argument(
        "-e",
        "--ensemblId",
        dest="ensemblId",
        action="store",
        required=True,
        help="Name of the column with ENSEMBL IDs.",
    )
    output = parser.add_argument_group(description="Output")
    output.add_argument(
        "-o",
        "--output",
        dest="output",
        action="store",
        required=True,
        help="Gene expression Annotation File with Gene Symbols.",
    )

    args = parser.parse_args()

    args.geneAnnot = os.path.abspath(args.geneAnnot)
    args.output = os.path.abspath(args.output)

    return args


def main():
    """
    Take a gene expression matrix and extract the column with ENSEMBL IDs. Then, translate ENSEMBL
    IDs into Gene_Symbol needed for the rest of the pipeline. Create a table with Unique
    Identifiers, ENSEMBL IDs, gene symbols, match scores, and selection resul, which is useful in
    cases of multiple matches.

    Arguments:
        :param species: Species to download information from mygene
        :type species: string

        :param geneAnnot: Gene Expression Annotation file with ENSEMBL IDs column
        :type geneAnnot: file

        :param ensemblId: Name of the column with ENSEMBL IDs
        :type ensemblId: string
    """

    args = getOptions()
    logger = logging.getLogger()
    sl.setLogger(logger)
    logger.info(
        u"""Importing data with following parameters: \
        \n\tSpecies: {0}\
        \n\tGene Annotation File: {1}\
        \n\tUnique ID column: {2}\
        \n\tENSEMBL ID Column: {3}""".format(
            args.species, args.geneAnnot, args.uniqId, args.ensemblId
        )
    )
    modules.checkForDuplicates(args.geneAnnot, args.uniqId)

    # Original Gene Expression Annotation Dataset with ENSEMBL IDs
    genesTable = pd.read_table(args.geneAnnot, delimiter="\t", header=0)

    # Find Gene Symbol
    mg = mygene.MyGeneInfo()
    genes = genesTable[args.ensemblId].tolist()
    genesTransformed = mg.querymany(
        genes,
        scopes="ensembl.gene",
        fields="symbol",
        species=args.species,
        verbose=False,
        returnall=True,
        as_dataframe=True,
        df_index=False,
    )
    genesTransformedTable = genesTransformed["out"]
    if genesTransformedTable.shape[0] != genesTransformed["missing"].shape[0]:
        genesTransformedTable.drop(labels=["_id"], axis=1, inplace=True)
        genesTransformedTable = genesTransformedTable[["query", "symbol", "_score"]]
        genesTransformedTable.columns = [args.ensemblId, "GeneSymbol", "Score"]

        # Merge Both datasets
        newGenesTable = pd.merge(genesTable, genesTransformedTable, on=args.ensemblId)

        # In case of duplicated, select the first one (High score)
        newGenesTable["Selected"] = "Yes"
        isDup = newGenesTable.duplicated(subset=args.ensemblId, keep="first")
        newGenesTable["Selected"][isDup] = "No"

        # Write table
        newGenesTable.to_csv(args.output, sep="\t", index=False)
    else:
        with open(args.output, 'w') as f:
            f.write("no matching result! Please check the selected species and input files.")

if __name__ == "__main__":
    main()
