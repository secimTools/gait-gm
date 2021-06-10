#!/usr/bin/env python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
#
# DESCRIPTION: This tool takes a column with gene and/or metabolites names and returns an
# annotation file with linked KEGG information (Name in KEGG, KEGG ID, ...)
#
# VERSION: 1.0
#######################################################################################
import os
import argparse
from argparse import RawDescriptionHelpFormatter
import pandas as pd
import gaitGM.keggPeaModules as modules


def getOptions():

    parser = argparse.ArgumentParser(
        description="kegg_anno", formatter_class=RawDescriptionHelpFormatter
    )

    # Tool Input
    tool = parser.add_argument_group(title="Tool Specific Inputs")
    tool.add_argument(
        "-s",
        "--species",
        dest="species",
        action="store",
        required=True,
        help="Specie to download.",
    )
    tool.add_argument(
        "-ga",
        "--geneAnnot",
        dest="geneAnnot",
        action="store",
        required=False,
        help="Gene Annotation File.",
    )
    tool.add_argument(
        "-gid",
        "--geneUniqId",
        dest="geneUniqId",
        action="store",
        required=False,
        help="Name of the column with gene unique Ids.",
    )
    tool.add_argument(
        "-gn",
        "--geneName",
        dest="geneName",
        action="store",
        required=False,
        help="Name of the column with genes names.",
    )
    tool.add_argument(
        "-ma",
        "--metAnnot",
        dest="metAnnot",
        action="store",
        required=False,
        help="Metabolite Annotation File.",
    )
    tool.add_argument(
        "-mid",
        "--metUniqId",
        dest="metUniqId",
        action="store",
        required=False,
        help="Name of the column with metabolite unique Ids.",
    )
    tool.add_argument(
        "-mn",
        "--metName",
        dest="metName",
        action="store",
        required=False,
        help="Name of the column with metabolite names.",
    )

    # Tool Output
    output = parser.add_argument_group(description="Output")
    output.add_argument(
        "-go",
        "--geneOut",
        dest="geneOut",
        action="store",
        required=False,
        help="Gene Output file name.",
    )
    output.add_argument(
        "-mo",
        "--metOut",
        dest="metOut",
        action="store",
        required=False,
        help="Metabolite Output file name.",
    )

    args = parser.parse_args()

    # Standardized paths
    if args.geneAnnot:
        args.geneAnnot = os.path.abspath(args.geneAnnot)
        args.geneOut = os.path.abspath(args.geneOut)
    if args.metAnnot:
        args.metAnnot = os.path.abspath(args.metAnnot)
        args.metOut = os.path.abspath(args.metOut)

    return args


def main(annotFile, uniqId, featureNames, parser, outputFile, featureType):

    """
    Analyze input genes or metabolites finding the corresponding name and its identifier in KEGG.
    Select the best match in similarity terms and discard the other possible ocurrences.
    Inform about ties to facilitate user correction.

    Arguments:
        :param annotFile: Annotation Dataset with gene names information
        :type annotFile: file

        :param uniqId: Name of the column with the gene Unique Identifiers
        :type uniqId: string

        :param featureNames: Name of the column with the gene names
        :type featureNames: string

        :param parser: Downloaded table with KEGG information about gene name and identifier
        :param parser: file

        :param outputFile: Name of the output file to write the information
        :type outputFile: file

        :param featureType: Specification of the feature to analyze, one of: 'Gene' or "Metabolite'
        :type featureType: string
    """

    output = open(outputFile, "w")
    output.write(
        uniqId +
        "\t" +
        featureNames +
        "\tFeature_Type\tMatched\tName_in_KEGG\tKEGG_ID\tSimilarity\tTie\tSelected\n"
    )

    annotFile = pd.read_table(annotFile, delimiter="\t", header=0)
    if "Selected" in annotFile.columns:
        annotFile = annotFile.loc[annotFile["Selected"] != "No"]

    annotFile = annotFile[[uniqId, featureNames]]

    # Output file -> UniqueID Feature_Name Matched Name_In_Kegg KEGG_ID Similarity Tie Selected
    modules.keggAnno(annotFile, parser, output, uniqId, featureNames, featureType)

    output.close()


if __name__ == "__main__":
    args = getOptions()

    if args.geneAnnot:
        modules.checkForDuplicates(args.geneAnnot, args.geneUniqId)
        gene2keggs = modules.downloadGeneParser(args.species)
        main(
            args.geneAnnot,
            args.geneUniqId,
            args.geneName,
            gene2keggs,
            args.geneOut,
            featureType="Gene",
        )
    if args.metAnnot:
        modules.checkForDuplicates(args.metAnnot, args.metUniqId)
        met2keggs = modules.downloadMetParser()
        main(
            args.metAnnot,
            args.metUniqId,
            args.metName,
            met2keggs,
            args.metOut,
            featureType="Metabolite",
        )
