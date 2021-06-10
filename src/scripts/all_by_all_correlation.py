#!/usr/bin/env python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
# DESCRIPTION: Perform a correlation analysis of two datasets.
#######################################################################################

import os
import logging
import warnings
import argparse
from argparse import RawDescriptionHelpFormatter
import matplotlib
import pandas as pd
from rpy2.robjects import pandas2ri
from rpy2.rinterface import RRuntimeWarning
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
import gaitGM.keggPeaModules as modules
from secimtools.dataManager import logger as sl
from importlib import resources as ires
matplotlib.use("Agg")


def getOptions():
    parser = argparse.ArgumentParser(
        description="allByAllCorr", formatter_class=RawDescriptionHelpFormatter
    )
    tool = parser.add_argument_group(description="Input")
    tool.add_argument(
        "-g",
        "--geneDataset",
        dest="geneDataset",
        action="store",
        required=True,
        help="Gene Expression dataset.",
    )
    tool.add_argument(
        "-gid",
        "--geneId",
        dest="geneId",
        action="store",
        required=True,
        help="Gene Unique ID column name.",
    )
    tool.add_argument(
        "-ga",
        "--geneAnnot",
        dest="geneAnnot",
        action="store",
        required=False,
        help="Gene Expression Annotation Dataset.",
    )
    tool.add_argument(
        "-gn",
        "--geneName",
        dest="geneName",
        action="store",
        required=False,
        help="Gene Expression Annotation Dataset column.",
    )
    tool.add_argument(
        "-m",
        "--metDataset",
        dest="metDataset",
        action="store",
        required=True,
        help="Metabolomic Datset.",
    )
    tool.add_argument(
        "-mid",
        "--metId",
        dest="metId",
        action="store",
        required=True,
        help="Metabolite Unique ID column name.",
    )
    tool.add_argument(
        "-ma",
        "--metAnnot",
        dest="metAnnot",
        action="store",
        required=False,
        help="Metabolomic Annotation Dataset.",
    )
    tool.add_argument(
        "-mn",
        "--metName",
        dest="metName",
        action="store",
        required=False,
        help="Metabolomics Annotation Dataset column.",
    )
    tool.add_argument(
        "-me",
        "--meth",
        dest="meth",
        action="store",
        required=True,
        help="Correlation coefficient to be computed.",
    )
    tool.add_argument(
        "-t",
        "--thres",
        dest="thres",
        action="store",
        required=True,
        help="Pvalue threshold for the output.",
    )
    output = parser.add_argument_group(description="Output")
    output.add_argument(
        "-o",
        "--output",
        dest="output",
        action="store",
        required=True,
        help="Output Table.",
    )
    output.add_argument(
        "-c",
        "--corMat",
        dest="corMat",
        action="store",
        required=True,
        help="Correlation Matrix.",
    )
    output.add_argument(
        "-f",
        "--fig",
        dest="fig",
        action="store",
        required=True,
        help="Output figure name for results [pdf].",
    )

    args = parser.parse_args()

    args.geneDataset = os.path.abspath(args.geneDataset)
    args.metDataset = os.path.abspath(args.metDataset)
    args.output = os.path.abspath(args.output)
    args.corMat = os.path.abspath(args.corMat)
    args.fig = os.path.abspath(args.fig)

    if args.geneAnnot:
        args.geneAnnot = os.path.abspath(args.geneAnnot)
    else:
        args.geneAnnot = ""
        args.geneAnnotName = ""
    if args.metAnnot:
        args.metAnnot = os.path.abspath(args.metAnnot)
    else:
        args.metAnnot = ""
        args.metAnnotName = ""
    if not args.thres:
        args.thres = ""

    return args


def main():
    """
    Perform a correlation analysis of a Gene Expression Dataset and a Metabolomic Dataset.

    Arguments:
        :param geneDataset metDataset: Gene expression/Metabolomics wide dataset, respectively.
        :type geneDataset metDataset: files

        :param geneId metId: Name of the Genes/metabolites unique identifier column, respectively.
        :type geneId metId: strings

        :param geneAnnot metAnnot: Gene Expression/Metabolomics Annotation Datasets, respectively.
        :type geneAnnot metAnnot: files

        :param geneAnnotName metAnnotName: Name of the column of the Annotation file that contains
        genes/metabolites names respectively.
        :type geneAnnotName metAnnotName: strings

        :param meth: Methodology for the correlation function. One of 'pearson', 'spearman' or
        'kendall'.
        :type meth: string

        :param thres: PValue Threshold to cut the correlations for the output table.
        :type thres: float

    Returns:
        :return output: Output table with the following information: Metabolite "\t" Gene "\t"
        Correlation "\t" pvalue

        :rtype output: file

        :return corMat: Correlation Matrix
        :rtype corMat: file

        :return fig: Network-like output figure
        :rtype fig: pdf
    """

    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    args = getOptions()
    logger = logging.getLogger()
    sl.setLogger(logger)
    logger.info(
        u"Importing data with the following parameters: "
        "\n\tGene Dataset:  {}"
        "\n\tGene UniqueID:  {}"
        "\n\tMet Dataset:{}"
        "\n\tMet UniqueID:  {}"
        "\n\tMethod:  {}"
        "\n\tThreshold:  {}".format(
            args.geneDataset,
            args.geneId,
            args.metDataset,
            args.metId,
            args.meth,
            args.thres,
        )
    )

    modules.checkForDuplicates(args.geneDataset, args.geneId)
    modules.checkForDuplicates(args.metDataset, args.metId)
    pandas2ri.activate()
    with ires.path("gaitGM.data", "all_by_all_correlation.R") as my_r_script_path:
        f = open(my_r_script_path, "r")
        rFile = f.read()
    allByAllCorrScript = STAP(rFile, "corr_main_func")
    # Prepare Gene Expression Data
    geneTable = pd.read_table(args.geneDataset, sep="\t", header=0)
    if args.geneAnnot:
        R_gene_df = modules.Ids2Names(
            geneTable, args.geneId, args.geneAnnot, args.geneName
        )
    else:

        geneTable = geneTable.set_index(args.geneId)
        R_gene_df = pandas2ri.py2rpy(geneTable)

    # Prepare Metabolomics Data
    metTable = pd.read_table(args.metDataset, sep="\t", header=0)
    if args.metAnnot:
        R_met_df = modules.Ids2Names(metTable, args.metId, args.metAnnot, args.metName)
    else:
        metTable = metTable.set_index(args.metId)
        R_met_df = pandas2ri.py2rpy(metTable)

    allByAllCorrScript.corr_main_func(
        x=R_gene_df,
        y=R_met_df,
        meth=args.meth,
        thres=args.thres,
        corrMatPath=args.corMat,
        outputPath=args.output,
        figurePath=args.fig,
    )


if __name__ == "__main__":
    main()
