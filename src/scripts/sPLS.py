#!/usr/bin/env python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
# DESCRIPTION: Take a column with the gene and/or metabolites names and return an annotation file
# with KEGG information (Name in KEGG, KEGG ID, ...)
#######################################################################################

import os
import logging
import argparse
from argparse import RawDescriptionHelpFormatter
import matplotlib
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
import gaitGM.keggPeaModules as modules
from secimtools.dataManager import logger as sl
from importlib import resources as ires
matplotlib.use("Agg")


def getOptions():
    parser = argparse.ArgumentParser(
        description="sPLS", formatter_class=RawDescriptionHelpFormatter
    )
    genes = parser.add_argument_group(description="Gene Expression")
    genes.add_argument(
        "-g",
        "--geneDataset",
        dest="geneDataset",
        action="store",
        required=True,
        help="Gene Expression dataset.",
    )
    genes.add_argument(
        "-gid",
        "--geneId",
        dest="geneId",
        action="store",
        required=True,
        help="Gene Unique ID column name.",
    )
    genes.add_argument(
        "-ga",
        "--geneAnno",
        dest="geneAnno",
        action="store",
        required=False,
        help="Gene Expression Annotation File\
                      (Only if user want gene names in output).",
    )
    genes.add_argument(
        "-gn",
        "--geneName",
        dest="geneName",
        action="store",
        required=False,
        help="Gene Names column name in geneAnno File\
                      (Only if user want gene names in output).",
    )
    genes.add_argument(
        "-k",
        "-keepX",
        dest="keepX",
        action="store",
        required=True,
        help="Number of genes to keep in each component.",
    )
    genes.add_argument(
        "-t",
        "-thres",
        dest="thres",
        action="store",
        required=True,
        help="Threshold to cut the sPLS output file.",
    )
    genes.add_argument(
        "-go",
        "--geneOption",
        dest="geneOption",
        action="store",
        required=True,
        help="One of: all, geneList, path, pana.",
    )
    genes.add_argument(
        "-gl",
        "--geneList",
        dest="geneList",
        action="store",
        required=False,
        help="Relevant genes to make a sublist\
                      (Only required in geneList option).",
    )
    genes.add_argument(
        "-gkp",
        "--geneKeggPath",
        dest="geneKeggPath",
        action="store",
        required=False,
        help="Gene Expression Add KEGG Pathway Info File\
                      (Only required in path option).",
    )
    genes.add_argument(
        "-gka",
        "--geneKeggAnno",
        dest="geneKeggAnno",
        action="store",
        required=False,
        help="Gene Expression Add KEGG Annotation Info File\
                      (Only required in pana option).",
    )
    genes.add_argument(
        "-gkn",
        "--geneKeggName",
        dest="geneKeggName",
        action="store",
        required=False,
        help="Name of the column in geneKeggAnno that has Gene Symbols\
                      (Only required in pana option).",
    )
    genes.add_argument(
        "-p2g",
        "--path2genes",
        dest="path2genes",
        action="store",
        required=False,
        help="Pathway2Genes File, from Add KEGG Pathway Info Tool\
                      (Only required in pana option).",
    )
    genes.add_argument(
        "-p2n",
        "--path2names",
        dest="path2names",
        action="store",
        required=False,
        help="PathId2PathNames File, from Add KEGG Pathway Info Tool\
                      (Only required in pana option and if desired).",
    )
    genes.add_argument(
        "-cu",
        "--cutoff",
        dest="cutoff",
        action="store",
        required=False,
        help="Variability cut-off value\
                      Default: 0.2 (Only required in pana option).",
    )
    genes.add_argument(
        "-f",
        "--facSel",
        dest="facSel",
        action="store",
        required=False,
        help="criterion to select components. One of 'accum', 'single', 'abs.val' or 'rel.abs'\
                      Default: 'single'. (Only required in pana option).",
    )

    metabolites = parser.add_argument_group(description="Metabolites")
    metabolites.add_argument(
        "-m",
        "--metDataset",
        dest="metDataset",
        action="store",
        required=True,
        help="Metabolomic Datset.",
    )
    metabolites.add_argument(
        "-mid",
        "--metId",
        dest="metId",
        action="store",
        required=True,
        help="Metabolite Unique ID column name.",
    )
    metabolites.add_argument(
        "-ma",
        "--metAnno",
        dest="metAnno",
        action="store",
        required=False,
        help="Metabolomics Annotation File \
             (Only if user want metabolite names in output).",
    )
    metabolites.add_argument(
        "-mn",
        "--metName",
        dest="metName",
        action="store",
        required=False,
        help="Metabolite Names column name \
             (Only if user want metabolite names in output, or in metOption=generic or both).",
    )
    metabolites.add_argument(
        "-mo",
        "--metOption",
        dest="metOption",
        action="store",
        required=True,
        help="One of generic, mmc or both.",
    )
    metabolites.add_argument(
        "-mka",
        "--metKeggAnno",
        dest="metKeggAnno",
        action="store",
        required=False,
        help="Metabolomic Add KEGG Annotation Info File. (Only required if metOption != mmc).",
    )
    metabolites.add_argument(
        "-mkp",
        "--metKeggPath",
        dest="metKeggPath",
        action="store",
        required=False,
        help="Metabolomic Add KEGG Pathway Info File\
                      (Only if geneOption = path).",
    )
    metabolites.add_argument(
        "-d",
        "--design",
        dest="design",
        action="store",
        required=False,
        help="Design file (Only required in mmc or both option).",
    )

    mmc = parser.add_argument_group(description="MMC")
    mmc.add_argument(
        "-c",
        "--correlation",
        dest="correlation",
        default="pearson",
        choices=("pearson", "kendall", "spearman"),
        required=False,
        help=(
            "Compute correlation coefficients using either "
            "'pearson' (standard correlation coefficient), "
            "'kendall' (Kendall Tau correlation coefficient), or "
            "'spearman' (Spearman rank correlation)."
        ),
    )
    mmc.add_argument(
        "-sl",
        "--sigmaLow",
        dest="sigmaLow",
        type=float,
        required=False,
        default=0.05,
        help="Low value of sigma (Default: 0.05).",
    )
    mmc.add_argument(
        "-sh",
        "--sigmaHigh",
        dest="sigmaHigh",
        type=float,
        required=False,
        default=0.50,
        help="High value of sigma (Default: 0.50).",
    )
    mmc.add_argument(
        "-sn",
        "--sigmaNum",
        dest="sigmaNum",
        type=int,
        required=False,
        default=451,
        help="Number of values of sigma to search" " (Default: 451).",
    )

    plot = parser.add_argument_group(title="Plot options")
    plot.add_argument(
        "-pal",
        "--palette",
        dest="palette",
        action="store",
        required=False,
        default="diverging",
        help="Name of the palette to use.",
    )
    plot.add_argument(
        "-col",
        "--color",
        dest="color",
        action="store",
        required=False,
        default="Spectral_10",
        help="Name of a valid color scheme" " on the selected palette",
    )

    output = parser.add_argument_group(description="Output")
    output.add_argument(
        "-f1",
        "--figure1",
        dest="figure1",
        action="store",
        required=True,
        help="sPLS heatmaps",
    )
    output.add_argument(
        "-o1",
        "--splsOut",
        dest="splsOut",
        action="store",
        required=True,
        help="Output Table.",
    )
    output.add_argument(
        "-f2",
        "--figure2",
        dest="figure2",
        action="store",
        required=False,
        help="MMC Heatmaps (Only if MMC Option)",
    )
    output.add_argument(
        "-o2",
        "--mmcOut",
        dest="mmcOut",
        action="store",
        required=False,
        help="MMC Output TSV name (Only if MMC Option)",
    )
    output.add_argument(
        "-o3",
        "--panaOut",
        dest="panaOut",
        action="store",
        required=False,
        help="PANA Output TSV name (Only if PANA Option)",
    )

    args = parser.parse_args()

    args.geneDataset = os.path.abspath(args.geneDataset)
    if args.geneAnno:
        args.geneAnno = os.path.abspath(args.geneAnno)
    if args.geneOption == "geneList":
        args.geneList = os.path.abspath(args.geneList)
    elif args.geneOption == "path":
        args.geneKeggPath = os.path.abspath(args.geneKeggPath)
        args.metKeggPath = os.path.abspath(args.metKeggPath)
    elif args.geneOption == "pana":
        args.geneKeggAnno = os.path.abspath(args.geneKeggAnno)
        args.path2genes = os.path.abspath(args.path2genes)
        args.panaOut = os.path.abspath(args.panaOut)
        if args.path2names:
            args.path2names = os.path.abspath(args.path2names)

    args.metDataset = os.path.abspath(args.metDataset)
    if args.metAnno:
        args.metAnno = os.path.abspath(args.metAnno)
    if args.metOption != "generic":
        args.design = os.path.abspath(args.design)
        args.mmcOut = os.path.abspath(args.mmcOut)
        args.figure2 = os.path.abspath(args.figure2)
    elif args.metOption != "mmc":
        args.metKeggAnno = os.path.abspath(args.metKeggAnno)

    args.figure1 = os.path.abspath(args.figure1)
    args.splsOut = os.path.abspath(args.splsOut)

    return args


def main():
    """
    Performs a Sparse Partial Least Squares (sPLS) analysis over subsets of gene expression and
    metabolomic data. To perform this subsetting, three different methodologies can be used for the
    metabolites:
    - By generic metabolite (sphingomyelin, ...)
    - By MMC cluster
    - By generic metabolite and then by MMC cluster
    and four for the genes:
    - All the genes
    - Genes contained in a list with interesting genes for the analysis
    - Pathway related genes for an specific generic metabolite
    - Metagenes (PANA approach)

    The outputs depend on the inputs.

    Arguments:
        :param geneDataset metDataset: Gene expression/Metabolomics wide dataset, respectively.
        :type geneDataset metDataset: files

        :param geneId metId: Name of the Genes/metabolites unique identifier column, respectively.
        :type geneId metId: strings

        :param geneAnnot metAnnot: Gene Expression/Metabolomics Annotation Dataset.
        :type geneAnnot metAnnot: files

        :param geneAnnotName metAnnotName: annotation file column with gene/metabolite names.
        :type geneAnnotName metAnnotName: strings

        :param design: Design File
        :type design: file

        :param keepX: Number of genes to keep in the sPLS model
        :param keepX: integer

        :param geneOption metOption: Options for metabolite subsetting (one of 'generic', 'mmc' or
        'both') and for gene expression subsetting (one of 'all', 'geneList', 'path' or 'pana')
        :type geneOption metOption: strings

        :param geneKeggAnno metKeggAnno: KEGG Annotation files for gene expression and
        metabolomics, respectively. From Add KEGG Anno Info Tool
        :type geneKeggAnno metKeggAnno: files

        :param geneKeggPath metKeggPath: KEGG Pathway files for gene expression and metabolomics,
        respectively. From Add KEGG Pathway Info Tool
        :type geneKeggPath metKeggPath: files

        :param path2genes: Downloaded KEGG file with this information: pathway_ID "\t" geneKEGG_ID
        :type path2genes: file

    Returns:
        :return figure1: sPLS heatmaps
        :rtype figure1: pdf

        :return splsOut: sif-like correlation matrix including a column describing the comparison.
        :rtype splsOut: file

        :return figure2: MMC plots if mmc or both metabolite subsetting option is selected.
        :rtype figure2: pdf

        :return mmcOut: MMC Output table if mmc or both metabolite subsetting option is selected.
        :rtype mmcOut: file

        :return panaOut: Table describing genes that forms the metagenes (1/0)
        :rtype panaOut: file
    """
    args = getOptions()
    logger = logging.getLogger()
    sl.setLogger(logger)
    logger.info(
        u"Importing data with the following parameters: "
        "\n\tGene Dataset:  {}"
        "\n\tGene UniqueID:  {}"
        "\n\tGene Option: {}"
        "\n\tMetabolite Dataset:{}"
        "\n\tMetabolite UniqueID:  {}"
        "\n\tMetabolite Option: {}".format(
            args.geneDataset,
            args.geneId,
            args.geneOption,
            args.metDataset,
            args.metId,
            args.metOption,
        )
    )
    pandas2ri.activate()
    with ires.path("gaitGM.data", "sPLS.R") as my_r_script_path:
        f = open(my_r_script_path, "r")
        rFile = f.read()
    sPLSScript = STAP(rFile, "sPLS")
    rGeneData, rMetData, multipleNames, multipleNamesId = modules.prepareSPLSData(args)
    rData = []
    data_counter = 0
    for R_met_df in rMetData:
        R_gene_df = rGeneData[data_counter]
        rData.append(
            sPLSScript.sPLS(geneData=R_gene_df, metData=R_met_df, keepX=args.keepX)
        )
        if args.geneOption == "path":
            data_counter += 1
    if args.metOption == "both":
        sPLSScript.plotInPdf(
            splsObjects=rData, figurePath=args.figure1, multipleNames=multipleNamesId
        )
        # Correlation Matrix
        corMatrix = sPLSScript.corrMat(
            splsObjects=rData, multipleNames=multipleNamesId, threshold=args.thres
        )
        robjects.r["write.table"](
            corMatrix,
            file=args.splsOut,
            sep="\t",
            quote=False,
            row_names=False,
            col_names=True,
        )
    else:
        sPLSScript.plotInPdf(
            splsObjects=rData, figurePath=args.figure1, multipleNames=multipleNames
        )
        # Correlation Matrix
        corMatrix = sPLSScript.corrMat(
            splsObjects=rData, multipleNames=multipleNames, threshold=args.thres
        )
        robjects.r["write.table"](
            corMatrix,
            file=args.splsOut,
            sep="\t",
            quote=False,
            row_names=False,
            col_names=True,
        )


if __name__ == "__main__":
    main()
