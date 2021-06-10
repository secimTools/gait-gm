#!/usr/bin/env python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
#
# DESCRIPTION: Download the Kegg information of the selected species and associate gene expression
# and/or metabolomics data with pathways.
#
# VERSION: 1.0
#######################################################################################
import os
import argparse
from argparse import RawDescriptionHelpFormatter
import gaitGM.keggPeaModules as modules


def getOptions():
    parser = argparse.ArgumentParser(
        description="Kegg Downloader", formatter_class=RawDescriptionHelpFormatter
    )

    tool = parser.add_argument_group(title="Tool Specific Inputs")
    tool.add_argument(
        "-sp",
        "--species",
        dest="species",
        action="store",
        required=True,
        help="Species to download.",
    )
    tool.add_argument(
        "-gka",
        "--geneKeggAnnot",
        dest="geneKeggAnnot",
        action="store",
        required=False,
        help="Gene KEGG Annotation File.",
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
        "-gkid",
        "--geneKeggId",
        dest="geneKeggId",
        action="store",
        required=False,
        help="Name of the column with gene KEGG Identifiers.",
    )
    tool.add_argument(
        "-mka",
        "--metKeggAnnot",
        dest="metKeggAnnot",
        action="store",
        required=False,
        help="Metabolite KEGG Annotation File.",
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
    tool.add_argument(
        "-mkid",
        "--metKeggId",
        dest="metKeggId",
        action="store",
        required=False,
        help="Name of the column with Metabolite KEGG Identifiers.",
    )

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
    output.add_argument(
        "-kg2p",
        "--kgen2pathways",
        dest="kgen2pathways",
        action="store",
        required=False,
        help="Gene2Pathway file.",
    )
    output.add_argument(
        "-km2p",
        "--kmet2pathways",
        dest="kmet2pathways",
        action="store",
        required=False,
        help="Metabolite2Pathway file.",
    )
    output.add_argument(
        "-p",
        "--pathways",
        dest="pathways",
        action="store",
        required=True,
        help="PathwaysNames file.",
    )

    args = parser.parse_args()

    # Standardized paths
    if args.geneKeggAnnot:
        args.geneKeggAnnot = os.path.abspath(args.geneKeggAnnot)
        args.kgen2pathways = os.path.abspath(args.kgen2pathways)
        args.geneOut = os.path.abspath(args.geneOut)
    if args.metKeggAnnot:
        args.metKeggAnnot = os.path.abspath(args.metKeggAnnot)
        args.kmet2pathways = os.path.abspath(args.kmet2pathways)
        args.metOut = os.path.abspath(args.metOut)

    args.pathways = os.path.abspath(args.pathways)

    return args


def main():
    """
    Add KEGG Pathway information to each feature (gene/metabolite) that has been previously
    analyzed with add_kegg_anno_info.py.

    Arguments:
        :params geneKeggAnnot metKeggAnnot: Output file from add_kegg_anno_info tool for genes or
        metabolites respectively.
        :types geneKeggAnnot metKeggAnnot: files

        :params geneUniqId metUniqId: Name of the column with Unique identifiers for genes or
        metabolites respectively.
        :types geneUniqId metUniqId: strings

        :params geneName metName: Name of the column with the names of the genes or metabolites
        respectively.
        :types geneName metName: strings

        :params geneKeggId metKeggId: Name of the column with the KEGG Identifiers of the genes or
        metabolites respectively.
        :types geneKeggId metKeggId: strings

    Return:
        :returns geneOut metOut: Output tables with this structure: UniqueID FeatureName
        FeatureType Kegg_ID Pathway_ID Pathway_Name
        :rtypes geneOut metOut: files

        :returns kgen2pathway kmet2pathway: Table with this information: kegg_identifiers "\t"
        pathway_identifiers
        :rtypes kgen2pathway kmet2pathway: files

        :returns pathways: Table with this information: pathway_identifier "\t"
        Pathway_name"-"Specified_organism
        :rtypes pathways: file
    """
    args = getOptions()
    args = modules.downloadKeggInfo(args)
    # Add KEGG Pathway Info for genes
    if args.geneKeggAnnot:
        geneDict, geneList = modules.keggAnnot2list(
            args.geneKeggAnnot,
            args.geneUniqId,
            args.geneName,
            args.geneKeggId,
            featureType="Gene",
        )
        modules.add_path_info(
            geneDict,
            geneList,
            "Gene",
            args.kgen2pathways,
            args.pathways,
            args.species,
            args.geneOut,
        )

    # Add KEGG Pathway Info for Metabolites
    if args.metKeggAnnot:
        metDict, metList = modules.keggAnnot2list(
            args.metKeggAnnot,
            args.metUniqId,
            args.metName,
            args.metKeggId,
            featureType="Metabolite",
        )
        modules.add_path_info(
            metDict,
            metList,
            "Metabolite",
            args.kmet2pathways,
            args.pathways,
            args.species,
            args.metOut,
        )


if __name__ == "__main__":
    main()
