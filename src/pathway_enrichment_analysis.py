#!/usr/bin/env python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
#
# DESCRIPTION: Perform Pathway Enrichment Analysis
#
# VERSION: 1.0
#######################################################################################
import argparse
from argparse import RawDescriptionHelpFormatter
import os
import keggPeaModules as modules
import relations_enrichment as RE
import scipy.stats as st


def getOptions():
    parser = argparse.ArgumentParser(
        description="Pathway Enrichment", formatter_class=RawDescriptionHelpFormatter
    )
    standard = parser.add_argument_group(description="Input")
    standard.add_argument(
        "-k",
        "--kegg",
        dest="kegg_downloader",
        action="store",
        required=True,
        help="Kegg Downloader Output File.",
    )
    standard.add_argument(
        "-s",
        "--steps",
        dest="steps_away",
        action="store",
        required=True,
        help="Number of steps away from desired meatbolite.",
    )

    standard.add_argument(
        "-g",
        "--gene",
        dest="deaGeneDataset",
        action="store",
        required=True,
        help="Differential Expression Analysis Gene Dataset.",
    )
    standard.add_argument(
        "-gid",
        "--gene_id",
        dest="gene_id_col",
        action="store",
        required=True,
        help="DEA Gene Dataset Unique Identificator column.",
    )
    standard.add_argument(
        "-gf",
        "--gene_flag",
        dest="gene_flag_col",
        action="store",
        required=True,
        help="DEA Gene Dataset Flag column.",
    )

    standard.add_argument(
        "-m",
        "--met",
        dest="deaMetDataset",
        action="store",
        required=True,
        help="Differential Expression Analysis Metabolomic Dataset.",
    )
    standard.add_argument(
        "-mid",
        "--met_id",
        dest="met_id_col",
        action="store",
        required=True,
        help="DEA Metabolomic Dataset Unique Identificator column.",
    )
    standard.add_argument(
        "-mf",
        "--met_flag",
        dest="met_flag_col",
        action="store",
        required=True,
        help="DEA Metabolomic Dataset Flag column.",
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

    args = parser.parse_args()

    args.kegg_downloader = os.path.abspath(args.kegg_downloader)
    args.deaGeneDataset = os.path.abspath(args.deaGeneDataset)
    args.deaMetDataset = os.path.abspath(args.deaMetDataset)

    args.output = os.path.abspath(args.output)

    return args


def createDict(args):
    """
    Create a dictionary with the association between different pathways and the related features,
    genes and metabolites present in the study.

    Arguments:
        :param kegg_downloader: Table with kegg_downloader program output information
        :type kegg_downloader: file

        :params deaGeneDataset deaMetDataset: Tables with Differential Expression Analysis information for gene expression
            and metabolomics, respectively
        :types deaGeneDataset deaMetDataset: files

        :params gene_id_col, met_id_col, gene_flag_col, met_flag_col: Column names of unique identifiers and desired
            flag column of gene expression and metabolomics datasets, respectively
        :type gene_id_col, met_id_col, gene_flag_col, met_flag_col: strings

    Returns:
        :return output: Table with this structure: Pathway_Name Odds_Ratio	P_value
        :rtype output: file
    """

    print("# Step 1) Creating dictionary {pathway: [[Genes], [Metabolites]]}\n")

    path_feat = {}
    pathways = []
    pathNames_feat = {}
    pathwaysNames = {}

    kegg_downloader = open(args.kegg_downloader, "r")

    # Total of different pathways
    kegg_downloader.readline()
    for line in kegg_downloader:
        pathway = line.split("\t")[3].strip()
        pathwayName = line.split("\t")[4].strip()
        if pathway not in pathways:
            pathways.append(pathway)
            pathwaysNames[pathway] = pathwayName

    # Dictionary {pathway: [[Genes], [Metabolites]]}
    kegg_downloader.readline()
    for pathway in pathways:
        gene_list = []
        geneName_list = []
        met_list = []
        metName_list = []
        kegg_downloader.seek(0)
        for line in kegg_downloader:
            if pathway == line.split("\t")[3].strip():
                feature = line.split("\t")[2].strip()
                featureType = line.split("\t")[1].strip()
                featureName = line.split("\t")[0].strip()
                if featureType == "Gene":
                    gene_list.append(feature)
                    geneName_list.append(featureName)
                else:
                    met_list.append(feature)
                    metName_list.append(featureName)

        path_feat[pathway] = [gene_list, met_list]
        pathNames_feat[pathwaysNames[pathway]] = [geneName_list, metName_list]

    with open("1-path_feat.txt", "w") as step1_path_feat:
        step1_path_feat.write(str(path_feat))
    with open("1-pathNames_feat.txt", "w") as step1_pathNames_feat:
        step1_pathNames_feat.write(str(pathNames_feat))

    kegg_downloader.close()

    return (args, path_feat, pathNames_feat, pathwaysNames)


def pathwayFisher(genes1step, step):
    """
    Perform a Fisher Exact test.
    """
    print("# Step 4) Performing Fisher Exact test\n")
    numberOfItems = {}
    output = open("output" + str(step) + ".txt", "w")
    output.write(
        "Genes"
        + str(step)
        + "stepAway\tDEG_with_DEM\tNoDEG_with_DEM\tDEG_with_NoDEM\tNoDEG_with_NoDEM\tPvalue\tOddsRatio\n"
    )

    for metabolite, genes in genes1step.iteritems():
        try:
            numberOfItems[len(genes)].append(metabolite)
        except KeyError:
            numberOfItems[len(genes)] = [metabolite]

    for numberOfGenes, metabolites in numberOfItems.iteritems():
        DEM_DEG, NoDEM_DEG, DEM_NoDEG, NoDEM_NoDEG = 0, 0, 0, 0
        for metabolite in metabolites:
            if "*" in metabolite:
                for gene in genes1step[metabolite]:
                    if "*" in gene:
                        DEM_DEG += 1
                    else:
                        DEM_NoDEG += 1
            else:
                for gene in genes1step[metabolite]:
                    if "*" in gene:
                        NoDEM_DEG += 1
                    else:
                        NoDEM_NoDEG += 1
        oddsratio, pvalue = st.fisher_exact(
            [[DEM_DEG, DEM_NoDEG], [NoDEM_DEG, NoDEM_NoDEG]]
        )
        output.write(
            str(numberOfGenes)
            + "\t"
            + str(DEM_DEG)
            + "\t"
            + str(DEM_NoDEG)
            + "\t"
            + str(NoDEM_DEG)
            + "\t"
            + str(NoDEM_NoDEG)
            + "\t"
            + str(pvalue)
            + "\t"
            + str(oddsratio)
            + "\n"
        )

    return ()


def main():
    args = getOptions()
    args, path_feat, pathNames, pathwaysNames = createDict(args)
    steps = int(args.steps_away)
    for step in range(steps):
        step += 1
        print("\t#", step, "STEPS AWAY #")
        path_feat, genes1step = RE.relationsFinder(args, path_feat, pathwaysNames)
        with open("2-step" + str(step) + "-path_feat.txt", "w") as step2_path_feat:
            step2_path_feat.write(str(path_feat))
        with open("2-step" + str(step) + "-genes1step.txt", "w") as step2_genes1step:
            step2_genes1step.write(str(genes1step))
        stepsTable = RE.stepsFrequency(args, genes1step, pathwaysNames)
        with open("3-step" + str(step) + "cpdDict.txt", "w") as step3_cpdDict:
            step3_cpdDict.write(str(stepsTable))
        pathwayFisher(stepsTable, step)


if __name__ == "__main__":
    main()
