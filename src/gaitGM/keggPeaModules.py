#!/home/moskalenko/projects/software/bioinformatics/secim/GAIT-GM/_conda/1.0.0/bin/python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
# DESCRIPTION: gaitGM module with functions for KEGG PEA Tools
#######################################################################################

import re
import sys
import csv
import requests
import logging
import tempfile
from difflib import SequenceMatcher
import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as st
from rpy2 import robjects as robjects
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
from matplotlib.backends.backend_pdf import PdfPages
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
from secimtools.dataManager import logger as sl
from secimtools.dataManager.interface import wideToDesign
from secimtools.visualManager.manager_color import colorHandler
from secimtools.visualManager.manager_figure import figureHandler
from secimtools.visualManager.module_mmc import expansion, get_clustering
from importlib import resources as ires


def checkForDuplicates(dataset, uniqID):
    """
    Check for duplicated values in a dataset column. If found, terminate the program and provide
    information to stderr

    Arguments:
        :param dataset: Input dataset to check if uniqID column is unique
        :type dataset: file

        :param uniqID: Unique Identifier Column Name
        :type uniqID: string
    """

    table = pd.read_table(dataset, sep="\t", header=0)

    # In case of a Selected Yes/No Column:
    if "Selected" in table.columns:
        table = table.loc[table["Selected"] != "No"]

    ids = table[uniqID]
    duplicates = table[ids.isin(ids[ids.duplicated()])][uniqID]

    if not duplicates.empty:
        pd.Series.__unicode__ = pd.Series.to_string
        sys.stderr.write(
            "Duplicated IDs!\n\nRow_Number\tRepeated_ID\n" + str(duplicates)
        )
        sys.exit(1)


def downloadGeneParser(species):
    """
    Add KEGG Annotation Info (kegg_id, gene_name) from KEGG database.

    Arguments:
        :param species: species identifier in kegg to download
        :type species: string

    Returns:
        :return gene2keggsArray: Dictionary with KEGG information about gene name and identifier
        :rtype gene2keggsArray: Dictionary
    """

    gene2keggsArray = {}
    with requests.get("http://rest.kegg.jp/list/" + species) as genes:
        gene2keggs = genes.content.splitlines()
    for line in gene2keggs:
        line = line.decode('utf-8')
        geneId = line.split("\t")[0]
        geneNames = line.split("\t")[1]
        gene2keggsArray[geneId] = geneNames

    return gene2keggsArray


def downloadMetParser():
    """
    Download metabolite information (kegg_id, cpd_name) from KEGG database.

    Returns:
        :return met2keggsArray: Dictionary with KEGG information about compound name and identifier
        :rtype met2keggsArray: dictionary
    """

    met2keggsArray = {}
    with requests.get("http://rest.kegg.jp/list/compound") as metabolites:
        met2keggs = metabolites.content.splitlines()
    for line in met2keggs:
        line = line.decode('utf-8')
        cpdId = line.split("\t")[0]
        cpdNames = line.split("\t")[1]
        met2keggsArray[cpdId] = cpdNames

    return met2keggsArray


def keggAnno(feature_table, feature2keggs, featureOut, uniqIDNameCol, featureNameCol, featureType):
    """
    Takes the Name column of a Dataset and find KEGG-related information.
    It works for Gene Expression and Metabolomic data (depending on featureType parameter).
    It will return an Annotation Dataset with this information and list possible ties.

    Arguments:
        :param feature_table: Table that contains at least a column with a Unique Identifier and
        another one with the feature name.
        :type feature_table: pandas dataset

        :param feature2keggs: Dictionary created with KEGG Information: {KEGG_ID: FeatureName}
        :type feature2keggs: dictionary

        :param featureOut: Initializated output table with this information: UniqId FeatureName
        Feature_Type Matched Name_in_KEGG KEGG_ID Similarity Tie Selected
        :type featureOut: file

        :param uniqIDNameCol: Name of the column with Unique Identifiers.
        :type uniqIDNameCol: string

        :param featureNameCol: Name of the column with feature names.
        :type featureNameCol: string

        :param featureType: One of: 'Gene' or 'Metabolite'.
        :type featureType: string

    Returns:
        :return featureOut: Finished output table
        :rtype featureOut: File
    """

    tmpList_v1 = []
    tmpList_v2 = []
    emptyList = ["nan", "", "na"]

    for index, featureRow in feature_table.iterrows():
        uniqueID = str(featureRow[uniqIDNameCol]).strip()
        featureName = str(featureRow[featureNameCol]).strip()
        featureName = featureName.replace("\t", "_")
        if (featureType == "Metabolite") and (featureName not in emptyList):
            newMetNames = metaboliteModification(featureName)
        # To not match NaN to KEGG
        if featureName.lower() in emptyList:
            tmpList_v1.append(
                uniqueID +
                "\t" +
                featureName +
                "\t" +
                featureType +
                "\tNA\tNA\tNA\tNA\tNA\tNA\n"
            )
        else:
            featureDict_v1 = {}
            # metabolite dict: {keggMetName1: [similarity1, kegg_cpd1], keggMetName2: [similarity2,
            # kegg_cpd2], ...}
            for kegg_id, keggFeatureNames in feature2keggs.items():
                if featureType == "Gene":
                    keggGeneNames = keggFeatureNames.split(";")[0].split(",")
                    for keggGeneName in keggGeneNames:
                        if re.search(
                            ".*" + re.escape(str(featureName)) + ".*",
                            keggGeneName,
                            re.IGNORECASE,
                        ):
                            similarity = calculateSimilarity(
                                str(featureName).strip(), str(keggGeneName).strip()
                            )
                            featureDict_v1[str(keggGeneName).strip()] = [
                                similarity,
                                kegg_id,
                            ]
                else:
                    for newMetName in newMetNames:
                        if re.search(
                            ".*" + re.escape(newMetName) + ".*",
                            keggFeatureNames,
                            re.IGNORECASE,
                        ):
                            featureDict_v1 = add2Dictionary(
                                featureName,
                                newMetName,
                                kegg_id,
                                keggFeatureNames,
                                featureDict_v1,
                            )
            if featureDict_v1:
                # Solving ties and Sort by similarity
                sortedMetDict = sorted(featureDict_v1.values(), reverse=True)
                maximum_v1 = sortedMetDict[0]
                if len(featureDict_v1) > 1:
                    maximum_v2 = sortedMetDict[1]
                    if (maximum_v1[0] == maximum_v2[0]) or (
                        maximum_v1[0] - maximum_v2[0] <= 0.05 * maximum_v1[0]
                    ):
                        # Warning Message
                        print(
                            "Warning! There is a tie with",
                            featureName + ":",
                            list(featureDict_v1.keys())[
                                list(featureDict_v1.values()).index(maximum_v1)
                            ],
                            "selected,",
                            list(featureDict_v1.keys())[
                                list(featureDict_v1.values()).index(maximum_v2)
                            ],
                            "rejected.",
                        )
                        is_tie = "Yes"
                    else:
                        is_tie = "No"
                else:
                    is_tie = "No"
            # UniqueID FeatureName featureType Matched KEGG_Name Kegg_cpd Similarity Tie  elected
                tmpList_v1.append(
                    uniqueID +
                    "\t" +
                    featureName +
                    "\t" +
                    featureType +
                    "\t" +
                    "Yes" +
                    "\t" +
                    str(
                        list(featureDict_v1.keys())[
                            list(featureDict_v1.values()).index(maximum_v1)
                        ]
                    ) +
                    "\t" +
                    str(maximum_v1[1]) +
                    "\t" +
                    str(round(maximum_v1[0], 2)) +
                    "\t" +
                    is_tie +
                    "\tYes\n")
                featureDict_v2 = dict(featureDict_v1)
                del featureDict_v2[
                    list(featureDict_v2.keys())[
                        list(featureDict_v2.values()).index(maximum_v1)
                    ]
                ]
                for keggName in featureDict_v2:
                    if is_tie == "Yes":
                        if (
                            keggName == list(featureDict_v1.keys())[
                                list(featureDict_v1.values()).index(maximum_v2)
                            ]
                        ):
                            tmpList_v1.append(
                                uniqueID +
                                "\t" +
                                featureName +
                                "\t" +
                                featureType +
                                "\tYes\t" +
                                str(keggName.strip()) +
                                "\t" +
                                str(featureDict_v2[keggName][1]) +
                                "\t" +
                                str(round(featureDict_v2[keggName][0], 2)) +
                                "\t" +
                                is_tie +
                                "\tNo\n"
                            )
                        else:
                            tmpList_v2.append(
                                uniqueID +
                                "\t" +
                                featureName +
                                "\t" +
                                featureType +
                                "\tYes\t" +
                                str(keggName.strip()) +
                                "\t" +
                                str(featureDict_v2[keggName][1]) +
                                "\t" +
                                str(round(featureDict_v2[keggName][0], 2)) +
                                "\tNo\tNo\n"
                            )
                    else:
                        tmpList_v2.append(
                            uniqueID +
                            "\t" +
                            featureName +
                            "\t" +
                            featureType +
                            "\tYes\t" +
                            str(keggName.strip()) +
                            "\t" +
                            str(featureDict_v2[keggName][1]) +
                            "\t" +
                            str(round(featureDict_v2[keggName][0], 2)) +
                            "\t" +
                            is_tie +
                            "\tNo\n")
            else:
                tmpList_v1.append(
                    uniqueID +
                    "\t" +
                    featureName +
                    "\t" +
                    featureType +
                    "\tNo\tNA\tNA\tNA\tNA\tNA\n"
                )

    # Writing results - first selected metabolites, last non-selected metabolites.
    for line in tmpList_v1:
        featureOut.write(line)
    for line2 in tmpList_v2:
        featureOut.write(line2)

    return featureOut


def metaboliteModification(metabolite):
    """
    Modify input metabolite name in order to find its kegg identifier. It will be modified
    following this pipeline:
    - Remove common metabolite prefixes (d-, alpha-, ...)
    - Use long name if aminoacid abreviation provided (cys = cysteine, ...)
    - Use long name if metabolite abreviation provided (orn = ornithine, ...)
    - Change -ic acid ending by -ate
    - Remove common chemical words that can worsen the matching (sulfoxide, ...)
    - Use long name if lipid abreviation provided (pc = phosphatidylcholine, ...)
    - If everything fails, use input metabolite name

   Arguments:
        :param metabolite: Input metabolite name
        :type metabolite: string

    Returns:
        :return newMetNames: Modified metabolite names in a list.
        :rtype newMetNames: list
    """

    mainPrefixes = [
        re.compile(r"^n?(-|\s)?[0-9]?(-|\s)?.*yl(-|\s)"),
        re.compile(r"^(((l|d|dl|ld)|[0-9])(-|\s))?tert(-|\s)?"),
        re.compile(r"^[0-9]?(-|\s)?hydroxy(-|\s)?"),
        re.compile(r"^n?(-|\s)?[0-9]?(-|\s)?boc(-|\s)?"),
        re.compile(r"\(.*\)"),
        re.compile(r"^n?(-|\s)?[0-9]?(-|\s)?acetyl(-|\s)?"),
        re.compile(
            r"^((cis|trans)?(-|\s))?((s|r|\(s\)|\(r\))?(-|\s))?((([0-9]?)(,?))+?(-|\s))?(n(-|\s))?((d|l|dl|ld)(-|\s))?((alpha|beta|a|b)(-|\s))?(([0-9]?)(,?))+?(-|\s)?"
        ),
    ]
    # Add as many chemical words as desired
    chemWords = [
        "sulfoxide",
        "-sulfoxide",
        "sulfate",
        "-sulfate",
        "sulfonate",
        "-sulfonate",
    ]
    # Commonly abbreviated words, like aminoacids
    aminoacids = {
        "l-cysteine": ["c", "cys"],
        "l-aspartate": ["d", "asp"],
        "l-serine": ["s", "ser"],
        "l-glutamine": ["q", "gln"],
        "l-lysine": ["k", "lys"],
        "l-isoleucine": ["i", "ile"],
        "l-proline": ["p", "pro"],
        "l-threonine": ["t", "thr"],
        "l-phenylalanine": ["f", "phe"],
        "l-asparagine": ["n", "asn"],
        "glycine": ["g", "gly"],
        "l-histidine": ["h", "his"],
        "l-leucine": ["l", "leu"],
        "l-arginine": ["r", "arg"],
        "l-tryptophan": ["w", "trp"],
        "l-alanine": ["a", "ala"],
        "l-valine": ["v", "val"],
        "l-glutamate": ["e", "glu"],
        "l-tyrosine": ["y", "tyr"],
        "l-methionine": ["m", "met"],
    }
    abbreviations = {
        "citrate": "cit",
        "ornithine": "orn",
        "thyroxine": "thyr",
        "butoxycarbonyl": "boc",
    }
    lipids = {
        "sphingomyelin": "sm",
        "lysophosphatidylcholine": "lysopc",
        "phosphatidylcholine": "pc",
        "phosphatidylethanolamine": "pe",
        "lysophosphatidylethanolamine": "lysope",
    }

    metabolite = metabolite.lower()
    newMetNames = []
    newMetNamesNoPrefix = []
    # If metabolite is an abbreviation of an aminoacid
    if (
        (metabolite in dict(tuple(aminoacids.values())).keys()) or
        (metabolite in dict(tuple(aminoacids.values())).values()) or
        (metabolite in str(tuple(aminoacids.keys())))
    ):
        for aminoacid in aminoacids.keys():
            if (metabolite in aminoacids[aminoacid]) and (aminoacid not in newMetNames):
                newMetNames.append(aminoacid)
            elif (aminoacid.endswith(metabolite)) and (aminoacid not in newMetNames):
                newMetNames.append(aminoacid)
    # If metabolite is an abbreviation of another commonly abbreviated metabolites
    if any(completeName in metabolite for completeName in abbreviations.values()):
        for completeName in abbreviations.keys():
            if abbreviations[completeName] in metabolite:
                newMetNames.append(
                    re.sub(abbreviations[completeName], completeName, metabolite)
                )
                newMetNames.append(metabolite)
    # If metabolite is an acid, use the -ate nomenclature
    if (
        ("ic acid" in metabolite) or
        ("ic_acid" in metabolite) or
        ("icacid" in metabolite)
    ) and (re.sub(r"ic.?acid", "ate", metabolite) not in newMetNames):
        newMetNames.append(re.sub(r"ic.?acid", "ate", metabolite))
        newMetNames.append(re.sub(r"ic.?acid", "ic acid", metabolite))
    # If metabolite contains any chemical word, remove it
    if any(chemWord in metabolite for chemWord in chemWords):
        for chemWord in chemWords:
            if (chemWord in metabolite) and (
                re.sub(chemWord, "", metabolite) not in newMetNames
            ):
                newMetNames.append(re.sub(chemWord, "", metabolite))
                newMetNames.append(metabolite)
    # Search in pubchem for synonyms - skipped, but retained for the future
    #    if (pcp.get_synonyms(metabolite, 'name', 'compound')) and (not newMetNames):
    #        synonymsList = pcp.get_synonyms(metabolite, 'name', 'compound')
    #        for dictionarySynonyms in synonymsList:
    #            synonyms = dictionarySynonyms[u'Synonym']
    #            for synonym in synonyms:
    #                # similarity > 0.2 to increase the speed
    #                if (str(synonym) not in newMetNames) and
    #                   (SequenceMatcher(a=metabolite, b=synonym).ratio() > 0.2):
    #                    newMetNames.append(str(synonym))
    # If metabolite is an abbreviation of a lipid
    if any(lipid in metabolite for lipid in lipids.values()):
        for lipid in lipids.keys():
            if (metabolite.startswith(lipids[lipid])) and (lipid not in newMetNames):
                newMetNames.append(lipid)
    # If everything fails, take original name
    if not newMetNames:
        newMetNames.append(metabolite.strip())
    # Check for prefixes
    for newMetName in newMetNames:
        for mainPrefix in mainPrefixes:
            if (
                mainPrefix.search(newMetName) and
                (re.sub(mainPrefix, "", newMetName) not in newMetNamesNoPrefix) and
                (newMetName not in aminoacids)
            ):
                newMetNameNoPrefix = re.sub(mainPrefix, "", newMetName)
                newMetNames.insert(len(newMetNames), newMetNameNoPrefix)
                newMetNamesNoPrefix.append(newMetNameNoPrefix)
        if newMetName not in newMetNamesNoPrefix:
            newMetNamesNoPrefix.append(newMetName)

    return newMetNamesNoPrefix


def add2Dictionary(metabolite, newMetName, kegg_cpd, keggMetNames, metDict):
    """
    Compare input metabolite name with KEGG name.

    If newMetName (metabolite without prefix) == keggMetName:
        Calculate similarity between metabolite (oirginal name) and keggMetName

    It creates different dictionaries for each metabolite (one per match)

     Arguments:
        :param metabolite: Input metabolite name
        :type metabolite: string

        :param newMetName: Input metabolite name without prefix and chemical words.
        :type newMetName: string

        :param kegg_cpd: Kegg identifier of the metabolite
        :type kegg_cpd: string

        :param keggMetNames: "List" of Kegg metabolite names associated to a Kegg identifier
        :type keggMetNames: string

    Returns:
        :return metDict: Output dictionary with this structure: {Metabolite_Name_in_Kegg:
        [similarity, Kegg_compound_identifier]
        :rtype metDict: dictionary
    """

    keggMetNames = keggMetNames.split(";")
    for keggMetName in keggMetNames:
        if re.search(
            ".*" + re.escape(newMetName) + ".*", keggMetName.strip(), re.IGNORECASE
        ):
            similarity = calculateSimilarity(metabolite.strip(), keggMetName.strip())
            metDict[keggMetName.strip()] = [similarity, kegg_cpd]

    return metDict


def calculateSimilarity(featureName, keggName):
    """
    Compare two feature (gene/metabolite) names and return the similarity between them. If the only
    difference between the names is one of the mainPrefixes a similarity of 90% is returned.

     Arguments:
        :param metabolite: Input metabolite name
        :type metabolite: string

        :param keggName: Name to check the similarity with
        :type keggName: string

    Returns:
        :return similarity: Percentage of similarity between 2 input names.
        :rtype similarity: float
    """

    mainPrefixes = [
        "",
        "cis-",
        "trans-",
        "d-",
        "l-",
        "(s)-",
        "alpha-",
        "beta-",
        "alpha ",
        "beta ",
        "alpha-d-",
        "beta-d-",
        "alpha-l-",
        "beta-l-",
        "l-beta-",
        "l-alpha-",
        "d-beta-",
        "d-alpha-",
    ]

    if featureName == keggName:
        similarity = 1.0
    elif featureName.lower() == keggName.lower():
        similarity = 0.9
    elif (featureName.lower() in mainPrefixes) or (keggName.lower() in mainPrefixes):
        similarity = SequenceMatcher(a=featureName.lower(), b=keggName.lower()).ratio()
    elif keggName.lower().replace(featureName.lower(), "") in mainPrefixes:
        similarity = 0.9
    elif featureName.lower().replace(keggName.lower(), "") in mainPrefixes:
        similarity = 0.9
    else:
        similarity = SequenceMatcher(a=featureName, b=keggName).ratio()

    return similarity


def downloadKeggInfo(args):
    """
    Download necessary information from Kegg Database for parsing.

    Arguments:
        :param geneKeggAnnot: Gene to KEGG ID Link file
        :type geneKeggAnnot: file

        :param metKeggAnnot: Metabolite to KEGG ID Link file
        :type metKeggAnnot: file

    Returns:
        :return gen2kegg: kegg_gene_identifier "\t" Gene_Symbol ";" Gene_name
        :rtype gen2kegg: file

        :return kgen2pathway: kegg_gene_identifier "\t" pathway_identifier_for_gene
        :rtype kgen2pathway: file

        :return met2kegg: kegg_metabolite_identifier "\t" Metabolite_names_list_sep_;
        :rtype met2kegg: file

        :return kmet2pathway: kegg_metabolite_identifier "\t" pathway_identifier_for_metabolite
        :rtype kmet2pathway:similarity file

        :return pathways: pathway_identifier_for_gene "\t" Pathway_name "-" Specified_organism
        :rtype pathways: file
    """

    # GeneKeggID2PathwayID
    if args.geneKeggAnnot:
        geneKeggAnnot = requests.get(
            "http://rest.kegg.jp/link/" + args.species + "/pathway"
        )
        geneKeggAnnot_file = args.species + "_geneKeggAnnot"
        with open(geneKeggAnnot_file, 'w') as fh:
            fh.write(geneKeggAnnot.content.decode("utf-8"))
        args.kgen2pathways = geneKeggAnnot_file

    # MetaboliteKeggID2PathwayID
    if args.metKeggAnnot:
        metKeggAnnot = requests.get(
            "http://rest.kegg.jp/link/compound/pathway"
        )
        metKeggAnnot_file = "metKeggAnnot"
        with open(metKeggAnnot_file, 'w') as fh:
            fh.write(metKeggAnnot.content.decode("utf-8"))
        args.kmet2pathways = metKeggAnnot_file

    # PathwayID2PathwayNames
    if args.pathways:
        pathways_data = requests.get(
            "http://rest.kegg.jp/list/pathway/" + args.species
        )
        pathways_file = args.species + "_pathways"
        with open(pathways_file, 'w') as fh:
            fh.write(pathways_data.content.decode("utf-8"))
        args.pathways = pathways_file

    return args


def keggAnnot2list(keggAnnotFile, UniqueID, featureName, featureKeggId, featureType):
    """
    Create a dictionary that for the next function to find KEGG pathway information.
    Use when the input dataset has been parsed with metabolite_parser tool.

    Arguments:
        :param keggAnnotFile: Output file from Add KEGG Annotation Info Tool
        :type keggAnnotFile: file

        :param UniqueID: Name of the column with the unique identifiers.
        :type UniqueID: string

        :param featureName: Name of the column with feature names.
        :type featureName: string

        :param featureKeggId: Name of the column with KEGG Identifiers.
        :type featureKeggId: string

        :param featureType: One of 'Gene' or 'Metabolite'
        :type featureType: string

    Returns:
        :return featureDict: Dictionary with this information: (uniqueID + "\t" + featureName) =
        [nameInKegg, kegg_id]
        :rtype featureDict: dictionary

        :return featureList: List with the information of the features without kegg_id
        :return featureList: list
    """

    featureDict = {}
    featureList = []

    with open(keggAnnotFile, "r") as keggAnnot:
        header = keggAnnot.readline()
        header = header.strip().split("\t")
        keggIdCol = header.index(featureKeggId)
        uniqueIdCol = header.index(UniqueID)
        featureNameCol = header.index(featureName)
        for line in keggAnnot:
            uniqueID = str(line.split("\t")[uniqueIdCol])
            featureName = str(line.split("\t")[featureNameCol])
            if "Selected" in header:
                selectedCol = header.index("Selected")
                selected = line.split("\t")[selectedCol].strip()
                if selected == "Yes":
                    kegg_id = line.split("\t")[keggIdCol]
                    featureDict[uniqueID + "\t" + featureName] = kegg_id
                # For those that does not have kegg_id
                if selected == "NA":
                    featureList.append(
                        uniqueID +
                        "\t" +
                        featureName +
                        "\t" +
                        featureType +
                        "\t" +
                        "NA" +
                        "\t" +
                        "NA" +
                        "\t" +
                        "NA" +
                        "\n"
                    )
            else:
                kegg_id = line.split("\t")[keggIdCol]
                if kegg_id:
                    featureDict[uniqueID + "\t" + featureName] = kegg_id
                else:
                    featureList.append(
                        uniqueID +
                        "\t" +
                        featureName +
                        "\t" +
                        featureType +
                        "\t" +
                        "NA" +
                        "\t" +
                        "NA" +
                        "\t" +
                        "NA" + "\n"
                    )

    return (featureDict, featureList)


def add_path_info(
    featureDict,
    featureList,
    featureType,
    keggId2pathway,
    pathId2pathName,
    species,
    outputFile,
):
    """
    Find all pathways in KEGG related to a gene or a metabolite.

    Arguments:
        :param featureDict: Dictionary with this information: (uniqueID + "\t" + featureName) =
        [nameInKegg, kegg_id]
        :type featureDict: dictionary

        :param featureList: List with the information of the features without kegg_id
        :type featureList: list

        :param featureType: Type of the feature. One of 'Gene' or 'Metabolite'
        :type featureType: string

        :param keggId2pathway: Downloaded information from KEGG with the gene/metabolite KEGG Id
        and the Pathway ID
        :type keggId2pathway: file

        :param pathId2pathName: Downloaded information from KEGG with the Pathway ID and the
        Pathway Names
        :type pathId2pathName: file

        :param species: species identifier in kegg
        :type species: string

        :param outputFile: Output File Name to write the results.
        :type outputFile: string
    """

    output = open(outputFile, "w")
    output.write(
        "UniqueID\tFeature_Name\tFeature_Type\tKEGG_ID\tPathway_ID\tPathway_Name\n"
    )

    features = []
    for inputFeature, feature in featureDict.items():
        kegg_ids = [feature]
        # Step 1) Obtain path_id
        path_ids = parseWord(kegg_ids, keggId2pathway, 1, species)
        # Step 2) Obtain path_name
        paths = parseWord(path_ids, pathId2pathName, 0, species)
        if paths:
            for path in paths:
                features.append(inputFeature + "\t" + featureType + "\t" + path + "\n")
    features.extend(featureList)
    for feature in sorted(features, key=natural_keys):
        output.write(feature)

    output.close()


def parseWord(toParseList, parser, nameIndex, species):
    """
    Read a column from a KEGG file, find a given word and return a value from another column.
    Use to parse genes and metabolites into kegg identifiers and find related pathways via names.

    Arguments:
        :param toParseList: Input words list to parse
        :type toParseList: list

        :param parser: Kegg downloaded file
        :type parser: file

        :param nameIndex: index value of the column where the input word is in the parser file
        :type nameIndex: integer (0 or 1)

        :param species: species identifier in kegg
        :type species: string

    Returns:
        :return outputList: Input list with parsed value included. outputList_v2 if removing gene
        symbol needed
        :rtype outputList: list
    """

    parser = open(parser, "r")
    outputList = []
    outputList_v2 = []
    for toParseLine in toParseList:
        toParseLineList = toParseLine.split("\t")
        toParseWord = toParseLineList[len(toParseLineList) - 1]
        parser.seek(0)
        for parserLine in parser:
            names2findIn = parserLine.split("\t")[nameIndex].strip()
            if ";" in names2findIn:
                names2findIn = names2findIn.split(";")[0].strip()
            if findWholeWord(toParseWord)(names2findIn):
                parsedName = parserLine.split("\t")[abs(nameIndex - 1)].strip()
                # Remove " - species information" from pathway name
                if " - " in parsedName:
                    pathName = parsedName.split(" - ")[0]
                    outputList.append(toParseLine + "\t" + pathName)
                # Organism specific pathways
                elif "map" in parsedName:
                    pathId = parsedName.replace("map", species)
                    outputList.append(toParseLine + "\t" + pathId)
                else:
                    outputList.append(toParseLine + "\t" + parsedName)
        if not outputList:
            outputList.append(toParseLine + "\tNA")
    # In cases of ensembl-gene_symbol-kegg_id-path_id-path_name
    for outputLine in outputList:
        if len(outputLine.split("\t")) == 4:
            if outputLine.split("\t")[3] != "NA":
                outputLineList = outputLine.split("\t")
                del outputLineList[0]
                outputList_v2.append("\t".join(outputLineList))
            else:
                outputList.remove(outputLine)

    parser.close()

    return outputList_v2 if outputList_v2 else outputList


def findWholeWord(word):
    """
    Find one word inside another word.

    Arguments:
        :param word: word to check if is contained inside another word
        :type word: string

    Returns:
        :return is_contained: Returns "yes" if it's contained, "no" if not
        :rtype is_contained: boolean
    """

    return re.compile(
        r"\b(?![-])({0})(?![-])\b".format(word), flags=re.IGNORECASE
    ).search


def atoi(text):
    """
    Sort numbers in human readable order.
    """
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human readable order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    """
    return [atoi(c) for c in re.split("(\d+)", text)]


def fisherExactTest(args, path_feat):
    """
    Perform a complete Pathway enrichment analysis:
    1) Calculate all values of the contingency table
    2) Perform a Fisher Exact test
    3) Perform an FDR correction of p-values
    *********************
    * Fisher exact test *
    *********************
             | Molecules associated to pathway | Molecules not associated to pathway | Total
    -----------------------------------------------------------------------------------------
    DEG/M    |              a                  |                 b                   |  a+b
    -----------------------------------------------------------------------------------------
    No DEG/M |              c                  |                 d                   |  c+d
    -----------------------------------------------------------------------------------------
    Total    |             a+c                 |                b+d                  |   N
    -----------------------------------------------------------------------------------------

    N = Total of molecules
    a+b = Total 1 flags
    c+d = Total 0 flags

    Arguments:
        :params deaGeneDataset deaMetDataset: Tables with Differential Expression Analysis
        information for gene expression and metabolomics, respectively
        :types deaGeneDataset met_dataset: files

        :params gene_id_col, met_id_col, gene_flag_col, met_flag_col: Column names of unique
        identifiers and desired flag column of gene expression and metabolomics datasets,
        respectively
        :type gene_id_col, met_id_col, gene_flag_col, met_flag_col: strings

        :params alpha, method: alpha-value and method desired for the FDR correction
        :type alpha, method: strings

    Returns:
        :return output: Table with this structure: Pathway Name Odds_Ratio	P_value	FDR_Correction
        Flag_#
        :rtype output: file
    """

    # N_gene = sum(1 for line in open(args.deaGeneDataset)) - 1
    # N_met = sum(1 for line in open(args.deaMetDataset)) - 1
    # N = N_gene + N_met
    # a+b, c+d)
    metDeCounter = 0
    metNonDeCounter = 0
    with open(args.deaGeneDataset, "r") as geneDataset:
        header_gene = geneDataset.readline().split("\t")
        indice_gene_flag = header_gene.index(args.gene_flag_col)
        for line in geneDataset:
            flag = int(line.split("\t")[indice_gene_flag])
            if flag == 1:
                metDeCounter += 1
            elif flag == 0:
                metNonDeCounter += 1

    geneDeCounter = 0
    geneNonDeCounter = 0
    with open(args.deaMetDataset, "r") as metDataset:
        header_met = metDataset.readline().split("\t")
        indice_met_flag = header_met.index(args.met_flag_col)
        for line in metDataset:
            flag = int(line.split("\t")[indice_met_flag])
            if flag == 1:
                geneDeCounter += 1
            elif flag == 0:
                geneNonDeCounter += 1

    a_b = metDeCounter + geneDeCounter
    c_d = metNonDeCounter + geneNonDeCounter

    # a)
    geneDataset = open(args.deaGeneDataset, "r")
    metDataset = open(args.deaMetDataset, "r")
    PEAList = []
    pvalues = []
    indice_gene_id = header_gene.index(args.gene_id_col)
    indice_met_id = header_met.index(args.met_id_col)

    for pathway in path_feat.keys():
        a_gen = 0
        a_met = 0
        a_c = 0

        for value in path_feat[pathway]:
            geneDataset.seek(0)
            geneDataset.readline()
            for line in geneDataset:
                gene = line.split("\t")[indice_gene_id].replace('"', "")
                if value == gene:
                    a_gen += int(line.split("\t")[indice_gene_flag])
                    a_c += 1
                    break
            metDataset.seek(0)
            metDataset.readline()
            for line_v2 in metDataset:
                metabolite = line_v2.split("\t")[indice_met_id].replace('"', "")
                if value == metabolite:
                    a_met += int(line_v2.split("\t")[indice_met_flag])
                    a_c += 1
                    break

        a = a_gen + a_met

        # b), c), d)
        b = a_b - (a_gen + a_met)
        c = a_c - (a_gen + a_met)
        d = c_d - (a_c - (a_gen + a_met))

        oddsratio, pvalue = st.fisher_exact([[a, b], [c, d]])
        pvalues.append(pvalue)

        PEAList.append(pathway + "\t" + str(oddsratio) + "\t" + str(pvalue))

    geneDataset.close()
    metDataset.close()

    return PEAList


##########################
# All vs All Correlation #
##########################


def Ids2Names(dataset, Id, annot, annotName):
    """
    Change unique identifiers to feature names. It could be used for genes or for metabolites.

    Arguments:
        :param dataset: Wide dataset (Gene Expression/Metabolomics)
        :type dataset: pandas dataframe

        :param Id: Name of the column with Unique Identifiers.
        :type Id: string

        :param annot: Annotation File (For Gene Expression/Metabolomics). This file must contain at
        least 2 columns.
        :type annot: file

        :param annotName: Name of the column of the Annotation File that contains the Feature Name
        (Genes/Metabolites)
        :type annotName: file

    Returns:
        :return new_dataset: Wide dataset with feature names (genes/metabolites) instead of unique
        identifiers
        :rtype new_dataset: pandas dataset
    """

    annotTable = pd.read_table(annot, sep="\t", header=0)
    if "Selected" in annotTable.columns:
        annotTable = annotTable[annotTable["Selected"] == "Yes"]
        annotTable = annotTable.drop_duplicates(subset=Id, keep="first")

    for index, row in dataset.iterrows():
        ID = str(row[Id])
        name = str(annotTable.loc[annotTable[Id] == ID, annotName].item())
        dataset.loc[index, Id] = ID + ": " + name

    new_dataset = dataset.rename(columns={Id: annotName})
    new_dataset = new_dataset.set_index(annotName)

    return new_dataset


########
# sPLS #
########


def prepareSPLSData(args):
    """
    Perform subsetting of the data for the sPLS tool.

    Arguments:
        :param geneDataset metDataset: Gene expression and Metabolomics wide dataset, respectively.
        :type geneDataset metDataset: files

        :param geneOption metOption: Options for subsetting Gene Expression and Metabolomics
        datasets, respectively.
        :type geneOption metOption: string
    """

    args.geneDataset = pd.read_table(args.geneDataset, sep="\t", header=0)
    metTable = pd.read_table(args.metDataset, sep="\t", header=0)

    if args.metOption == "generic":
        metSubsetDict = prepareSPLSMetGenericData(args)
        if all(len(metList) < 3 for metList in metSubsetDict.values()):
            err_msg = """
                      There is no group of Metabolites with more than 3 identical KEGG ID.
                      Try changing the Metabolite subsetting option to 'By MMC Pattern'."
                      """
            sys.stderr.write(err_msg)
            sys.exit(2)
    elif args.metOption == "mmc":
        metSubsetDict = prepareSPLSMetMMCData(args)
        if all(len(metList) < 3 for metList in metSubsetDict.values()):
            err_msg = """
                      There is no block with more than 3 Metabolites.
                      Try changing the Metabolite subsetting option to 'By Metabolite Class'.
                      """
            sys.stderr.write(err_msg)
            sys.exit(2)
    elif args.metOption == "both":
        metSubsetDict = prepareSPLSMetBothData(args)
        if all(len(metList) < 3 for metList in metSubsetDict.values()):
            err_msg = """
                      Tere is no block with more than 3 Metabolites inside Metabolite Classes. Try
                      changing the Metabolite subsetting option to 'By MMC Pattern' or 'By
                      Metabolite Class'.
                       """
            sys.stderr.write(err_msg)
            sys.exit(2)

    metData = []
    rMetData = []
    rGeneData = []
    multipleNames, multipleNamesId = [], []
    for keggName, uniqueIDs in metSubsetDict.items():
        if 1 < len(uniqueIDs) < 3:
            print(
                "Metabolite",
                keggName.split(";")[0],
                "has only",
                len(uniqueIDs),
                "repetitions. At least 3 repetitions are needed to perform an sPLS.",
            )
        elif len(uniqueIDs) >= 3:
            multipleNames.append(keggName.split(";")[0])
            multipleNamesId.append(keggName)
            metSubTable = pd.DataFrame(columns=metTable.columns)
            for uniqueID in uniqueIDs:
                for index, row in metTable.iterrows():
                    if uniqueID == row[args.metId]:
                        metSubTable = metSubTable.append(
                            pd.Series(row, index=metTable.columns), ignore_index=True
                        )
            # Add Annotation
            if args.metAnno:
                metSubTable = Ids2Names(
                    metSubTable, args.metId, args.metAnno, args.metName
                )
            else:
                metSubTable = metSubTable.set_index(args.metId)
            metData.append(metSubTable)
            # convert to R dataframe
            with localconverter(robjects.default_converter + pandas2ri.converter):
                R_met_subdf = robjects.conversion.py2rpy(metSubTable)
            # Subset Gene Dataset
            if args.geneOption == "path":
                R_gene_df = prepareSPLSGenePathData(args, uniqueIDs)
                if R_gene_df == "":
                    print(keggName, "has no pathways associated.")
                else:
                    rMetData.append(R_met_subdf)
                    rGeneData.append(R_gene_df)
            else:
                rMetData.append(R_met_subdf)
    if args.geneOption == "all":
        R_gene_df = prepareSPLSGeneAllData(args)
        rGeneData.append(R_gene_df)
    elif args.geneOption == "geneList":
        R_gene_df = prepareSPLSGeneListData(args)
        rGeneData.append(R_gene_df)
    elif args.geneOption == "pana":
        R_gene_df = prepareSPLSGenePanaData(args)
        rGeneData.append(R_gene_df)

    return (rGeneData, rMetData, multipleNames, multipleNamesId)


def prepareSPLSMetGenericData(args):
    """
    Subsetting by Generic Metabolite in the sPLS Tool.

    Arguments:
        :param metKeggAnno: Metabolite output file from Add KEGG Annotation Info Tool.
        :type metKeggAnno: file

        :param metId: Name of the column with Metabolite Unique identifiers.
        :type metId: string

    Returns:
        :return metSubsetDict: Dictionary with this structure: keggName;keggCpdId =
        uniqueId;originalName
        :rtype metSubsetDict: dictionary
    """

    metKeggAnno = pd.read_table(args.metKeggAnno, sep="\t", header=0)
    metKeggAnno = metKeggAnno[metKeggAnno["Selected"] != "No"]
    metSubsetDict = {}

    for index, row in metKeggAnno.iterrows():
        keggName = str(row["Name_in_KEGG"])
        uniqueId = str(row[args.metId])
        keggCpdId = str(row["KEGG_ID"])
        if keggName != "nan":
            try:
                metSubsetDict[keggName + ";" + keggCpdId].add(uniqueId)
            except KeyError:
                metSubsetDict[keggName + ";" + keggCpdId] = set([uniqueId])

    return metSubsetDict


def prepareSPLSMetMMCData(args):
    """
    Subsetting by MMC clusters in the sPLS Tool.

    Arguments:
        :param metDataset: Metabolite wide dataset
        :type metDataset: file

        :param design: Design File, necessary for MMC program
        :type design: file

        :param metKeggAnno: Metabolite output file from Add KEGG Annotation Info Tool.
        :type metKeggAnno: file

        :param metId: Name of the column with Metabolite Unique identifiers.
        :type metId: string

        :param figure2: Full path for the MMC figure output.
        :type figure2: string

    Returns:
        :return metSubsetDict: Dictionary with this structure: keggName;keggCpdId =
        uniqueId;originalName
        :rtype metSubsetDict: dictionary
    """

    metSubsetDict = {}
    fh1List, fh2List, fh3List = [], [], []
    MMCResultTable = pd.DataFrame(
        columns=[args.metId, "Module", "Entry Index", "Average Degree", "Degree"]
    )

    args.input = args.metDataset
    args.uniqID = args.metId
    # MMC function
    mmc_table, fh1, fh2, fh3, MMCResultTable = mmc_main(args, "", MMCResultTable)
    mmc_tables = mmc_table.groupby("Module", as_index=False)
    mmc_dict = dict(iter(mmc_tables))
    for module, data in mmc_dict.items():
        metSubsetDict["Module_" + str(module)] = data[args.metId].tolist()

    # Create output from maps
    fh1List.append(fh1)
    fh2List.append(fh2)
    fh3List.append(fh3)
    mmcPlot(args.figure2, fh1List, fh2List, fh3List)
    mmcWriteOutput(args, MMCResultTable)

    return metSubsetDict


def prepareSPLSMetBothData(args):
    """
    Subsetting by Generic Metabolite and MMC clusters in the sPLS Tool.

    Arguments:
        :param metDataset: Metabolite wide dataset
        :type metDataset: file

        :param design: Design File, necessary for MMC program
        :type design: file

        :param metKeggAnno: Metabolite output file from Add KEGG Annotation Info Tool.
        :type metKeggAnno: file

        :param metId: Name of the column with Metabolite Unique identifiers.
        :type metId: string

        :param figure2: Full path for the MMC figure output.
        :type figure2: string

    Returns:
        :return metSubsetDict: Dictionary with this structure: keggName;keggCpdId =
        uniqueId;originalName
        :rtype metSubsetDict: dictionary
    """
    metDataset = pd.read_table(args.metDataset, sep="\t", header=0)
    genericMetsDict = prepareSPLSMetGenericData(args)
    metSubsetDict = {}
    fh1List, fh2List, fh3List = [], [], []
    MMCResultTable = pd.DataFrame(
        columns=[
            args.metId,
            "Module",
            "Entry Index",
            "Average Degree",
            "Degree",
            "Subset",
        ]
    )

    for keggName, uniqueIds in genericMetsDict.items():
        if 1 < len(uniqueIds) < 3:
            print(
                "Metabolite",
                keggName.split(";")[0],
                "has only",
                len(uniqueIds),
                "repetitions. At least 3 repetitions are needed to perform an sPLS.",
            )
        elif len(uniqueIds) >= 3:
            genericMetsList = []
            for uniqueId in uniqueIds:
                genericMetsList.append(uniqueId.split(";")[0])
            metSubsetTable = metDataset[metDataset[args.metId].isin(genericMetsList)]
            # MMC function needs file, not pandas dataframe. Create temporary file
            temp = tempfile.NamedTemporaryFile()
            metSubsetTable.to_csv(temp.name, sep="\t", index=False)
            args.input = temp.name
            args.uniqID = args.metId
            mmc_table, fh1, fh2, fh3, MMCResultTable = mmc_main(
                args, keggName.split(";")[0], MMCResultTable
            )
            fh1List.append(fh1)
            fh2List.append(fh2)
            fh3List.append(fh3)
            mmc_tables = mmc_table.groupby("Module", as_index=False)
            mmc_dict = dict(iter(mmc_tables))
            for module, data in mmc_dict.items():
                metSubsetDict[keggName.split(";")[0] + "_Module_" + str(module)] = data[
                    args.metId
                ].tolist()
            temp.close()

    # Create output from maps
    mmcPlot(args.figure2, fh1List, fh2List, fh3List)
    mmcWriteOutput(args, MMCResultTable)

    return metSubsetDict


def prepareSPLSGeneAllData(args):
    """
    Subsetting by All Genes (no subseting) in the sPLS Tool.

    Arguments:
        :param geneDataset: Gene Expression wide dataset
        :type geneDataset: file

        :param geneAnno: Gene Expression Annotation file.
        :type geneAnno: file

        :param geneId: Name of the column with Gene Expression Unique identifiers.
        :type geneId: string

        :param geneName: Name of the column in the geneAnno file that contains genes names.
        :type geneName: string

    Returns:
        :return R_gene_df: Wide dataset (annotated or not) in R format.
        :rtype R_gene_df: R object
    """

    geneTable = args.geneDataset

    # Add Annotation
    if args.geneAnno:
        geneTable = Ids2Names(geneTable, args.geneId, args.geneAnno, args.geneName)
    else:
        geneTable = geneTable.set_index(args.geneId)

    # convert to R dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
        R_gene_df = robjects.conversion.py2rpy(geneTable)
    return R_gene_df


def prepareSPLSGeneListData(args):
    """
    Subsetting by geneList in the sPLS Tool.

    Arguments:
        :param geneDataset: Gene Expression wide dataset
        :type geneDataset: file

        :param geneList: File with relevant genes (one per row) for the study.
        :type geneList: file

        :param geneAnno: Gene Expression Annotation file.
        :type geneAnno: file

        :param geneId: Name of the column with Gene Expression Unique identifiers.
        :type geneId: string

        :param geneName: Name of the column in the geneAnno file that contains genes names.
        :type geneName: string

    Returns:
        :return R_gene_df: Wide dataset (annotated or not) in R format.
        :rtype R_gene_df: R object
    """

    with open(args.geneList) as geneListFile:
        geneList = geneListFile.read().splitlines()
    geneTable_subset = args.geneDataset[args.geneDataset[args.geneId].isin(geneList)]

    # Add Annotation
    if args.geneAnno:
        geneTable_subset = Ids2Names(
            geneTable_subset, args.geneId, args.geneAnno, args.geneName
        )
    else:
        geneTable_subset = geneTable_subset.set_index(args.geneId)

    # convert to R dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
        R_gene_df = robjects.conversion.py2rpy(geneTable_subset)

    return R_gene_df


def prepareSPLSGenePathData(args, uniqueIDs):
    """
    Subsetting by metabolite's pathways in the sPLS Tool.

    Arguments:
        :param geneDataset: Gene Expression wide dataset
        :type geneDataset: file

        :param geneKeggPath metKeggPath: Output files from Add KEGG Pathway Info Tool
        :type geneKeggPath metKeggPath: files

        :param geneAnno: Gene Expression Annotation file.
        :type geneAnno: file

        :param geneId: Name of the column with Gene Expression Unique identifiers.
        :type geneId: string

        :param geneName: Name of the column in the geneAnno file that contains genes names.
        :type geneName: string

        :param keggName: metabolite KEGG Name to extract the pathway info.
        :type keggName: string

    Returns:
        :return R_gene_df: Wide dataset (annotated or not) in R format.
        :rtype R_gene_df: R object
    """

    geneKeggPath = pd.read_table(args.geneKeggPath, header=0, sep="\t")
    metKeggPath = pd.read_table(args.metKeggPath, header=0, sep="\t")

    pathsMet = (
        metKeggPath.loc[(metKeggPath[args.metId].isin(uniqueIDs)), "Pathway_Name"]
        .drop_duplicates()
        .values.tolist()
    )
    # Avoid NA:
    pathsMet = [path for path in pathsMet if str(path) != "nan"]
    if "Metabolic pathways" in pathsMet:
        pathsMet.remove("Metabolic pathways")
    if not pathsMet:
        R_gene_df = ""
    else:
        genesInPathsMet = (
            geneKeggPath.loc[(geneKeggPath["Pathway_Name"].isin(pathsMet)), args.geneId]
            .drop_duplicates()
            .values.tolist()
        )
        geneTable_subset = args.geneDataset.loc[
            (args.geneDataset[args.geneId].isin(genesInPathsMet))
        ]
        # Add Annotation
        if args.geneAnno:
            geneTable_subset = Ids2Names(
                geneTable_subset, args.geneId, args.geneAnno, args.geneName
            )
        else:
            geneTable_subset = geneTable_subset.set_index(args.geneId)
        # convert to R dataframe
        with localconverter(robjects.default_converter + pandas2ri.converter):
            R_gene_df = robjects.conversion.py2rpy(geneTable_subset)

    return R_gene_df


def prepareSPLSGenePanaData(args):
    """
    Subsetting by metagenes (PANA) in the sPLS Tool.

    Arguments:
        :param geneDataset: Gene Expression wide dataset
        :type geneDataset: file

        :param geneAnno: Gene Expression Annotation file.
        :type geneAnno: file

        :param geneKeggAnno: Output file from Add KEGG Annotation Info Tool for genes.
        :type geneKeggAnno: file

        :param geneId: Name of the column with Gene Expression Unique identifiers.
        :type geneId: string

        :param geneName: Name of the column in the geneAnno file that contains genes names.
        :type geneName: string

        :param path2genes: Downloaded file form KEGG with this information: path_id "\t" gene_id
        :type path2genes: file

    Returns:
        :return R_gene_df: Wide dataset (annotated or not) in R format.
        :rtype R_gene_df: R object
    """
    pandas2ri.activate()

    with ires.path("gait-gm.data", "PCA2GO.2.R") as my_r_script_path:
        f = open(my_r_script_path, "r")
        rFile = f.read()
    PANAScript = STAP(rFile, "PCA2GO")

    # Prepare input_file
    anno = pd.read_table(args.geneKeggAnno, sep="\t", header=0)
    anno = anno[anno["Selected"] != "No"]
    anno = anno.drop_duplicates(subset=args.geneId, keep="first")
    anno = anno[[args.geneId, "KEGG_ID"]]
    input_file = pd.merge(anno, args.geneDataset, on=args.geneId)
    # Only keep KEGG_ID
    input_file = input_file.drop(labels=args.geneId, axis=1)
    input_file = input_file.drop_duplicates(subset="KEGG_ID", keep="first")
    # Remove genes without KEGG_ID
    input_file = input_file.dropna()
    input_file = input_file.set_index("KEGG_ID")
    # convert to R dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
        R_input_file = robjects.conversion.py2rpy(input_file)

    # Prepare genes2pathway
    pathway2genes = pd.read_table(
        args.path2genes, sep="\t", header=None, names=["pathId", "geneId"]
    )
    genes2pathway = pathway2genes[["geneId", "pathId"]]
    # convert to R dataframe
    with localconverter(robjects.default_converter + pandas2ri.converter):
        R_genes2pathway = robjects.conversion.py2rpy(genes2pathway)

    # Run PANA script
    panaOutput = PANAScript.PCA2GO(
        X=R_input_file,
        annotation=R_genes2pathway,
        var_cutoff=args.cutoff,
        fac_sel=args.facSel,
    )
    with localconverter(robjects.default_converter + pandas2ri.converter):
        panaOutputTable = robjects.conversion.rpy2py(panaOutput[1])

    # Add Annotation
    if args.path2names:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            gene_df = robjects.conversion.rpy2py(panaOutput[0])
        gene_df.set_index("metagene_name", inplace=True)
        path2name = pd.read_table(
            args.path2names, sep="\t", names=["pathId", "pathName"]
        )
        path2nameDict = pd.Series(
            path2name.pathName.values, index=path2name.pathId
        ).to_dict()
        pathNames = []
        for pathway, row in gene_df.iterrows():
            pathNames.append(reduce_path_name(pathway, path2nameDict))

        # Gene "Wide" Dataset
        gene_df.insert(loc=0, column="pathNames", value=pathNames)
        gene_df = gene_df.set_index("pathNames")
        # convert to R dataframe
        with localconverter(robjects.default_converter + pandas2ri.converter):
            R_gene_df = robjects.conversion.py2rpy(gene_df)

        # panaOutputTable
        pathNames.insert(0, "KEGG_ID")
        panaOutputTable.columns = pathNames
        panaOutputTable["KEGG_ID"] = panaOutputTable.KEGG_ID.astype(str)   # AMM and Zihao
        if args.geneKeggAnno:
            panaOutputTable = Ids2Names(
                panaOutputTable, "KEGG_ID", args.geneKeggAnno, args.geneKeggName
            )
    else:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            gene_df = robjects.conversion.rpy2py(panaOutput[0])
        gene_df.set_index("metagene_name", inplace=True)

        with localconverter(robjects.default_converter + pandas2ri.converter):
            R_gene_df = robjects.conversion.py2rpy(gene_df)

    # Write PANA Output table
    panaOutputTable = panaOutputTable.astype(str)
    panaOutputTable.to_csv(args.panaOut, sep="\t", header=True, index=True)

    return R_gene_df


def reduce_path_name(pathway, path2nameDict):
    """
    Take a pathway name from KEGG and remove spaces and vowels to save space. This is very useful
    when using PANA and annotation file for metagenes for plots.

    Arguments:
        :param pathway: metagene ID: "PathwayID_panaIndex"
        :type pathway: string

        :param path2nameDict: Dictionary with KEGG Information linking pathway IDs with pathway
        names
        :type path2nameDict: Dictionary

    Returns:
        :return newPathName: Name of the pathway capitalized and without spaces and vowels.
        :rtype newPathName: string
    """

    # pathway = Pathway Name - organism_panaIndex
    pathway = str(pathway)
    pathId = pathway.split("_")[0]
    panaIndex = pathway.split("_")[1]
    pathName = path2nameDict[pathId]
    # Remove - Organism
    pathNameList = pathName.split(" - ")
    if len(pathNameList) > 2:
        newPathName = pathNameList[0] + "," + pathNameList[1]
    else:
        newPathName = pathNameList[0]

    # Capitalize, Remove vowels and spaces
    vowels = ["a", "e", "i", "o", "u"]
    newPathName = "".join([i for i in newPathName if i not in vowels])
    newPathName = newPathName.title()
    newPathName = newPathName.replace(" ", "")
    newPathName = pathway.split(":")[1] + ": " + newPathName + "_" + panaIndex

    return newPathName


# The following 4 functions have been modified from SECIM Tools modulated_modularity_clustering.py


def mmc_main(args, outputName, MMCResultTable):
    """
    Run a modulated modularity clustering (MMC) analysis
    """
    logger = logging.getLogger()
    sl.setLogger(logger)
    # Import data through the SECIMTools interface
    dat = wideToDesign(
        wide=args.input, design=args.design, uniqID=args.uniqID, logger=logger
    )

    # If there is no variance in a row, the correlations cannot be computed.
    dat.wide["variance"] = dat.wide.apply(lambda x: ((x - x.mean()).sum() ** 2), axis=1)
    dat.wide = dat.wide[dat.wide["variance"] != 0.0]
    dat.wide.drop("variance", axis=1, inplace=True)
    logger.info("Table arranged")

    # Compute the matrix of correlation coefficients.
    C = dat.wide.T.corr(method=args.correlation).values
    logger.info("Correlated")

    # For now, ignore the possibility that a variable
    # will have negligible variation.
    mask = np.ones(dat.wide.shape[0], dtype=bool)

    # Count the number of variables not excluded from the clustering.
    # p = np.count_nonzero(mask)

    # Consider all values of tuning parameter sigma in this array.
    sigmas, step = np.linspace(
        args.sigmaLow, args.sigmaHigh, num=args.sigmaNum, retstep=True
    )

    # Compute the clustering for each of the several values of sigma.
    # Each sigma corresponds to a different affinity matrix,
    # so the modularity matrix is also different for each sigma.
    # The goal is to the clustering whose modularity is greatest
    # across all joint (sigma, partition) pairs.
    # In practice, we will look for an approximation of this global optimum.
    # exit()
    logger.info("Begin clustering")
    clustering, sigma, m = get_clustering(C, sigmas)

    # Report a summary of the results of the technical analysis.
    logger.info("After partition refinement:")
    logger.info("Sigma: {0}".format(sigma))
    logger.info("Number of clusters: {0}".format(clustering.max() + 1))
    logger.info("Modulated modularity: {0}".format(m))

    # Run the nontechnical analysis using the data frame and the less nerdy
    # of the outputs from the technical analysis.
    global palette
    palette = colorHandler(pal=args.palette, col=args.color)
    mmc_table = mmc_nontechnical_analysis(
        args, dat.wide, mask, C, clustering, outputName, MMCResultTable
    )
    logger.info("Script Complete!")

    return mmc_table


def mmc_nontechnical_analysis(args, df, mask, C, clustering, outputName, MMCResultTable):
    """
    Re-order values for user convenience based on the results of the technical analysis.
    """

    # Get the map from the name to the original row index.
    all_row_names = df.index.values
    row_index_map = {s: i for i, s in enumerate(all_row_names)}

    # If some variables are uninformative for clustering,
    # the correlation matrix and the cluster vector will have smaller
    # dimensions than the number of rows in the original data frame.
    remaining_row_names = df[mask].index.values

    # Count the variables included in the clustering.
    p = clustering.shape[0]

    # Count the clusters.
    k = clustering.max() + 1

    # To sort the modules and to sort the variables within the modules,
    # we want to use absolute values of correlations.
    C_abs = np.abs(C)

    # For each cluster, get its indices and its submatrix of C_abs.
    selections = []
    submatrices = []
    degrees = np.zeros(p, dtype=float)
    for i in range(k):
        selection = np.flatnonzero(clustering == i)
        selections.append(selection)
        submatrix = C_abs[np.ix_(selection, selection)]
        submatrices.append(submatrix)
        if selection.size > 1:
            denom = selection.size - 1
            degrees[selection] = (submatrix.sum(axis=0) - 1) / denom

    # Modules should be reordered according to decreasing "average degree".
    cluster_sizes = []
    average_degrees = []
    for selection in selections:
        cluster_sizes.append(selection.size)
        average_degrees.append(degrees[selection].mean())

    module_to_cluster = np.argsort(average_degrees)[::-1]
    cluster_to_module = {v: k for k, v in enumerate(module_to_cluster)}

    triples = [(cluster_to_module[clustering[i]], -degrees[i], i) for i in range(p)]

    _a, _b, new_to_old_idx = zip(*sorted(triples))

    # Make a csv file if requested.
    if args.metOption == "both":
        header = (
            args.metId,
            "Module",
            "Entry Index",
            "Average Degree",
            "Degree",
            "Subset",
        )
    else:
        header = (args.metId, "Module", "Entry Index", "Average Degree", "Degree")
    temp = tempfile.NamedTemporaryFile()
    with open(temp.name, "w") as fout:
        writer = csv.writer(
            fout, "excel-tab"
        )  # problematic; need to switch to tsv file!
        writer.writerow(header)
        for old_i in new_to_old_idx:
            name = remaining_row_names[old_i]
            cluster = clustering[old_i]
            if args.metOption == "both":
                row = (
                    name,
                    cluster_to_module[cluster] + 1,
                    row_index_map[name] + 1,
                    average_degrees[cluster],
                    degrees[old_i],
                    outputName,
                )
            else:
                row = (
                    name,
                    cluster_to_module[cluster] + 1,
                    row_index_map[name] + 1,
                    average_degrees[cluster],
                    degrees[old_i],
                )
            writer.writerow(row)
    mmc_table = pd.read_table(temp, sep="\t", header=0)

    # Append to global mmc_table
    MMCResultTable = pd.concat([MMCResultTable, mmc_table], ignore_index=True)

    # Create Output
    fh1 = figureHandler(proj="2d")
    fh2 = figureHandler(proj="2d")
    fh3 = figureHandler(proj="2d")

    # Prepare to create the sorted heatmaps. (fh2)
    C_sorted = C[np.ix_(new_to_old_idx, new_to_old_idx)]
    clustering_new = clustering[np.ix_(new_to_old_idx)]

    # Draw the third heatmap (smoothed).
    # Make a smoothed correlation array. (fh3)
    S = expansion(clustering_new)
    block_mask = S.dot(S.T)
    denom = np.outer(S.sum(axis=0), S.sum(axis=0))
    small = S.T.dot(C_sorted).dot(S) / denom
    C_all_smoothed = S.dot(small).dot(S.T)
    C_smoothed = C_all_smoothed * (1 - block_mask) + C_sorted * block_mask

    # Getting list of names for heatmaps 2 and 3
    hpnames = [remaining_row_names[old_i] for old_i in new_to_old_idx]

    newNames = []
    sortNewNames = []
    # Traslating if annotation file provided
    if args.metAnno:
        annotTable = pd.read_table(args.metAnno, sep="\t", header=0)

        for noSortName in remaining_row_names:
            newName = annotTable.loc[
                annotTable[args.metId] == noSortName, args.metName
            ].item()
            newNames.append(noSortName + ": " + newName)
        for sortName in hpnames:
            newName = annotTable.loc[
                annotTable[args.metId] == sortName, args.metName
            ].item()
            sortNewNames.append(noSortName + ": " + newName)
    else:
        newNames = remaining_row_names
        sortNewNames = hpnames

    # Plot using something like http://stackoverflow.com/questions/15988413/
    # Drawing heatmaps
    # Draw first heatmap [C]
    plotHeatmap(C, fh1.ax[0], cmap=palette.mpl_colormap, xlbls=newNames, ylbls=newNames)

    # Draw second heatmap [C_sorted](reordered according to the clustering).
    plotHeatmap(
        C_sorted,
        fh2.ax[0],
        cmap=palette.mpl_colormap,
        xlbls=sortNewNames,
        ylbls=sortNewNames,
    )

    # Draw the heatmap [C_smoothed](smoothed version of C_sorted)
    plotHeatmap(
        C_smoothed,
        fh3.ax[0],
        cmap=palette.mpl_colormap,
        xlbls=sortNewNames,
        ylbls=sortNewNames,
    )

    if args.metOption == "both":
        fh1.formatAxis(xTitle="sampleID", figTitle="Correlations for: " + outputName)
        fh2.formatAxis(
            xTitle="sampleID", figTitle="Re-Ordered correlations for: " + outputName
        )
        fh3.formatAxis(
            xTitle="sampleID", figTitle="Smoothed correlations for: " + outputName
        )
    else:
        fh1.formatAxis(xTitle="sampleID", figTitle="Correlations")
        fh2.formatAxis(xTitle="sampleID", figTitle="Re-Ordered correlations")
        fh3.formatAxis(xTitle="sampleID", figTitle="Smoothed correlations")

    return (mmc_table, fh1, fh2, fh3, MMCResultTable)


def plotHeatmap(data, ax, cmap=False, xlbls=True, ylbls=True):
    """
    Create and plot a heatmap or a hcheatmap (hiearchical cluster heat map)

    :Arguments:
        :type data: pandas data frame.
        :param data: Data that is going to be used to plot

        :type ax: matplotlib axis
        :param ax: axis to be drawn on.
    :Returns:
        :return ax: matplotlib axis.
        :rtype ax: axis with the drawn figure.
    """

    # Create a custom colormap for the heatmap values
    # You can have two possible ways to create a pallete.
    # 1) sns.diverging_palette (a seaborn function)
    # 2) palettable.colorbrewer.diverging.[color].mpl_colors

    if not cmap:
        colors = colorHandler(pal="diverging", col="Spectral_10")
        cmap = colors.mpl_colormap

    # Draw the full plot
    sns.set(font_scale=0.1)
    sns.heatmap(
        data, linewidths=0.0, ax=ax, cmap=cmap, xticklabels=xlbls, yticklabels=ylbls
    )


def mmcPlot(figurePath, fh1List, fh2List, fh3List):
    """
    Create a PDF file and write all MMC plots to it.

    Arguments:
        :param figurePath: Full path for the MMC output figure.
        :type figurePath: string

        :params fh#List: Lists with first, second and third MMC heatmaps, respectively.
        :types fh#List: lists
    """

    # Create output from maps
    with PdfPages(figurePath) as pdf:
        for i in range(0, len(fh1List)):
            fh1List[i].addToPdf(pdf)
            fh2List[i].addToPdf(pdf)
            fh3List[i].addToPdf(pdf)


def mmcWriteOutput(args, MMCTable):
    """
    Write MMC Output table to a file. It's possible to change IDs for Names in an annotation file.

    Arguments:
        :param MMCTable: Output file from MMC Analysis.
        :type MMCTable: pandas dataframe

        :param metId: Name of the column with Metabolite Unique identifiers.
        :type metId: string

        :param metAnno: Metabolite Annotation File
        :type metAnno: file

        :param metName: Name of the column in the metAnno file that contains metabolite names.
        :type metName: string

        :param mmcOut: Full path to the MMC output table file.
        :type mmcOut: string
    """

    # Add Annotation
    if args.metAnno:
        MMCTable = Ids2Names(MMCTable, args.metId, args.metAnno, args.metName)
    else:
        MMCTable = MMCTable.set_index(args.metId)

    MMCTable.to_csv(args.mmcOut, sep="\t", header=True, index=True)
