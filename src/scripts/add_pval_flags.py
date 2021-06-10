#!/usr/bin/env python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
# DESCRIPTION: Add flags to a Differential Expression Analysis Dataset.
# Output shows whether the feature belongs (1) or not (0) to a certain flag.
#######################################################################################

import os
import logging
import argparse
import linecache
from argparse import RawDescriptionHelpFormatter
import numpy as np
from numpy import genfromtxt
from secimtools.dataManager import logger as sl
import gaitGM.keggPeaModules as modules


def getOptions():
    parser = argparse.ArgumentParser(
        description="Add Pval Flags", formatter_class=RawDescriptionHelpFormatter
    )
    tool = parser.add_argument_group(description="Tool Input")
    tool.add_argument(
        "-de",
        "--deaDataset",
        dest="deaDataset",
        action="store",
        required=True,
        help="Differential Expression Analysis Datset.",
    )
    tool.add_argument(
        "-id",
        "--uniqID",
        dest="uniqID",
        action="store",
        required=True,
        help="Name of the column with unique identifiers.",
    )
    tool.add_argument(
        "-p",
        "--pvalue",
        dest="pvalue",
        action="store",
        required=True,
        help="Name of the column with P-values.",
    )
    tool.add_argument(
        "-t",
        "--thres",
        dest="thresholds",
        action="store",
        required=True,
        help="P-value thresholds.",
    )
    output = parser.add_argument_group(description="Output")
    output.add_argument(
        "-o",
        "--output",
        dest="output",
        action="store",
        required=True,
        help="Output file name.",
    )
    output.add_argument(
        "-fl",
        "--flags",
        dest="flags",
        action="store",
        required=True,
        help="Flags file name.",
    )
    args = parser.parse_args()

    args.deaDataset = os.path.abspath(args.deaDataset)
    args.output = os.path.abspath(args.output)
    args.flags = os.path.abspath(args.flags)

    return args


def main():
    """
    Add binary flags (0/1) to a differential expression dataset depending on p-value thresholds.

    Arguments:
        :param deaDataset: Matrix with Differential Expression Analysis information
        :type deaDataset: file

        :param pvalue: Name of the column with the p-value information
        :type pvalue: string

        :param uniqid: Name of the column with the unique identifier
        :type uniqid: string

        :param thresholds: Desired flag thresholds. Must be separed with ",", no spaces allowed.
        :type thresholds: string

    Returns:
        :return output: Table with input and added correspondent flags columns
        :rtype output: file

        :return flags: Table with only the correspondent flags columns
        :rtype flags: file
    """
    args = getOptions()
    logger = logging.getLogger()
    sl.setLogger(logger)
    logger.info(
        u"""Importing data with following parameters: \
        \n\tDEA Dataset: {0}\
        \n\tUnique ID: {1}\
        \n\tPvalues: {2}\
        \n\tThresholds: {3}""".format(
            args.deaDataset, args.uniqID, args.pvalue, args.thresholds
        )
    )

    modules.checkForDuplicates(args.deaDataset, args.uniqID)

    output = open(args.output, "w")
    flags = open(args.flags, "w")

    with open(args.deaDataset, "r") as data:
        header = data.readline().strip().split("\t")

    thresholds = args.thresholds.split(",")

    header_list = []
    for word in header:
        if word == "":
            output.write("NA")
            header_list.append("NA")
        elif header.index(word) == len(header) - 1:
            word = word.replace('"', "")
            output.write(word)
            header_list.append(word)
        else:
            word = word.replace('"', "")
            output.write(word + "\t")
            header_list.append(word)

    flags.write(str(args.uniqID))
    for threshold in thresholds:
        flags.write("\tFlag_" + threshold)
        output.write("\tFlag_" + threshold)
        header_list.append("\tFlag_" + threshold)

    flags.write("\n")
    output.write("\n")
    # Get P value column from a DEA dataset
    deaTable = genfromtxt(
        args.deaDataset,
        delimiter="\t",
        usecols=header_list.index(args.pvalue),
        dtype=None,
    )
    deaTable = np.delete(deaTable, 0)

    # Add 1/0 if smaller/greater than threshold
    i = 2
    for pvalue in deaTable:
        line = linecache.getline(args.deaDataset, i).strip()
        pvalue = float(pvalue.strip())
        flags.write(line.split("\t")[header_list.index(args.uniqID)])
        output.write(line)
        for threshold in thresholds:
            if pvalue <= float(threshold):
                flags.write("\t1")
                output.write("\t1")
            else:
                flags.write("\t0")
                output.write("\t0")
        flags.write("\n")
        output.write("\n")
        i += 1

    return args


if __name__ == "__main__":
    main()
