#!/usr/bin/env python
######################################################################################
# AUTHOR: Francisco Huertas <f.huertas@ufl.edu>
# CONTRIBUTORS: Alison Morse <ammorse@ufl.edu>, Oleksandr Moskalenko <om@rc.ufl.edu>
# DESCRIPTION: Take a wide dataset and append a unique identifier to it.
#######################################################################################

import os
import logging
import argparse
from argparse import RawDescriptionHelpFormatter
import gaitGM.keggPeaModules as modules
from secimtools.dataManager import logger as sl


def getOptions():
    parser = argparse.ArgumentParser(
        description="split_wide_dataset", formatter_class=RawDescriptionHelpFormatter
    )
    standard = parser.add_argument_group(description="Required Input")
    standard.add_argument(
        "-i",
        "--input",
        dest="input",
        action="store",
        required=True,
        help="Input dataset in wide format.",
    )
    standard.add_argument(
        "-id",
        "--ID",
        dest="uniqID",
        action="store",
        required=False,
        help="Name of the column with unique identifiers.",
    )
    tool = parser.add_argument_group(title="Tool Specific Inputs")
    tool.add_argument(
        "-s",
        "--samples",
        dest="samples",
        action="store",
        required=True,
        help="Select sample columns.",
    )
    tool.add_argument(
        "-p",
        "--prefix",
        dest="prefix",
        action="store",
        required=False,
        help="Prefix to add to the new unique ID.",
    )
    tool.add_argument(
        "-p2",
        "--prefix2",
        dest="prefix2",
        action="store",
        required=False,
        help="Prefix to add to the old unique ID\
                      (only if the Unique ID is numeric).",
    )
    output = parser.add_argument_group(description="Output")
    output.add_argument(
        "-w",
        "--wide",
        dest="wide",
        action="store",
        required=True,
        help="Wide dataset file.",
    )
    standard.add_argument(
        "-d",
        "--design",
        dest="design",
        action="store",
        required=True,
        help="Design file.",
    )
    output.add_argument(
        "-a",
        "--annot",
        dest="annot",
        action="store",
        required=True,
        help="Annotation file.",
    )
    args = parser.parse_args()

    args.input = os.path.abspath(args.input)
    args.wide = os.path.abspath(args.wide)
    args.design = os.path.abspath(args.design)
    args.annot = os.path.abspath(args.annot)

    return args


def main():
    """
    Creates a Wide Dataset, a Design File, and an Annotation File.

    Arguments:
        :param input: Input dataset without Wide Format
        :type input: file

        :param uniqID: Unique Identifier Column Name (optional)
        :type uniqID: string

        :param samples: Position of the samples columns in the input dataset, separated by commas
        :type samples: string

        :param prefix: Prefix to add to the new Unique ID (optional)
        :type prefix: string

        :param prefix2: Prefix to add to the old Unique ID (optional). Necessary if provided unique
            ID is numeric.
        :type prefix2: string

    Returns:
        :return wide: Input dataset in Wide Format (With Unique ID and samples in columns)
        :rtype wide: file

        :return design: Design File template for samples
        :rtype design: file

        :return annot: Annotation File (With Unique ID and annotations in columns)
        :rtype annot: file
    """
    args = getOptions()
    logger = logging.getLogger()
    sl.setLogger(logger)
    logger.info(
        u"""Importing data with following parameters: \
        \n\tInput: {0}\
        \n\tSamples: {1}""".format(
            args.input, args.samples
        )
    )

    if args.uniqID:
        modules.checkForDuplicates(args.input, args.uniqID)

    wide = open(args.wide, "w")
    design = open(args.design, "w")
    annot = open(args.annot, "w")

    header_v2 = []
    rowNum = 1
    with open(args.input, "r") as inputDat:
        header = inputDat.readline().strip().split("\t")
        # Design File
        design.write("sampleID\t\n")
        for sample in args.samples.split(","):
            design.write(header[int(sample.strip()) - 1] + "\t\n")
            header_v2.append(header[int(sample.strip()) - 1])
        noSamples = list(set(header) - set(header_v2))
        # Annotation File & Wide Dataset Headers
        if args.uniqID:
            noSamples.remove(args.uniqID)
            wide.write(str(args.uniqID) + "\t" + "\t".join(header_v2) + "\n")
            if noSamples:
                annot.write(str(args.uniqID) + "\t" + "\t".join(noSamples) + "\n")
            else:
                annot.write(str(args.uniqID) + "\n")
        else:
            wide.write("UniqueID\t" + "\t".join(header_v2) + "\n")
            if noSamples:
                annot.write("UniqueID\t" + "\t".join(noSamples) + "\n")
            else:
                annot.write("UniqueID\n")
        # Annotation File & Wide Dataset Filling
        for row in inputDat:
            if args.uniqID:
                idNcol = header.index(args.uniqID)
                uniqueID = row.split("\t")[idNcol]
                if args.prefix2:
                    uniqueID = args.prefix2 + "_" + uniqueID
            else:
                uniqueID = str(args.prefix) + "_" + str(rowNum)
            data = ""
            for ncol in args.samples.split(","):
                data = str(data) + "\t" + str(row.strip().split("\t")[int(ncol) - 1]).strip()
            annot.write(uniqueID)
            for noSample in noSamples:
                ncol = header.index(noSample)
                annot.write("\t" + row.split("\t")[ncol])
            annot.write("\n")
            wide.write(uniqueID + data + "\n")
            rowNum += 1

    wide.close()
    design.close()
    annot.close()

    return args


if __name__ == "__main__":
    main()
