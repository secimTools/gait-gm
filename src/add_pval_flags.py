#!/usr/bin/env python
######################################################################################
# DATE: 2018/April/06
# 
# MODULE: add_flags.py
#
# VERSION: 1.0
# 
# AUTHOR: Francisco Huertas (f.huertas@ufl.edu)
#
# DESCRIPTION: The tool adds flags to a Differential Expression Analysis Dataset.
# Output says if the feature belongs (1) or not (0) to a certain flag.
#
#######################################################################################
# Import built-in packages
import os
import logging
import argparse
import linecache
from argparse import RawDescriptionHelpFormatter

# Import add-on packages
import numpy as np
from numpy import genfromtxt

# Import local data libraries
import keggPeaModules as modules
from secimtools.dataManager import logger as sl

def getOptions():

    parser = argparse.ArgumentParser(description="Add Pval Flags",
	   formatter_class=RawDescriptionHelpFormatter)

    # Tool input
    tool = parser.add_argument_group(description="Tool Input")
    tool.add_argument('-de',"--deaDataset", dest="deaDataset", action='store', 
	   required=True, help="Differential Expression Analysis Datset.")
    tool.add_argument('-id',"--uniqID", dest="uniqID", action='store', 
        required=True, help="Name of the column with unique identifiers.")
    tool.add_argument('-p',"--pvalue", dest="pvalue", action='store', 
	   required=True, help="Name of the column with P-values.")
    tool.add_argument('-t',"--thres", dest="thresholds", action='store', 
	   required=True, help="P-value thresholds.")    
    # Tool output
    output = parser.add_argument_group(description="Output")
    output.add_argument('-o',"--output", dest="output", action="store", 
            required=True, help="Output file name.")
    output.add_argument('-fl',"--flags", dest="flags", action="store", 
            required=True, help="Flags file name.")
    args = parser.parse_args()

    # Standardized paths    
    args.deaDataset = os.path.abspath(args.deaDataset)
    args.output = os.path.abspath(args.output)
    args.flags = os.path.abspath(args.flags)
    
    return(args)
    
def main(args):
    
    """
    Function that add binary flags (0/1) depending on p-value thresholds.
    
    Arguments:
        :param deaDataset: Matrix with Differential Expression Analysis information
        :type deaDataset: file
        
        :param pvalue: Name of the column with the p-value information
        :type pvalue: string
        
        :param uniqid: Name of the column with the unique identifier
        :type uniqid: string
        
        :param thresholds: Desired thresholds to decide the flag. Must be separed with ",", no spaces allowed.
        :type thresholds: string
    
    Returns:
        :return output: Table with the same information as the input, adding correspondent flags columns
        :rtype output: file
        
        :return flags: Table with only the correspondent flags columns
        :rtype flags: file
    """

    output = open(args.output, "w")
    flags = open(args.flags, "w")    
    
    with open(args.deaDataset, "r") as data:
        header = data.readline().strip()
    
    thresholds = args.thresholds.split(",") 

    # Solving empty header and don't add "\t" at the end of the header
    header_v1 = header.strip().split('\t')
    header_v2 = []
    for word in header_v1:
        if word == "":
            output.write("NA")
            header_v2.append("NA")
        elif header_v1.index(word) == len(header_v1) - 1:
            word = word.replace('"', '')
            output.write(word)
            header_v2.append(word)
        else:
            word = word.replace('"', '')
            output.write(word + "\t")
            header_v2.append(word)
    
    flags.write(str(args.uniqID))    
    # Add thresholds to header
    for threshold in thresholds:
        flags.write("\tFlag_" + threshold)
        output.write("\tFlag_" + threshold)
        header_v2.append("\tFlag_" + threshold)
        
    flags.write("\n")
    output.write("\n")
    # Get P value column from DEA dataset
    deaTable = genfromtxt(args.deaDataset, delimiter='\t', usecols=header_v2.index(args.pvalue), dtype=None)
    deaTable = np.delete(deaTable, 0)
    
    # Add 1/0 if smaller/greater than threshold
    i = 2   
    for pvalue in deaTable:
        line = linecache.getline(args.deaDataset, i).strip()
        pvalue = float(pvalue.strip())
        flags.write(line.split("\t")[header_v2.index(args.uniqID)])
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
    
    return(args)    
        
if __name__ == '__main__':
    # Command line options
    args = getOptions()
    
    # Setting up logger
    logger = logging.getLogger()
    sl.setLogger(logger)
    logger.info(u"""Importing data with following parameters: \
        \n\tDEA Dataset: {0}\
        \n\tUnique ID: {1}\
        \n\tPvalues: {2}\
        \n\tThresholds: {3}"""
        .format(args.deaDataset, args.uniqID, args.pvalue, args.thresholds))
        
    # Check for duplicates
    modules.checkForDuplicates(args.deaDataset, args.uniqID)
    
    # Calling main scripts
    main(args)