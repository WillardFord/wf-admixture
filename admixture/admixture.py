#!/usr/bin/env python
"""
Command-line script to perform admixture algorithm of individuals given Plink input

Similar to ADMIXTURE
"""

import argparse
from . import utils as utils
import os
import sys
import time

import scipy as sp
import numpy as np

def main():
    parser = argparse.ArgumentParser(
        prog="",
        description=""
    )

    # Inputs
    parser.add_argument("bed",\
                        help="A .bed file as described by Plink.",
                        metavar="BED", type = str)
    
    parser.add_argument("k", help="An integer representing the number of populations "\
                        "to divide your dataset into.",
                        metavar="K", type = int)
    
    # Output
    parser.add_argument("-o", "--output", help="PREFIX will send your outputs to "\
                        "PREFIX.Q and PREFIX.P." \
                        "Default: stdout", metavar="PREFIX", type = str, required=False)
    
    parser.add_argument("-m", "--metrics", help="Used with -o. Will generate a runtime "\
                        "metrics file at PREFIX.metrics" \
                        "Default: No metrics file.", action='store_true', required=False)

    # Ideas:
    #   What types of metrics should I store?

    # Parse args
    args = parser.parse_args()

    # Check for valid inputs:
    K = args.k
    bed_file = args.bed
    metrics = args.metrics

    if K < 1:
        utils.ERROR("K must be greater than 0.")

    if not bed_file.endswith(".bed"):
        utils.ERROR(f"{bed_file} must have '.bed' suffix.")

    # Check for existence of .bim, .fam, and .bed files.
    input_path = os.path.join(os.path.dirname(bed_file),os.path.basename(bed_file)[:-4])
    bim_file = input_path + ".bim"
    fam_file = input_path + ".fam"
    
    if not os.path.isfile(bed_file):
        utils.ERROR(f"{bed_file} does not exist.")

    if not os.path.isfile(bim_file):
        utils.ERROR(f"{bim_file} does not exist.")

    if not os.path.isfile(fam_file):
        utils.ERROR(f"{fam_file} does not exist.")

    print("Congratulations! You loaded v0.0.0 of wf-admixture")

    # Read in input files
    # How do I store the input information

    # Set up output file
    if args.output == None:
        outf = sys.stdout
    else: 
        outf = open(args.output, "w")

    if metrics:
        pass


if __name__ == '__main__':
    main()
