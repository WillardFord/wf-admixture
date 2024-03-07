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
                        "PREFIX.Q and PREFIX.P and create a PREFIX.metrics file. " \
                        "Default: stdout", metavar="PREFIX", type = str, required=False)

    # Ideas:
    #   

    # Parse args
    args = parser.parse_args()

    # Set up output file
    if args.output == None:
        metrics = True
        outf = sys.stdout
    else: 
        outf = open(args.output, "w")
        metrics = False

    print("Congratulations! You loaded v0.0.0")

if __name__ == '__main__':
    main()
