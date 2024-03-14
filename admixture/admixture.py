#!/usr/bin/env python
"""
Command-line script to perform admixture algorithm of individuals given Plink input

Similar to ADMIXTURE
"""

import argparse
from . import utils as utils
import os
from collections import namedtuple
import random

from datetime import datetime

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
                        "Default: stdout", metavar="PREFIX", type = str, required=True)
    
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

    print("Congratulations! You loaded v0.0.1 of wf-admixture")

    # Read in input files
    """
    Scratch Notes, please disregard. 
    Relevant documentation is in function descriptions and code is more easily 
        understood from reading it directly.
    bim file:
        Necessary Information:
            1. J: Number of snps
        Potentially Useful Information:
            1. SNP identifying information for output 
                (if we filter SNPs we'll need this)
                a. Chromosome code
                b. variant identifier
            2. SNP filtering information
                a. position in morgans or centimorgans
                b. base-pair coordinates
        Present Not Useful Information:
            1. Allele variant information
    fam file:
        Necessary Information:
            1. I: Number of individuals
        Potentially Useful Information:
            (To filter out individuals we'll need this information)
            1. Family information to filter out individuals from the same family
                a. Family ID
                b. within Family ID of father and mother
            2. Unique ID's to identify individuals if we filter some out
                a. Within-family ID with family ID
        Present Not Useful Information:
            1. Sex labels -- given equal proportions of each sex from each population,
                    sex information shouldn't add anything to our model
            2. Phenotype values -- not used here
    bed file:
        Necessary Information:
            1. Genotype information for every snp of each individual
    """
    #datetime.time(15, 8, 24, 78915)
    print(datetime.now())

    snps: list[utils.SNP] = utils.readBIM(bim_file)
    samples: list[utils.SAMPLE] = utils.readFAM(fam_file)

    I = len(samples)
    J = len(snps)
    genotypes :np.ndarray= utils.readBED(bed_file, I, J)

    print()
    print("Loaded input files")
    print("\tNum Samples:", I)
    print("\tNum Variants:", J)

    """
    Scratch Notes, please disregard. 
    Relevant documentation is in function descriptions and code is more easily 
        understood from reading it directly.
    admixture_proportions: Q[i,k] = population k's contribution to individual i
        Constraints:
            1. Each row should sum to 1 
                i.e. the contribution of all populations to an individual i sum to 1
            2. For all x in Q, x >= 0
        Init Possibilities: TODO (Start with 1, try others later)
            1. Uniform across rows
            2. Randomly assign value of 1 to a single population k for each individual.
            3. Randomly assign values summing to 1 to each population k for each individual.
    """
    # Init Option 3
    admixture_proportions :np.ndarray = np.ones((I, K)) / (10*(K-1))
    random.seed(11)
    for i in range(I):
        val = int(random.random()*K)
        admixture_proportions[i][val] = 0.9

    allele_frequencies :np.ndarray = utils.initializeAlleleFrequencies(I, J, K, 
                                        genotypes, admixture_proportions)

    print()
    print("Begin Linear Optimization")
    print(datetime.now())
    print()

    # Begin Linear Optimization
    epsilon = 1e-4 # Should probably be as low as 10e-4
    print("Epsilon:", epsilon)


    iterations = 0
    oldll = utils.logLiklihood(I, J, K, genotypes, admixture_proportions, 
                                allele_frequencies)
    print("\tstarting ll:", oldll)
    
    while True:
        admixture_proportions, allele_frequencies = utils.frappeEM(I, J, K, 
                                genotypes, admixture_proportions, 
                                allele_frequencies)
        #print("Admixture Proportions:\t", admixture_proportions)
        #print("Allele Frequencies:\t", allele_frequencies)
        iterations += 1
        print(f"Iteration: {iterations}")
        ll = utils.logLiklihood(I, J, K, genotypes, admixture_proportions, 
                                allele_frequencies)
        print("\tnew ll:", ll)
        dif =  ll - oldll
        print("\tdif:",  dif)
        if dif < 0:
            utils.ERROR("Log Liklihood went down!!")
        if dif < epsilon: break
        oldll = ll

        """
        if iterations % 2 == 0: admixture_proportions = utils.updateQ(I, J, K, genotypes, 
                                                                admixture_proportions, allele_frequencies)
        else: allele_frequencies = utils.updateF(I, J, K, genotypes, 
                                           admixture_proportions, allele_frequencies)
        num_iterations += 1
        if endCondition(admixture_proportions, allele_frequencies):
            break
        """

    # Set up output file
    prefix = args.output
    outQ = prefix+".Q"
    outF = prefix+".F"

    print()
    print(datetime.now())
    print("Writing to:")
    print("\t"+ outQ)
    print("\t"+ outF)
    print()

    np.savetxt(outQ, admixture_proportions)
    np.savetxt(outF, allele_frequencies)

    if metrics:
        pass

    print(datetime.now())

    # Write outputs

if __name__ == '__main__':
    main()
