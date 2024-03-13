"""
Utilities for admixture
"""

import sys
from collections import namedtuple
import numpy as np
import math

# Named Tuples
SNP = namedtuple('BIM', ['ID','chromosome','locCM', 'locBP'])
SAMPLE = namedtuple('SAMPLE', ['ID', 'famID', 'fatherID', 'motherID'])

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def ERROR(msg:str, E=None) -> None:
    """
    Print an error message and die

    :param msg: Error message to print
    """
    sys.stderr.write(bcolors.FAIL + "[ERROR]: " + bcolors.ENDC + "{msg}\n".format(msg=msg) )
    if E: raise E
    sys.exit(1)

# TODO add method description
def readFile(filepath, delim='\t'): 
    with open(filepath, 'r') as f: 
        for line in f:
            yield tuple(line.split(delim))

def readBIM(bim_file):
    """
    Read in a .bim file and returns a tuple. Kills process if error while reading.

    :param bim_file: Existing bim_file path
    :return: array of named tuples, SNP, with fields: ID, chromosome, locCM, and locBP
            ID: variant identifier
            chromosome: Chromosome code
            loc_CM: position in morgans or centimorgans
            loc_BP: base-pair coordinates
    """
    snps :list[SNP] = []
    try:
        for values in readFile(bim_file):
            snps.append(SNP(values[1], values[0], values[2], values[3]))
    except Exception as E:
        message = "Failed to read in bim file. Are you sure it's in the correct format?"
        ERROR(message, E)
    return snps

def readFAM(fam_file):
    """
    Read in a .fam file and returns a tuple. Kills process if error while reading.

    :param fam_file: Existing fam_file path
    :return: array of named tuples, SAMPLE, with fields: ID, famID, fatherID, and motherID
            ID: within family ID of sample
            famID: unique identifier for family
            fatherID: within family ID of father
            motherID: within family ID of mother
    """
    # TODO maybe change within family ID's to unique IDs to make filtering 
    # easy by appending with family ID
    samples :list[SAMPLE] = []
    try:
        for values in readFile(fam_file):
            samples.append(SAMPLE(values[1], values[0], values[2], values[3]))
    except Exception as E:
        message = "Failed to read in fam file. Are you sure it's in the correct format?"
        ERROR(message, E)
    return samples

def readBED(bed_file, num_samples: int, num_variants:int) -> np.ndarray:
    """
    bed file:
        Necessary Information:
            1. Genotype information for every snp of each individual
    """
    """
    Read in a .bed file and returns a 2d numpy array of samples by variants
        representing minor allele counts. Kills process if error while reading.

    :param bed_file: Existing bed_file path
    :param num_samples: Number of samples in input
    :param num_variants: Number of variants in input
    :return: 2d numpy array of samples by variants representing minor allele counts.

    """
    genotypes = np.zeros((num_samples, num_variants))
    buffer_size = math.ceil(num_samples/4)
    bed_header =b"\x6c\x1b\x01"
    try:
        # TODO turn these into isolated functions
        with open(bed_file, "rb") as f:
            header = f.read(3)
            if header != bed_header: 
                ERROR("bed file missing proper header.")
            for j in range(num_variants):
                snp_genotypes = ""
                for _ in range(buffer_size): # Can we read in an entire buffer_size at once?
                    buffer = f.read(1)
                    # Convert a single byte into a binary string
                    snp_genotypes += "{0:b}".format(int.from_bytes(buffer))
                for i in range(num_samples):
                    first_allele = int(snp_genotypes[2*i])
                    second_allele = int(snp_genotypes[2*i+1])
                    if first_allele < second_allele:
                        # TODO: 01 indicates a missing genotype
                        pass
                    minor_allele_count = 2 - first_allele - second_allele
                    genotypes[i][j] = minor_allele_count
    except Exception as E:
        message = "Failed to read in bed file. Are you sure it's in the correct format?"
        ERROR(message, E)
    return genotypes


def initializeAlleleFrequencies(num_samples:int, num_snps:int, num_populations: int,
                                genotypes:np.ndarray, 
                                admixture_proportions: np.ndarray) -> np.ndarray:
    """
    TODO add comment
    TODO make into subroutines
    allele_freqs: F[k,j] = allele 1 frequency (maf by input standard) of snp j in pop k
    Constraints:
        1. For all x in F, 0 <= x <= 1
    Initialization should be completely given by Q init: 
        Alleles from random individuals contribute with different weights
        --> F[k,j] = weighted sum of snp j across all individuals
            i.e. F[k,j] = sum over all i of given_freq(i,j)*Q(i,k)/(2*Q(i,k)) 
            We only divide by the proportion of an indidividual that contributed 
                to that population. (Times 2 chromosomes)
    """
    allele_freqs = np.zeros((num_populations, num_snps))
    for k in range(num_populations):
        for j in range(num_snps):
            count_samples = 0
            for i in range(num_samples):
                count_samples += admixture_proportions[i][k]
                allele_freqs[k][j] += (admixture_proportions[i][k] * genotypes[i][j])
            if count_samples > 0:
                allele_freqs[k][j] /= (count_samples*2)
    return allele_freqs

def updateQ(I, J, K, G:np.ndarray, Q:np.ndarray, F:np.ndarray):
    new_Q = Q.copy()

    for i in range(I):
        for k in range(K):
            firstDifferentialFirstTerm = 0
            firstDifferentialSecondTerm = 0
            secondDifferentialFirstTerm = 0
            secondDifferentialSecondTerm = 0
            firstTermDenominator = 0
            secondTermDenominator = 0
            for j in range(J):
                firstDifferentialFirstTerm += G[i][j] * F[k][j]
                firstDifferentialSecondTerm += (2-G[i][j])*(1-F[k][j])
                
                for m in range(K):
                    firstTermDenominator += Q[i][m] * F[m][j]
                    secondTermDenominator += Q[i][m]*(1-F[m][j])




def updateF(I, J, K, G, Q, F):
    firstDifferential = 0
    secondDifferential = 0