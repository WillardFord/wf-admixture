import random
import sys
import math
import random

"""
Usage:
python generateFiles.py samples variants populations
"""
def main():
    bed_file = "toy.bed"
    bim_file = "toy.bim"
    fam_file = "toy.fam"

    samples, variants = sys.argv[1:]
    num_samples = int(samples)
    num_variants = int(variants)

    # toy.fam
    # We only need unique ids for each user, nothing else we use in v0.0
    f = open(fam_file, 'w')
    for i in range(num_samples):
        f.write(f"{i}\t{i}\t0\t0\t0\t0\n")
    f.close()

    # toy.bim
    chrom = 1
    max_range = 240_000_000
    f = open(bim_file, 'w')
    for i in range(num_variants):
        location = random.randrange(1, max_range)
        f.write(f"{chrom}\t{i}\t0\t{location}\tA\tC\n")
    f.close()

    # TODO set up populations
    # TODO add admixture

    # toy.bed
    header = b"\x6c\x1b\x01"
    block_size = math.ceil(num_samples/4)
    f = open(bed_file, 'wb')
    f.write(header)
    random.seed(1)
    for i in range(num_variants):
        # Write two bits per for each sample
        buffer = ""
        for j in range(block_size*4):
            if j < num_samples:
                for _ in range(2):
                    if random.random() > 0.5: buffer += '1'
                    else: buffer += '0'
            else: 
                buffer += '00'
            if len(buffer) == 8:
                s=int(buffer, 2).to_bytes()
                f.write(s)
                buffer = ""
    f.close()

if __name__ == "__main__":
    main()