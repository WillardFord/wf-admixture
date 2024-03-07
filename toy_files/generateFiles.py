import random
import sys
import math

"""
Usage:
python generateFiles.py samples variants populations
"""
def main():
    bed_file = "toy.bed"
    bim_file = "toy.bim"
    fam_file = "toy.fam"

    samples, variants, populations = sys.argv[1:]
    num_samples = int(samples)
    num_variants = int(variants)
    num_populations = int(populations)

    # toy.fam
    # We only need unique ids for each user, nothing else we use
    f = open(fam_file, 'w')
    for i in range(num_samples):
        f.write(f"{i}\t{i}0\t0\t0\t0\n")
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
    prefix = b"\x6c\x1b\x01"
    num_blocks = math.ceil(num_samples/4)
    f = open(bed_file, 'wb')
    f.write(prefix)
    for i in range(num_variants):
        # Write two bits per iteration
        for j in range(num_blocks * 4):
            sample_id = j//2 + 1
            if sample_id > num_samples:
                f.write(b"\x00\x00")
                continue
            f.write(b"\x01\x01")
    f.close()

if __name__ == "__main__":
    main()