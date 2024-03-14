# wf-admixture

This is a demonstration project for CSE284. It implements a smaller, simpler version of the admixture algorithm. See [ADMIXTURE](https://dalexander.github.io/admixture/index.html) page for more details. Refer to the final-project-files directory for presentation materials and summarization document.

I plan to implement several improvements listed in [Alexander 2011](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-246)


## Install instructions

Installation requires [NumPy](https://numpy.org) and [argparse](https://docs.python.org/3/library/argparse.html).

Navigate to the directory in which you would like to download wf-admixture and install:

```
git clone https://github.com/WillardFord/wf-admixture.git
cd wf-admixture
python setup.py install
```
Note: To install locally (i.e. without root access) use the `--user` tag when calling `setup.py`. (You might then have to make sure the local installation directory is your PATH.)

If the install was successful, typing `wf-admixture -h` should show a useful message.

## Basic usage

The basic usage of `wf-admixture` requires a `.bed`, `.bim`, and `.fam` files in the same directory.

Note: This is not a .bed file from UCSC genome browser. Refer to the [plink file documentation](https://www.cog-genomics.org/plink/1.9/formats#bed) for information about input file types.


To run `wf-admixture` on a small test example (using files in this repo):
```
wf-admixture -bed toy_files/toy.bed -k 3 -o toy_files/toy
```

Output files will then be stored at as numpy text arrays at the following location.
```
toy_files/toy.Q
toy_files/toy.F
```

To compare to output to `ADMIXTURE`, run:
```
admixture toy_files/toy.bed 3
```

## wf-admixture options

There are several optional inputs to `wf-admixture`

* `-q [int]`, `--threads`: Number of threads to use. Default will not use multithreading. Only recommended for large inputs.

* `-t [float]`, `--threshold`: Threshold between subsequent log liklihoods to indicate completion. Smaller values will be more accurate but take significantly more time. Default is 1.

* `-v`, `--verbose`: Indicates whether to output log liklihood information for each iteration.


## Real Data Example

Benchmark wf-admixture against ADMIXTURE using 1000 genomes data

Available in `analysis/benchmark.ipynb` and `analysis/benchmark/visualization.ipynb`

## File format

The output file formats are the same as [ADMIXTURE](https://dalexander.github.io/admixture/admixture-manual.pdf), a whitespace seperated table. Any continuing analysis should be interchangable.

## Testing

UNDER CONSTRUCTION: unit tests not yet implemented

To run tests:
```
# Run command line tests
sh tests/cmdline_tests.sh

# Run unit tests
python -m pytest --cov=
```

## Methodology

The wf-admixture tool implements a linear programming algorithm based off of Admixture. We are trying to estimate the contribution of population, k, to each individual, i. Let this be represented as an I x K matrix Q. In doing so we must calculate the minor allele frequency, j, of each population. Let this be represented as a K x J matrix F. 

1. Prune input SNPs in LD with each other using centimorgan distance.
2. Assume individuals are a independent unions of random gametes. This gives rise to a Hardy-Weinberg equilibrium of genotype at each SNP for each individual based on Q and F.
3. Calculate Log Liklihood equation from Hardy-Weinberg equilibrium.
4. Use linear programming techniques to iteratively calculate MLE from Log Liklihood. Primarily Expectation-Maximization algorithm.

## Next Steps

Given additional time I'd also like to implement the following improvements:

1. Use a Block Relaxation algorithm described in ADMIXTURE as opposed to EM algorithm which should be much faster when combined with updating only Q or F on a single step which allows us to use Taylor approximation given convexity assumption.
2. Add acceleration to linear programming method to accelerate our arrival at the optimum.
3. Add some level of confidence or error rate in these results.

## Sources

https://journals.plos.org/plosgenetics/article%3Fid=10.1371%2Fjournal.pgen.1003925#references

https://dalexander.github.io/admixture/admixture-manual.pdf

https://genome.cshlp.org/content/19/9/1655

https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-246

## Contributors

This repository was generated by Willard Ford, with inspiration from the [CSE 185 Example Repository](https://github.com/gymreklab/cse185-demo-project#readme) and the work of my fellow students.

Group 27 for the purposes of CSE 284.

Please submit a pull request with any corrections or suggestions.
