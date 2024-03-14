# wf-admixture

This is a demonstration project for CSE284. It implements a smaller, simpler version of the admixture algorithm. See [ADMIXTURE](https://dalexander.github.io/admixture/index.html) page for more details. Refer to the final-project-files directory for presentation materials and summarization document.

I plan to implement several improvements listed in [Alexander 2011](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-246)


## Install instructions

Installation requires [NumPy](https://numpy.org), [SciPy](https://scipy.org), and (PuLP)[https://pypi.org/project/PuLP/].

To install PuLP use the below commands:
```
pip install pulp
# Test installation and complete setup
pulptest
```

Navigate to the directory in which you would like to download wf-admixture and install:

```
git clone https://github.com/WillardFord/wf-admixture.git
cd wf-admixture
python setup.py install
```
Note: To install locally (i.e. without root access) use the `--user` tag when calling `setup.py`.

If the install was successful, typing `wf-admixture --help` should show a useful message.

## Basic usage

The basic usage of `wf-admixture` requires a `.bed`, `.bim`, and `.fam` files in the same directory.

Note: This is not a .bed file from UCSC genome browser. Refer to the [plink file documentation](https://www.cog-genomics.org/plink/1.9/formats#bed) for information about input file types.

```
wf-admixture your-favorite.bed k [-o PREFIX]
```

To run `wf-admixture` on a small test example (using files in this repo):
```
wf-admixture toy_files/toy.bed 3 -o toy_files/toy
```

To see the output:
```
cat toy_files/toy.Q
cat toy_files/toy.F
```

To compare to output of `ADMIXTURE`, run:
```

```

## wf-admixture options

There are 2 required inputs to `wf-admixture`, a reference ___TODO

* `-o PREFIX`, `--output PREFIX`: Write ancestry fractions to `PREFIX.Q` and allele frequencies to `PREFIX.P`. By default outputs are written to stdout.

* `-m`, `--metrics`: Used with -o. Will generate a runtime metrics file at `PREFIX.metrics`.

## Real Data Example

TODO:
Benchmark wf-admixture against ADMIXTURE using real data.

```
CURDIR=chr21Example
mkdir $CURDIR
cd $CURDIR

VCF_DOWNLOAD=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
VCF=1000G_ALL.chr21.vcf.gz

SAMPLEINFO=igsr_samples.tsv

PREFIX=chr21-test
VALID_SAMPLES=admixed_sample_ids.txt

curl $VCF_DOWNLOAD -o $VCF


echo "" > $VALID_SAMPLES

cat $SAMPLEINFO | grep "1000 Genomes phase 3 release" | awk -F"\t" '($4=="CEU")' | awk '{print $1 "\t" $1}' >> $VALID_SAMPLES

cat $SAMPLEINFO | grep "1000 Genomes phase 3 release" | awk -F"\t" '($4=="PEL")' | awk '{print $1 "\t" $1}' >> $VALID_SAMPLES

cat $SAMPLEINFO | grep "1000 Genomes phase 3 release" | awk -F"\t" '($4=="GWD")' | awk '{print $1 "\t" $1}' >> $VALID_SAMPLES

cat $SAMPLEINFO | grep "1000 Genomes phase 3 release" | awk -F"\t" '($4=="ASW")' | awk '{print $1 "\t" $1}' >> $VALID_SAMPLES

cat $SAMPLEINFO | grep "1000 Genomes phase 3 release" | awk -F"\t" '($4=="PUR")' | awk '{print $1 "\t" $1}' >> $VALID_SAMPLES

plink --vcf $VCF --keep $VALID_SAMPLES --double-id --maf 0.01 --out $PREFIX --make-bed

plink --bed $PREFIX.bed --bim $PREFIX.bim --fam $PREFIX.fam --chr 21 --indep-pairwise 50 10 0.1

plink --bed $PREFIX.bed --bim $PREFIX.bim --fam $PREFIX.fam --chr 21 --extract plink.prune.in --make-bed --out $PREFIX.pruned

plink --bed $PREFIX.pruned.bed --bim $PREFIX.pruned.bim --fam $PREFIX.pruned.fam --thin-count 1000 --make-bed --out $PREFIX.pruned.thinned

cd ..
wf-admixture $CURDIR/$PREFIX.pruned.thinned.bed 3 -o $CURDIR/$PREFIX.pruned.thinned.wf-admixture

---

TINY_SAMPLES=admixed_sample_ids.tiny.txt
head -n 100 $VALID_SAMPLES | tail -n 5 > $TINY_SAMPLES

plink --bed $PREFIX.pruned.bed --bim $PREFIX.pruned.bim --fam $PREFIX.pruned.fam --thin-count 3 --make-bed --out $PREFIX.pruned.tiny

plink --vcf $VCF --keep $TINY_SAMPLES --double-id --maf 0.01 --out $PREFIX.tiny --make-bed
plink --bed $PREFIX.tiny.bed --bim $PREFIX.tiny.bim --fam $PREFIX.tiny.fam --chr 21 --indep-pairwise 50 10 0.1
plink --bed $PREFIX.tiny.bed --bim $PREFIX.tiny.bim --fam $PREFIX.tiny.fam --extract plink.prune.in --thin-count 3 --make-bed --out $PREFIX.tiny.pruned.thinned

cd ..
wf-admixture $CURDIR/$PREFIX.tiny.pruned.thinned.bed 3 -o $CURDIR/$PREFIX.tiny.pruned.thinned.wf-admixture

```

## File format

The output file formats are the same as [ADMIXTURE](https://dalexander.github.io/admixture/admixture-manual.pdf), a tab seperated table. Any continuing analysis should be interchangable.

## Testing

UNDER CONSTRUCTION: unit tests not yet implemented

To run tests:
```
# Run command line tests
sh tests/cmdline_tests.sh

# Run unit tests
python -m pytest --cov=.
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
