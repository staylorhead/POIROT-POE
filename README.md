# POIROT: Parent-of-Origin Inference using Robust Omnibus Test

## Introduction
POIROT is a tool developed for identifying genetic variants harboring parent-of-origin effects in unrelated samples that considers multiple quantitative phentoypes simulateneously. The method is based on a robust test for equality of phenotypic covariance matrices between heterozygotes and homozygotes at a given locus. POIROT can handle both normally-distributed and non-normal continuous phenotypes. The method also adjusts for the effects of important covariates. For more information on this method, please see the following pre-print:

>[*POIROT: A powerful test for parent-of-origin effects in unrelated samples leveraging multiple phenotypes.* bioRxiv.](https://www.biorxiv.org/)

POIROT is implemented as a series of R functions which can be loaded and executed in an R environment and are described in this manual.

For questions or issues related to these R routines, please contact Taylor Head (<taylor.fischer@emory.edu>).

### The R Environment
R is a widely-used, free and open source software environment for statistical computing and graphics. The most recent version of R can be downloaded from the 
[Comprehensive R Archive Network (CRAN)](http://cran.r-project.org/)
CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions. For complete details on how to compile, install, and manage R and R packages, please refer to the manual [R Installation and Administration](http://cran.r-project.org/doc/manuals/r-release/R-admin.html).

## Input Files

Three tab-delimited text files are required as input for running POIROT. 

### 1. Phenotype File

* A file containing quantitative phenotype values as in `./ExampleData/phenotypes.txt`, with one subject per row and one phenotype per column.
* The first row of the file should contain phenotype labels.
* No missing data is allowed. 
* It is assumed all columns are phenotypes to be analyzed, and rows (subjects) should already be sorted to match the order of rows (subjects) in the other two input files. 

### 2. Covariate File

* A file containing quantitative or categorical covariate values as in `./ExampleData/covariates.txt`, with one subject per row and one covariate per column. 
* The first row of the file should contain covariate labels.
* No missing data is allowed. 
* It is assumed all phenotypes should be adjusted for each of the covariates in the file, and rows (subjects) should already be sorted to match the order of rows (subjects) in the other two input files. 

### 3. Individual-Level Genotype Data

* A file containing individual-level genotype data as in `./ExampleData/variants.txt`, with one subject per row and one SNP per column.
* The first row of the file should contain variant names.
* No missing data is allowed. 
* Rows (subjects) should already be sorted to match the order of rows (subjects) in the other two input files. 
* Variants must be coded by minor allele count (0/1/2).

## Example Analysis

For the example analysis, we have provided sample files containing simulated data on 1000 unrelated subjects. We have included data on 5 phenotypes we are interested in testing for parent-of-origin effects. While the method is scalable to hundreds of thousands of SNPs, we have included data on two SNPs for this toy analysis.

### 1. Load R Functions and Input Files

```
# Set working directory to that housing the POIROT functions R script and input files

source("POIROT-functions.R")

PHENO <- read.delim("phenotypes.txt")
GENO <- read.delim("variants.txt")
COVAR <- read.delim("covariates.txt")
```

### 2. Covariate Adjustment

```
PHENO_ADJ <- extract_residuals(PHENO,COVAR)
```

### 3. Perform Test

```
# the following assumes all variants in GENO will be tested for parent-of-origin effects
out <- data.frame(t(sapply(1:ncol(GENO),
                           FUN=do_POIROT_by_snp,
                           phenodat=PHENO_ADJ,
                           genodat=GENO)))
out$variant <- colnames(GENO)
head(out) # data frame of final results
```
