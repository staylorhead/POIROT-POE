## POIROT: Parent-of-Origin Inference using Robust Omnibus Test

### Introduction
In the pre-print available [here], we created a test for identifying genetic variants harboring parent-of-origin effects in unrelated samples that considers multiple quantitative phentoypes simulateneously. The method is based on a robust test for equality of phenotypic covariance matrices between heterozygotes and homozygotes at a given locus. POIROT can handle both normally-distributed and non-normal continuous phenotypes. The method also adjusts for the effects of important covariates.

POIROT is implemented as an R function, which can be loaded and executed in an R environment.

### The R Environment
R is a widely-used, free and open source software environment for statistical computing and graphics. The most recent version of R can be downloaded from the 
[Comprehensive R Archive Network (CRAN)](http://cran.r-project.org/)
CRAN provides precompiled binary versions of R for Windows, MacOS, and select Linux distributions. For complete details on how to compile, install, and manage R and R packages, please refer to the manual [R Installation and Administration](http://cran.r-project.org/doc/manuals/r-release/R-admin.html).

### Required R Packages for POIROT

There are no dependencies required for this method at the moment.

### Running POIROT
For the example analysis, we have provided sample files containing simulated data on 1000 unrelated subjects. We have included data on 5 phenotypes we are interested in testing for parent-of-origin effects. While the method is scalable to hundreds of thousands of SNPs, we have included data on two SNPs for ease.
