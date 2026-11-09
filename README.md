# MAGICMAP

Mixtures Aggregating Genes Into Coordinated Modules And Pathways 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17566967.svg)](https://doi.org/10.5281/zenodo.17566967)

## Description

MAGICMAP is designed to compare genetic effect sizes (with standard errors) of lead variants from a GWAS of a target trait to GWAS of a comparitor trait to identify distinct genetic components. This is done by fitting a regression mixture model where each component is a York regression with the intercept fixed to zero. The York regression accounts for uncertainty in the effect sizes, including potentially correlated sampling error (e.g. due to sample overlap between the two GWAS). Correlated sampling error can be estimated based on LD Score regression.

Details of the model are provided in the [vignette](/vignettes/MAGICMAP_model_description.pdf).

## Installation

This R package can be installed using the devtools package:

    devtools::install_github("rkwalters/MAGICMAP")

It depends on:

* [SCOUTJOY](https://github.com/aelliott08/SCOUTJOY)
* [yorktools](https://github.com/rkwalters/yorktools)
* mvtnorm
* Ckmeans.1d.dp

### Rust engine branch

The `dev-jp-rust` branch adds an optional compiled engine for the EM fit
(`magicmap(..., engine = "rust")`). The default engine is still the R code, and
the rust engine is a like-for-like port of it, so the two reach the same
solutions. Building this branch from source requires a Rust toolchain with
`cargo` and `rustc` on `PATH` or in `~/.cargo/bin` (see https://rustup.rs/).
Rust is needed to build the package but not to run it.

    devtools::install_github("carbocation/MAGICMAP", ref = "dev-jp-rust")

`benchmarks/benchmark-rust-em.R` times the two engines against each other.

## Usage

### Inputs

For lead variants from GWAS of the target phenotype:

* GWAS effect sizes (betas) and corresponding standard errors for the target trait
* GWAS effect sizes (betas) and corresponding standard errors for the comparator trait

From [LD Score regression](https://github.com/bulik/ldsc) analysis of genetic correlation between the two GWAS:

* Covariance intercept (`gcov_int`)
* Intercepts for GWAS of the target and comparator traits

If the two GWAS do not have overlapping samples, then the LDSC information can be omitted if necessary (e.g. if full genome-wide summary statistics are unavailble for estimation.)


### Running MAGICMAP

Documentation is accessible within R by running `?magicmap`

Briefly, if `dat` is a data frame that includes columns:

* `lead_SNP`: lead variant for the locus from GWAS of the target trait
* `BETA.x`: estimated effect sizes of the named SNP from GWAS of the target trait 
* `SE.x`: standard error of the effect sizes of the named SNP from GWAS of the target trait 
* `BETA.y`: estimated effect sizes of the named SNP from GWAS of the comparison trait 
* `SE.y`: standard error of the effect sizes of named lead SNP from GWAS of the comparison trait

Then MAGICMAP can be run using:

    require("MAGICMAP")
    fit <- magicmap(data = dat, "BETA.x","SE.x","BETA.y","SE.y",
                    ids="lead_SNP",
                    k=2:6,
                    CovIntercept=gcov_int, TargetXIntercept = x_int, ComparatorYIntercept = y_int)

where `k` specifies the number of mixture components to consider, and `gcov_int`, `x_int`, and `y_int` are the intercept estimates from LDSC.


### Outputs

`magicmap()` outputs an R object that contains:

* SCOUTJOY's global test of whether then effect sizes are consistent with having a single homogeneous relationship (i.e. whether k=1 can be rejected)
* Fit statistics for each considered value of k (number of components), including the log likelihood, AIC, and BIC
* For each considered value of k, the estimated slope and prior for each component, and the posterior probabilities and best fit effect sizes for each input variant

The best fit model can be visualized by running `plot.magicmap()` on the output from the `magicmap()` function. See the documentation for details on customizing the plot.


