# Matrix eQTL: Ultra Fast eQTL Analysis via Large Matrix Operations

Matrix eQTL is designed for fast eQTL analysis on large datasets.
Matrix eQTL can test for association between genotype
and gene expression using linear regression
with either additive or ANOVA genotype effects.
The models can include covariates to account for factors
as population stratification, gender, and clinical variables.
It also supports models with heteroscedastic and/or correlated errors,
false discovery rate estimation and
separate treatment of local (cis) and distant (trans) eQTLs.
For more details see
[Shabalin (2012)](https://academic.oup.com/bioinformatics/article/28/10/1353/213326).

# Key features

* Designed for eQTL analysis of large datasets.
* Performs testing for all or only local transcript-SNP pairs.
* Ultra-fast, no loss of precision.
* Equally fast for models with covariates.
* Supports:
    * Linear additive and ANOVA models.
    * Testing for the effect of genotype-covariate interaction.
    * Covariates to account for sex, population structure, surrogate variables, etc.
    * Limited support for correlated and heteroskedastic errors.
    * Correction for multiple testing using FDR [![External link](website/external.png "External link")](http://en.wikipedia.org/wiki/False_discovery_rate)
    * Separate p-value thresholds and FDR control for local and distant tests ([more info](website/runit.md#cis "Local and distant tests"))
* Convenient R package on [CRAN](https://CRAN.R-project.org/package=MatrixEQTL) and GitHub.

## Installation

### Install CRAN Version

To install the
[CRAN version](https://CRAN.R-project.org/package=MatrixEQTL)
of `MatrixEQTL`, run

```
install.packages("MatrixEQTL")
```

### Install GitHub Version

To install `MatrixEQTL` directly from GitHub, run

```
if(!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("andreyshabalin/MatrixEQTL")
```

## Getting started

The best way to start using Matrix eQTL is to first run the sample code on the provided toy dataset.
The sample code and data are part of the package and are also available ![here](data).
The package manual contains the sample code at the bottom of the help page for Matrix eQTL main function available via ?Matrix_eQTL_main command.
