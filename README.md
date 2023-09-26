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
The sample code and data are part of the package and are also available [here](data).
The package manual contains the sample code at the bottom of the help page for the Matrix eQTL main function.
Access it via `?Matrix_eQTL_main` command.

The sample code performs eQTL analysis of a toy data set consisting of three files:
[genotype](data/SNP.txt), [expression](data/GE.txt), and [covariates](data/Covariates.txt).
For every gene-SNP pair it runs linear regression analysis accounting for the set of covariates.

Let's go over the sample code line by line.

```r

# First step is to load the package:

library("MatrixEQTL")

# The toy data set files are stored with the package at the following location.

base.dir = find.package("MatrixEQTL")

# Then we set the parameters such as selected linear model and
# names of genotype and expression data files.

useModel = modelLINEAR # modelANOVA or modelLINEAR or modelLINEAR_CROSS
SNP_file_name = paste0(base.dir, "/data/SNP.txt")
expression_file_name = paste0(base.dir, "/data/GE.txt")

# A separate file may be provided with extra covariates.
# In case of no covariates set the variable covariates_file_name to character().

covariates_file_name = paste0(base.dir, "/data/Covariates.txt")
output_file_name = tempfile()

# The p-value threshold determines which gene-SNP associations are
# saved in the output file output_file_name.
# Note that for larger datasets the threshold should be lower.
# Setting the threshold to a high value for a large dataset may
# cause excessively large output files.

pvOutputThreshold = 1e-2

# Finally, define the covariance matrix for the error term.
# This parameter is rarely used.
# If the covariance matrix is a multiple of identity, set it to numeric().

errorCovariance = numeric()

# The next section of the sample code contains three very similar parts
# loading the files with genotype, gene expression, and covariates.
# In each part one can set the file delimiter
# (i.e. tabulation "\t", comma ",", or space " "),
# the string representation for missing values,
# the number of rows with column labels, and
# the number of columns with row labels.
# Finally, one can change the number of the variables
# in a slice for the file reading procedure (do not change if not sure).

snps = SlicedData$new()
snps$fileDelimiter = "\t"      # the TAB character
snps$fileOmitCharacters = "NA" # denote missing values
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name )

gene = SlicedData$new()
gene$fileDelimiter = "\t"      # the TAB character
gene$fileOmitCharacters = "NA" # denote missing values
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in pieces of 2,000 rows
gene$LoadFile( expression_file_name )

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      # the TAB character
cvrt$fileOmitCharacters = "NA" # denote missing values
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
cvrt$fileSliceSize = 2000      # read file in pieces of 2,000 rows
cvrt$LoadFile( covariates_file_name )

# Finally, the main Matrix eQTL function is called:

me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)

```

Each significant gene-SNP association is recorded in a separate line in the output file and in the returned object `me`.
In case of cis/trans eQTL analysis described below, two output files are produced, one with cis-eQTLs, another only with trans.
Every record contains a SNP name, a transcript name, estimate of the effect size, t- or F-statistic, p-value, and FDR.


## Cis- and trans- eQTL analysis

Matrix eQTL can distinguish local (cis-) and distant (trans-) eQTLs and
perform separate correction for multiple comparisons for those groups.

<!--- The <a href="R.html#cis">second sample code</a> shows how to run such analysis. --->

The main Matrix eQTL function `Matrix_eQTL_main` requires several extra parameters for local/distant analysis:

* `pvOutputThreshold.cis` &ndash; p-value threshold for local eQTLs.
* `output_file_name.cis` &ndash; detected local eQTLs are saved in this file.
* `cisDist` &ndash; maximum distance at which gene-SNP pair is considered local.</li>
* `snpspos` &ndash; data frame with information about SNP locations, must have 3 columns - SNP name, chromosome, and position.
    * See [sample SNP location file](data/snpsloc.txt).

| snpid  | chr  |	pos  |
| ------ | ---- | ------ |
| Snp_01 | chr1 | 721289 |
| Snp_02 | chr1 | 752565 | 
| Snp_03 | chr1 | 777121 |
|  ...   | ...  |  ...   |

* `genepos` &ndash; data frame with information about gene locations, must have 4 columns - the name, chromosome, and positions of the left and right ends.
    * See [sample gene location file](data/geneloc.txt)

| geneid   | chr  | left   | right  |
| -------- | ---  | ------ | ------ |
| Gene_01  | chr1 | 721289 | 731289 |
| Gene_02  | chr1 | 752565 | 762565 |
| Gene_03  | chr1 | 777121 | 787121 |
|  ...     | ...  |  ...   |  ...   |

<!--- For more information see Matrix eQTL reference manual via command `Matrix_eQTL_main` in R or click <code><a href="html/Matrix_eQTL_main.html">Matrix_eQTL_main</a></code>.</p> --->







































