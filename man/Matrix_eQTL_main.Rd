\name{Matrix_eQTL_main}
\alias{Matrix_eQTL_main}
\alias{Matrix_eQTL_engine}
\title{
    Main function for fast eQTL analysis in MatrixEQTL package
}
\description{
    \code{Matrix_eQTL_engine} function tests association of
    every row of the \code{snps} dataset with every row of the
    \code{gene} dataset using a linear regression model
    defined by the \code{useModel} parameter (see below).

    The testing procedure accounts for extra 
    covariates specified by the \code{cvrt} parameter.
    
    The \code{errorCovariance} parameter can be set to the 
    error variance-covariance matrix to account 
    for heteroskedastic and/or correlated errors.
    
    Associations significant at \code{pvOutputThreshold} 
    (\code{pvOutputThreshold.cis}) levels are saved to 
    \code{output_file_name} (\code{output_file_name.cis}), 
    with corresponding estimates of effect size (slope coefficient),
    test statistics, p-values, and q-values (false discovery rate).

    Matrix eQTL can perform separate analysis for
    local (cis) and distant (trans) eQTLs.
    For such analysis one has to set the cis-analysis specific
    parameters \code{pvOutputThreshold.cis > 0},
    \code{cisDist}, \code{snpspos} and 
    \code{genepos} in the call of \code{Matrix_eQTL_main} function.
    A gene-SNP pair is considered local if the 
    distance between them is less or equal to \code{cisDist}.
    The genomic location of genes and SNPs is defined by 
    the data frames \code{snpspos} and \code{genepos}.
    Depending on p-value thresholds \code{pvOutputThreshold} and
    \code{pvOutputThreshold.cis} Matrix eQTL runs in 
    one of three different modes:
    \itemize{
        \item Set \code{pvOutputThreshold > 0} and 
            \code{pvOutputThreshold.cis = 0} (or use \code{Matrix_eQTL_engine})
            to perform eQTL analysis without using gene/SNP locations.
            Associations significant at the \code{pvOutputThreshold} level
            are be recorded in \code{output_file_name}
            and in the returned object.
        \item Set \code{pvOutputThreshold = 0} and
            \code{pvOutputThreshold.cis > 0} to perform eQTL analysis for
            local gene-SNP pairs only. Local associations significant at
            \code{pvOutputThreshold.cis} level will be recorded in
            \code{output_file_name.cis} and in the returned object.
        \item Set \code{pvOutputThreshold > 0} and
            \code{pvOutputThreshold.cis > 0} to perform eQTL analysis
            with separate p-value thresholds for local and distant eQTLs.
            Distant and local associations significant at corresponding
            thresholds are recorded in \code{output_file_name} and
            \code{output_file_name.cis} respectively and in the returned object.
            In this case the false discovery rate is calculated
            separately for these two sets of eQTLs.
    }
    
    \code{Matrix_eQTL_engine} is a wrapper for \code{Matrix_eQTL_main}
    for eQTL analysis without regard to gene/SNP location and provided
    for compatibility with the previous versions of the package.

    The parameter \code{pvalue.hist} allows to record information sufficient
    to create a histogram or QQ-plot of all the p-values
    (see \code{\link[=plot.MatrixEQTL]{plot}}).
}
\usage{
Matrix_eQTL_main(
        snps, 
        gene, 
        cvrt = SlicedData$new(), 
        output_file_name = "", 
        pvOutputThreshold = 1e-5,
        useModel = modelLINEAR, 
        errorCovariance = numeric(), 
        verbose = TRUE, 
        output_file_name.cis = "", 
        pvOutputThreshold.cis = 0,
        snpspos = NULL, 
        genepos = NULL,
        cisDist = 1e6,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)

Matrix_eQTL_engine(
        snps, 
        gene, 
        cvrt = SlicedData$new(), 
        output_file_name, 
        pvOutputThreshold = 1e-5, 
        useModel = modelLINEAR, 
        errorCovariance = numeric(), 
        verbose = TRUE,
        pvalue.hist = FALSE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)
}
\arguments{
    \item{snps}{
        \code{\linkS4class{SlicedData}} object with genotype information. 
        Can be real-valued for linear models and
        must take at most 3 distinct values for ANOVA unless
        the number of ANOVA categories is set to
        a higher number (see \code{useModel} parameter).
    }
    \item{gene}{
        \code{\linkS4class{SlicedData}} object with gene expression information.
        Must have columns matching those of \code{snps}.
    }
    \item{cvrt}{
        \code{\linkS4class{SlicedData}} object with additional covariates. 
        Can be an empty \code{SlicedData} object in case of no covariates.
        The constant is always included in the model and
        would cause an error if included in \code{cvrt}.
        The order of columns must match those in \code{snps} and \code{gene}.
    }
    \item{output_file_name}{
        \code{character}, \code{connection}, or \code{NULL}.
        If not \code{NULL}, significant associations are saved to
        this file (all significant associations if \code{pvOutputThreshold=0}
        or only distant if \code{pvOutputThreshold>0}).
        If the file with this name exists, it is overwritten.
    }
    \item{output_file_name.cis}{
        \code{character}, \code{connection}, or \code{NULL}.
        If not \code{NULL}, significant local associations
        are saved to this file.
        If the file with this name exists, it is overwritten.
    }
    \item{pvOutputThreshold}{
        \code{numeric}. Significance threshold for all/distant tests.
    }
    \item{pvOutputThreshold.cis}{ 
        \code{numeric}. Same as \code{pvOutputThreshold}, but for local eQTLs.
    }
    \item{useModel}{
        \code{integer}. 
            Eigher \code{modelLINEAR}, 
            \code{modelANOVA}, or
            \code{modelLINEAR_CROSS}.
        \enumerate{
            \item Set \code{useModel = \link{modelLINEAR}} to model
                the effect of the genotype as additive linear and
                test for its significance using t-statistic.
            \item Set \code{useModel = \link{modelANOVA}} to treat genotype
                as a categorical variables and use ANOVA model and
                test for its significance using F-test.
                The default number of ANOVA categories is 3.
                Set otherwise like this:
                \code{options(MatrixEQTL.ANOVA.categories=4)}.
            \item Set \code{useModel = \link{modelLINEAR_CROSS}} to add
                a new term to the model equal to the product of genotype and
                the last covariate; the significance of this term is
                then tested using t-statistic.
        }
    }
    \item{errorCovariance}{
        \code{numeric}. The error covariance matrix.
        Use \code{numeric()} for homoskedastic independent errors. 
    }
    \item{verbose}{
        \code{logical}. Set to \code{TRUE} to display
        more detailed report about the progress.
    }
    \item{snpspos}{
        \code{data.frame} object with information about SNP locations,
        must have 3 columns - SNP name, chromosome, and position, like this:
        \tabular{ccc}{
            snpid   \tab chr    \tab pos     \cr
            Snp_01  \tab 1      \tab 721289  \cr
            Snp_02  \tab 1      \tab 752565  \cr
            \ldots  \tab \ldots \tab \ldots  \cr
        }
    }
    \item{genepos}{
        \code{data.frame} with information about transcript locations,
        must have 4 columns - the name, chromosome, and 
        positions of the left and right ends, like this:
        \tabular{cccc}{
            geneid  \tab  chr    \tab left   \tab right  \cr
            Gene_01 \tab   1     \tab 721289 \tab 731289 \cr
            Gene_02 \tab   1     \tab 752565 \tab 762565 \cr
            \ldots  \tab \ldots  \tab \ldots \tab \ldots \cr
        }
    }
    \item{cisDist}{
        \code{numeric}. SNP-gene pairs within this distance are 
        considered local. The distance is measured from the
        nearest end of the gene. SNPs within a gene are always considered local.
    }
    \item{pvalue.hist}{
        \code{logical}, \code{numerical}, or \code{"qqplot"}.
        Defines whether and how the distribution of p-values is recorded
        in the returned object.
        If \code{pvalue.hist = FALSE}, the information is not recorded
        and the analysis is performed faster. 
        Alternatively, set \code{pvalue.hist = "qqplot"} to record information
        sufficient to create a QQ-plot of the p-values
        (use \code{\link[=plot.MatrixEQTL]{plot}} on the returned
        object to create the plot).
        To record information for a histogram set \code{pvalue.hist}
        to the desired number of bins of equal size.
        Finally, \code{pvalue.hist} can also be set to a
        custom set of bin edges.
    }
    \item{min.pv.by.genesnp}{
        \code{logical}. Set \code{min.pv.by.genesnp = TRUE} to record
        the minimum p-value for each SNP and each gene in the returned object.
        The minimum p-values are recorded even if if they are above the
        corresponding thresholds of \code{pvOutputThreshold} and
        \code{pvOutputThreshold.cis}.
        The analysis runs faster when the parameter is set to \code{FALSE}.
    }
    \item{noFDRsaveMemory}{
        \code{logical}. Set \code{noFDRsaveMemory = TRUE} to save
        significant gene-SNP pairs directly to the output files,
        reduce memory footprint and skip FDR calculation.
        The eQTLs are not recorded in the returned object
        if \code{noFDRsaveMemory = TRUE}.
    }
}
\details{
    Note that the columns of 
    \code{gene}, \code{snps}, and \code{cvrt} must match.
    If they do not match in the input files,
    use \code{ColumnSubsample} method to subset and/or reorder them.
}
\value{
    The detected eQTLs are saved in \code{output_file_name}
    and/or \code{output_file_name.cis} if they are not \code{NULL}.
    The method also returns a list with a summary of the performed analysis.
    \item{param}{
        Keeps all input parameters and also records
        the number of degrees of freedom for the full model.
    }
    \item{time.in.sec}{
        Time difference between the start and 
        the end of the analysis (in seconds).
    }
    \item{all}{
        Information about all detected eQTLs.
    }
    \item{cis}{
        Information about detected local eQTLs.
    }
    \item{trans}{
        Information about detected distant eQTLs.
    }
    The elements \code{all}, \code{cis}, and \code{trans}
    may contain the following components
    \describe{
        \item{\code{ntests}}{
            Total number of tests performed. This is used for FDR calculation.
        }
        \item{\code{eqtls}}{
            Data frame with recorded significant associations.
            Not available if \code{noFDRsaveMemory=FALSE}
        }
        \item{\code{neqtls}}{
            Number of significant associations recorded.
        }
        \item{\code{hist.bins}}{
            Histogram bins used for recording p-value distribution.
            See \code{pvalue.hist} parameter.
        }
        \item{\code{hist.counts}}{
            Number of p-value that fell in each histogram bin.
            See \code{pvalue.hist} parameter.
        }
        \item{\code{min.pv.snps}}{
            Vector with the best p-value for each SNP.
            See \code{min.pv.by.genesnp} parameter.
        }
        \item{\code{min.pv.gene}}{
            Vector with the best p-value for each gene.
            See \code{min.pv.by.genesnp} parameter.
        }
    }
}
\references{
    The package website:
    \url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\seealso{
    The code below is the sample code for eQTL analysis NOT using gene/SNP locations.    
    
    See \code{\link{MatrixEQTL_cis_code}} for sample code for
    eQTL analysis that separates local and distant tests.
}
\examples{
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste0(base.dir, "/data/SNP.txt");

# Gene expression file name
expression_file_name = paste0(base.dir, "/data/GE.txt");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste0(base.dir, "/data/Covariates.txt");

# Output file name
output_file_name = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold = 1e-2;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0){
    cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

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
    noFDRsaveMemory = FALSE);
        
unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)
}
