###########################################################################
## generic function summary() and new methods for different classes
###########################################################################

# for prepGOglm class, define method summary.prepGOglm():

##' Summerize DE test results (from the \code{\link{prepare}}
##' function)
##'
##' @title Summerize Data Preparation (DE Test)
##'
##' @param x An object of class \code{prepGOglm}
##'
##' @param ... Other parameters (for future use)
##'
##' @return Some descriptive summaries based on the results from the
##' \code{\link{prepare}} function
##'
##' @method summary prepGOglm
##' @rdname summary.prepGOglm
##' @export
##'
##' @seealso \code{\link{summary}}, \code{\link{summary.goglm}}
##'
##' @author Gu Mi \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
##'
##' @examples
##' ## Load the datasets into R session:
##' data(ProsCan_DE)
##' DE_data <- ProsCan_DE
##' data(ProsCan_Length)
##' Length_data <- ProsCan_Length
##'
##' ## Prepare a data frame to be passed to goglm():
##' gene_table <- prepare(DE_data, Length_data, trans.p = "d.log", trans.l = TRUE)
##'
##' ## For a summary of the DE test results:
##' summary(gene_table)
##'
summary.prepGOglm <- function(x, ...){
    n.gene <- dim(x)[1]
    na.dep <- sum(is.na(x[,1]))
    na.len <- sum(is.na(x[,2]))
    cat("-------------------------------------------------------------- \n")
    cat("| Total number of genes under study is", n.gene, "\n")
    cat("| There are", na.dep, "genes without DE test p-values",
        "(", round(na.dep/n.gene*100,2), "% ) \n")
    cat("| There are", na.len, "genes without length information",
        "(", round(na.len/n.gene*100,2), "% ) \n")
    cat("| Quantiles of significance statistics are",
        round(quantile(x[,1], na.rm=TRUE),2), "\n")
    cat("| Quantiles of gene lengths (in bp) are",
        round(quantile(x[,2], na.rm=TRUE),2), "\n")
    cat("-------------------------------------------------------------- \n")
    cat("| Please make sure that significance statistics and gene lengths do not require any further transformations. \n")
}
