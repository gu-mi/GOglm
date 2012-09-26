###########################################################################
## generic function summary() and new methods for different classes
###########################################################################

# generic function summary():

##' A generic function for summerizing data preparation and GOglm results
##'
##' @title Summerize Data Preparation and GOglm Results
##'
##' @param x An object of class \code{prepGOglm} or \code{goglm}
##'
##' @param ... Other parameters (for future use)
##'
##' @return Some descriptive summaries depending on the object passed
##' to \code{x}
##'
##' @S3method summary prepGOglm
##'
##' @S3method summary goglm
##'
##' @seealso \code{\link{summary.prepGOglm}}, \code{\link{summary.goglm}}
##' @author Gu Mi \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
##'
##' @export
##'
##' @examples
##' ## Load the datasets into R session:
##' data(ProsCan_DE)
##' DE.data <- ProsCan_DE
##' data(ProsCan_Length)
##' Length.data <- ProsCan_Length
##'
##' ## Prepare a data frame to be passed to goglm():
##' gene_table <- prepare(DE.data, Length.data, trans.p = "d.log", trans.l = TRUE)
##'
##' ## For a summary of the DE test results:
##' summary(gene_table)
##'
##' ## For illustration, only consider a subset of genes:
##' gene_data <- gene_table[1:100,1:2]
##'
##' ## Prepare the "category-to-genes" list:
##' library(goseq)
##' gene2cats <- getgo(rownames(gene_data), "hg18", "ensGene")
##' cat2genes <- revMap(gene2cats)
##'
##' ## Run goglm():
##' res <- goglm(gene_data, cat2genes, n=5)
##' names(res)  # "GOID"   "over.p" "anno"   "rank"
##'
##' ## For a summary of the GOglm results:
##' summary(res)
##'
summary <- function(x, ...){
    UseMethod("summary")
}
