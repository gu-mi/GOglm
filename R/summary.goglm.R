###########################################################################
## generic function summary() and new methods for different classes
##########################################################################

# for goglm class, define method summary.goglm():

##' Summerize GOglm results (from the \code{\link{goglm}} function)
##'
##' @title Summerize GOglm Results
##'
##' @param x An object of class \code{goglm}
##'
##' @param en.cut \emph{P}-value cut-off for declaring enriched
##' categories (default = 0.05)
##'
##' @param ... Other parameters (for future use)
##'
##' @method summary goglm
##' @rdname summary.goglm
##' @export
##'
##' @return Some descriptive summaries based on the \code{goglm} result
##'
##' @seealso \code{\link{summary.prepGOglm}}, \code{\link{summary}}
##'
##' @author  Gu Mi \email{mig@@stat.oregonstate.edu}, Yanming Di
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
summary.goglm <- function(x, en.cut = 0.05, ...){
    n.cat <- length(x$GOID)
    n.anno <- x$anno
    n.en.cat <- sum(x$over.p < en.cut)
    cat("-------------------------------------------------------------- \n")
    cat("| Total number of categories under study is", n.cat, "\n")
    cat("| The number of genes annotated to these categories ranges from",
        min(n.anno), "to", max(n.anno), "\n")
    cat("| Under", en.cut,"cut-off, the number of enriched categories is",
        n.en.cat, "\n")
    cat("-------------------------------------------------------------- \n")
}
