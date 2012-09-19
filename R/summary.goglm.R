##' Print Summary of a GOglm Object. Generally it is the
##' output from \code{\link{goglm}}.
##'
##' @title Print Summary of a GOglm Object
##' @param obj An object of class \code{goglm}. Generally it is the
##' output from \code{\link{goglm}}.
##' @param en.pval \emph{P}-value cut-off for declaring enriched categories.
##' @param print.level Control the amount of messages printed: 0 for
##' summary information only; 1 (default) for adding a ranking list of
##' GO terms; 2 for adding details for the enriched GO terms (from the
##' \code{GO.db} package).
##' @param ... Other parameters (for future use).
##' @return Summaries of the dataset, a complete list of GO ranking
##' list, and detailed information for the enriched categories
##' (depending on the value of the \code{print.level} argument).
##' @references Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
##' "Length bias correction in Gene Ontology enrichment analysis using
##' logistic regression", PLOS ONE, in press.
##' @import GO.db
##' @author Gu Mi \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
##' @seealso \code{\link{goglm}}.
##' @export
##' @examples
##' ## See the example for ProsCancer.
##'
summary.goglm <- function(obj, en.pval = 0.05, print.level = 1, ...){

    ## some basic summaries...
    cat("We test a total of", length(obj$GOID),"GO categories", "\n")
    cat("-----------------------------------------------------", fill=TRUE)
    cat("The number of enriched categories is (p-value cut-off at",
        en.pval,")", sum(obj$over.p < en.pval), "\n")
    cat("-----------------------------------------------------", fill=TRUE)
    cat("The number of annotated genes ranges from",
        min(obj$anno), "to", max(obj$anno), "\n")
    cat("-----------------------------------------------------", fill=TRUE)
    cat("The mean gene length for these categories ranges from",
        min(obj$mean.len), "to", max(obj$mean.len), "\n")
    cat("-----------------------------------------------------", fill=TRUE)

    ## a complete list of GO terms...
    rank.list <- cbind(over.pvals = obj$over.p,
                       mean.length = obj$mean.len,
                       No.anno = obj$anno,
                       ranking = obj$rank)
    rownames(rank.list) <- obj$GOID

    if (print.level == 1){
        cat("The GO ranking list is printed below:", "\n")
        cat("-----------------------------------------------------",
            fill=TRUE)
        return(rank.list)
    }

    if (print.level == 2){
        cat("The GO ranking list is printed below:", "\n")
        cat("-----------------------------------------------------",
            fill=TRUE)
        print(rank.list)
        cat("-----------------------------------------------------",
            fill=TRUE)
        cat("More info. on the", sum(obj$over.p < en.pval),
            "enriched GO categories:", "\n")
        m <- sum(obj$over.p < en.pval)
        tops <- obj$GOID
        # GOTERM is in the GO.db package; print out GO term details:
        # concise version:
        for (i in 1:m) {
            cat("|--------------------",i,"--------------------|","\n")
            print(GOID(GOTERM[[as.character(tops[i])]]))
            print(Term(GOTERM[[as.character(tops[i])]]))
            print(Ontology(GOTERM[[as.character(tops[i])]]))
        }
    }
}
