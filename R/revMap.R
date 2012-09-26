##' Reverse a list of genes-to-categories to a list of categories-to-genes
##'
##' This function (\code{reversemapping}) is contained in the
##' \code{goseq} package but not exported. We use this function to
##' facilitate reverse-mapping of the gene-to-category list.
##'
##' @title Reverse Mapping
##'
##' @param map Output from the \code{getgo} function in \code{goseq}
##'
##' @return A list with the reversed \code{gene2cats}, i.e., entry
##' names are GO terms, and elements are corresponding gene names
##'
##' @note This function is written by Matthew D. Young and is used in
##' the \code{goseq} package. See \code{goseq:::reversemapping}.
##'
##' @author Matthew D. Young \email{myoung@@wehi.edu.au}, Gu Mi
##' \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
##'
##' @export
##'

revMap <- function (map)
{
    tmp <- unlist(map, use.names = FALSE)
    names(tmp) <- rep(names(map), times = as.numeric(summary(map)[,
        1]))
    return(split(names(tmp), as.vector(tmp)))
}
