###########################################################################
## generic function plot() and new methods for different classes
###########################################################################

# for prepGOglm class, define method plot.prepGOglm():

##' Plot DE test results (from the \code{\link{prepare}} function)
##'
##' @title Plot Two Histograms For Significance Statistics and Gene Lengths
##'
##' @param x An object of class \code{prepGOglm}
##'
##' @param ... Other parameters (for future use)
##'
##' @return A panel of two histograms for significance statistics and
##' (transformed) gene lengths, respectively.
##'
##' @method plot prepGOglm
##'
##' @author Gu Mi \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
##'
##' @seealso \code{\link{plot}}
##'
##' @examples
##' ## Load the datasets into R session:
##' data(ProsCan_DE)
##' DE.data <- ProsCan_DE
##' data(ProsCan_Length)
##' Length.data <- ProsCan_Length
##'
##' ## Prepare a data frame of class "prepGOglm":
##' gene_table <- prepare(DE.data, Length.data, trans.p = "d.log", trans.l = TRUE)
##'
##' ## Call the generic plot() function:
##' plot(gene_table)
##'
plot.prepGOglm <- function(x, ...){
    par(mfrow=c(1,2))
    hist(x[,1], main="Histogram of Siginificance Statistics",
         xlab="Significance Statistics")
    hist(x[,2], main="Historgram of Gene Lengths",
         xlab="Gene Length (in bp)")
}
