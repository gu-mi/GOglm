##' Gene Ontology (GO) Enrichment Analysis Using Logistic Regression
##'
##' This is the main function that implements the GOglm method for GO
##' enrichment analysis using logistic regression. Users need to
##' specify three arguments, which will be illustrated in the Argument
##' and Examples sections below. A DE test output with DE
##' \emph{p}-values and gene length information, and a
##' category-to-gene mapping list, are required to implement
##' \code{goglm}. In general, the DE test output is obtained by the
##' \code{\link{prepare}} function, and the mapping list can be
##' obtained by reverse-mapping the results from the \code{getgo}
##' function in the \code{goseq}.
##'
##' @title Implement GOglm method for GO enrichment analysis
##'
##' @param gene_data Output from the \code{\link{prepare}}
##' function. It contains valid gene identifiers as row names. Two
##' columns are (1) (transformed) DE test \emph{p}-values (significance
##' statistics) and (2) (transformed) gene lengths.
##'
##' @param cat2genes A list. Entry names are GO terms, and elements
##' are corresponding gene names. This mapping is obtained by the
##' \code{getgo} and \code{revMap} functions.
##'
##' @param n If a category has fewer than \emph{n} genes annotated,
##' then this cagtegory will be excluded in the final GO ranking list.
##'
##' @return An object of class \code{goglm} to be passed to
##' \code{summary} for more readable results. See Examples below.
##'
##' @references Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
##' "Length bias correction in Gene Ontology enrichment analysis using
##' logistic regression", PLOS ONE, in press.
##'
##' @author Gu Mi \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
##'
##' @seealso \code{\link{summary.goglm}} which summarizes GOglm
##' results and produces more readable outputs.
##'
##' @importFrom goseq getgo
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
##' ## For more readable outputs:
##' output <- cbind(res$over.p, res$anno, res$rank)
##' rownames(output) <- unfactor(res$GOID)
##' colnames(output) <- c("over.p", "n.anno", "rank")
##' head(output)
##'
##' ## For a summary of the GOglm results:
##' summary(res)
##'

goglm <- function(gene_data, cat2genes, n=5){
    n.cat <- length(cat2genes)
    n.gene <- length(gene2cats)
    over.p <- numeric(n.cat)
    mean.len <- numeric(n.cat)

    for (i in 1:n.cat){

        # count how many genes are in the i^th category
        ng.incat <- dim(gene_data[cat2genes[[i]], ])[1]

        # add a column to gene_data "y" to indicate category presence
        g_data_ext <- as.data.frame(cbind(gene_data, y = as.numeric(rownames(gene_data) %in% cat2genes[[i]])))

        glm.fit <- glm(y ~ Sig.stat + Length, family = "quasibinomial",
                       data = g_data_ext, control = list(maxit = 50))
        sumry <- summary(glm.fit)
        over.p[i] <- sign(coefficients(glm.fit)[[2]])*
            sumry$coef[ ,"Pr(>|t|)"][[2]]
    }

    over.de.p <- over.p
    over.de.p[over.de.p < 0] = 1
    anno <- as.vector(sapply(cat2genes, length))
    GO.info <- data.frame(GOID = names(cat2genes), over.p = over.de.p,
                          anno = anno)
    rownames(GO.info) <- unfactor(GO.info$GOID)
    sub.GO.info <- GO.info[GO.info$anno > n, ]
    ord.GO.info <- sub.GO.info[order(sub.GO.info$over.p), ]
    ord.GO.info$rank <- seq(1, dim(ord.GO.info)[1])
    class(ord.GO.info) <- "goglm"
    return(invisible(ord.GO.info))
}

