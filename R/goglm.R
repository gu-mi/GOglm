##' Gene Ontology (GO) Enrichment Analysis Using Logistic Regression
##'
##' This is the main function that implements the GOglm method
##' for GO enrichment analysis. It takes a data frame with two
##' columns: one for DE testing raw \emph{p}-values and the other for
##' median transcript lengths (either transformed or not). The row
##' names of the data frame must be valid gene identifiers (in the
##' Example below, we use Ensembl gene IDs, e.g. ENSG00000127954),
##' supported by the \code{getgo} function in the \code{goseq}
##' package.
##'
##' To obtain a mapping between category names and gene names, we use
##' the \code{reversemapping} (internal) function in
##' \code{goseq}. Also the \code{nullp} function in \code{goseq} is
##' used to facilitate constructing the category-to-gene list.
##'
##' @title Implement GOglm method for the GO enrichment analysis
##' @param data An extended data frame from a differential expression
##' (DE) test output. Rows are gene identifiers and columns include at
##' least the following quantities: "PValue" obtained from testing DE
##' for each gene; "length" of median transcript length for each gene
##' (in bp).
##' @param trans.p How to transform \emph{p}-values into significance
##' statistics. Value of "\code{log}" gives -log(PValue) and value of
##' "\code{d.log}" gives log(1-log(PValue)).
##' @param trans.l Whether to log-transform length (default is
##' \code{TRUE}).
##' @param genome A string identifying the genome that genes refer to
##' (cf. the \code{goseq} function in the \code{goseq} package).
##' @param id A string identifying the gene identifier used by genes
##' (cf. the \code{goseq} function in the \code{goseq} package).
##' @param cutoff \emph{P}-value cut-off for declaring DE genes
##' (default = 0.05).
##' @param min.cat minimum number of genes within a category (default
##' = 10).
##' @note In the data example (see \code{ProsCancer}), we use a subset
##' of genes for illustrations. When all genes are included in the
##' analysis, it takes longer time to get the results.
##' @return \code{goglm} returns a list with 5 entries:
##' \item{GOID}{Gene Ontology category names.}
##' \item{over.p}{\emph{P}-values for the enrichment tests.}
##' \item{mean.len}{Mean gene lengths for all categories.}
##' \item{anno}{Number of genes annotated to each category.}
##' \item{rank}{Ranking of all categories being studied.}
##' @references Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
##' "Length bias correction in Gene Ontology enrichment analysis using
##' logistic regression", PLOS ONE, in press.
##' @import goseq
##' @author  Gu Mi \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
##' @seealso \code{\link{summary.goglm}} which summarizes GOglm results and produces more readable outputs.
##' @export
##' @examples
##' ## See the example for ProsCancer.

goglm <- function(data, trans.p = c("log", "d.log"), trans.l = TRUE,
                  genome, id, cutoff = 0.05, min.cat = 10){

    trans <- match.arg(trans.p)
    sig.stat <- switch(trans,
                       log = -log(data$PValue),
                       d.log = log(1-log(data$PValue)))
    data$sig.stat <- sig.stat
    if (trans.l){
        data$log.gl <- log(data$length)
    }
    # relate gene IDs to GO IDs
    gene2cats <- getgo(rownames(data), genome, id)
    # exclude genes without GO terms
    # names of genes with annotation
    gnames <- unique(names(gene2cats)[!is.na(names(gene2cats))])
    # without annotation
    gnames.wo <- rownames(data[!rownames(data) %in% gnames, ])
    # subset data further: work on "dataset"
    dataset <- data[gnames, ]
    dataset$de.ind <- ifelse(dataset$PValue < cutoff, 1, 0)

    cat2genes <- goseq:::reversemapping(gene2cats)
    genes <- dataset$de.ind
    names(genes) <- rownames(dataset)
    n.gene <- length(genes)
    pwf <- invisible(nullp(genes, genome, id))
    cat2genenum <- relist(flesh=match(unlist(cat2genes),rownames(pwf)),
                          skeleton=cat2genes)
    n.cat <- length(cat2genenum)
    all.GOID <- names(cat2genenum) # a list of GO ID

    # main GOglm function
    over.p <- numeric(n.cat)
    mean.len <- numeric(n.cat)
    for (i in 1:n.cat){
        ones <- dim(dataset[cat2genenum[[i]], ])[1]
        y <- c(rep(1,ones), rep(0,(n.gene-ones)))
        sig.stat <- c(dataset[cat2genenum[[i]],"sig.stat"],
               dataset[-cat2genenum[[i]],"sig.stat"])
        length <- c(dataset[cat2genenum[[i]],"log.gl"],
            dataset[-cat2genenum[[i]],"log.gl"])
        glm.fit <- glm(y ~ sig.stat + length, family=quasibinomial,
                 control = list(maxit = 50))
        sumry <- summary(glm.fit)
        over.p[i] <- sign(coefficients(glm.fit)[[2]])*
            sumry$coef[ ,"Pr(>|t|)"][[2]]
        mean.len[i] <- round(mean(dataset[cat2genenum[[i]],"length"],
                                  na.rm=TRUE))
    }
    over.de.p <- over.p
    over.de.p[over.de.p < 0] = 1
    anno <- as.vector(sapply(cat2genenum, length))
    GO.info <- data.frame(GOID = all.GOID, over.p = over.de.p,
                          mean.len = mean.len, anno = anno)
    head(GO.info)
    rownames(GO.info) <- unfactor(GO.info$GOID)
    GOglm = GO.info[GO.info$anno > min.cat, ]
    ord.GOglm = GOglm[order(GOglm$over.p), ]
    ord.GOglm$rank = seq(1,dim(ord.GOglm)[1])
    class(ord.GOglm) <- "goglm"
    return(invisible(ord.GOglm))
}
