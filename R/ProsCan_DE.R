##' The Prostate Cancer Data Example: DE Test Result
##'
##' We use this data example to show one of the expected inputs from
##' end-users to prepare an data frame object for GO enrichment
##' analysis using \code{GOglm}.
##'
##' The prostate cancer data consist of seven samples:
##' three from mock treated prostate cancer cells and four from
##' treated cancer cells. The data originally consisted of 49605 genes
##' annotated using Ensembl gene ID (ensGene) and NCBI Build 36.3
##' (hg18). Descriptions of data preparations can be found in the
##' additional file of Young \emph{et al.} (2010).
##'
##' We performed DE tests on the original dataset, so the data frame
##' provided here is an formatted version that the
##' \code{\link{prepare}} function can further take as the first
##' argument. The Example below shows the format of the required
##' input.
##'
##' In this illustrative example, we used \code{edgeR} (Robinson
##' \emph{et al.} (2010)) with a common dispersion estimate to obtain
##' DE test \emph{p}-values, which will be further transformed into
##' significance statistics by the \code{\link{prepare}}
##' function. Other DE testing methods based on the negative binomial
##' (NB) model for RNA-Seq data can also be adopted, such as the
##' tagwise or trend options in \code{edgeR}, the \code{NBPSeq} and
##' the \code{DESeq} approaches. All of these methods use the same
##' exact NB test for assessing DE, but differ in how they estimate
##' the dispersion parameter as a function of the mean frequency. We
##' use this example to illustrate the DE test results expected from
##' the end-users, from whatever DE testing procedures.
##'
##' @name ProsCan_DE
##'
##' @docType data
##'
##' @references Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
##' "Length bias correction in Gene Ontology enrichment analysis using
##' logistic regression", PLOS ONE, 7(10): e46128.
##'
##' Li H, Lovci M, Kwon Y, Rosenfeld M, Fu X, et al. (2008)
##' "Determination of tag density required for digital transcriptome
##' analysis: application to an androgen-sensitive prostate cancer
##' model", Proc Natl Acad Sci U S A 105: 20179-20184.
##'
##' Young M, Wakefield M, Smyth G, Oshlack A (2010) "Gene ontology
##' analysis for RNA-seq: accounting for selection bias", Genome Biol
##' 11: R14.
##'
##' Robinson M, McCarthy D, Smyth G (2010) "edgeR: a Bioconductor
##' package for differential expression analysis of digital gene
##' expression data", Bioinformatics 26: 139-140.
##'
##' @source \url{http://www.ncbi.nlm.nih.gov/pubmed/19088194}
##'
##' @keywords datasets
##'
##' @examples
##' ## Load the dataset into R session:
##' data(ProsCan_DE)
##' DE.data <- ProsCan_DE
##'
##' ## Another dataset from this package:
##' data(ProsCan_Length)
##' Length.data <- ProsCan_Length
##'
##' ## Prepare a data frame to be passed to goglm():
##' gene_table <- prepare(DE.data, Length.data, trans.p = "d.log", trans.l = TRUE)
##' ## Check first 10 rows of the data frame:
##' gene_table[1:10,1:2]
##'
##' ## We can call the summary() function:
##' summary(gene_table)
##'
##' ## We can also call the plot() function:
##' plot(gene_table)
##'
NULL
