##' The Prostate Cancer Data Example: Gene Length Information
##'
##' Similar to \code{\link{ProsCan_DE}}, we include in this dataset
##' each gene's length to be passed as the second argument to the
##' \code{\link{prepare}} function.
##'
##' A gene's length is defined as the median length of all its
##' corresponding mature transcripts. In the \code{goseq}  package,
##' length data are obtained from the UCSC genome browser for each
##' combination of \code{genome} and \code{id}. In \code{GOglm}, we
##' assume that end-users already have the length information
##' available, so that a column of gene lengths is properly matched to
##' the gene names (row names) in the data frame. For an illustration
##' on accessing gene length information, please see the Example
##' section of \code{\link{goglm}}.
##'
##' @name ProsCan_Length
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
##' data(ProsCan_Length)
##' Length_data <- ProsCan_Length
##'
##' ## Another dataset from this package:
##' data(ProsCan_DE)
##' DE_data <- ProsCan_DE
##'
##' ## Prepare a data frame to be passed to goglm():
##' gene_table <- prepare(DE_data, Length_data, trans.p = "d.log", trans.l = TRUE)
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
