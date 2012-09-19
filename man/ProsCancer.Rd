\docType{data}
\name{ProsCancer}
\alias{ProsCancer}
\title{The Prostate Cancer Data Example}
\format{A 22743 by 9 data frame. Each row has a valid Ensembl gene ID
as its name which will be mapped to GO categories. Initial DE testing
outputs (log fold changes, \emph{p}-values, etc.) and gene length
information are also included in this data frame.}
\source{
  \url{http://www.ncbi.nlm.nih.gov/pubmed/19088194}
}
\description{
  The prostate cancer data consist of seven samples: three
  from mock treated prostate cancer cells and four from
  treated cancer cells. The data originally consisted of
  49605 genes annotated using Ensembl gene ID (ensGene) and
  NCBI Build 36.3 (hg18). Descriptions of data preparations
  can be found in the additional file of Young \emph{et
  al.} (2010).
}
\details{
  We performed DE tests on the original dataset, so the
  data frame provided here is an extended formatted version
  that \code{\link{goglm}} can further take as the first
  argument. See the Examples below.

  Essentially a data frame is required to be passed as the
  first argument to \code{\link{goglm}}. This data frame
  must have valid gene identifiers as row names to be
  mapped to GO terms. Two columns, DE testing
  \emph{p}-values and gene lengths, must be available for
  each gene. Missing values are allowed and will be coded
  as \code{NA}.

  In this illustrative example, we used \code{edgeR}
  (Robinson \emph{et al.} (2010)) with a common dispersion
  estimate to obtain DE test \emph{p}-values, which will be
  further transformed into significance statistics. Other
  DE testing methods based on the negative binomial model
  for RNA-Seq data can also be adopted, such as the tagwise
  or trend options in \code{edgeR}, the \code{NBPSeq} and
  the \code{DESeq} approaches. All of these methods use the
  same exact NB test for assessing DE, but differ in how
  they estimate the dispersion parameter as a function of
  the mean frequency.

  A gene's length is defined as the median length of all
  its corresponding mature transcripts. In the \code{goseq}
  package, length data are obtained from the UCSC genome
  browser for each combination of \code{genome} and
  \code{id}. In \code{GOglm}, we assume that end-users
  already have the length information available, so that a
  column of gene lengths is properly matched to the gene
  names (row names) in a data frame.

  We used the \code{org.Hs.egGO2ALLEGS} R object in the
  Bioconductor annotation package \code{org.Hs.eg.db} to
  obtain mappings between a given GO category and all its
  annotated Entrez Gene identifiers. Gene length
  information is also accessible from the Ensembl Project
  online.
}
\examples{
## Load the prostate cancer dataset:
data(ProsCancer)
head(ProsCancer)

## Make sure row names of dataset are valid gene identifiers:
rownames(ProsCancer) <- as.character(ProsCancer$X)

## For reproducible results:
set.seed(999)

## Sample a subset of data (500 genes indices):
idx <- sample(1:dim(ProsCancer)[1], 500)

## Essentially we only need two columns: p-values and lengths:
psc.subset <- ProsCancer[idx, c("PValue", "length")]
head(psc.subset)

## Run GOglm:
result <- goglm(data = psc.subset, trans.p = "d.log", genome = "hg18",
id = "ensGene")
## check list entry names:
names(result)
# "GOID"  "over.p"  "mean.len"  "anno"  "rank"

## Print GOglm results with default messages:
summary(result, en.pval = 0.05, print.level = 1)
## ... or with all messages:
summary(result, en.pval = 0.05, print.level = 2)
}
\references{
  Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
  "Length bias correction in Gene Ontology enrichment
  analysis using logistic regression", PLOS ONE, in press.

  Li H, Lovci M, Kwon Y, Rosenfeld M, Fu X, et al. (2008)
  "Determination of tag density required for digital
  transcriptome analysis: application to an
  androgen-sensitive prostate cancer model", Proc Natl Acad
  Sci U S A 105: 20179-20184.

  Young M, Wakefield M, Smyth G, Oshlack A (2010) "Gene
  ontology analysis for RNA-seq: accounting for selection
  bias", Genome Biol 11: R14.

  Robinson M, McCarthy D, Smyth G (2010) "edgeR: a
  Bioconductor package for differential expression analysis
  of digital gene expression data", Bioinformatics 26:
  139-140.
}
\keyword{datasets}

