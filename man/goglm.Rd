\name{goglm}
\alias{goglm}
\title{Implement the GOglm method for GO enrichment analysis}
\usage{
  goglm(gene_data, cat2genes, n = 5)
}
\arguments{
  \item{gene_data}{Output from the \code{\link{prepare}}
  function. It contains valid gene identifiers as row
  names. Two columns are (1) (transformed) DE test
  \emph{p}-values (significance statistics) and (2)
  (transformed) gene lengths.}

  \item{cat2genes}{A list. Entry names are GO terms, and
  elements are corresponding gene names. This mapping is
  obtained by the \code{getgo} and \code{revMap}
  functions.}

  \item{n}{If a category has fewer than \emph{n} genes
  annotated, then this cagtegory will be excluded in the
  final GO ranking list.}
}
\value{
  An object of class \code{goglm} to be passed to
  \code{summary} for more readable results. See Examples
  below.
}
\description{
  Gene Ontology (GO) Enrichment Analysis Using Logistic
  Regression
}
\details{
  This is the main function that implements the GOglm
  method for GO enrichment analysis using logistic
  regression. Users need to specify three arguments, which
  will be illustrated in the Argument and Examples sections
  below. A DE test output with DE \emph{p}-values and gene
  length information, and a category-to-gene mapping list,
  are required to implement \code{goglm}. In general, the
  DE test output is obtained by the \code{\link{prepare}}
  function, and the mapping list can be obtained by
  reverse-mapping the results from the \code{getgo}
  function in the \code{goseq}.
}
\examples{
## Load the datasets into R session:
data(ProsCan_DE)
DE_data <- ProsCan_DE
data(ProsCan_Length)
Length_data <- ProsCan_Length

## Prepare a data frame to be passed to goglm():
gene_table <- prepare(DE_data, Length_data, trans.p = "d.log", trans.l = TRUE)

## For illustration, only consider a subset of genes:
gene_data <- gene_table[1:100,1:2]

## Prepare the "category-to-genes" list:
library(goseq)
gene2cats <- getgo(rownames(gene_data), "hg18", "ensGene")
cat2genes <- revMap(gene2cats)

## Run goglm():
res <- goglm(gene_data, cat2genes, n=5)
names(res)  # "GOID"   "over.p" "anno"   "rank"

## For more readable outputs:
output <- cbind(res$over.p, res$anno, res$rank)
rownames(output) <- unfactor(res$GOID)
colnames(output) <- c("over.p", "n.anno", "rank")
head(output)

## For a summary of the GOglm results:
summary(res)
}
\author{
  Gu Mi \email{mig@stat.oregonstate.edu}, Yanming Di
  \email{diy@stat.oregonstate.edu}
}
\references{
  Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
  "Length bias correction in Gene Ontology enrichment
  analysis using logistic regression", PLOS ONE, 7(10):
  e46128
}
\seealso{
  \code{\link{summary.goglm}} which summarizes GOglm
  results and produces more readable outputs.
}

