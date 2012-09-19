\name{goglm}
\alias{goglm}
\title{Implement GOglm method for the GO enrichment analysis}
\usage{
  goglm(data, trans.p = c("log", "d.log"), trans.l = TRUE,
    genome, id, cutoff = 0.05, min.cat = 10)
}
\arguments{
  \item{data}{An extended data frame from a differential
  expression (DE) test output. Rows are gene identifiers
  and columns include at least the following quantities:
  "PValue" obtained from testing DE for each gene; "length"
  of median transcript length for each gene (in bp).}

  \item{trans.p}{How to transform \emph{p}-values into
  significance statistics. Value of "\code{log}" gives
  -log(PValue) and value of "\code{d.log}" gives
  log(1-log(PValue)).}

  \item{trans.l}{Whether to log-transform length (default
  is \code{TRUE}).}

  \item{genome}{A string identifying the genome that genes
  refer to (cf. the \code{goseq} function in the
  \code{goseq} package).}

  \item{id}{A string identifying the gene identifier used
  by genes (cf. the \code{goseq} function in the
  \code{goseq} package).}

  \item{cutoff}{\emph{P}-value cut-off for declaring DE
  genes (default = 0.05).}

  \item{min.cat}{minimum number of genes within a category
  (default = 10).}
}
\value{
  \code{goglm} returns a list with 5 entries:
  \item{GOID}{Gene Ontology category names.}
  \item{over.p}{\emph{P}-values for the enrichment tests.}
  \item{mean.len}{Mean gene lengths for all categories.}
  \item{anno}{Number of genes annotated to each category.}
  \item{rank}{Ranking of all categories being studied.}
}
\description{
  Gene Ontology (GO) Enrichment Analysis Using Logistic
  Regression
}
\details{
  This is the main function that implements the GOglm
  method for GO enrichment analysis. It takes a data frame
  with two columns: one for DE testing raw \emph{p}-values
  and the other for median transcript lengths (either
  transformed or not). The row names of the data frame must
  be valid gene identifiers (in the Example below, we use
  Ensembl gene IDs, e.g. ENSG00000127954), supported by the
  \code{getgo} function in the \code{goseq} package.

  To obtain a mapping between category names and gene
  names, we use the \code{reversemapping} (internal)
  function in \code{goseq}. Also the \code{nullp} function
  in \code{goseq} is used to facilitate constructing the
  category-to-gene list.
}
\note{
  In the data example (see \code{ProsCancer}), we use a
  subset of genes for illustrations. When all genes are
  included in the analysis, it takes longer time to get the
  results.
}
\examples{
## See the example for ProsCancer.
}
\author{
  Gu Mi \email{mig@stat.oregonstate.edu}, Yanming Di
  \email{diy@stat.oregonstate.edu}
}
\references{
  Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
  "Length bias correction in Gene Ontology enrichment
  analysis using logistic regression", PLOS ONE, in press.
}
\seealso{
  \code{\link{summary.goglm}} which summarizes GOglm
  results and produces more readable outputs.
}

