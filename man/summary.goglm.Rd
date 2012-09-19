\name{summary.goglm}
\alias{summary.goglm}
\title{Print Summary of a GOglm Object}
\usage{
  summary.goglm(obj, en.pval = 0.05, print.level = 1, ...)
}
\arguments{
  \item{obj}{An object of class \code{goglm}. Generally it
  is the output from \code{\link{goglm}}.}

  \item{en.pval}{\emph{P}-value cut-off for declaring
  enriched categories.}

  \item{print.level}{Control the amount of messages
  printed: 0 for summary information only; 1 (default) for
  adding a ranking list of GO terms; 2 for adding details
  for the enriched GO terms (from the \code{GO.db}
  package).}

  \item{...}{Other parameters (for future use).}
}
\value{
  Summaries of the dataset, a complete list of GO ranking
  list, and detailed information for the enriched
  categories (depending on the value of the
  \code{print.level} argument).
}
\description{
  Print Summary of a GOglm Object. Generally it is the
  output from \code{\link{goglm}}.
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
  \code{\link{goglm}}.
}

