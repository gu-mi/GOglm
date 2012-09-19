\docType{package}
\name{GOglm-package}
\alias{GOglm}
\alias{GOglm-package}
\title{Length Bias Correction in Gene Ontology Enrichment Analysis Using
Logistic Regression}
\description{
  This \code{GOglm} package is a beta version under
  development. The \code{goglm} function implements the
  GOglm method discussed in Mi et al. (PLOS ONE, in press).
  The package includes a summarized RNA-Seq data example
  (the prostate cancer dataset) for methodological
  illustrations.
}
\details{
  \tabular{ll}{ Package: \tab GOglm\cr Type: \tab
  Package\cr Version: \tab 0.1.0\cr Depends: \tab goseq,
  GO.db, org.Hs.eg.db, R (>= 2.10)\cr LazyData: \tab yes\cr
  LazyLoad: \tab yes\cr License: \tab GPL-2\cr URL: \tab
  http://people.oregonstate.edu/~mig\cr Collate: \tab
  'goglm.R' 'summary.goglm.R'\cr Built: \tab R 2.15.1; ;
  2012-09-13 17:19:48 UTC; unix\cr }
}
\examples{
## Load package into R session:
library(GOglm)
}
\author{
  Gu Mi <mig@stat.oregonstate.edu>, Yanming Di
  <diy@stat.oregonstate.edu>.

  Maintainer: Gu Mi <mig@stat.oregonstate.edu>
}
\references{
  Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
  "Length bias correction in Gene Ontology enrichment
  analysis using logistic regression", PLOS ONE, in press.
}
\keyword{package}

