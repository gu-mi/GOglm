\name{summary.prepGOglm}
\alias{summary.prepGOglm}
\title{Summerize Data Preparation (DE Test)}
\usage{
  \method{summary}{prepGOglm} (x, ...)
}
\arguments{
  \item{x}{An object of class \code{prepGOglm}}

  \item{...}{Other parameters (for future use)}
}
\value{
  Some descriptive summaries based on the results from the
  \code{\link{prepare}} function
}
\description{
  Summerize DE test results (from the \code{\link{prepare}}
  function)
}
\examples{
## Load the datasets into R session:
data(ProsCan_DE)
DE.data <- ProsCan_DE
data(ProsCan_Length)
Length.data <- ProsCan_Length

## Prepare a data frame to be passed to goglm():
gene_table <- prepare(DE.data, Length.data, trans.p = "d.log", trans.l = TRUE)

## For a summary of the DE test results:
summary(gene_table)
}
\author{
  Gu Mi \email{mig@stat.oregonstate.edu}, Yanming Di
  \email{diy@stat.oregonstate.edu}
}
\seealso{
  \code{\link{summary}}, \code{\link{summary.goglm}}
}

