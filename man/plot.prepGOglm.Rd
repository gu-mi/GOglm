\name{plot.prepGOglm}
\alias{plot.prepGOglm}
\title{Plot Two Histograms For Significance Statistics and Gene Lengths}
\usage{
  \method{plot}{prepGOglm} (x, ...)
}
\arguments{
  \item{x}{An object of class \code{prepGOglm}}

  \item{...}{Other parameters (for future use)}
}
\value{
  A panel of two histograms for significance statistics and
  (transformed) gene lengths, respectively.
}
\description{
  Plot DE test results (from the \code{\link{prepare}}
  function)
}
\examples{
## Load the datasets into R session:
data(ProsCan_DE)
DE_data <- ProsCan_DE
data(ProsCan_Length)
Length_data <- ProsCan_Length

## Prepare a data frame of class "prepGOglm":
gene_table <- prepare(DE_data, Length_data, trans.p = "d.log", trans.l = TRUE)

## Call the generic plot() function:
plot(gene_table)
}
\author{
  Gu Mi \email{mig@stat.oregonstate.edu}, Yanming Di
  \email{diy@stat.oregonstate.edu}
}
\seealso{
  \code{\link{plot}}
}

