\name{plot}
\alias{plot}
\title{Plot Data Preparation (DE Test) Results}
\usage{
  plot(x, ...)
}
\arguments{
  \item{x}{An object of class \code{prepGOglm}. Other
  plotting features for new classes will be added as
  package develops}

  \item{...}{Other parameters (for future use)}
}
\value{
  Plots depending on the class of the first argument. For
  an object of class \code{prepGOglm}, it returns a panel
  of two histograms for significance statistics and
  (transformed) gene lengths, respectively.
}
\description{
  A generic function for summerizing data preparation
  results
}
\examples{
## Load the datasets into R session:
data(ProsCan_DE)
DE.data <- ProsCan_DE
data(ProsCan_Length)
Length.data <- ProsCan_Length

## Prepare a data frame of class "prepGOglm":
gene_table <- prepare(DE.data, Length.data, trans.p = "d.log", trans.l = TRUE)

## Call the generic plot() function:
plot(gene_table)
}
\author{
  Gu Mi \email{mig@stat.oregonstate.edu}, Yanming Di
  \email{diy@stat.oregonstate.edu}
}
\seealso{
  \code{\link{plot.prepGOglm}}
}

