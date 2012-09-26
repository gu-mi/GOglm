\name{prepare}
\alias{prepare}
\title{Prepare a Data Frame for using GOglm.}
\usage{
  prepare(DE.data, Length.data,
    trans.p = c("n.log", "d.log"), trans.l = FALSE)
}
\arguments{
  \item{DE.data}{A data frame with valid gene identifiers
  as row names and one column for the untransformed DE test
  \emph{p}-values.}

  \item{Length.data}{A data frame with valid gene
  identifiers as row names and one column for the
  untransformed gene lengths (in bp). Note that the gene
  identifiers must be in the same order as the row names of
  \code{DE.data}, or an error message will appear.}

  \item{trans.p}{How to transform DE test \emph{p}-values
  to get significance statistics. Users can use "n.log" for
  \code{-log(p)}, or use "d.log" for \code{log(1-log(p))}.
  A panel of histograms can be obtained to evaluate the
  ranges of the transformed quantities. See Examples
  below.}

  \item{trans.l}{A logical value indicating whether to
  transform gene lengths. Default is \code{TRUE}.}
}
\value{
  An object of class \code{prepGOglm} with gene identifiers
  as row names and two columns. The first column contains
  (transformed) DE test \emph{p}-values (significance
  statistics), and the second column contains (transformed)
  gene lengths. End-users have the choice of applying
  different transformations to either DE test
  \emph{p}-values or gene lengths, or both.
}
\description{
  This function incorporates DE test results and gene
  length information provided by the end-users. It
  generates a data frame ready to be used by other
  functions in the \code{GOglm} package.
}
\details{
  This function takes two data frames: one is the output
  from a DE test, and the other contains gene length
  information. The end-users should confirm that these two
  data frames have the same gene identifiers as row names.
  Otherwise an error message will show up.
}
\examples{
## Load the datasets into R session:
data(ProsCan_DE)
DE.data <- ProsCan_DE
data(ProsCan_Length)
Length.data <- ProsCan_Length

## Prepare a data frame to be passed to goglm():
gene_table <- prepare(DE.data, Length.data, trans.p = "d.log", trans.l = TRUE)

## Check first 10 rows of the data frame:
gene_table[1:10,1:2]

## We can call the summary() function:
summary(gene_table)

## We can also call the plot() function:
plot(gene_table)
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

