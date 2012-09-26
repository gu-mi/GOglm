\name{unfactor}
\alias{unfactor}
\title{Purge factors}
\usage{
  unfactor(var)
}
\arguments{
  \item{var}{The variable from which you want the factors
  removed}
}
\value{
  The variable with all factors converted to characters or
  numbers
}
\description{
  Removes all factors from a variable in a sensible way
}
\details{
  As factors are their own type, to remove factors we must
  convert each level into another type. This is currently
  done using "typeless" behaviour: a factor is converted to
  a numeric vector if this can be done without inducing
  NAs, otherwise it is coerced using as.character.
  Currently supported types are: /codefactor,
  /codedata.frame and /codelist.
}
\note{
  This function is written by Matthew D. Young and is used
  in the \code{goseq} package.
}
\author{
  Matthew D. Young \email{myoung@wehi.edu.au}, Gu Mi
  \email{mig@stat.oregonstate.edu}, Yanming Di
  \email{diy@stat.oregonstate.edu}
}

