##' Removes all factors from a variable in a sensible way
##'
##' As factors are their own type, to remove factors we must convert
##' each level into another type. This is currently done using
##' "typeless" behaviour: a factor is converted to a numeric vector if
##' this can be done without inducing NAs, otherwise it is coerced
##' using as.character. Currently supported types are: /codefactor,
##' /codedata.frame and /codelist.
##'
##' @title Purge factors
##'
##' @param var The variable from which you want the factors removed
##'
##' @return The variable with all factors converted to characters or
##' numbers
##'
##' @note This function is written by Matthew D. Young and is used in
##' the \code{goseq} package.
##'
##' @author Matthew D. Young \email{myoung@@wehi.edu.au}, Gu Mi
##' \email{mig@@stat.oregonstate.edu}, Yanming Di
##' \email{diy@@stat.oregonstate.edu}
##'
##' @export
##'

unfactor <- function (var)
{
    if (is.factor(var)) {
        tmp = names(var)
        tmpopt = getOption("warn")
        options(warn = -1)
        out = as.numeric(levels(var))[as.integer(var)]
        options(warn = tmpopt)
        if (any(is.na(out)) & any(is.na(out) != is.na(var))) {
            out = as.character(levels(var))[as.integer(var)]
        }
        names(out) = tmp
    }
    else if (is.data.frame(var)) {
        out = var
        for (i in 1:dim(var)[2]) {
            out[, i] = unfactor(var[, i])
        }
    }
    else if (is.list(var)) {
        out = lapply(var, unfactor)
    }
    else {
        out = var
    }
    return(out)
}
