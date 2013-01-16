################################################################################
## Copyright (C) 2013 Gu Mi <mig@stat.oregonstate.edu>
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA
################################################################################

# Package Documentation
# 
# Author: Gu Mi
################################################################################

##' Length Bias Correction in Gene Ontology Enrichment Analysis Using
##' Logistic Regression
##'
##' This GOglm package is a beta version under development. The goglm
##' function implements the GOglm method discussed in Mi et al. (PLOS
##' ONE, 7(10): e46128).  The package includes a summarized RNA-Seq data
##' example (the prostate cancer dataset) for methodological
##' illustrations.
##'
##'
##' @name GOglm-package
##' @aliases GOglm-package GOglm
##' @docType package
##' @author Gu Mi <mig@@stat.oregonstate.edu>, Yanming Di
##' <diy@@stat.oregonstate.edu>.
##' Maintainer: Gu Mi <http://people.oregonstate.edu/~mig>
##' 
##' @references Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012)
##' "Length bias correction in Gene Ontology enrichment analysis using
##' logistic regression", 7(10): e46128.
##' @keywords package
##' 
NULL

.onLoad <- function(libname, pkgname){
  message <- "The GOglm package is experimental and under development: 
  Function syntax may have minor changes in future versions.
  See https://github.com/gu-mi/GOglm for more information. \n"
  packageStartupMessage(message)
}  
