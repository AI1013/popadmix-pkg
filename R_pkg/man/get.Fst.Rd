% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Fst-function.R
\name{get.Fst}
\alias{get.Fst}
\title{Calculate Fst.}
\usage{
get.Fst(pop.1, pop.2)
}
\arguments{
\item{pop.1}{Numeric vector representing allele frequencies in the same locus for populations 1.}

\item{pop.2}{Numeric vector representing allele frequencies in the same locus for populations 2.}
}
\value{
If all inputs are positive and logical, then the output
  will be an numeric vector. If there are missing values in the
  inputs, the output will be NA with a warning.
}
\description{
\code{get.Fst} returns the Fst between two populations.
}
\details{
library(mspath)

The input pop.1 and pop.2 represent allele frequencies in the same locus for populations 1
and 2 respectively. Ht is the expected heterozygosity of the overall population, and Hs is
the mean expected heterozygosity across subpopulations.
}
\examples{
pop.1 <- c(0.3,0.7)
pop.2 <- c(0.4,0.6)
get.Fst(pop.1, pop.2)
}
\references{
Boca, S. M., & Rosenberg, N. A. (2011). Mathematical properties of Fst between admixed populations 
and their parental source populations. Theoretical population biology, 80(3), 208-216.
}
\keyword{Fst}
