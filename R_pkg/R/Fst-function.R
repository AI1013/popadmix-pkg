#' library(mspath)
#' @title Calculate Fst.

#' @description \code{get.Fst} returns the Fst between two populations.

#' @details The input pop.1 and pop.2 represent allele frequencies in the same locus for populations 1
#' and 2 respectively. Ht is the expected heterozygosity of the overall population, and Hs is
#' the mean expected heterozygosity across subpopulations.

#' @param pop.1 Numeric vector.
#' @param pop.2 Numeric vector.
#' @inheritParams get.het.
#' @return If all inputs are positive and logical, then the output
#'   will be an numeric vector. If there are missing values in the
#'   inputs, the output will be NA with a warning.
#' @export get.Fst
#' @keywords Fst
#' @references Boca, S. M., & Rosenberg, N. A. (2011). Mathematical properties of Fst between admixed populations 
#' and their parental source populations. Theoretical population biology, 80(3), 208-216.
#' @examples
#' pop.1 <- 0.3
#' pop.2 <- 0.4
#' get.Fst(pop.1, pop.2)

get.Fst <- function(pop.1, pop.2) {
  if (any(is.na(pop.1) == TRUE)) warning("There should not be NA in the pop.1")
  if (any(is.na(pop.2) == TRUE)) warning("There should not be NA in the pop.2")
  meta.pop <- (pop.1 + pop.2)/2
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pop.1), get.het(pop.2)))
  Fst <- (Ht - Hs)/Ht
  Fst
}

#' @title Calculate Fst: between a founder population and the admixed population.

#' @description \code{get.Fst.admix} returns the Fst between a founder population and the admixed population for the K=2 founder population case.

#' @details The input pop.1 and pop.2 represent allele frequencies in the same locus for populations 1
#' and 2 respectively, and gamma is the ancestry fraction corresponding to population 1. Ht is
#' the expected heterozygosity of the overall population, and Hs is the mean expected
#' heterozygosity across subpopulations.

#' @param pop.1 Numeric vector.
#' @param pop.2 Numeric vector.
#' @param gamma Numeric vector.
#' @inheritParams get.het.
#' @return If all inputs are positive and logical, then the output
#'   will be an numeric vector. If there are missing values in the
#'   inputs, the output will be NA with a warning.
#' @export get.Fst.admix
#' @keywords Fst
#' @references Boca, S. M., & Rosenberg, N. A. (2011). Mathematical properties of Fst between admixed populations 
#' and their parental source populations. Theoretical population biology, 80(3), 208-216.
#' @examples
#' pop.1 <- 0.3
#' pop.2 <- 0.4
#' gamma <- 0.25
#' get.Fst.admix(pop.1, pop.2, gamma)

get.Fst.admix <- function(pop.1, pop.2, gamma) {
  if (any(is.na(pop.1) == TRUE)) warning("There should not be NA in pop.1")
  if (any(is.na(pop.2) == TRUE)) warning("There should not be NA in pop.2")
  if (any(is.na(gamma) == TRUE)) warning("There should not be NA in gamma")
  pop.3 <- gamma * pop.1 + (1 - gamma) * pop.2
  meta.pop <- (pop.1 + pop.3)/2
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pop.1), get.het(pop.3)))
  Fst <- (Ht - Hs)/Ht
  Fst
}

#' @title Calculate Fst for pop.1 and admixed population.

#' @description \code{get.Fst.admix.general} returns Fst for pop.1 and admixed population for the general case of K founding populations.

#' @details The input pops is matrix with K rows and I columns (sums on rows is 1), and
#' gamma are the ancestry fractions corresponding to the K founding populations
#' (add up to 1). Ht is the expected heterozygosity of the overall population,
#' and Hs is the mean expected heterozygosity across subpopulations. If the sums
#' on rows in the matrix is not 1 or the sum of vectors of gamma is not 1, then
#' cannot use this function

#' @param pops Matirx.
#' @param gamma Numeric vectors.
#' @inheritParams get.het.
#' @return If all inputs are positive and logical, then the output
#'   will be an numeric vector.If there are missing values in the
#'   inputs, the output will be NA with a warning.
#' @export get.Fst.admix.general
#' @keywords Fst
#' @references Boca, S. M., & Rosenberg, N. A. (2011). Mathematical properties of Fst between admixed populations 
#' and their parental source populations. Theoretical population biology, 80(3), 208-216.
#' @examples
#' pops <- matrix(nrow = 2, ncol = 10, c(0.059, 0.153, 0.132, 0.101, 0.08, 0.055, 0.106, 0.145,
#' 0.048, 0.121, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000), byrow=TRUE)
#' gamma <- c(0.3, 0.7)
#' get.Fst.admix.general(pops, gamma)

get.Fst.admix.general <- function(pops, gamma) {
  if (!(rowSums(pops) == 1)) stop("The sum of the allele frequency for each population should add up to 1")
  if (!(sum(gamma) == 1)) stop("The sum of gamma should add up to 1")
  if (any(is.na(pops) == TRUE)) warning("There should not be NA in pops")
  if (any(is.na(gamma) == TRUE)) warning("There should not be NA in gamma")
  ## get admixed population
  admix.pop <- t(pops) %*% gamma

  ## get allele frequencies of 'meta-population'
  meta.pop <- (pops[1, ] + admix.pop)/2

  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pops[1, ]), get.het(admix.pop)))
  Fst <- (Ht - Hs)/Ht
  Fst
}


