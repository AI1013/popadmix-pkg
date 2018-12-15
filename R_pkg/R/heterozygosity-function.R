#' @title Calculate the heterozygosity.

#' @description \code{get.het} returns the heterozygosity of a population at a locus.

#' @details The input is the vector of allele ferquencies (should add up to 1).
#'   If there are missing values in the inputs or the sum of the inputs don't
#'   equal to 1, then cannot use this function.

#' @param allele.freq Numeric vectors.
#' @return If all inputs are positive and logical, then the output
#'   will be an numeric vector. If there are missing values in the
#'   inputs, the output will be NA with a warning.
#' @export get.het
#' @keywords heterozygosity
#' @references Boca, S. M., & Rosenberg, N. A. (2011). Mathematical properties of Fst between admixed populations 
#' and their parental source populations. Theoretical population biology, 80(3), 208-216.
#' @examples
#' c <- c(0.2, 0.3, 0.5)
#' get.het(c)

get.het <- function(allele.freq) {
  
  if (any(is.na(allele.freq) == TRUE)) warning("There should not be NA in allele.freq")
  1 - sum(allele.freq^2)
}

#' @title Calculate the heterozygosity of the admixed population
#' @description \code{get.Hadm.admix.general} returns the heterozygosity of the admixed population.
#' @details The input is a matrix with K rows and I columns (sums on rows is 1), vector of gammas (add up to 1).
#' @param pops A matrix that contains the allele frequencies for each population, for each row, the sum should be 1.
#' @param gamma Numeric vectors indictating the proportion of each subpopulation in the mixed population, and the sum should be 1.
#' @return If there are missing values in the inputs, the output will be NA with a warning. If the sums of rows in martix is not
#'   1 or the sum of gamma is not 1, it will stop.
#' @export get.Hadm.admix.general
#' @examples
#' pops <- matrix(nrow = 2, ncol = 10, c(0.059, 0.153, 0.132, 0.101, 0.08, 0.055, 0.106, 0.145, 0.048, 
#'                0.121, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000), byrow=TRUE)
#' gamma <- c(0.3, 0.7)
#' get.Hadm.admix.general(pops, gamma)

get.Hadm.admix.general <- function(pops, gamma){
  if (!(rowSums(pops) == 1)) stop("The sum of the allele frequency for each population should add up to 1")
  if (!(sum(gamma) == 1)) stop("The sum of gamma should add up to 1")
  if (any(is.na(pops) == TRUE)) warning("There should not be NA in pops")
  if (any(is.na(gamma) == TRUE)) warning("There should not be NA in gamma")
  ##get admixed population
  admix.pop <- t(pops) %*% gamma
  ##calculate the heterozygosity of the admixed population
  Hadm <- get.het(admix.pop)
  Hadm
}
