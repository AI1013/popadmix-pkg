#' @title Create the allele frequency data we use in this package.
#' @description \code{get_allele_freq} returns the allele frequency data.
#' @export get_allele_freq
#' @references Wang, S., Ray, N., Rojas, W., Parra, M. V., Bedoya, G., Gallo, C., ... & Camrena, B. (2008). 
#' Geographic patterns of genome admixture in Latin American Mestizos. PLoS genetics, 4(3), e1000037.
#' @examples
#' get_allele_freq()

get_allele_freq <- function(){
  dataset <- read.table(system.file("extdata", "Genotypes_Wang_etal_Mestizo.txt", package = "popadmix"), header = FALSE, 
                        sep = "", stringsAsFactors = TRUE, fill = TRUE)
  cn <- dataset[1,]
  cn <- cn[!is.na(cn)]
  dataset <- dataset[-2,]
  dataset <- dataset[2:nrow(dataset), 6:ncol(dataset)]
  names(dataset)[2:ncol(dataset)] <- cn
  dataset[dataset == -9] <- NA
  
  pop <- unique(dataset[,1])
  pop <- as.character(pop)
  rn <- pop
  rn[2] <- "AMERICA"
  allele_frequencies <- list()
  
  for (n in 2:ncol(dataset)) {
    loci <- unique(dataset[n])
    loci <- loci[!is.na(loci)]
    allele_matrix <- NULL
    
    for (popname in pop) {
      subpop <- subset(dataset[n], dataset[1] == popname)
      subpop <- subpop[!is.na(subpop)]
      sum_number <- length(subpop)
      H <- NULL
      for (ele in loci) {
        number <- length(unlist(subset(subpop, subpop == ele)))
        freq <- number / sum_number
        H <- c(H, freq)
      }
      allele_matrix <- rbind(allele_matrix, H)
    }
    rownames(allele_matrix) <- rn
    colnames(allele_matrix) <- loci
    allele_matrix<- as.matrix(allele_matrix)
    allele_frequencies[[n-1]] <- allele_matrix
  }
  
  names(allele_frequencies) <- cn
  allele_frequencies
  
}