## ----echo=TRUE, message=FALSE, warning=FALSE, fig.cap=TRUE---------------
knitr::opts_chunk$set(echo = TRUE, comment = NA)

## ------------------------------------------------------------------------
library(knitr)
library(kableExtra)

## ------------------------------------------------------------------------
get.het <- function(allele.freq)
{
1 - sum(allele.freq^2)
}

## ------------------------------------------------------------------------
sum(0.2, 0.3, 0.5)

## ------------------------------------------------------------------------
is.na(c(0.2, 0.3, 0.5))

## ------------------------------------------------------------------------
get.het(c(0.2, 0.3, 0.5))

## ------------------------------------------------------------------------
get.Fst <- function(pop.1, pop.2)
{
  meta.pop <- (pop.1 + pop.2)/2
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pop.1), get.het(pop.2)))
  Fst <- (Ht - Hs)/Ht
  Fst
}

## ------------------------------------------------------------------------
is.na(c(0.3, 0.4))

## ------------------------------------------------------------------------
get.Fst(0.3, 0.4)

## ------------------------------------------------------------------------
get.Fst.admix <- function(pop.1, pop.2, gamma)
{
  pop.3 <- gamma*pop.1 + (1 - gamma)*pop.2
  meta.pop <- (pop.1 + pop.3)/2
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pop.1), get.het(pop.3)))
  Fst <- (Ht - Hs)/Ht
  Fst
}

## ------------------------------------------------------------------------
is.na(c(0.3, 0.4, 0.25))

## ------------------------------------------------------------------------
get.Fst.admix(0.3, 0.4, 0.25)

## ------------------------------------------------------------------------
get.Fst.admix.general <- function(pops, gamma)
{
  ##get admixed population
  admix.pop <- t(pops) %*% gamma
  
  ##get allele frequencies of "meta-population"
  meta.pop <- (pops[1, ] + admix.pop)/2
  
  Ht <- get.het(meta.pop)
  Hs <- mean(c(get.het(pops[1, ]), get.het(admix.pop)))
  Fst <- (Ht - Hs)/Ht
  Fst
}

## ----results='asis'------------------------------------------------------
pops <- matrix(nrow = 2, ncol = 10, c(0.059, 0.153, 0.132, 0.101, 0.08, 0.055, 0.106, 0.145, 0.048, 0.121, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000), byrow=TRUE)

## ----echo=TRUE-----------------------------------------------------------
class(pops)
rownames(pops) <- c("pop1", "pop2")
colnames(pops) <- paste("Allele",1:ncol(pops))
pops %>% 
    kable() %>%
    kable_styling()

## ----echo=TRUE-----------------------------------------------------------
rowSums(pops)

## ----echo=TRUE-----------------------------------------------------------
gamma <- c(0.3, 0.7)
sum(gamma)

## ------------------------------------------------------------------------
is.na(pops)
is.na(gamma)

## ----echo=TRUE-----------------------------------------------------------
get.Fst.admix.general(pops, gamma)

