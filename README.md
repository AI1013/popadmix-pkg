# popadmix-pkg
A package to calculate heterozygosity, Fst and make plots about them.

## Overview
This package can calculate the *Fst* and heterozygosity value in different cases, and make the triangle plots about them. The main functions are `get.Fst.admix.general`, `get.Hadm.admix.general` and `plot_HorF`; please see the corresponding help files for details, as well as the companion [paper](https://www.sciencedirect.com/science/article/pii/S0040580911000463):

Boca, S. M., & Rosenberg, N. A. (2011). Mathematical properties of Fst between admixed populations and their parental source populations. Theoretical population biology, 80(3), 208-216.

## Installation
To install the package in R, first install the devtools package, and then use the commands

```
library(devtools)
install_github('AI1013/popadmix-pkg/R_pkg', build_vignettes = TRUE)
```
