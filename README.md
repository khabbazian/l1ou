
### Tools for detecting past changes in the expected mean trait values and studying trait evolution from comparative data

The l1ou package provides functions to study trait evolution from comparative data and detect past changes in the expected mean trait values, as well as convergent evolution. It uses the Ornstein-Uhlenbeck process along a phylogenetic tree, which can model a changing adaptive landscape over time and over lineages. 
<!--Detection of evolutionary shifts in trait evolution from extant taxa is motivated by the study of convergent evolution, or to correlate shifts in traits with habitat changes or with changes in other phenotypes.-->
Shifts can be detected from multiple traits, assuming that all traits shifted along the same lineages. Estimation is very fast thanks to lasso techniques, and the user can choose from multiple information criteria for model selection, including a phylognetic BIC. 
Citation: 

- M. Khabbazian, R. Kriebel, K. Rohe, and C&eacute;cile An&eacute;.
Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models.
Methods in Ecology and Evolution, 7(7):811–824.
[doi:10.1111/2041-210X.12534](http://dx.doi.org/10.1111/2041-210X.12534)

#### [l1ou Reference manual](https://github.com/khabbazian/l1ou/blob/master/l1ou.pdf)

### Install using the devtools package

from within R:
```r
install.packages("devtools")
library(devtools)
install_github("khabbazian/l1ou")
```
Windows users will first need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

### Install without the devtools package

To resolve dependencies, first install the following packages from CRAN, then the knitr package.
From within R:
```r
install.packages(c("igraph", "phylolm", "lars", "grplasso", "magic", "genlasso", "Rcpp"))
install.packages("knitr")
```
Now in the shell, with asterisks to be replaced with the correct version number:
```shell
git clone https://github.com/khabbazian/l1ou.git 
R CMD build l1ou 
R -e 'install.packages("l1ou_*.**.tar.gz")'
```

### Version notes 

major changes are indicated below.

- v1.40:
  * intercept correctly handled after noise-whitening
  (results may change for variables with a mean far from 0)
  * bug fix in the function calculating the square-root (and inverse) of the
  phylogenetic covariance matrix.
- v1.25: the penalty term in the AICc score now considers the intercept as a free variable.
  The change only affects the final value of the AICc score.
- v1.23: "estimate\_convergent\_regimes" function accepts multiple traits. 
- v1.22: 
	* the scores returned by "estimate\_shift\_configuration” function are now for the non-normalized, original data.
	* "estimate\_shift\_configuration" function also accepts multiple traits with missing values. 

