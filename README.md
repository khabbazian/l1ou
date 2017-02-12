
### Tools for detecting past changes in the expected mean trait values and studying trait evolution from comparative data
The l1ou package provides functions to study trait evolution from comparative data and detect past changes in the expected mean trait values, as well as convergent evolution. It uses the Ornstein-Uhlenbeck process along a phylogenetic tree, which can model a changing adaptive landscape over time and over lineages. 
<!--Detection of evolutionary shifts in trait evolution from extant taxa is motivated by the study of convergent evolution, or to correlate shifts in traits with habitat changes or with changes in other phenotypes.-->
Shifts can be detected from multiple traits, assuming that all traits shifted along the same lineages. Estimation is very fast thanks to lasso techniques, and the user can choose from multiple information criteria for model selection, including a phylognetic BIC. 
Citation: 

- M. Khabbazian, R. Kriebel, K. Rohe, and C&eacute;cile An&eacute;. "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models". Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534  

#### [l1ou Reference manual](http://www.columbia.edu/~mhk2154/pdfs/l1ou.pdf)

### Version notes 
  Starting with version v1.22, the scores returned by "estimate\_shift\_configuration” function 
  are for the non-normalized, original data.  

  Starting with version v1.25, the penalty term in AICc score considers the intercept as a free variable. The change only affects the final value of the AICc score.


### Install using the devtools package.
```
R> install.packages("devtools")
R> install_github("khabbazian/l1ou")
```
Windows users will first need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

### Install without the devtools package.
To resolve the dependencies, first install the following packages from CRAN
```
R> install.packages(c("igraph", "phylolm", "lars", "grplasso", "magic", "genlasso", "Rcpp"))
```
Install the knitr package
```
R> install.packages("knitr")
```
Now in bash
```
$> git clone https://github.com/khabbazian/l1ou.git 
$> R CMD build l1ou 
$> R -e 'install.packages("l1ou_*.**.tar.gz")'
```
Replace the asterisks with the correct version number.


