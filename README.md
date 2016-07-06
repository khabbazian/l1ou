
### Tools for detecting past changes in the expected mean trait values and studying trait evolution from comparative data
The l1ou package provides functions to study trait evolution from comparative data and detect past changes in the expected mean trait values, as well as convergent evolution. It uses the Ornstein-Uhlenbeck process along a phylogenetic tree, which can model a changing adaptive landscape over time and over lineages. 
<!--Detection of evolutionary shifts in trait evolution from extant taxa is motivated by the study of convergent evolution, or to correlate shifts in traits with habitat changes or with changes in other phenotypes.-->
Shifts can be detected from multiple traits, assuming that all traits shifted along the same lineages. Estimation is very fast thanks to lasso techniques, and the user can choose from multiple information criteria for model selection, including a phylognetic BIC. 
Citation: 

- M. Khabbazian, R. Kriebel, K. Rohe, and C&eacute;cile An&eacute;. "Fast and accurate detection of evolutionary shifts in Ornstein-Uhlenbeck models". Methods in Ecology and Evolution. doi:10.1111/2041-210X.12534  

### Version notes 
  Starting with version v1.22, the scores returned by "estimate\_shift\_configuration” function 
  are for the non-normalized, original data.  


### Install using the devtools package.
```
install.packages("devtools")
require(devtools)
install_github("khabbazian/l1ou")
require(l1ou)
```
Windows users will first need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/).

#### [l1ou Reference manual](http://homepages.cae.wisc.edu/~khabbazian/pdfs/l1ou.pdf)
