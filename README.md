# phyloLSnoDist
This package performs phylogenetic inference under the least squares framework but without the use of traditional evolutionary distances, as described in [``Phylogenetic least squares estimation without genetic distances''](https://arxiv.org/abs/2311.12717).

## Preliminaries
First, two packages need to be installed: 

  - `robustDist` 
  - `phyloLSnoDist`
  

Installation of `robustDist` can be accomplished with the following:
```
install.packages("robustDist", repos=c("http://R-Forge.R-project.org",
                                       "http://cran.at.r-project.org"),dependencies=TRUE))
```
You may also need to include the option `INSTALL_opts = c('--no-lock').`

Then, the easiest way to install `phyloLSnoDist` will be to use `install_github` from the `devtools` package:

```
install.packages("devtools") # if not already installed
library(devtools)
install_github("peterbchi/phyloLSnoDist")
```

