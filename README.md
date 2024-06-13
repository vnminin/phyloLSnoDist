# phyloLSnoDist
This package performs phylogenetic inference under the least squares framework but without the use of traditional evolutionary distances, as described in [``Phylogenetic least squares estimation without genetic distances''](https://arxiv.org/abs/2311.12717).

## Installation
First, two packages need to be installed: 

  - `robustDist` 
  - `phyloLSnoDist`
  

Installation of `robustDist` can be accomplished with the following:
```
install.packages("robustDist", repos=c("http://R-Forge.R-project.org",
                                       "http://cran.at.r-project.org"),dependencies=TRUE)
```
You may also need to include the option `INSTALL_opts = c('--no-lock').`

Then, the easiest way to install `phyloLSnoDist` will be to use `install_github` from the `devtools` package:

```
install.packages("remotes")
remotes::install_github("vnminin/phyloLSnoDist")
```


## Simple Example
First we must load the `phyloLSnoDist` package, which should also load its dependencies:
```
library(phyloLSnoDist)
```

Using the `rtree` function from the `ape` package, we generate a 4-taxon tree which we will use for demonstration:
```
set.seed(100)
my_tree <- rtree(4)
plot(my_tree, type='u')
```

With this tree, let us generate some DNA nucleotide sequence data, and convert it to `phyDat` format for usage by the `phylo.ls.nodist` function:
```
my_DNA <- simSeq(my_tree, 2000)
my_DNA_pD <- as.phyDat(my_DNA)
```

Then, we run our `phylo.ls.nodist` function to run phylogenetic inference with our new loss function. Since there are only three unrooted topologies, we do an exhaustive search by setting `search.all = TRUE`:
```
nodist_tree <- phylo.ls.nodist(my_DNA_pD, search.all = TRUE)
plot(nodist_tree, type='u')
```


## Extended Vignette
A longer vignette, showing the above code with its intended output, along with replication of select simulation experiments from the manuscript, is available in the `vignette` directory. The output pdf is `vignette.pdf` along with associated R Markdown file, `vignette.Rmd`.

## Other Notes
The `analysis` directory of this repository contains scratchwork, older versions of code, and unpublished simulation experiments. All final code for the manuscript is in the `manuscript` directory. 
