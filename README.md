
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CRANE v1.0.0

Constrained Random Alteration of Network Edges (CRANE), a computational
method for sampling networks with fixed node strengths,

## Installation

You can install the released version of CRANE from
[GitHub](https://github.com/PadiLab/CRANE) with:

``` r
install_github("PadiLab/CRANE")
```

## Publication

Coming soon…

## Example

We have included network edge list from angiogenic ovarian cancer (ang)
and non-angiogenic ovarian cancer (nonAng). The following are some
example use:

``` r
library(CRANE)
#> 
#> 
head(ang)
#>       TF Gene    weight
#> 1    AHR A1CF -1.090067
#> 2     AR A1CF -1.296286
#> 3 ARID3A A1CF  3.347115
#> 4   ARNT A1CF  2.394743
#> 5  BRCA1 A1CF  3.084647
#> 6  CREB1 A1CF -1.384014

# Running CRANE for bipartite network
newElist=crane.bipartite(ang,alpha=0.3)
#> [1] "Converting Edgelist to Adj Matrix"
#> [1] "NORMAL PERTURBATION"
#> [1] "Applying Alpha = 0.3"
#> [1] "Constructing Iteratuvely Perturbed Network"
#> [1] 100
#> [1] "Sorting Nodes"

head(newElist)
#>     from   to      weight
#> 1    AHR A1CF -1.46548313
#> 2     AR A1CF -1.46621086
#> 3 ARID3A A1CF  3.30632912
#> 4   ARNT A1CF  2.43699996
#> 5  BRCA1 A1CF  3.57989571
#> 6  CREB1 A1CF -0.06841034
```
