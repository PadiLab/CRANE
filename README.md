
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CRANE v1.0.0

Constrained Random Alteration of Network Edges (CRANE) is a
computational method for sampling networks with fixed node strengths.
This package also include CRANE integraction with ALPACA.

## Installation

You can install the released version of CRANE from
[GitHub](https://github.com/PadiLab/CRANE) with:

``` r
library(devtools)
install_github("PadiLab/CRANE")
```

## Publication

Coming soon…

## Examples

We have included network edge list from angiogenic ovarian cancer (ang)
and non-angiogenic ovarian cancer (nonAng). The following are some
example use:

1.  CRANE for network perturbation

<!-- end list -->

``` r
library(CRANE)

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
#> [1] "Applying Alpha = 0.3"
#> [1] "Constructing Iteratuvely Perturbed Network"
#> [1] 100
#> [1] "Sorting Nodes"

head(newElist)
#>     from   to     weight
#> 1    AHR A1CF -0.4863448
#> 2     AR A1CF -2.3563477
#> 3 ARID3A A1CF  2.9850825
#> 4   ARNT A1CF  2.3343717
#> 5  BRCA1 A1CF  2.9642476
#> 6  CREB1 A1CF -1.9202671
```

2.  CRANE integration with ALPACA

<!-- end list -->

``` r
library(ALPACA)
library(CRANE)

# Create Input Matrix nonAng as the baseline.
input=cbind(nonAng,ang[,3])

# Run ALPACA
alp=alpaca(input,NULL,verbose = F)
#> [1] "Detecting communities in control network..."
#> Weights detected. Building condor object with weighted edges.
#> [1] "modularity of projected graph 0.130474098565802"
#> [1] "Q = 0.138287700459528"
#> [1] "Q = 0.138287700459528"
#> [1] "Computing differential modularity matrix..."
#> [1] "Computing differential modules..."
#> [1] "Merging 12401 communities"
#> [1] 1
#> [1] 2
#> [1] 3
#> [1] 4
#> [1] 5
#> [1] 6
#> [1] "Merging 73 communities"
#> [1] 1
#> [1] 2
#> [1] 3
#> [1] 4
#> [1] "Merging 16 communities"
#> [1] 1
#> [1] "Computing node scores..."

# Apply Crane
alpListObject=alpaca.crane(input, alp, isParallel = T)
#> [1] "Converting to Adjecency Matrix..."
#> [1] "Generating CRANE Networks, this may take a long time (upto 10-20 minutes for Networks with 19000+ nodes) ..."
#> [1] "Creating Test Distribution..."
#> [1] "Detecting communities in control network..."
#> Weights detected. Building condor object with weighted edges.
#> [1] "modularity of projected graph 0.130474098565802"
#> [1] "Q = 0.138287700459528"
#> [1] "Q = 0.138287700459528"
#> [1] "Comuting Differential Scores..."
#> [1] "Computing Differential Score 20 %"
#> [1] "Computing Differential Score 40 %"
#> [1] "Computing Differential Score 60 %"
#> [1] "Computing Differential Score 80 %"
#> [1] "Computing Differential Score 100 %"

# TF Results
head(alpListObject$TF)
#>        Membership Differential Modularity Score       Pvalue t-statistic
#> ARID3A          1                    0.11078839 2.354490e-42  -134.03163
#> ARNT            1                    0.08121685 1.000000e+00   136.38618
#> BRCA1           2                    0.06560763 1.000000e+00    25.99591
#> ELF5            1                    0.03373797 1.000000e+00    82.42068
#> ETS1            2                    0.05777723 3.955381e-31   -54.75762
#> FEV             2                    0.05693467 1.000000e+00    83.51181

# Gene Results
head(alpListObject$Gene)
#>        Membership Differential Modularity Score       Pvalue t-statistic
#> A1CF            4                  0.0014013382 9.999995e-01    6.175483
#> A2M             1                  0.0007451022 1.127484e-05   -5.042341
#> A4GALT          6                  0.0055156253 1.982201e-13  -12.413065
#> A4GNT           2                  0.0003932609 1.974168e-03   -3.131707
#> AAAS            1                  0.0003320119 1.000000e+00    9.451461
#> AACS            9                  0.0027239872 1.545997e-02   -2.268509
```
