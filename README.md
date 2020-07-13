
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

**Generating Ensembles of Gene Regulatory Networks to Assess Robustness
of Disease Modules**  
James T Lim, Chen Chen, Adam D Grant, Megha Padi  
bioRxiv 2020.07.12.198747; doi:
<https://doi.org/10.1101/2020.07.12.198747>  
url: <https://www.biorxiv.org/content/10.1101/2020.07.12.198747v1>

## Examples

We have included network edge list from angiogenic ovarian cancer (ang)
and non-angiogenic ovarian cancer (nonAng). The following are some
example use:

1.  CRANE for network perturbation

<!-- end list -->

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
#> [1] "Applying Alpha = 0.3"
#> [1] "Constructing Iteratuvely Perturbed Network"
#> [1] 100
#> [1] "Sorting Nodes"

head(newElist)
#>     from   to     weight
#> 1    AHR A1CF -0.7996633
#> 2     AR A1CF -0.5410527
#> 3 ARID3A A1CF  2.9311514
#> 4   ARNT A1CF  2.3026895
#> 5  BRCA1 A1CF  2.8180388
#> 6  CREB1 A1CF -1.0686323
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
#> ARID3A          1                    0.11078839 7.100957e-43  -139.69642
#> ARNT            1                    0.08121685 1.000000e+00    98.68718
#> BRCA1           2                    0.06560763 1.000000e+00    23.98381
#> ELF5            1                    0.03373797 1.000000e+00    98.33036
#> ETS1            2                    0.05777723 3.636443e-30   -50.68560
#> FEV             2                    0.05693467 1.000000e+00    87.91174

# Gene Results
head(alpListObject$Gene)
#>        Membership Differential Modularity Score       Pvalue t-statistic
#> A1CF            4                  0.0014013382 9.999974e-01   5.5703422
#> A2M             1                  0.0007451022 6.637922e-06  -5.2329474
#> A4GALT          6                  0.0055156253 1.529687e-16 -16.4454533
#> A4GNT           2                  0.0003932609 4.630008e-02  -1.7392323
#> AAAS            1                  0.0003320119 1.000000e+00   9.1197478
#> AACS            9                  0.0027239872 5.605418e-01   0.1536929
```
