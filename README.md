
<!-- README.md is generated from README.Rmd. Please edit that file -->

# chantrics

<!-- badges: start -->

![R-CMD-check](https://github.com/tbruckbauer/chantrics/workflows/R-CMD-check/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/tbruckbauer/chantrics/branch/master/graph/badge.svg?token=EQQ177WODM)](https://codecov.io/gh/tbruckbauer/chantrics)
![pkgdown](https://github.com/tbruckbauer/chantrics/workflows/pkgdown/badge.svg)
[![CodeFactor](https://www.codefactor.io/repository/github/tbruckbauer/chantrics/badge?s=8292aacb6a947f6280b202c4f29e21e4510dce53)](https://www.codefactor.io/repository/github/tbruckbauer/chantrics)
<!-- badges: end -->

`chantrics` applies the Chandler-Bate loglikelihood adjustment (Chandler
and Bate 2007) implemented in the
[chandwich](https://cran.r-project.org/package=chandwich) package to
different models frequently used in basic Econometrics applications.
`adj_loglik()` adjusts the parameter covariance matrix of the models to
incorporate clustered data, and can mitigate for model misspecification
by wrapping `chandwich::adjust_loglik` for the supported models.

The returned model of class `chantrics` can be plugged into standard
model evaluation and model comparison methods, for example `summary()`,
`confint()` and `anova`, and a hypothesis test framework provided by
`alrtest()`.

## Installation

<!--You can install the released version of chantrics from [CRAN](https://CRAN.R-project.org) with:


```r
# install.packages("chantrics")
```
-->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
## Remove the # if "devtools" is not installed yet.
# install.packages("devtools")
devtools::install_github("tbruckbauer/chantrics")
```

## Usage

This example is using the misspecified count data regression from
Chapter 5.1 in the Object-Oriented Computation of Sandwich Estimators
vignette from the *sandwich* package (Zeileis 2006). First, data from a
negative binomial model is generated, and then a Poisson model is fit,
which is clearly misspecified.

``` r
library(chantrics)

set.seed(123)
x <- rnorm(250)
y <- rnbinom(250, mu = exp(1 + x), size = 1)

## Fit the Poisson glm model, which is not correctly specified
fm_pois <- glm(y~x +I(x^2), family = poisson)
lmtest::coeftest(fm_pois)
#> 
#> z test of coefficients:
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.063268   0.041357 25.7094  < 2e-16 ***
#> x            0.996072   0.053534 18.6062  < 2e-16 ***
#> I(x^2)      -0.049124   0.023146 -2.1223  0.03381 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# The I(x^2) term is spuriously significant.
```

In order to mitigate the misspecification, the loglikelihood adjustment
is applied to the model object using `adj_loglik(model_object)`. If
clustered data is available, a vector or factor indicating from which
cluster the observation originates can be passed into the function using
`cluster`. If it is not supplied, it is assumed that each observation
originates from its own cluster and is independent. More information on
this can be found in the [clustering
vignette](https://chantrics.theobruckbauer.eu/articles/clustering-vignette.html).

``` r
## Apply the loglikelihood adjustment to the model
fm_pois_adj <- adj_loglik(fm_pois)
lmtest::coeftest(fm_pois_adj)
#> 
#> z test of coefficients:
#> 
#>              Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  1.063268   0.083776 12.6918   <2e-16 ***
#> x            0.996072   0.105217  9.4668   <2e-16 ***
#> I(x^2)      -0.049124   0.036284 -1.3539   0.1758    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# The I(x^2) term is no longer significant.
```

After this, some of the functionality is demonstrated. More information
on the different methods can be found in the [Introducing chantrics
vignette](https://chantrics.theobruckbauer.eu/articles/chantrics-vignette.html).

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-chanbate07" class="csl-entry">

Chandler, Richard, and Steven Bate. 2007. “Inference for Clustered Data
Using the Independence Loglikelihood.” *Biometrika* 94 (1): 167–83.
<https://doi.org/doi:10.1093/biomet/asm015>.

</div>

<div id="ref-zeileis06" class="csl-entry">

Zeileis, Achim. 2006. “Object-Oriented Computation of Sandwich
Estimators.” *Journal of Statistical Software, Articles* 16 (9): 1–16.
<https://doi.org/10.18637/jss.v016.i09>.

</div>

</div>
