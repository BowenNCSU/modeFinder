
<!-- README.md is generated from README.Rmd. Please edit that file -->

# modeFinder

<!-- badges: start -->

<!-- badges: end -->

The goal of modeFinder is to compute the univariate mode estiamtes via
simple empirical density estimates, kernel density estimates or
Bernstein polynomials.

## Mode estimation via Bernstein polynomials

The R function of mode estimation via Bernstein polynomials is based on
the draft of:

Liu, Bowen, and Sujit K. Ghosh. *“On empirical estimation of mode based
on weakly dependent samples.”* Computational Statistics & Data Analysis
152 (2020): 107046.

## Simple empirical mode estimation

It is natural to estimate the population distribution function by the
empirical distribution function, \(F_n\) given by

\[
F_n = \frac{1}{n} \sum^n_{i = 1} I \{X_i \leq x\}
\]

and for any \(X^*\) such that \(X^{(i)} < X^* < X^{(i + 1)}\),
\(i = 1, 2, \dots, n\), the slope of the empirical distribution function
can be approximated by

\[
\frac{ F_n (X^{(i + 1)}) - F_n (X^{(i)})}{X^{(i + 1)} - X^{(i)}}.
\]

Consequently, we can give a simple mode estimator, \(u^{(\hat{k})}\)
based on the empirical distribution function by

\[
\hat{k} = \underset{0 \leq i \leq m^* - 1}{argmax} \frac{ F_n (u^{(i + 1)}) - F_n (u^{(i)})}{u^{(i + 1)} - u^{(i)}}.
\]

where
\(\{u^{(1)}, u^{(2)}, \dots, u^{(m^*)}\} = \{ \frac{1}{m^*}, \frac{2}{m^*}, \dots, \frac{m^*}{m^*}\}\).
However, the simple empirical mode estimation does not collaborate with
the possible smoothness characteristics of the underlying distribution
function.

## Installation

<!-- You can install the released version of modeFinder from [CRAN](https://CRAN.R-project.org) with: -->

<!-- And -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("BowenNCSU/modeFinder")
```

## Example

This is a basic example which shows you how to solve a mode estimation
problem:

``` r
library(modeFinder)
#> Loading required package: Rcpp
n = 1e6
set.seed(1)
data = rbeta(n, 2, 5)
## Return Bernstein polynomials based mode estimate
emp_mode(data)
#> 
#> maximum point: 0.200715, probability around maximum: 0.038835, probability of left boundary: 0.00354, probability of right boundary: 0
#> The mode estimate is chosen as the maximum point of Bernstein polynomials density estimate.
#> $mode
#> [1] 0.2007149
#> 
#> $optimal_alpha
#> [1] 0.35
#> 
#> $`optimal_p-value`
#> [1] 0.9963385
## Return kenerl density mode estimate
emp_mode(data, smooth_option = "kernel")
#> $mode
#> [1] 1.952236e-312
## Return simple empirical mode estimate
emp_mode(data, smooth = FALSE)
#> The mode estimate is choosen as the mid point of interval with maximum difference among knots equally-spaced intervals.
#> $`number of knots`
#> [1] 1000000
#> 
#> $mode
#> [1] 0.2081525
```

### Note:

macOS requires Xcode developer tools otherwise *xcrun: error* flags. Run
the following: *xcode-select –install*.
