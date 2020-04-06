
<!-- README.md is generated from README.Rmd. Please edit that file -->

# IntegrativeHTEcf

The goal of *IntegrativeHETcf* is to implement integrative analyses for
the heterogenous treatment effect and confounding functions combining a
randomized trial and a confounded real-world evidence study with
unmeasured confounders.

Two datasets

  - The randomized trial contains observations on (A,X,Y), where the
    treatment assignment A is randomized.

  - The real-world evidence study contains observations on (A,X,Y),
    where the treatment assignment A may be confounded and there exist
    unmeasured confounders.

## Installation with `devtools`:

``` r
devtools::install_github("shuyang1987/IntegrativeHTEcf")
```

### Main Paper: coming soon

The reference paper will come soon.

### Usage

IntHTEcf( A, X, X.hte, X.cf, Y, S,
nboots=50)

### Arguments

| Argument |                                                                                                                                                                  |
| -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A        | is a vector of treatment ( (n+m) x 1), where n is the RT sample size and m is the RW sample size                                                                 |
| X        | is a matrix of covariates without intercept ( (n+m) x p)                                                                                                         |
| X.hte    | is a matrix of covariates for the heterogeneous treatment effect function without intercept ( (n+m) x p1).                                                       |
| X.cf     | is a matrix of covariates for the confounding function without intercept ( (n+m) x p2).                                                                          |
| Y        | is a vector of outcome ( (n+m) x 1)                                                                                                                              |
| S        | is a vector of the binary indicator of belonging to the randomized trial (RT) sample; i.e., 1 if the unit belongs to the RT sample, and 0 otherwise ( (n+m) x 1) |
| nboots   | is the number of bootstrap samples                                                                                                                               |

### Value

|             |                                                                                                   |
| ----------- | ------------------------------------------------------------------------------------------------- |
| est.meta    | the HTE estimator based on a meta analysis                                                        |
| est.rct     | the HTE estimator using the semiparametric efficient estimation based only on RCT data            |
| est.int     | the HTE estimator using the semiparametric efficient estimation based on the combined RCT and RWD |
| att.meta    | the ATT estimator based on a meta analysis                                                        |
| att.rct     | the ATT estimator based only on RCT data                                                          |
| att.int     | the ATT estimator based on the combined RCT and RWD                                               |
| ve.meta     | variance estimates of est.meta                                                                    |
| ve.rct      | variance estimates of est.rct                                                                     |
| ve.int      | variance estimates of est.int                                                                     |
| ve.att.meta | variance estimates of att.meta                                                                    |
| ve.att.rct  | variance estimates of att.rct                                                                     |
| ve.att.int  | variance estimates of att.int                                                                     |

## Example

This is an example for illustration.

``` r

library(mgcv)
#> Loading required package: nlme
#> This is mgcv 1.8-26. For overview type 'help("mgcv-package")'.
library(MASS)
library(rootSolve)
#> Warning: package 'rootSolve' was built under R version 3.5.3
set.seed(1)
n=500
m=1000
## RCT: A,X,U
A1 <- rbinom(n,1,0.5)
X1 <- runif(n,-1,1)
U1 <- rnorm(n,1,1)
## RWE: A,X,U
A2 <- rbinom(m,1,0.5)
Sigma0 <- Sigma1 <- matrix(1,2,2)
Sigma0[1,2] <- Sigma0[2,1] <- -0.5
Sigma1[1,2] <- Sigma1[2,1] <-  0.5
XU0 <- mvrnorm(n = m, c(1,0), Sigma0)
XU1 <- mvrnorm(n = m, c(1,0), Sigma1)
X2 <- XU0[,1]*(1-A2)+XU1[,1]*A2
U2 <- XU0[,2]*(1-A2)+XU1[,2]*A2
## Combined data
A <- c(A1,A2)
X <- c(X1,X2)
U <- c(U1,U2)
S <- c(rep(1,n),rep(0,m))
Y <- 1+X+0.5*X^2+0.5*U+0.5*rnorm(n+m)+A*(1+0.75*X^2)
X.hte<-cbind(X,X^2)
X.cf <-cbind(X)
## True parameter values
psi<-c(1,0,0.75)
att.true<-2.5
X<-as.matrix(X)
X.hte<-as.matrix(X.hte)
X.cf<-as.matrix(X.cf)
IntegrativeHTEcf::IntHTEcf( A, X, X.hte, X.cf, Y, S, nboots=50)
#> $est.meta
#>      psi0      psi1      psi2 
#> 0.9464932 0.3981970 0.5516552 
#> 
#> $att.meta
#> [1] 2.447991
#> 
#> $est.rct
#>       psi0       psi1       psi2 
#>  1.0153667 -0.1781974  0.7627508 
#> 
#> $att.rct
#> [1] 1.276543
#> 
#> $est.int
#>       psi0       psi1       psi2       phi0       phi1 
#>  0.9802591 -0.1580110  0.7821850 -0.5007907  0.5754986 
#> 
#> $att.int
#> [1] 1.248274
#> 
#> $ve.meta
#>       psi0       psi1       psi2 
#> 0.01092300 0.03917656 0.02634711 
#> 
#> $ve.att.meta
#> [1] 0.006923964
#> 
#> $ve.rct
#>       psi0       psi1       psi2 
#> 0.01196633 0.01370939 0.05931311 
#> 
#> $ve.att.rct
#> [1] 0.003533183
#> 
#> $ve.int
#>         psi0         psi1         psi2         phi0         phi1 
#> 0.0038744174 0.0128593576 0.0008165552 0.0061280205 0.0170160074 
#> 
#> $ve.att.int
#> [1] 0.003828164
```
