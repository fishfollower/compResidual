# Composition residuals 
Git page for `compResidual` R-package containing methods to compute residuals for observed compositions. The computed residuals will be independent, standardized, and normally distributed if the model is correctly describing the observations. The distributions currently implemented corresponds to common choices in age or length based assessment models (Multinomial, Dirichlet, Dirichlet-multinomial, multivariate-normal, and Logistic-normal).

The package is intended to be mainly useful in the situation where a purely fixed effects model is producing predictions/estimates corresponding to a set of observations, then those can be given and the proper residuals can be produced. 

Can be installed by typing: 

```R
devtools::install_github("fishfollower/compResidual/compResidual")
```

(To ensure only installing in 64-bit some windows users may need to install with:

```R
devtools::install_github("fishfollower/compResidual/compResidual", INSTALL_opts=c("--no-multiarch"))
```
)

## A quick example:

Imagine that a model assumes that the composition observations $(X)_{i,j}$ are following a multinomial distribution. $X$ is a matrix where each column is an independent multinomially distributed vector. Corresponding to the matrix $X$ a matrix with estimated probability vectors are collected in a matrix $P$. Then the residuals can be computed like:

```R
library(compResidual)

X<-matrix(c(
  9,    1,    5,    9,    9,    3,    6,    3,    3,     4,
  1,    5,    7,    1,    3,    6,    8,    9,    6,     5,
  3,    7,    4,    2,    1,    8,    4,    5,    1,     2,
  4,    2,    3,    5,    3,    1,    0,    1,    0,     2,
  1,    3,    3,    0,    5,    5,    4,    3,    1,     0,
  7,    7,    3,    8,    4,    2,    3,    4,   14,    12
), nrow=6, ncol=10, byrow=TRUE) 

P<-matrix(c(
  0.32, 0.08, 0.16, 0.24, 0.32, 0.20, 0.20, 0.16, 0.16, 0.16,
  0.16, 0.16, 0.24, 0.20, 0.12, 0.16, 0.32, 0.28, 0.20, 0.20,
  0.12, 0.24, 0.20, 0.12, 0.04, 0.24, 0.16, 0.12, 0.04, 0.08,
  0.04, 0.16, 0.16, 0.12, 0.04, 0.08, 0.00, 0.08, 0.04, 0.20,
  0.12, 0.08, 0.12, 0.00, 0.12, 0.12, 0.12, 0.08, 0.04, 0.04,
  0.24, 0.28, 0.12, 0.32, 0.36, 0.20, 0.20, 0.28, 0.52, 0.32
), nrow=6, ncol=10, byrow=TRUE) 

res<-resMulti(X,P)

plot(res)
```

<p align="center">
  <img src="figs/fig1.png?raw=true">
</p>

Instead of providing the probability matrix it is also allowed to provide the predicted observations (N*P).  

