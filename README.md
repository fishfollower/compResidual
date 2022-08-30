# Composition residuals 
Git page for `compResidual` R-package containing methods to compute residuals for observed compositions. The computed residuals will be independent, standardized, and normally distributed if the model is correctly describing the observations. The distributions currently implemented corresponds to common choices in age or length based assessment models (Multinomial, Dirichlet, Dirichlet-multinomial, multivariate-normal, and Logistic-normal).

The package is intended to be mainly useful in the situation where a purely fixed effects model is producing predictions/estimates corresponding to a set of observations, then those can be given and the proper residuals can be produced. 

Can be installed by typing: 

```R
TMB:::install.contrib("https://github.com/vtrijoulet/OSA_multivariate_dists/archive/main.zip")
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


## Principle

### One-step-ahead (OSA) quantile residuals for continuous non-normal independent observations

If the observations $(x_1,\ldots,x_n)$ are continuous, univariate, and independent, but originating from a distribution of which is not a normal distribution, but has cumulative distribution function (cdf) $F_x$, then independent and normally distributed residuals can be obtained via transformation. First, transforming the observations via the cdf will lead to quantities which follow a uniform distribution $u_i = F_x(x_i)$. This can be explained by the fact that $u_i \in (0,1)$ and the cdf of $u$ is $F_u(u) = P(F_x(X)\unicode{x003C} u) = P(X\unicode{x003C} F^{-1}_x(u)) = F_x(F^{-1}_x(u)) = u$, which is the distribution function for the uniform distribution. Secondly, transforming these uniformly distributed quantities $u_1,...,u_k$ by the inverse cdf of the standard normal distribution $\Phi^{-1}$ will lead to residuals which follow a standard normal distribution if the model is correct. This can be seen by calculating the cdf for the transformed $P(\Phi^{-1}(U)\unicode{x003C} r) = P(U\unicode{x003C} \Phi(r)) = \Phi(r)$ , so the wanted cdf. Collectively these quantile residuals are defined simply as: $r_i =\Phi^{-1}(F_x(x_i))$. The model is defining the cdf of the observations, so if the model is incorrect, then the residuals will deviate systematically from a standard normal distribution. 

### OSA randomized quantile residuals for discrete independent observations 

If the observations are discrete, but still univariate and independent, then the distribution function $F_x$ is a step function. In this discrete case, the transformation by the distribution function needs an additional step. The probability mass at a given value $x_i$ needs to be transformed onto the interval from $F(x_i-\epsilon)$ to $F(x_i)$. The transformation into uniform(0,1) distributed quantities can be achieved by sampling from the uniform distribution for that interval, so $u_i \sim U(F(x_i-\epsilon), F(x_i))$. The final step to get standard normal residuals is the same transformation by the inverse cdf of the standard normal distribution $r_i = \Phi^{-1}(u_i)$. These randomized quantile residuals will again have the desired properties (independence and normality).

### OSA residuals to remove dependence in dependent observations

The OSA residual of the i'th observation is computed using one of the methods described above depending on the observations being continuous or discrete (quantile or randomized quantile residuals), but instead of using the cdf of the observation in isolation, the cdf of the predicted distribution of the i'th prediction conditioned on all previous observations is used. This allows the resulting residuals to become independent standard normal if the model is correct. 

More details are available in the following paper \url{...}

### Additional data examples (Gulf of Maine Haddock)

Additional data for Gulf of Maine Haddock used in the paper is supplied in the file <a href="GOMhaddock.RData">GOMhaddock.RData</a>. 

```R
# example for fleet 2
load("GOMhaddock.RData") # load df 
obs<-GOMhaddock[GOMhaddock$Fleet==2, grep("^obsP",colnames(GOMhaddock))] # extract observations fleet 2
obs<-round(obs*GOMhaddock[GOMhaddock$Fleet==2, "ESS"]) # multiply by effective sample size and round  
pred<-GOMhaddock[GOMhaddock$Fleet==2, grep("^predP",colnames(GOMhaddock))] # extract predictions fleet 2

library(compResidual) # load library
res<-resMulti(t(obs), t(pred)) # calculate residuals 
plot(res) 
```

<p align="center">
  <img src="figs/fig2.png?raw=true">
</p>
