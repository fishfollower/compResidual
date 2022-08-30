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

### Additional data examples (Gulf of Maine Haddock)

| Fleet | Year | ESS | Obs. P1 | Obs. P2 | Obs. P3 | Obs. P4 | Obs. P5 | Obs. P6 | Obs. P7 | Obs. P8 | Obs. P9 | 
|  :--- |:--- |:--- |:--- |:--- |:--- |:--- |:--- |:--- |:--- |:--- |:--- | 
| 1 | 1977 | 60 | 0.01531181 | 0.6782595 | 0.02046705 | 0.1411919 | 0.07101912 | 0.0728273 | 0 | 0 | 0.0009233255 |
 | 1 | 1978 | 60 | 0 | 0.1093663 | 0.6688071 | 0.05031961 | 0.1059514 | 0.06091474 | 0.003093897 | 0 | 0.001546948 |
 | 1 | 1979 | 60 | 0 | 0.02647209 | 0.2201156 | 0.6203044 | 0.0720214 | 0.03898045 | 0.01781851 | 0.004287456 | 0 |
 | 1 | 1980 | 60 | 0 | 0.2598913 | 0.03057718 | 0.2220003 | 0.4014686 | 0.04223821 | 0.02808048 | 0.008165663 | 0.007578205 |
 | 1 | 1981 | 60 | 0.0005604633 | 0.428274 | 0.1925859 | 0.0783848 | 0.09154234 | 0.1454803 | 0.02460701 | 0.03133257 | 0.007232646 |
 | 1 | 1982 | 60 | 0.007888932 | 0.1610484 | 0.3942909 | 0.1610743 | 0.02610614 | 0.07811081 | 0.1239133 | 0.02787077 | 0.01969638 |
 | 1 | 1983 | 60 | 0.003092164 | 0.003550262 | 0.2394995 | 0.2795259 | 0.2265583 | 0.04254588 | 0.0724368 | 0.09966502 | 0.03312623 |
 | 1 | 1984 | 60 | 0.0007342144 | 0.05445423 | 0.03053108 | 0.3658835 | 0.1570607 | 0.2233235 | 0.03805678 | 0.03964758 | 0.09030837 |
 | 1 | 1985 | 60 | 0.0006824386 | 0.02289961 | 0.2650895 | 0.06513497 | 0.270094 | 0.1152563 | 0.1835002 | 0.03594177 | 0.04140127 |
 | 1 | 1986 | 60 | 0.004496497 | 0.01129353 | 0.1918854 | 0.3751961 | 0.08501516 | 0.1192095 | 0.09034822 | 0.1071839 | 0.01537175 |
 | 1 | 1987 | 60 | 0 | 0.058556 | 0.09863559 | 0.3015918 | 0.1387152 | 0.09778283 | 0.1617396 | 0.09607732 | 0.04690165 |
 | 1 | 1988 | 20 | 0.001845018 | 0.003075031 | 0.07626076 | 0.07564576 | 0.3370234 | 0.3419434 | 0.04674047 | 0.09225092 | 0.02521525 |
 | 1 | 1989 | 20 | 0.01072797 | 0.1777778 | 0.02681992 | 0.3249042 | 0.1478927 | 0.183908 | 0.1149425 | 0.006130268 | 0.006896552 |
 | 1 | 1990 | 20 | 0.03019845 | 0.008628128 | 0.6173425 | 0.007333909 | 0.124245 | 0.07592752 | 0.1186368 | 0.01768766 | 0 |
 | 1 | 1991 | 20 | 0.0190184 | 0.04417178 | 0.1 | 0.3595092 | 0.1742331 | 0.1711656 | 0.07730061 | 0.03558282 | 0.0190184 |
 | 1 | 1992 | 20 | 0.01058201 | 0.07701352 | 0.5549677 | 0.2145797 | 0.1122869 | 0.01293357 | 0.006466784 | 0 | 0.0111699 |
 | 1 | 1993 | 60 | 0.03318386 | 0.1802691 | 0.3255605 | 0.206278 | 0.08878924 | 0.09865471 | 0.04125561 | 0.01524664 | 0.01076233 |
 | 1 | 1994 | 60 | 0.05963303 | 0.2174312 | 0.4082569 | 0.1247706 | 0.03119266 | 0.08440367 | 0.05229358 | 0.01559633 | 0.006422018 |
 | 1 | 1995 | 60 | 0.01004838 | 0.2653517 | 0.3368068 | 0.2817268 | 0.03796055 | 0.02344622 | 0.01749163 | 0.01600298 | 0.01116487 |
 | 1 | 1996 | 60 | 0.01134982 | 0.0952574 | 0.5249291 | 0.2290231 | 0.0664775 | 0.01661938 | 0.02877989 | 0.02269964 | 0.004864208 |
 | 1 | 1997 | 60 | 0.003070809 | 0.01318642 | 0.3013006 | 0.4638728 | 0.1627529 | 0.03414017 | 0.01246387 | 0.005057803 | 0.004154624 |
 | 1 | 1998 | 60 | 0.01243301 | 0.05101822 | 0.05380493 | 0.2844587 | 0.4132905 | 0.1129689 | 0.03729904 | 0.01843516 | 0.01629153 |
 | 1 | 1999 | 60 | 0.01596867 | 0.01144923 | 0.1190118 | 0.1982525 | 0.2916541 | 0.2084965 | 0.1159988 | 0.02139199 | 0.01777644 |
 | 1 | 2000 | 60 | 0.004232058 | 0.1209663 | 0.1165579 | 0.1883266 | 0.1147946 | 0.2265914 | 0.1271381 | 0.05607477 | 0.04531829 |
 | 1 | 2001 | 60 | 0.0004122012 | 0.04053311 | 0.3230283 | 0.1835669 | 0.1330036 | 0.1199505 | 0.1108821 | 0.05550976 | 0.03311349 |
 | 1 | 2002 | 60 | 0.0005682625 | 0.003409575 | 0.03949425 | 0.3911067 | 0.1663589 | 0.1568405 | 0.04560307 | 0.1000142 | 0.09660463 |
 | 1 | 2003 | 140 | 0.0001213298 | 0.01310362 | 0.008371754 | 0.06563941 | 0.6150206 | 0.1098034 | 0.07643776 | 0.02620723 | 0.08529483 |
 | 1 | 2004 | 140 | 0.003031834 | 0.004926731 | 0.02337039 | 0.04181405 | 0.09310258 | 0.6526023 | 0.07465892 | 0.04257201 | 0.06392117 |
 | 1 | 2005 | 140 | 0.0002017553 | 0.05074145 | 0.008372844 | 0.04993443 | 0.08564511 | 0.139312 | 0.539191 | 0.05417129 | 0.07243014 |
 | 1 | 2006 | 140 | 0.001666239 | 0.002947962 | 0.1663676 | 0.008972058 | 0.06870033 | 0.09189951 | 0.1072802 | 0.4738529 | 0.07831325 |
 | 1 | 2007 | 140 | 0.008162411 | 0.02710339 | 0.01883633 | 0.3501465 | 0.01192968 | 0.05713688 | 0.04531185 | 0.09187945 | 0.3894935 |
 | 1 | 2008 | 140 | 0.0001244555 | 0.02115744 | 0.06695706 | 0.02414437 | 0.511761 | 0.00609832 | 0.05289359 | 0.03696329 | 0.2799004 |
 | 1 | 2009 | 140 | 0.0001690046 | 0.006253169 | 0.0704749 | 0.08720635 | 0.02551969 | 0.4985635 | 0.008957242 | 0.05408146 | 0.2487747 |
 | 1 | 2010 | 140 | 0.004573439 | 0.01284081 | 0.02638522 | 0.06772208 | 0.09076517 | 0.03324538 | 0.5189094 | 0.006332454 | 0.239226 |
 | 1 | 2011 | 140 | 0.01509078 | 0.04715869 | 0.0108465 | 0.01108229 | 0.1200189 | 0.1110587 | 0.03819854 | 0.4263145 | 0.2202311 |
 | 1 | 2012 | 140 | 0.00536441 | 0.2389937 | 0.0963744 | 0.03884573 | 0.0227525 | 0.1250462 | 0.06992231 | 0.02552719 | 0.3771735 |
 | 1 | 2013 | 140 | 0.03674322 | 0.06945024 | 0.4942241 | 0.08796103 | 0.03729993 | 0.01475296 | 0.08058455 | 0.04272791 | 0.1362561 |
 | 1 | 2014 | 140 | 0.04247267 | 0.1426619 | 0.1145921 | 0.5550883 | 0.0366905 | 0.01408747 | 0.005361648 | 0.03206476 | 0.05698066 |
 | 1 | 2015 | 140 | 0.004388186 | 0.2724895 | 0.2860759 | 0.09518987 | 0.295443 | 0.0178903 | 0.008945148 | 0.001350211 | 0.01822785 |
 | 1 | 2016 | 140 | 0.002700449 | 0.01219173 | 0.5296057 | 0.1796593 | 0.04991859 | 0.2105953 | 0.005718597 | 0.001628212 | 0.007982209 |




More details are available in the following paper \url{...}




