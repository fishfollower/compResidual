##' Simulate Dirichlet observations 
##' @param n Number of simulations   
##' @param alpha concentration parameter
##' @return Matrix of Dirichlet observations 
##' @importFrom stats rgamma
##' @details This function simulates Dirichlet observations.
##' @useDynLib compResidual
##' @export

rdirichlet <- function(n,alpha){
  rd1 <- function(alpha){
    x <- rgamma(length(alpha), alpha)
    p <- x/sum(x)
    p
  }
  replicate(n,rd1(alpha))
}


##' Simulate Dirichlet-multinomial observations 
##' @param n Number of simulations 
##' @param N Sum of observations  
##' @param alpha concentration parameter
##' @return Matrix of Dirichlet-multinomial observations 
##' @details This function simulates Dirichlet-multinomial observations.
##' @useDynLib compResidual
##' @export

rdirM <- function(n,N,alpha){
  rd1 <- function(alpha){
    x <- rgamma(length(alpha), alpha)
    p <- x/sum(x)
    p
  }
  replicate(n, as.vector(stats::rmultinom(1, N, rd1(alpha))))
}


##' Simulate logistic-normal observations 
##' @param mu Vector of mean values 
##' @param Sigma Covariance matrix  
##' @param mult Multiplicative logistic-normal if TRUE, otherwise additive
##' @return Matrix of logistic-normal observations
##' @details This function simulates logistic-normal observations.
##' @useDynLib compResidual
##' @export
rlogistN = function(mu, Sigma, mult=FALSE){
  x = MASS::mvrnorm(1,mu,Sigma)
  if(mult) p = exp(x)/cumprod(1+exp(x))
  else p = exp(x)/(1+sum(exp(x)))
  return(c(p, 1-sum(p)))
}

##' Transform logistic-normal to normal observations
##' @param p Vector of mean values 
##' @param mult Multiplicative logistic-normal if TRUE, otherwise additive
##' @return Transformed observations
##' @details This function transforms logistic-normal observations to normally distributed observations with dimension (length(p)-1).
##' @useDynLib compResidual
##' @export
logistictransf = function(p, mult = FALSE){
  ps = p[-length(p)]
  if(mult) ys = log(ps/(1-cumsum(ps))) #multiplicative
  else ys = log(ps) - log(1-sum(ps)) #additive
  return(ys)
}
