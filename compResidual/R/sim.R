##' Simulate Dirichlet observations 
##' @param n Number of simulations   
##' @param alpha concentration parameter
##' @return Matrix of Dirichlet observations 
##' @details The model ...
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
##' @details The model ...
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
