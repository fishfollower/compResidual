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
