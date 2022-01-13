##' Residual multinomial model 
##' @param obs An observation assumed to be multinomial   
##' @param pred predicted observations
##' @return One step ahead randomized quantile residuals 
##' @details The model ...
##' @useDynLib compResidual
##' @export
##' @examples
resMulti <- function(obs, pred){
  ret<-(obs-pred)/pred
  ret
}
