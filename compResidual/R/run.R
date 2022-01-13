##' Residual multinomial model 
##' @param obs An observation assumed to be multinomial   
##' @param pred predicted observations
##' @return One step ahead randomized quantile residuals 
##' @details The model ...
##' @useDynLib compResidual
##' @export
##' @examples
##' o<-rmultinom(100, 25, c(.2,.2,.1,.1,.1,.3))
##' p<-matrix(rep(25*c(.2,.2,.1,.1,.1,.3), 100), nrow=6)
##' res<-resMulti(o,p)
##' plot(as.vector(res))

resMulti <- function(obs, pred){
  dat<-list()
  dat$code<-1 # multinomial    
  dat$dim<-nrow(obs)  
  dat$obs<-as.vector(obs)  
  dat$pred<-as.vector(pred)
  dat$idx<-seq(1,length(obs), by=nrow(obs))-1
  param<-list(dummy=0)
  obj <- TMB::MakeADFun(dat, param, DLL="compResidual", silent=TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  res <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", discrete=TRUE, method="cdf", trace=FALSE)
  res <- matrix(res$residual, nrow=dat$dim)
  res
}
