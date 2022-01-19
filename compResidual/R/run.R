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
  use <- 1:nrow(res)%%dat$dim!=0 # no residual for last group
  res <- matrix(res$residual[use], nrow=(dat$dim-1))
  res
}


##' Residual Dirichlet model 
##' @param obs An observation (proportions) assumed to be Dirichlet distributed   
##' @param alpha concentration parameter
##' @return One step ahead randomized quantile residuals 
##' @details The model ...
##' @useDynLib compResidual
##' @export
##' @examples
##' a<-matrix(rep(1000*c(.2,.2,.1,.1,.1,.3), 100), nrow=6)
##' o<-rdirichlet(100,a[,1])
##' res<-resDir(o,a)
##' plot(as.vector(res))

resDir <- function(obs, alpha){
  if(sum(apply(obs, 2, sum))!=ncol(obs)) stop("Dirichlet observations should be proportions")
  dat<-list()
  dat$code <- 2 # Dirichlet   
  dat$dim <-nrow(obs)  
  dat$obs<-as.vector(obs)
  dat$alpha<-as.vector(alpha)
  dat$idx<-seq(1,length(obs), by=nrow(obs))-1
  param<-list(dummy=0)
  obj <- TMB::MakeADFun(dat, param, DLL="compResidual", silent=TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  res <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", method="cdf", trace=FALSE)
  use <- 1:nrow(res)%%dat$dim!=0 # no residual for last proportion
  res <- matrix(res$residual[use], nrow=(dat$dim-1))
  res
}
