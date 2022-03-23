##' Residuals for multivariate-normal observations 
##' @param obs Matrix of observations assumed to be multivariate-normal distributed   
##' @param mu Matrix of predicted observations
##' @param S Covariance matrix, either one matrix to use for columns of obs or a list of covariance matrices for each column
##' @param ... Additional arguments to pass to TMB::oneStepPredict
##' @return One step ahead quantile residuals 
##' @details The model estimates the one step ahead quantile residuals for observations that are multivariate-normal distributed.
##'  
##' The matrices should be created such as each column is one sample from a multivariate-normal distribution, the columns are assumed independent (i.e. the compositions are the rows). 
##' @useDynLib compResidual
##' @export
##' @examples
##' mu <- matrix(0,ncol=100, nrow=6)
##' sigma <- diag(nrow(mu))+0.5
##' o <- replicate(100, MASS::mvrnorm(1,mu[,1],sigma))
##' res <- resmvnorm(o, mu, sigma)
##' plot(res)

resmvnorm <- function(obs, mu, S, ...){
  obs<-as.matrix(obs)
  res <- c()
  for (k in 1:ncol(obs)){ # for each column (in case S is a list)
    dat<-list()
    dat$code <- 0 # Multivariate-normal
    dat$obs<-as.vector(obs[,k])
    dat$mu<-as.vector(mu[,k])
    if (is.matrix(S)) dat$S <- as.matrix(S) else dat$S <- as.matrix(S[[k]])
    param<-list(dummy=0)
    obj <- TMB::MakeADFun(dat, param, DLL="compResidual", silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr)
    tmp <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", trace=FALSE, ...)
    res <- cbind(res, tmp$residual)
  }
  class(res)<-"cres"
  res
}


##' Residuals for multinomial observations 
##' @param obs Matrix of observations assumed to be multinomial distributed  
##' @param pred Matrix of predicted observations or probabilities
##' @param ... Additional arguments to pass to TMB::oneStepPredict
##' @return One step ahead randomized quantile residuals 
##' @details The model estimates the one step ahead randomized quantile residuals for observations that are multinomial distributed.
##'  
##' The matrices should be created such as each column is one sample from a multinomial distribution, the columns are assumed independent (i.e. the compositions are the rows). 
##' @useDynLib compResidual
##' @export
##' @examples
##' o<-rmultinom(100, 25, c(.2,.2,.1,.1,.1,.3))
##' p<-matrix(rep(25*c(.2,.2,.1,.1,.1,.3), 100), nrow=6)
##' res<-resMulti(o,p)
##' plot(res)

resMulti <- function(obs, pred, ...){
  obs<-as.matrix(obs)
  dat<-list()
  dat$code<-1 # multinomial    
  dat$dim<-nrow(obs)  
  dat$obs<-as.vector(obs)  
  dat$pred<-as.vector(pred)
  dat$idx<-seq(1,length(obs), by=nrow(obs))-1
  param<-list(dummy=0)
  obj <- TMB::MakeADFun(dat, param, DLL="compResidual", silent=TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  res <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", discrete=TRUE, method="cdf", trace=FALSE, ...)
  use <- 1:nrow(res)%%dat$dim!=0 # no residual for last group
  res <- matrix(res$residual[use], nrow=(dat$dim-1))
  class(res)<-"cres"
  res
}


##' Residuals for Dirichlet observations 
##' @param obs Matrix of observations (proportions) assumed to be Dirichlet distributed   
##' @param alpha Matrix of concentration parameters
##' @param ... Additional arguments to pass to TMB::oneStepPredict
##' @return One step ahead quantile residuals 
##' @details The model estimates the one step ahead quantile residuals for observations that are Dirichlet distributed.
##'  
##' The matrices should be created such as each column is one sample from a Dirichlet distribution, the columns are assumed independent (i.e. the compositions are the rows). 
##' @useDynLib compResidual
##' @export
##' @examples
##' a<-matrix(rep(1000*c(.2,.2,.1,.1,.1,.3), 100), nrow=6)
##' o<-rdirichlet(100,a[,1])
##' res<-resDir(o,a)
##' plot(res)

resDir <- function(obs, alpha, ...){
  obs<-as.matrix(obs)
  if(!all.equal(apply(obs, 2, sum), rep(1,ncol(obs)))) stop("Dirichlet observations should be proportions, so sum to 1")
  dat<-list()
  dat$code <- 2 # Dirichlet   
  dat$dim <-nrow(obs)  
  dat$obs<-as.vector(obs)
  dat$alpha<-as.vector(alpha)
  dat$idx<-seq(1,length(obs), by=nrow(obs))-1
  param<-list(dummy=0)
  obj <- TMB::MakeADFun(dat, param, DLL="compResidual", silent=TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  res <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", method="cdf", trace=FALSE, ...)
  use <- 1:nrow(res)%%dat$dim!=0 # no residual for last proportion
  res <- matrix(res$residual[use], nrow=(dat$dim-1))
  class(res)<-"cres"
  res
}


##' Residuals for Dirichlet-multinomial observations 
##' @param obs Matrix of observations assumed to be Dirichlet-multinomial distributed   
##' @param alpha Matrix of concentration parameters
##' @param ... Additional arguments to pass to TMB::oneStepPredict
##' @return One step ahead randomized quantile residuals 
##' @details The model estimates the one step ahead randomized quantile residuals for observations that are Dirichlet-multinomial distributed.
##'  
##' The matrices should be created such as each column is one sample from a Dirichlet-multinomial distribution, the columns are assumed independent (i.e. the compositions are the rows). 
##' @useDynLib compResidual
##' @export
##' @examples
##' a<-matrix(rep(1000*c(.2,.2,.1,.1,.1,.3), 100), nrow=6)
##' o<-rdirM(100,1000, a[,1])
##' res<-resDirM(o,a)
##' plot(res)

resDirM <- function(obs, alpha, ...){
  obs<-as.matrix(obs)
  dat<-list()
  dat$code <- 3 # Dirichlet-multinomial   
  dat$dim <-nrow(obs)  
  dat$obs<-as.vector(obs)
  dat$alpha<-as.vector(alpha)
  dat$idx<-seq(1,length(obs), by=nrow(obs))-1
  param<-list(dummy=0)
  obj <- TMB::MakeADFun(dat, param, DLL="compResidual", silent=TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  res <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", method="cdf", discrete=TRUE, trace=FALSE, ...)
  use <- 1:nrow(res)%%dat$dim!=0 # no residual for last group
  res <- matrix(res$residual[use], nrow=(dat$dim-1))
  class(res)<-"cres"
  res
}


##' Residuals for logistic-normal observations 
##' @param obs Matrix of observations (proportions) assumed to be logistic-normal distributed   
##' @param mu Matrix of predicted observations
##' @param S Covariance matrix, either one matrix to use for columns of obs or a list of covariance matrices for each column
##' @param do_mult 0 if additive logistic-normal, 1 if multiplicative logistic-normal
##' @param ... Additional arguments to pass to TMB::oneStepPredict
##' @return One step ahead quantile residuals 
##' @details The model estimates the one step ahead quantile residuals for observations that are logistic-normal distributed.
##'  
##' The matrices should be created such as each column is one sample from a logistic-normal distribution, the columns are assumed independent (i.e. the compositions are the rows). 
##' @useDynLib compResidual
##' @export
##' @examples
##' mu <- matrix(0,ncol=100, nrow=5)
##' sigma <- diag(nrow(mu))
##' do_mult <- F
##' o <- sapply(1:100, function(x) rlogistN(mu[,1], sigma, do_mult))
##' res <- reslogistN(o, mu, sigma, do_mult)
##' plot(res)

reslogistN <- function(obs, mu, S, do_mult, ...){
  obs<-as.matrix(obs)
  res <- c()
  for (k in 1:ncol(obs)){ # for each column (in case S is a list)
    dat<-list()
    dat$code <- 4 # Logistic-normal
    dat$do_mult <- as.integer(do_mult) # 0 = Additive logistic-normal, 1 = Multiplicative logistic-normal
    dat$obs<-as.vector(obs[,k])
    dat$mu<-as.vector(mu[,k])
    if (is.matrix(S)) dat$S <- as.matrix(S) else dat$S <- as.matrix(S[[k]])
    param<-list(dummy=0)
    obj <- TMB::MakeADFun(dat, param, DLL="compResidual", silent=TRUE)
    opt <- nlminb(obj$par, obj$fn, obj$gr)
    suppressWarnings(tmp <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", method="oneStepGeneric", trace=FALSE, ...))
    use <- 1:(length(dat$obs)-1) # no residual for last group
    res <- cbind(res, tmp$residual[use])
  }
  class(res)<-"cres"
  res
}
