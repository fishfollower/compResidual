##' Residual multinomial model 
##' @param obs An observation matrix assumed to be multinomial   
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
##' @param obs An observation matrix (proportions) assumed to be Dirichlet distributed   
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
  res <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", method="cdf", trace=FALSE)
  use <- 1:nrow(res)%%dat$dim!=0 # no residual for last proportion
  res <- matrix(res$residual[use], nrow=(dat$dim-1))
  res
}


##' Residual Dirichlet-multinomial model 
##' @param obs An observation matrix (proportions) assumed to be Dirichlet-multinomial distributed   
##' @param alpha concentration parameter
##' @return One step ahead randomized quantile residuals 
##' @details The model ...
##' @useDynLib compResidual
##' @export
##' @examples
##' a<-matrix(rep(1000*c(.2,.2,.1,.1,.1,.3), 100), nrow=6)
##' o<-rdirM(100,1000, a[,1])
##' res<-resDirM(o,a)
##' plot(as.vector(res))

resDirM <- function(obs, alpha){
  dat<-list()
  dat$code <- 3 # Dirichlet-multinomial   
  dat$dim <-nrow(obs)  
  dat$obs<-as.vector(obs)
  dat$alpha<-as.vector(alpha)
  dat$idx<-seq(1,length(obs), by=nrow(obs))-1
  param<-list(dummy=0)
  obj <- TMB::MakeADFun(dat, param, DLL="compResidual", silent=TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  res <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", method="cdf", discrete=TRUE, trace=FALSE)
  use <- 1:nrow(res)%%dat$dim!=0 # no residual for last group
  res <- matrix(res$residual[use], nrow=(dat$dim-1))
  res
}


##' Residual Logistic-normal model 
##' @param obs An observation matrix (proportions) assumed to be Logistic-normal distributed   
##' @param mu predicted observations
##' @param S covariance matrix, either one matrix to use for columns of obs or a list of covariance matrices for each column
##' @param do_mult 0 if additive logistic-normal, 1 if multiplicative logistic-normal
##' @return One step ahead randomized quantile residuals 
##' @details The model ...
##' @useDynLib compResidual
##' @export
##' @examples
##' mu <- matrix(0,ncol=100, nrow=5)
##' sigma <- diag(nrow(mu))
##' do_mult <- F
##' o <- sapply(1:100, function(x) rlogistN(mu[,1], sigma, do_mult))
##' res <- reslogistN(o, mu, sigma, do_mult)
##' plot(as.vector(res))

reslogistN <- function(obs, mu, S, do_mult){
  res2 <- c()
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
    suppressWarnings(res <- TMB::oneStepPredict(obj, observation.name="obs", data.term.indicator="keep", method="oneStepGeneric", trace=FALSE))
    use <- 1:(length(dat$obs)-1) # no residual for last group
    res2 <- cbind(res2, res$residual[use])
  }
  res2
}
