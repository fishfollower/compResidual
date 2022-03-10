##' Plot residual object 
##' @method plot cres
##' @param  x residual object as returned from one of the residual functions.
##' @param  ... extra arguments
##' @importFrom graphics par
##' @details gives a 4 plot overview of the residuals.  
##' @export
plot.cres<-function(x, ...){
  op <- par(mfrow=c(2,2), mar=c(4,4,2,2))
  col <- rep(1:nrow(x), ncol(x))
  # bubble plot
  add_legend <- function(x, cex.text=1, ...) {
    # opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
    #             mar=c(0, 0, 0, 0), new=TRUE)
    # on.exit(par(opar))
    #plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    zscale <- pretty(x,min.n=4)
    uu<-par("usr")
    yy<-rep(uu[3]+.03*(uu[4]-uu[3]), length(zscale))
    xx<-seq(uu[1]+.10*(uu[2]-uu[1]),uu[1]+.4*(uu[2]-uu[1]), length=length(zscale))
    text(xx,yy,labels=zscale, cex=cex.text)
    colb <- ifelse(zscale<0, rgb(1, 0, 0, alpha=.5), rgb(0, 0, 1, alpha=.5))
    bs<-1
    if("bubblescale"%in%names(list(...))) bs <- list(...)$bubblescale
    points(xx,yy,cex=sqrt(abs(zscale))/max(sqrt(abs(zscale)), na.rm=TRUE)*5*bs, pch=19, col=colb)
  }
  sample <- rep(1:ncol(x), each=nrow(x))
  composition <- rep(1:nrow(x), ncol(x))
  stockassessment::plotby(sample, composition, x, bubblescale = 0.3, xlab="")
  add_legend(x, cex.text=0.8, bubblescale = 0.3)
  # ACF
  acf(as.vector(x), main="")  
  # Q-Q plot
  qqnorm(x, col=col, main="")
  abline(0,1)
  #abline(v = apply(x,1,mean), col = col, lwd = 2)
  legend("topleft", col=col, legend=1:nrow(x), pch=1, bty="n", cex=0.8, ncol=2)
  # plot vector
  plot(as.vector(x), col=col, ylab="residuals", xlab="")
  par(op)
}