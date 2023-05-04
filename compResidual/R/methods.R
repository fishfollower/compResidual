##' Plot residual object 
##' @method plot cres
##' @param  x Residual object as returned from one of the residual functions
##' @param pick_one Number of the plot to provide if one of the 4 plots printed by default needs to be extracted 
##' @param maxLag Maximum lag to use for ACF plot
##' @param  ... extra arguments
##' @importFrom graphics par plot legend abline lines points text
##' @importFrom stats qqnorm
##' @details The function produces a 4 plot overview of the residuals by default or one specific plot if the 'pick_one' argument is provided.
##' 
##' Plot 1: Bubble plot of the residuals
##' Plot 2: Autocorrelation and cross-correlation plot for all dimensions of the residual matrix (rows, columns and diagonal) with 95% confidence interval
##' Plot 3: Q-Q plot of the residuals, each color is a row of the residual matrix (i.e. compositional group)
##' Plot 4: Simple plot of the residuals, each color is a row of the residual matrix (i.e. compositional group)   
##' @export
plot.cres<-function(x, pick_one=NULL, maxLag=NULL, ...){
  if (missing(pick_one)){
    op <- par(mfrow=c(2,2), mar=c(4,4,2,2)) 
    } else {
      if (!pick_one %in% 1:4) stop("pick_one should be a number between 1 and 4")
      op <- par(mfrow=c(1,1))
    }
  if (!is.null(rownames(x))) yname <- rownames(x) else yname=1:nrow(x)
  if (!is.null(colnames(x))) xname <- colnames(x) else xname=1:ncol(x)
  col <- rep(1:nrow(x), ncol(x))
  # bubble plot
  if (pick_one==1 || missing(pick_one)) {
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
    plotby(sample, composition, x, bubblescale = 0.3, xlab="", yaxt="n", xaxt="n",...)
    xTicks <- pretty(1:ncol(x))
    xTicks <- xTicks[xTicks!=0]
    yTicks <- pretty(1:nrow(x))
    yTicks <- yTicks[yTicks!=0]
    axis(1, at=xTicks, labels = xname[xTicks])
    axis(2, at=yTicks, labels = yname[yTicks])
    add_legend(x, cex.text=0.8, bubblescale = 0.3)
  }
  # ACF
  if (pick_one==2 || missing(pick_one)) {
    acfr <- acf_res(x, "row")
    acfc <- acf_res(x, "column")
    acfd <- acf_res(x, "diagonal")
    ci <- 1.96/sqrt(length(x))
    low <- as.numeric(min(c(acfr$acf, acfc$acf, acfd$acf)))
    ylow <- max(abs(c(ci,low)))
    lcol=c("#7fc97f", "#beaed4", "#fdc086")
    if (is.null(maxLag)) maxLag <- min(c(10, max(length(acfr$acf), length(acfc$acf), length(acfc$acf))-1))
    plot(y=as.vector(acfr$acf), x=as.vector(acfr$lag), type = "h", ylim=c(-ylow*1.2,1), xlim=c(0, maxLag), col=lcol[1], lwd=2, ylab="ACF", xlab="lag", ...)
    lines(y=as.vector(acfc$acf), x=as.vector(acfc$lag)+0.05, type = "h", col=lcol[2], lwd=2)
    lines(y=as.vector(acfd$acf), x=as.vector(acfd$lag)+0.1, type = "h", col=lcol[3], lwd=2)
    abline(h=0)
    abline(h=c(-ci,ci), lty=2)
    legend("topright", col=lcol, legend=c("row", "column", "diagonal"), bty="n", lwd=2)
    #acf(as.vector(x), main="") 
  }
  # Q-Q plot
  if (pick_one==3 || missing(pick_one)) {
    qqnorm(x, col=col, main="")
    abline(0,1)
    #abline(v = apply(x,1,mean), col = col, lwd = 2)
    legend("topleft", col=col, legend=yname, pch=1, bty="n", cex=0.8, ncol=2)
  }
  # plot vector
  if (pick_one==4 || missing(pick_one)) plot(as.vector(x), col=col, ylab="residuals", xlab="")
  par(op)
}
