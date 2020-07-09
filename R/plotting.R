

#' @name plotmse.cw
#' @export
plotmse.cw <- function(dat, set, resMSE,
                       trendline=TRUE, uncert = TRUE, med = TRUE,
                       hcrs=NA, ylim = NULL){
    if(any(!is.na(hcrs))){
        resMSEnew <- vector("list",length(hcrs))
        for(i in 1:length(hcrs)){
            hcri <- hcrs[i]
            resMSEnew[[i]] <- resMSE[[hcri]]
        }
        resMSE <- resMSEnew
        set$hcr <- hcrs
    }
    ## summary (median, 95% CIs)
    res <- sumMSE(resMSE)
    ## vars
    nms <- length(resMSE)
    xhist <- seq(1,dat$ny,1)
    xsim <- seq(dat$ny,dat$ny+set$nysim,1)
    idxsim <- (dat$ny):(dat$ny+set$nysim)
    idxhist <- 1:dat$ny
    xall <- 1:(dat$ny+set$nysim)
    if(is.null(ylim)) ylim <- c(0.8,1.2) * range(lapply(res,function(x) x$CW))
    cols <- rainbow(nms)
    ## historic
    i <- 1 ## historic pattern the same between mss
    llhist <- res[[i]]$CW[1,idxhist]
    ulhist <- res[[i]]$CW[3,idxhist]
    medhist <- res[[i]]$CW[2,idxhist]
    msy <- resMSE[[1]][[1]]$refs$MSY ## median(resMSE[[1]][[1]]$refdist$MSY) ##
##    tmp <- quantile(resMSE[[1]][[1]]$refdist$MSY,probs = c(0.025,0.975))
    llmsy <- NA##tmp[1]
    ulmsy <- NA##tmp[2]
    ## plot
    plot(xhist, medhist,
         xlim = c(0, dat$ny+set$nysim),
         ylim = ylim,
         ty='n', lwd=2,
         ylab="Catch", xlab = "")
    polygon(x = c((dat$ny-(set$nyhist-1)):dat$ny,rev((dat$ny-(set$nyhist-1)):dat$ny)),
            y = c(rep(ylim[1]-1e4,set$nyhist),rep(1.5*ylim[2],set$nyhist)),
            border=NA, col="grey95")
    polygon(x = c(xhist,rev(xhist)), y = c(llhist,rev(ulhist)),
            border=NA, col=rgb(t(col2rgb("grey30"))/255,alpha=0.2))
    lines(xhist, medhist, lwd=2)
    polygon(x = c(-2,rep(1.2*max(xall),2),-2), y = c(rep(llmsy,2),rep(ulmsy,2)),
            border=NA, col=rgb(t(col2rgb("grey40"))/255,alpha=0.2))
    abline(h=msy,lty=2)
    abline(h=0,lty=2)
    ## projection
    if(uncert){
        for(i in 1:nms){
            llsim <- res[[i]]$CW[1,idxsim]
            ulsim <- res[[i]]$CW[3,idxsim]
            polygon(x = c(xsim,rev(xsim)), y = c(llsim,rev(ulsim)),
                    border=NA, col=rgb(t(col2rgb(cols[i]))/255,alpha=0.2))
        }
    }
    for(i in 1:nms){
        medsim <- res[[i]]$CW[2,idxsim]
        if(med) lines(xsim, medsim, lwd=2, col=cols[i])
        if(is.numeric(trendline)){
            for(j in 1:length(trendline))
                lines(xsim, resMSE[[i]][[trendline[j]]]$CW[idxsim], col=cols[i])
        }else if(trendline)
            lines(xsim, resMSE[[i]][[1]]$CW[idxsim], col=cols[i])
    }
    abline(v=dat$ny, col="grey60",lwd=2)
    abline(v=max(which(dat$Fvals==0)), col="grey60",lwd=2,lty=2)
    legend("topleft", legend=set$hcr,
           col=cols, bty="n", lwd=2,lty=1)
    title("Catch")
    box()
}


#' @name plotmse.b
#' @export
plotmse.b <- function(dat, set, resMSE,
                      trendline=TRUE, uncert = TRUE, med = TRUE,
                      hcrs=NA, ylim = NULL){
    if(any(!is.na(hcrs))){
        resMSEnew <- vector("list",length(hcrs))
        for(i in 1:length(hcrs)){
            hcri <- hcrs[i]
            resMSEnew[[i]] <- resMSE[[hcri]]
        }
        resMSE <- resMSEnew
        set$hcr <- hcrs
    }
    ## summary (median, 95% CIs)
    res <- sumMSE(resMSE)
    ## vars
    nms <- length(resMSE)
    xhist <- seq(1,dat$ny,1)
    xsim <- seq(dat$ny,dat$ny+set$nysim,1)
    idxsim <- (dat$ny):(dat$ny+set$nysim)
    idxhist <- 1:dat$ny
    xall <- 1:(dat$ny+set$nysim)
    if(is.null(ylim)) ylim <- c(0.8,1.2) * range(lapply(res,function(x) x$TSB))
    cols <- rainbow(nms)
    ## historic
    i <- 1 ## historic pattern the same between mss
    llhist <- res[[i]]$TSB[1,idxhist]
    ulhist <- res[[i]]$TSB[3,idxhist]
    medhist <- res[[i]]$TSB[2,idxhist]
    bmsy <- resMSE[[1]][[1]]$refs$Bmsy ##median(resMSE[[1]][[1]]$refdist$Bmsy) ##
##    tmp <- quantile(resMSE[[1]][[1]]$refdist$Bmsy,probs = c(0.025,0.975))
    llbmsy <- NA ##tmp[1]
    ulbmsy <- NA ##tmp[2]
    ## plot
    plot(xhist, medhist,
         xlim = c(0, dat$ny+set$nysim),
         ylim = ylim,
         ty='n', lwd=2,
         ylab="Total biomass", xlab = "")
    polygon(x = c((dat$ny-(set$nyhist-1)):dat$ny,rev((dat$ny-(set$nyhist-1)):dat$ny)),
            y = c(rep(ylim[1]-1e4,set$nyhist),rep(1.5*ylim[2],set$nyhist)),
            border=NA, col="grey95")
    polygon(x = c(xhist,rev(xhist)), y = c(llhist,rev(ulhist)),
            border=NA, col=rgb(t(col2rgb("grey30"))/255,alpha=0.2))
    lines(xhist, medhist, lwd=2)
    polygon(x = c(-2,rep(1.2*max(xall),2),-2), y = c(rep(llbmsy,2),rep(ulbmsy,2)),
            border=NA, col=rgb(t(col2rgb("grey40"))/255,alpha=0.2))
    abline(h=bmsy,lty=2)
    abline(h=resMSE[[1]][[1]]$refs$B0,lty=2)
    ## projection
    if(uncert){
        for(i in 1:nms){
            llsim <- res[[i]]$TSB[1,idxsim]
            ulsim <- res[[i]]$TSB[3,idxsim]
            polygon(x = c(xsim,rev(xsim)), y = c(llsim,rev(ulsim)),
                    border=NA, col=rgb(t(col2rgb(cols[i]))/255,alpha=0.2))
        }
    }
    for(i in 1:nms){
        medsim <- res[[i]]$TSB[2,idxsim]
        if(med) lines(xsim, medsim, lwd=2, col=cols[i])
        if(is.numeric(trendline)){
            for(j in 1:length(trendline))
                lines(xsim, resMSE[[i]][[trendline[j]]]$TSB[idxsim], col=cols[i])
        }else if(trendline)
            lines(xsim, resMSE[[i]][[1]]$TSB[idxsim], col=cols[i])
    }
    abline(v=dat$ny, col="grey60",lwd=2)
    abline(v=max(which(dat$Fvals==0)), col="grey60",lwd=2,lty=2)
    legend("topleft", legend=set$hcr,
           col=cols, bty="n", lwd=2,lty=1)
    title("Biomass")
    box()
}


#' @name plotmse.f
#' @export
plotmse.f <- function(dat, set, resMSE,
                      trendline=TRUE, uncert = TRUE, med = TRUE,
                      hcrs=NA, ylim=NULL){
    if(any(!is.na(hcrs))){
        resMSEnew <- vector("list",length(hcrs))
        for(i in 1:length(hcrs)){
            hcri <- hcrs[i]
            resMSEnew[[i]] <- resMSE[[hcri]]
        }
        resMSE <- resMSEnew
        set$hcr <- hcrs
    }
    ## summary (median, 95% CIs)
    res <- sumMSE(resMSE)

    ## vars
    nms <- length(resMSE)
    xhist <- seq(1,dat$ny,1)
    xsim <- seq(dat$ny,dat$ny+set$nysim,1)
    idxsim <- (dat$ny):(dat$ny+set$nysim)
    idxhist <- 1:dat$ny
    xall <- 1:(dat$ny+set$nysim)
    if(is.null(ylim)) ylim <- c(0.8,1.2) * range(lapply(res,function(x) x$FM))
    cols <- rainbow(nms)
    ## historic
    i <- 1 ## historic pattern the same between mss
    llhist <- res[[i]]$FM[1,idxhist]
    ulhist <- res[[i]]$FM[3,idxhist]
    medhist <- res[[i]]$FM[2,idxhist]
    fmsy <- resMSE[[1]][[1]]$refs$Fmsy ## median(resMSE[[1]][[1]]$refdist$Fmsy) ##
##    tmp <- quantile(resMSE[[1]][[1]]$refdist$Fmsy,probs = c(0.025,0.975))
    llfmsy <- NA ## tmp[1]
    ulfmsy <- NA ##tmp[2]
    ## plot
    plot(xhist, medhist,
         xlim = c(0, dat$ny+set$nysim),
         ylim = ylim,
         ty='n', lwd=2,
         ylab="Fishing mortality", xlab = "")
    polygon(x = c((dat$ny-(set$nyhist-1)):dat$ny,rev((dat$ny-(set$nyhist-1)):dat$ny)),
            y = c(rep(ylim[1]-10,set$nyhist),rep(1.5*ylim[2],set$nyhist)),
            border=NA, col="grey95")
    polygon(x = c(xhist,rev(xhist)), y = c(llhist,rev(ulhist)),
            border=NA, col=rgb(t(col2rgb("grey30"))/255,alpha=0.2))
    lines(xhist, medhist, lwd=2)
    polygon(x = c(-2,rep(1.2*max(xall),2),-2), y = c(rep(llfmsy,2),rep(ulfmsy,2)),
            border=NA, col=rgb(t(col2rgb("grey40"))/255,alpha=0.2))
    abline(h=fmsy,lty=2)
    abline(h=0,lty=2)
    ## projection
    if(uncert){
        for(i in 1:nms){
            llsim <- res[[i]]$FM[1,idxsim]
            ulsim <- res[[i]]$FM[3,idxsim]
            polygon(x = c(xsim,rev(xsim)), y = c(llsim,rev(ulsim)),
                    border=NA, col=rgb(t(col2rgb(cols[i]))/255,alpha=0.2))
        }
    }
    for(i in 1:nms){
        medsim <- res[[i]]$FM[2,idxsim]
        if(med) lines(xsim, medsim, lwd=2, col=cols[i])
        if(is.numeric(trendline)){
            for(j in 1:length(trendline))
                lines(xsim, resMSE[[i]][[trendline[j]]]$FM[idxsim], col=cols[i])
        }else if(trendline)
            lines(xsim, resMSE[[i]][[1]]$FM[idxsim], col=cols[i])
    }
    abline(v=dat$ny, col="grey60",lwd=2)
    abline(v=max(which(dat$Fvals==0)), col="grey60",lwd=2,lty=2)
    legend("topleft", legend=set$hcr,
           col=cols, bty="n", lwd=2,lty=1)
    title("Fishing mortality")
    box()
}



#' @name plotmse.tradeoff
#' @export
plotmse.tradeoff <- function(mets, metrics = c("avRelCatch","PBBlim"), hcrs=NA){
    if(any(!is.na(hcrs))){
        metsnew <- vector("list",length(hcrs))
        for(i in 1:length(hcrs)){
            hcri <- hcrs[i]
            mets[[i]] <- mets[[hcri]]
        }
        mets <- metsnew
    }
    ## vars
    nms <- length(mets)
    allmets <- rownames(mets[[1]])
    xvar <- metrics[2]
    yvar <- metrics[1]
    idx <- which(xvar == allmets)
    idy <- which(yvar == allmets)
    xlim <- c(0.8,1.2) * range(lapply(mets, function(x) x[allmets == xvar]),na.rm=TRUE)
    ylim <- c(0.8,1.2) * range(lapply(mets, function(x) x[allmets == yvar]),na.rm=TRUE)
    cols <- rainbow(nms)
    ## metric specific settings
    if(length(c(grep("converged",xvar))>1)) xlim <- c(0,100)
    if(length(c(grep("converged",yvar))>1)) ylim <- c(0,100)
    ## plot
    plot(xlim, ylim,
         xlim = xlim,
         ylim = ylim,
         ty='n', lwd=2,
         xlab=xvar,
         ylab=yvar)
    ## metric specific lines
    if(length(c(grep("RelCatch",xvar),grep("BBmsy",xvar))>1)) abline(v=1,lty=2)
    if(length(c(grep("RelCatch",yvar),grep("BBmsy",xvar))>1)) abline(h=1,lty=2)
    if(length(grep("PBBlim",xvar)==1)) abline(v=0.05,lty=2)
    if(length(grep("PBBlim",yvar)==1)) abline(h=0.05,lty=2)
    ## MS
    for(i in 1:nms){
        xs <- mets[[i]][idx,]
        ys <- mets[[i]][idy,]
        segments(x0 = xs[1], y0 = ys[2], x1 = xs[3], y1 = ys[2], col=cols[i], lwd=1.5)
        segments(x0 = xs[2], y0 = ys[1], x1 = xs[2], y1 = ys[3], col=cols[i], lwd=1.5)
        points(xs[2],ys[2], col=cols[i], pch=16, cex=1.5)
    }
    legend("topright", legend=set$hcr,
           col=cols, bty="n", lwd=2,lty=1)
    box()
}


#' @name plotmse.quant
#' @export
plotmse.quant <- function(dat, set, resMSE, hcrs=NULL,
                          quants = c("bpbmsy.sd","fmfmsy.sd","cp.sd"),
                          hline=NULL){

    if(is.null(hcrs)){
        hcr <- set$hcr
        spicthcr <- which("spict" == unlist(lapply(strsplit(hcr,"-"), "[[", 1)))
    }else spicthcr <- hcrs



    nspict <- length(spicthcr)
    nquants <- length(quants)
    quant1 <- quant2 <- quant3 <- quant4 <- list()
    for(i in 1:nspict){
        ind <- spicthcr[i]
        quant1[[i]] <- apply(do.call(rbind,lapply(resMSE[[ind]], function(x) x$tacs[[quants[1]]])), 2, mean, na.rm=TRUE)
        if(nquants > 1) quant2[[i]] <- apply(do.call(rbind,lapply(resMSE[[ind]], function(x) x$tacs[[quants[2]]])), 2, mean, na.rm=TRUE)
        if(nquants > 2) quant3[[i]] <- apply(do.call(rbind,lapply(resMSE[[ind]], function(x) x$tacs[[quants[3]]])), 2, mean, na.rm=TRUE)
        if(nquants > 3) quant4[[i]] <- apply(do.call(rbind,lapply(resMSE[[ind]], function(x) x$tacs[[quants[4]]])), 2, mean, na.rm=TRUE)
    }

    cols <- rainbow(nspict)
    years <- 1:set$nysim

    opar <- par(mfrow=c(length(quants),1),mar=c(0,5,2,2), oma = c(5,0,0,0))
    ## BpBmsy
    xaxt <- ifelse(nquants==1,"s","n")
    xlab <- ifelse(nquants==1,"Time","")
    plot(years, quant1[[1]], ty='n',
         ylim=range(0,unlist(quant1)),
         ylab="", xlab = xlab, xaxt=xaxt)
    abline(h=0,lty=2)
    if(!is.null(hline)) abline(h=hline)
    for(i in 1:nspict) lines(years, quant1[[i]], col = cols[i], lwd=1.5)
    mtext(quants[1],2,3)
    ## legend
    legend("topright", legend=set$hcr[spicthcr],
           col=cols, bty="n", lwd=2,lty=1)
    box()
    ## FmFmsy
    xaxt <- ifelse(nquants==2,"s","n")
    xlab <- ifelse(nquants==2,"Time","")
    if(nquants > 1){
        plot(years, quant2[[1]], ty='n',
             ylim=range(0,unlist(quant2)),
             ylab="", xlab = xlab, xaxt=xaxt)
    abline(h=0,lty=2)
    if(!is.null(hline)) abline(h=hline)
    for(i in 1:nspict) lines(years, quant2[[i]], col = cols[i], lwd=1.5)
    mtext(quants[2],2,3)
    box()
    }
    ## Cp
    xaxt <- ifelse(nquants==3,"s","n")
    xlab <- ifelse(nquants==3,"Time","")
    if(nquants > 2){
        plot(years, quant3[[1]], ty='n',
             ylim=range(0,unlist(quant3)),
             ylab="", xlab = xlab, xaxt=xaxt)
    abline(h=0,lty=2)
    if(!is.null(hline)) abline(h=hline)
    for(i in 1:nspict) lines(years, quant3[[i]], col = cols[i], lwd=1.5)
    mtext(quants[3],2,3)
    }

    box()

}
