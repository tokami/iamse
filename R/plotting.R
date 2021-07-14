

#' @name plotiamse.cw
#' @export
plotiamse.cw <- function(dat, set, resMSE,
                       trendline=TRUE, uncert = TRUE, med = TRUE,
                       hcrs=NA, ylim = NULL, plot.legend = TRUE,
                       ylab="Catch", xlab = ""){
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
    res <- summary.mse(resMSE)
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
    if(!is.null(dat$ref$MSY)) msy <- dat$ref$MSY
    llmsy <- NA##tmp[1]
    ulmsy <- NA##tmp[2]
    ## plot
    plot(xhist, medhist,
         xlim = c(0, dat$ny+set$nysim),
         ylim = ylim,
         ty='n', lwd=2,
         ylab=ylab, xlab = xlab)
    polygon(x = c((dat$ny-(dat$nyC-1)):dat$ny,rev((dat$ny-(dat$nyC-1)):dat$ny)),
            y = c(rep(ylim[1]-1e4,dat$nyC),rep(1.5*ylim[2],dat$nyC)),
            border=NA, col="grey95")
    polygon(x = c(xhist,rev(xhist)), y = c(llhist,rev(ulhist)),
            border=NA, col=rgb(t(col2rgb("grey30"))/255,alpha=0.2))
    lines(xhist, medhist, lwd=2)
    polygon(x = c(-2,rep(1.2*max(xall),2),-2), y = c(rep(llmsy,2),rep(ulmsy,2)),
            border=NA, col=rgb(t(col2rgb("grey40"))/255,alpha=0.2))
    if(!is.null(dat$ref$MSY)) abline(h=msy,lty=2)
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
                lines(xsim, apply(resMSE[[i]][[trendline[j]]]$CW,1,sum)[idxsim], col=cols[i])
        }else if(trendline)
            lines(xsim, apply(resMSE[[i]][[1]]$CW,1,sum)[idxsim], col=cols[i])
    }
    abline(v=dat$ny, col="grey60",lwd=2)
    abline(v=max(which(dat$FM==0)), col="grey60",lwd=2,lty=2)
    if(plot.legend) legend("topright", legend=set$hcr,
           col=cols, bty="n", lwd=2,lty=1)
    title("Catch")
    box()
}


#' @name plotiamse.b
#' @export
plotiamse.b <- function(dat, set, resMSE,
                        yearly = TRUE,
                      trendline=TRUE, uncert = TRUE, med = TRUE,
                      hcrs=NA, ylim = NULL, plot.legend = TRUE,
                      ylab="Total biomass", xlab = ""){
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
    res <- summary.mse(resMSE, yearly = yearly)

    ## vars
    ns <- dat$ns
    ny <- dat$ny
    nysim <- set$nysim
    nms <- length(resMSE)
    if(yearly){
        xhist <- seq(1,ny,1)
        xsim <- seq(ny,ny+nysim,1)
        xall <- 1:(ny+nysim)
    }else{
        xhist <- seq(0,ny-1/ns, 1/ns)
        xsim <- seq(ny-1/ns,ny+nysim-1/ns,1/ns)
        xall <- seq(0,ny+nysim-1/ns,1/ns)
    }
    idxhist <- seq(length(xhist))
    idxsim <- tail(idxhist,1) + c(0,seq(length(xsim)-1))
    seaFM <- dat$FM[dat$ny,] / max(dat$FM[dat$ny,])
    ## if all FM == 0
    if(any(is.na(seaFM))){
        seaFM <- 1
    }

    if(is.null(ylim)) ylim <- c(0.8,1.2) * range(lapply(res,function(x) x$TSBfinal))
    cols <- rainbow(nms)
    ## historical
    i <- 1 ## historic pattern the same between mss
    llhist <- res[[i]]$TSBfinal[1,idxhist]
    ulhist <- res[[i]]$TSBfinal[3,idxhist]
    medhist <- res[[i]]$TSBfinal[2,idxhist]

    ## reference
    if(!is.null(dat$ref$Bmsy)){
        bmsy <- dat$ref$Bmsy
        if(yearly){
            ## HERE:
            tmp <- aggregate(list(FM = fmsy * seaFM), by = list(years = dat$yvec), sum)
            bmsy <- c(tmp$FM, rep(tail(tmp$FM, 1), set$nysim))
        }else{
            ## HERE:
            bmsy <- c(bmsy * seaFM, rep(tail(fmsy * seaFM, dat$ns), set$nysim))
        }
    }
    llbmsy <- NA ##tmp[1]
    ulbmsy <- NA ##tmp[2]
    ## plot
    plot(xhist, medhist,
         xlim = c(0, dat$ny+set$nysim),
         ylim = ylim,
         ty='n', lwd=2,
         ylab=ylab, xlab = xlab)
    polygon(x = c((dat$ny-(dat$nyC-1)):dat$ny,rev((dat$ny-(dat$nyC-1)):dat$ny)),
            y = c(rep(ylim[1]-1e4,dat$nyC),rep(1.5*ylim[2],dat$nyC)),
            border=NA, col="grey95")
    polygon(x = c(xhist,rev(xhist)), y = c(llhist,rev(ulhist)),
            border=NA, col=rgb(t(col2rgb("grey30"))/255,alpha=0.2))
    lines(xhist, medhist, lwd=2)
    polygon(x = c(-2,rep(1.2*max(xall),2),-2), y = c(rep(llbmsy,2),rep(ulbmsy,2)),
            border=NA, col=rgb(t(col2rgb("grey40"))/255,alpha=0.2))
    if(!is.null(dat$ref$Bmsy)) abline(h=bmsy,lty=2)
    if(!is.null(dat$ref$B0)) abline(h=dat$ref$B0,lty=2)
    ## projection
    if(uncert){
        for(i in 1:nms){
            llsim <- res[[i]]$TSBfinal[1,idxsim]
            ulsim <- res[[i]]$TSBfinal[3,idxsim]
            polygon(x = c(xsim,rev(xsim)), y = c(llsim,rev(ulsim)),
                    border=NA, col=rgb(t(col2rgb(cols[i]))/255,alpha=0.2))
        }
    }
    for(i in 1:nms){
        medsim <- res[[i]]$TSBfinal[2,idxsim]
        if(med) lines(xsim, medsim, lwd=2, col=cols[i])
        if(is.numeric(trendline)){
            for(j in 1:length(trendline))
                lines(xsim, resMSE[[i]][[trendline[j]]]$TSBfinal[idxsim], col=cols[i])
        }else if(trendline)
            lines(xsim, resMSE[[i]][[1]]$TSBfinal[idxsim], col=cols[i])
    }
    abline(v=dat$ny, col="grey60",lwd=2)
    abline(v=max(which(dat$FM==0)), col="grey60",lwd=2,lty=2)
    if(plot.legend) legend("topright", legend=set$hcr,
                           col=cols, bty="n", lwd=2,lty=1)
    title("Biomass")
    box()
}


#' @name plotiamse.f
#' @export
plotiamse.f <- function(dat, set, resMSE,
                        yearly = TRUE,
                        trendline=TRUE, uncert = TRUE, med = TRUE,
                        hcrs=NA, ylim=NULL, plot.legend = TRUE,
                        ylab="Fishing mortality", xlab = ""){
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
    res <- summary.mse(resMSE, yearly = yearly)

    ## vars
    ns <- dat$ns
    ny <- dat$ny
    nysim <- set$nysim
    nms <- length(resMSE)
    if(yearly){
        xhist <- seq(1,ny,1)
        xsim <- seq(ny,ny+nysim,1)
        xall <- 1:(ny+nysim)
    }else{
        xhist <- seq(0,ny-1/ns, 1/ns)
        xsim <- seq(ny-1/ns,ny+nysim-1/ns,1/ns)
        xall <- seq(0,ny+nysim-1/ns,1/ns)
    }
    idxhist <- seq(length(xhist))
    idxsim <- tail(idxhist,1) + c(0,seq(length(xsim)-1))
    seaFM <- dat$FM[dat$ny,] / max(dat$FM[dat$ny,])
    ## if all FM == 0
    if(any(is.na(seaFM))){
        seaFM <- 1
    }

    if(is.null(ylim)) ylim <- c(0.8,1.2) * range(lapply(res,function(x) x$FM))
    cols <- rainbow(nms)
    ## historical
    i <- 1 ## historic pattern the same between mss
    llhist <- res[[i]]$FM[1,idxhist]
    ulhist <- res[[i]]$FM[3,idxhist]
    medhist <- res[[i]]$FM[2,idxhist]
    if(!is.null(dat$ref$Fmsy)){
        fmsy <- dat$ref$Fmsy
        if(yearly){
            tmp <- aggregate(list(FM = fmsy * seaFM), by = list(years = dat$yvec), sum)
            fmsy <- c(tmp$FM, rep(tail(tmp$FM, 1), set$nysim))
        }else{
            fmsy <- c(fmsy * seaFM, rep(tail(fmsy * seaFM, dat$ns), set$nysim))
        }
    }
    llfmsy <- NA ## tmp[1]
    ulfmsy <- NA ## tmp[2]
    ## plot
    plot(xhist, medhist,
         xlim = c(0, dat$ny+set$nysim),
         ylim = ylim,
         ty='n', lwd=2,
         ylab=ylab, xlab = xlab)
    ## data period
    polygon(x = c((dat$ny-(dat$nyC-1)):dat$ny,rev((dat$ny-(dat$nyC-1)):dat$ny)),
            y = c(rep(ylim[1]-10,dat$nyC),rep(1.5*ylim[2],dat$nyC)),
            border=NA, col="grey95")
    ## historical
    polygon(x = c(xhist,rev(xhist)), y = c(llhist, rev(ulhist)),
            border=NA, col=rgb(t(col2rgb("grey30"))/255,alpha=0.2))
    lines(xhist, medhist, lwd=2)
    ## reference level
    if(!is.null(dat$ref$Fmsy)) {
        ## CHECK: include seaFM
        ## polygon(x = c(-2,rep(1.2*max(xall),2),-2),
        ##         y = c(rep(llfmsy,2),rep(ulfmsy,2)),
        ##         border=NA, col=rgb(t(col2rgb("grey40"))/255,alpha=0.2))
        lines(xall, fmsy, lty=2)
    }
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
                lines(xsim, apply(resMSE[[i]][[trendline[j]]]$FM,1,sum)[idxsim], col=cols[i])
        }else if(trendline)
            lines(xsim, apply(resMSE[[i]][[1]]$FM,1,sum)[idxsim], col=cols[i])
    }
    abline(v=dat$ny, col="grey60",lwd=2)
    abline(v=max(which(dat$FM==0)), col="grey60",lwd=2,lty=2)
    if(plot.legend) legend("topright", legend=set$hcr,
           col=cols, bty="n", lwd=2,lty=1)
    title("Fishing mortality")
    box()

}



#' @name plotiamse.tradeoff
#' @export
plotiamse.tradeoff <- function(mets, metrics = c("avRelCatch","PBBlim"),
                             hcrs=NA, plot.legend = TRUE){
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
    if(plot.legend) legend("topright", legend=set$hcr,
           col=cols, bty="n", lwd=2,lty=1)
    box()
}


#' @name plotiamse.quant
#' @export
plotiamse.quant <- function(dat, set, resMSE, hcrs=NULL,
                          quants = c("bpbmsy.sd","fmfmsy.sd","cp.sd"),
                          hline=NULL, plot.legend = TRUE){

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


#' @name plotiamse.prod
#' @export
plotiamse.prod <- function(prod, cols = NULL){
    meds <- prod$meds
    means <- prod$means
    lo <- prod$lo
    up <- prod$up
    alltv <- length(prod$meds)
    if(is.null(cols))
        cols <- rep(c("darkred","dodgerblue","darkgreen","darkorange","purple","gray","black","goldenrod"),100)
    plot(meds[[1]]$TSB, meds[[1]]$SP, ty = 'n',
         ylim = range(sapply(meds,function(x) x$SP),
                      sapply(lo,function(x) x$SP),
                      sapply(up,function(x) x$SP), na.rm =TRUE),
         xlim = range(sapply(meds,function(x) x$TSB),
                      sapply(lo,function(x) x$TSB),
                      sapply(up,function(x) x$TSB), na.rm =TRUE),
         xlab = "TSB", ylab = "SP")
    for(i in 1:alltv){
        if(i <= 3){
            polygon(c(lo[[i]]$TSB, rev(up[[i]]$TSB)), c(lo[[i]]$SP, rev(up[[i]]$SP)), border = NA,
                    col = rgb(t(col2rgb(cols[i])/255), alpha = 0.2))
        }
        lines(meds[[i]]$TSB, meds[[i]]$SP, col=cols[i], lwd=2)
    }
}
