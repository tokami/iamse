#' @name est.ref.levels
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @export
#'
est.ref.levels <- function(dat, set=NULL, fvec = seq(0,5,0.1),
                   ncores=parallel::detectCores()-1,
                   ref = c("Fmsy","Bmsy","MSY","ESBmsy","SSBmsy","B0"),
                   plot = FALSE){

    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    asmax <- dat$asmax
    amax <- dat$amax

    if(is.null(set)) set <- check.set()
    ## Remove variability
    set$noiseF <- c(0,0,0)
    set$noiseR <- c(0,0,0)
    set$noiseM <- c(0,0,0)
    set$noiseH <- c(0,0,0)
    set$noiseW <- c(0,0,0)
    set$noiseR0 <- c(0,0,0)
    set$noiseMat <- c(0,0,0)
    set$noiseSel <- c(0,0,0)
    set$noiseImp <- c(0,0,0)
    errs <- vector("list", 7)
    errs$eF <- gen.noise(nyref, set$noiseF[1], set$noiseF[2], set$noiseF[3])
    errs$eR <- gen.noise(nyref, set$noiseR[1], set$noiseR[2], set$noiseR[3])
    errs$eM <- gen.noise(nyref, set$noiseM[1], set$noiseM[2], set$noiseM[3])
    errs$eH <- gen.noise(nyref, set$noiseH[1], set$noiseH[2], set$noiseH[3])
    errs$eW <- gen.noise(nyref, set$noiseW[1], set$noiseW[2], set$noiseW[3])
    errs$eR0 <- gen.noise(nyref, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
    errs$eMat <- gen.noise(nyref, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
    errs$eSel <- gen.noise(nyref, set$noiseSel[1], set$noiseSel[2], set$noiseSel[3])
    errs$eImp <- gen.noise(nyref, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
    setx <- c(set, errs)

    if(any(ref %in% c("Fmsy","Bmsy","MSY"))){
        res0 <- parallel::mclapply(as.list(fvec),
                                   function(x){
                                       with(simpop(log(x), dat, setx, out = 0),
                                            c(x,
                                              head(tail(TSB,2),1),
                                              head(tail(SP,2),1),
                                              head(tail(CW,2),1),
                                              head(tail(ESB,2),1),
                                              head(tail(SSB,2),1)
                                              )
                                            )
                                   },
                                   mc.cores = ncores)

        ## refs
        res <- do.call(cbind, res0)
        msy <- res[4, which.max(res[3,])]
        fmsy <- res[1, which.max(res[3,])]
        bmsy <- res[2, which.max(res[3,])]
        binf <- res[2, which.max(res[2,])]
        esbmsy <- res[5, which.max(res[3,])]
        ssbmsy <- res[6, which.max(res[3,])]

        if(plot){
            plot(fvec, res[3,], ty='b')
        }
    }


    ## B0
    if(any(ref == "B0")){
        unfished <- simpop(log(1e-20), dat, setx, out = 0)
        b0 <- tail(unfished$TSB,1)
    }


    refs <- list()
    if(any(ref == "Fmsy")) refs$Fmsy = fmsy
    if(any(ref == "MSY")) refs$MSY = msy
    if(any(ref == "Bmsy")) refs$Bmsy = bmsy
    if(any(ref == "ESBmsy")) refs$ESBmsy = esbmsy
    if(any(ref == "SSBmsy")) refs$SSBmsy = ssbmsy
    if(any(ref == "B0")) refs$B0 = b0

    dat$ref <- refs

    ## return
    return(dat)
}


#' @name est.ref.levels.stochastic
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @export
#'
est.ref.levels.stochastic <- function(dat, set=NULL, fmax = 10,
                        ncores = parallel::detectCores()-1,
                        ref = c("Fmsy","Bmsy","MSY","B0","ESBmsy","SSBmsy","ESB0"),
                        plot = FALSE, get.final = FALSE){

    ## Checks
    if(is.null(set)) set <- check.set()
    dist <- NULL
    if(!(set$refMethod %in% c("mean","median"))){
        stop("'set$refMethod' not known! Has to be 'mean' or 'median'!")
    }

    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    amax <- dat$amax + 1
    asmax <- amax * ns
    nyref <- set$refYears
    nrep <- set$refN
    nyrefmsy <- set$refYearsMSY
    tvflag <- FALSE

    ## Time-variant processes
    ## natural mortality
    mtv <- length(unique(as.numeric(dat$M)))  ## CHECK: if M is matrix
    ms <- unique(dat$M)
    mind <- match(dat$M, ms)
    if(length(dat$Msel) > 1){
        msel <- dat$Msel[!duplicated(dat$Msel)]
        mseltv <- length(msel)
    }else{
        msel <- dat$Msel[1]
        mseltv <- 1
    }
    if(mseltv > 1 && mseltv != mtv) stop("Both natural mortality over time (dat$M) and over age (dat$Msel) are time-variant, but do not have the same dimensions. This is not yet implemented, please let both vary equally or keep one of them constant.")
    alltv <- max(c(mtv, mseltv))
    ## selectivity
    if(length(dat$sel) > 1){
        sel <- dat$sel[!duplicated(dat$sel)]
        seltv <- length(sel)
    }else{
        sel <- dat$sel[1]
        seltv <- 1
    }
    if(seltv > 1 && alltv > 1 && seltv != alltv) stop("Both gear selectivity (dat$sel) and natural mortality (dat$M or dat$Msel) are time-variant, but do not have the same dimensions. This is not yet implemented, please let both vary equally or keep one of them constant.")
    alltv <- max(c(alltv,seltv))

    ##
    refall <- c("Fmsy","MSY","Bmsy","ESBmsy","SSBmsy","B0","ESB0")
    ##

    ## errors (have to be re-used for estimation of Bmsy)
    errs <- vector("list", nrep)
    for(i in 1:nrep){
        errs[[i]] <- list() ##vector("list", 9)
        errs[[i]]$eF <- gen.noise(nyref, set$noiseF[1], set$noiseF[2], set$noiseF[3])
##        errs[[i]]$eF <- rep(1.0, nyref)
        errs[[i]]$eR <- gen.noise(nyref, set$noiseR[1], set$noiseR[2], set$noiseR[3])
        errs[[i]]$eM <- gen.noise(nyref, set$noiseM[1], set$noiseM[2], set$noiseM[3])
        errs[[i]]$eH <- gen.noise(nyref, set$noiseH[1], set$noiseH[2], set$noiseH[3])
        errs[[i]]$eW <- gen.noise(nyref, set$noiseW[1], set$noiseW[2], set$noiseW[3])
        errs[[i]]$eR0 <- gen.noise(nyref, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
        errs[[i]]$eMat <- gen.noise(nyref, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
        errs[[i]]$eSel <- gen.noise(nyref, set$noiseSel[1], set$noiseSel[2], set$noiseSel[3])
        errs[[i]]$eImp <- gen.noise(nyref, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
    }
    ##
    datx <- dat

    ##
    ## datx$yvec <- rep(1:nyref, each = ns)
    ## datx$svec <- rep(1:ns, each = nyref)
    ## datx$s1vec <- seq(1, nyref * ns, ns)
    ## datx$as2a <- rep(1:amax, each = ns)
    ## datx$as2s <- rep(1:ns, amax)
    ## datx$inds <- seq(1,asmax,ns)

    ## ## HERE: debugging
    ## browser()

    ## x = 1
    ## setx <- c(set, errs[[x]])
    ## setx$refMethod <- "median"
    ## setx$refYearsMSY <- 10
    ## f <- 1
    ## str(datx,2)
    ## datx$weight <- dat$weight[,40]
    ## datx$weightF <- dat$weightF[,40]
    ## datx$M

    ## ## datx$h
    ## ## str(datx)

    ## f <- 0.7
    ## simpop(log(f), datx, setx, out=0)



    ## For now
    alltv <- 1


    if(any(ref %in% c("Fmsy","Bmsy","MSY","ESBmsy","SSBmsy"))){

        ## Fmsy
        res <- parallel::mclapply(as.list(1:nrep), function(x){
            setx <- c(set, errs[[x]])
            tmp <- rep(NA, alltv)
            for(i in 1:alltv){
                ## datx$M <- as.matrix(ms[mtv[i],])
                ## datx$Msel <- msels[mseltv[i]]
                ## datx$sel <- sels[seltv[i]]
                ## setx$tvm <- 1
                ## setx$tvmsel <- 1
                ## setx$tvsel <- 1
                ## datx$M <- rep(dat$M[i], nyref)
                ## ind <- (i-1)*ns+1
                ## datx$Ms <- rep(dat$Ms[ind:(ind+ns)], nyref)
                ## ind2 <- ifelse(ntv2 > 1, i, 1)
                ## datx$Msels <- msels[ind2]
                ## datx$Msel <- lapply(datx$Msels, rowMeans)
                ## TODO: for now overall average ref levels for time-varying pars
                if(get.final){
                    if(inherits(dat$weight,"matrix")) datx$weight <- dat$weight[,ncol(dat$weight)]
                    if(inherits(dat$weightF,"matrix")) datx$weightF <- dat$weightF[,ncol(dat$weightF)]
                    datx$M <- matrix(datx$M[nrow(datx$M),], 1, dat$ns)
                    datx$h <- tail(datx$h,1)
                    datx$R0 <- tail(datx$R0,1)
                }else{
                    if(inherits(dat$weight,"matrix")) datx$weight <- rowMeans(dat$weight)
                    if(inherits(dat$weightF,"matrix")) datx$weightF <- rowMeans(dat$weightF)
                    datx$M <- matrix(colMeans(datx$M), 1, dat$ns)
                    datx$h <- mean(datx$h)
                    datx$R0 <- mean(datx$R0)
                }
                opt <- optimise(function(x) unlist(simpop(x, datx, setx, out=1)),
                                log(c(0.001,fmax)), maximum = TRUE)
                ## opt <- optimise(function(x){
                ##     datx$FM <- matrix(exp(x)/datx$ns, dat$ny, dat$ns)
                ##     datx$M <- matrix(datx$M[1,], dat$ny, dat$ns)
                ##     initpop(datx, setx, out=5)
                ## },
                ##                 log(c(0.001,fmax)), maximum = TRUE)
                tmp[i] <- exp(opt$maximum)
            }
            return(tmp)
        }, mc.cores = ncores)
        fmsys <- do.call(rbind, res)


        ## MSY and Biomass reference points
        res <- parallel::mclapply(as.list(1:nrep), function(x){
            setx <- c(set, errs[[x]])
            tmp <- vector("list", alltv)
            for(i in 1:alltv){
                ## datx$M <- rep(dat$M[i], nyref)
                ## ind <- (i-1)*ns+1
                ## datx$Ms <- rep(dat$Ms[ind:(ind+ns)], nyref)
                ## ind2 <- ifelse(ntv2 > 1, i, 1)
                ## datx$Msels <- msels[ind2]
                ## datx$Msel <- lapply(datx$Msels, rowMeans)
                setx$tvm <- 1
                setx$tvmsel <- 1
                setx$tvsel <- 1
                ## TODO: for now overall average ref levels for time-varying pars
                if(get.final){
                    if(inherits(dat$weight,"matrix")) datx$weight <- dat$weight[,ncol(dat$weight)]
                    if(inherits(dat$weightF,"matrix")) datx$weightF <- dat$weightF[,ncol(dat$weightF)]
                    datx$M <- matrix(datx$M[nrow(datx$M),], 1, dat$ns)
                    datx$h <- tail(datx$h,1)
                    datx$R0 <- tail(datx$R0,1)
                }else{
                    if(inherits(dat$weight,"matrix")) datx$weight <- rowMeans(dat$weight)
                    if(inherits(dat$weightF,"matrix")) datx$weightF <- rowMeans(dat$weightF)
                    datx$M <- matrix(colMeans(datx$M), 1, dat$ns)
                    datx$h <- mean(datx$h)
                    datx$R0 <- mean(datx$R0)
                }
                tmp0 <- simpop(log(fmsys[x,i]), datx, setx, out=0)
                if(set$refMethod == "mean"){
                    tmp[[i]] <- c(mean(tail(tmp0$CW,nyrefmsy)), mean(tail(tmp0$TSB,nyrefmsy)),
                                  mean(tail(tmp0$ESB,nyrefmsy)), mean(tail(tmp0$SSB,nyrefmsy)))
                }else if(set$refMethod == "median"){
                    tmp[[i]] <- c(median(tail(tmp0$CW,nyrefmsy)), median(tail(tmp0$TSB,nyrefmsy)),
                                  median(tail(tmp0$ESB,nyrefmsy)), median(tail(tmp0$SSB,nyrefmsy)))
                }
            }
            return(tmp)
        }, mc.cores = ncores)

        ## sort in list by reference point with matrix (nrep, ntv)
        brefs <- vector("list", 4) ## c("MSY","Bmsy","ESBmsy","SSBmsy")
        for(i in 1:4){
            brefs[[i]] <- do.call(rbind,
                                  lapply(as.list(1:nrep),
                                         function(x) sapply(1:alltv, function(j) res[[x]][[j]][[i]])))
        }

    }

    ## B0
    if(any(ref %in% c("B0","ESB0"))){

        res <- parallel::mclapply(as.list(1:nrep), function(x){
            setx <- c(set, errs[[x]])
            tmp <- vector("list", alltv)
            for(i in 1:alltv){
                ## datx$M <- rep(dat$M[i], nyref)
                ## ind <- (i-1)*ns+1
                ## datx$Ms <- rep(dat$Ms[ind:(ind+ns)], nyref)
                ## ind2 <- ifelse(ntv2 > 1, i, 1)
                ## datx$Msels <- msels[ind2]
                ## datx$Msel <- lapply(datx$Msels, rowMeans)
                setx$tvm <- 1
                setx$tvmsel <- 1
                setx$tvsel <- 1
                ## TODO: for now overall average ref levels for time-varying pars
                if(get.final){
                    if(inherits(dat$weight,"matrix")) datx$weight <- dat$weight[,ncol(dat$weight)]
                    if(inherits(dat$weightF,"matrix")) datx$weightF <- dat$weightF[,ncol(dat$weightF)]
                    datx$M <- matrix(datx$M[nrow(datx$M),], 1, dat$ns)
                    datx$h <- tail(datx$h,1)
                    datx$R0 <- tail(datx$R0,1)
                }else{
                    if(inherits(dat$weight,"matrix")) datx$weight <- rowMeans(dat$weight)
                    if(inherits(dat$weightF,"matrix")) datx$weightF <- rowMeans(dat$weightF)
                    datx$M <- matrix(colMeans(datx$M), 1, dat$ns)
                    datx$h <- mean(datx$h)
                    datx$R0 <- mean(datx$R0)
                }
                tmp0 <- simpop(log(1e-20), datx, setx, out=0)
                if(set$refMethod == "mean"){
                    tmp[[i]] <- c(mean(tail(tmp0$TSB,nyrefmsy)),
                                mean(tail(tmp0$ESB,nyrefmsy)))
                }else if(set$refMethod == "median"){
                    tmp[[i]] <- c(median(tail(tmp0$TSB,nyrefmsy)),
                                median(tail(tmp0$ESB,nyrefmsy)))
                }
            }
            return(tmp)
        }, mc.cores = ncores)

        ## b0s <- do.call(rbind, res)

        ## sort in list by reference point with matrix (nrep, ntv)
        b0s <- vector("list", 2) ## b("B0","ESB0")
        for(i in 1:2){
            b0s[[i]] <- do.call(rbind,
                                  lapply(as.list(1:nrep),
                                         function(x) sapply(1:alltv, function(j) res[[x]][[j]][[i]])))
        }

    }


    ## all refs in one list
    refs <- c(list(fmsys), brefs, b0s)

    ## remove runs where long-term SP is smaller or equal to 0
    for(i in 1:alltv){
        ind <- which(brefs[[1]][,i] <= 0)
        if(length(ind) > 0){
            for(j in 1:length(refs)) refs[[j]][ind,i] <- NA ## refs[[j]][-ind,i]
        }
    }

    ## overall refs
    if(set$refMethod == "mean"){
        meds <- lapply(refs, function(x){
            tmp <- rep(NA, alltv)
            for(i in 1:alltv) tmp[i] <- mean(x[,i], na.rm = TRUE)
            return(tmp)
        })
    }else if(set$refMethod == "median"){
        meds <- lapply(refs, function(x){
            tmp <- rep(NA, alltv)
            for(i in 1:alltv) tmp[i] <- median(x[,i], na.rm = TRUE)
            return(tmp)
        })
    }

    ind <- which(refall %in% ref)
    refdist <- refs[ind]
    names(refdist) <- refall[ind]
    dat$refdist <- refdist

    refmed <- matrix(NA, ncol=length(ref), nrow=nt)
    for(i in 1:length(ind)){
        refmed[,i] <- meds[[ind[[i]]]][mind]
    }
    colnames(refmed) <- refall[ind]
    dat$ref <- as.data.frame(refmed)

    if(plot){
        if(alltv < 4){
            cols <- c("dodgerblue2","darkorange","darkgreen","purple")
            nr <- ceiling(length(refdist)/3)
            par(mfrow=c(nr,3))
            for(i in 1:length(refdist)){
                hist(refdist[[i]][,1], main = names(refdist)[i],
                     breaks=20, freq = TRUE, xlim = range(refdist[[i]],na.rm=TRUE),
                     xlab = "", col = rgb(t(col2rgb(cols[1]))/255,alpha=0.4))
                if(alltv > 1){
                    for(j in 2:alltv) hist(refdist[[i]][,j],
                                         breaks=20, freq = TRUE,
                                         add = TRUE, col = rgb(t(col2rgb(cols[j]))/255,alpha=0.4))
                }
                ## abline(v=mean(dist[,i]), lty=1, lwd=1.5, col=4)
                ## abline(v=median(dist[,i]), lty=2, lwd=1.5, col=4)
                ## if(i == 1) legend("topright", legend = c("mean","median"),
                ##                   col=4, lty=c(1,2),lwd=1.5)
            }
        }else{
            nr <- floor(length(refdist)/2)
            par(mfrow=c(nr,2))
            for(i in 1:length(refdist)){
                plot(refmed[,i], main = colnames(refmed)[i],
                     ty = 'l', lwd=1.5, xlab = "Time", ylab = colnames(refmed)[i])
            }
        }
    }

    ## return
    return(dat)
}


#' @name sbr
#' @export
sbr <- function(FM, dat, out=0){
    ## some variables
    amax <- dat$amax + 1
    ## initialize starting values
    NAA <- rep(0, amax)
    NAA[1] <- 1
    for(a in 2:amax)
        NAA[a] <- NAA[a-1] * exp(-(dat$M[a-1] + FM * dat$sel[a-1]))

    SSBPR <- sum(NAA * dat$weight * dat$mat * dat$fecun)
    YPR <- sum(baranov(FM*dat$sel, dat$M, NAA) * dat$weightFs)

    if(out == 0) res <- SSBPR else  res <- YPR

    return(res)
}
