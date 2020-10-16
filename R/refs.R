#' @name estRef
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @export
#'
estRef <- function(dat, set=NULL, fvec = seq(0,5,0.1),
                   ncores=parallel::detectCores()-1,
                   ref = c("Fmsy","Bmsy","MSY","ESBmsy","SSBmsy","B0"),
                   plot = FALSE){

    ny <- dat$ny
    ns <- dat$nseasons
    nt <- ny * ns

    if(is.null(set)) set <- checkSet()
    ## Remove variability
    set$noiseF <- c(0,0,0)
    set$noiseR <- c(0,0,0)
    set$noiseM <- c(0,0,0)
    set$noiseH <- c(0,0,0)
    set$noiseR0 <- c(0,0,0)
    set$noiseMat <- c(0,0,0)
    set$noiseImp <- c(0,0,0)
    errs <- vector("list", 7)
    errs$eF <- genNoise(nyref, set$noiseF[1], set$noiseF[2], set$noiseF[3])
    errs$eR <- genNoise(nyref, set$noiseR[1], set$noiseR[2], set$noiseR[3])
    errs$eM <- genNoise(nyref, set$noiseM[1], set$noiseM[2], set$noiseM[3])
    errs$eH <- genNoise(nyref, set$noiseH[1], set$noiseH[2], set$noiseH[3])
    errs$eR0 <- genNoise(nyref, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
    errs$eMat <- genNoise(nyref, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
    errs$eImp <- genNoise(nyref, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
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


#' @name estRefStoch
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @export
#'
estRefStoch <- function(dat, set=NULL,
                        ncores = parallel::detectCores()-1,
                        ref = c("Fmsy","Bmsy","MSY","B0","ESBmsy","SSBmsy"),
                        plot = FALSE){

    ny <- dat$ny
    ns <- dat$nseasons
    nt <- ny * ns
    nyref <- set$refYears
    nrep <- set$refN
    nyrefmsy <- set$refYearsMSY

    if(is.null(set)) set <- checkSet()
    dist <- NULL
    if(!(set$refMethod %in% c("mean","median"))){
        stop("'set$refMethod' not known! Has to be 'mean' or 'median'!")
    }

    ## errors (have to be re-used for estimation of Bmsy)
    errs <- vector("list", nrep)
    for(i in 1:nrep){
        errs[[i]] <- vector("list", 7)
        errs[[i]]$eF <- genNoise(nyref, set$noiseF[1], set$noiseF[2], set$noiseF[3])
        errs[[i]]$eR <- genNoise(nyref, set$noiseR[1], set$noiseR[2], set$noiseR[3])
        errs[[i]]$eM <- genNoise(nyref, set$noiseM[1], set$noiseM[2], set$noiseM[3])
        errs[[i]]$eH <- genNoise(nyref, set$noiseH[1], set$noiseH[2], set$noiseH[3])
        errs[[i]]$eR0 <- genNoise(nyref, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
        errs[[i]]$eMat <- genNoise(nyref, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
        errs[[i]]$eImp <- genNoise(nyref, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
    }


    if(any(ref %in% c("Fmsy","Bmsy","MSY","ESBmsy","SSBmsy"))){

        ## Fmsy
        res <- parallel::mclapply(as.list(1:nrep), function(x){
            setx <- c(set, errs[[x]])
            opt <- optimise(function(x) unlist(simpop(x, dat, setx, out=1)),
                            log(c(0.001,10)), maximum = TRUE)
            c(exp(opt$maximum))
        }, mc.cores = ncores)
        dist <- do.call(rbind, res)

        ## Other ref points
        res <- parallel::mclapply(as.list(1:nrep), function(x){
            setx <- c(set, errs[[x]])
            tmp <- simpop(log(dist[x,1]), dat, setx, out=0)
            if(set$refMethod == "mean"){
                c(mean(tail(tmp$CW,nyrefmsy)), mean(tail(tmp$TSB,nyrefmsy)),
                  mean(tail(tmp$ESB,nyrefmsy)), mean(tail(tmp$SSB,nyrefmsy)))
            }else if(set$refMethod == "median"){
                c(median(tail(tmp$CW,nyrefmsy)), median(tail(tmp$TSB,nyrefmsy)),
                  median(tail(tmp$ESB,nyrefmsy)), median(tail(tmp$SSB,nyrefmsy)))
            }
        }, mc.cores = ncores)
        dist <- cbind(dist, do.call(rbind, res))
        colnames(dist) <- c("Fmsy","MSY","Bmsy","ESBmsy","SSBmsy")

    }

    ## B0
    if(any(ref == "B0")){
        b0 <- sapply(1:nrep, function(x){
            setx <- c(set, errs[[x]])
            if(set$refMethod == "mean"){
                mean(tail(simpop(log(1e-20), dat, setx, out = 0)$TSB, nyrefmsy))
            }else if(set$refMethod == "median"){
                median(tail(simpop(log(1e-20), dat, setx, out = 0)$TSB, nyrefmsy))
            }
        })
        dist <- cbind(dist, b0)
        colnames(dist) <- c(colnames(dist)[-ncol(dist)], "B0")
    }

    ## remove runs where long-term SP is smaller or equal to 0
    nami <- colnames(dist)
    if("MSY" %in% nami){
        dist2 <- dist[dist[,2] > 0,]
    }else{
        dist2 <- dist
    }
    if(set$refMethod == "mean"){
        meds <- apply(dist2, 2, mean)
    }else if(set$refMethod == "median"){
        meds <- apply(dist2, 2, median)
    }

    refs <- list()
    if(any(ref == "Fmsy")) refs$Fmsy = as.numeric(meds[which(nami == "Fmsy")])
    if(any(ref == "MSY")) refs$MSY = as.numeric(meds[which(nami == "MSY")])
    if(any(ref == "Bmsy")) refs$Bmsy = as.numeric(meds[which(nami == "Bmsy")])
    if(any(ref == "ESBmsy")) refs$ESBmsy = as.numeric(meds[which(nami == "ESBmsy")])
    if(any(ref == "SSBmsy")) refs$SSBmsy = as.numeric(meds[which(nami == "SSBmsy")])
    if(any(ref == "B0")) refs$B0 = as.numeric(meds[which(nami == "B0")])

    dat$ref <- refs
    dat$refdist <- dist

    if(plot){
        nr <- floor(ncol(dist)/2)
        par(mfrow=c(nr,2))
        for(i in 1:ncol(dist)){
            hist(dist[,i], main = colnames(dist)[i], breaks=50)
            abline(v=mean(dist[,i]), lty=1, lwd=1.5, col=4)
            abline(v=median(dist[,i]), lty=2, lwd=1.5, col=4)
            if(i == 1) legend("topright", legend = c("mean","median"),
                              col=4, lty=c(1,2),lwd=1.5)
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
