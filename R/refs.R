## estimate reference levels for OM



#' @name simpopR
#' @export
simpopR <- function(FM, dat, set){
    years <- set$refYears
    ## errors
    eF <- rnorm(years, 0, set$sigmaF) - set$sigmaF^2/2
    eR <- genDevs(years, set$sigmaR, set$rhoR)
    eM <- rnorm(years, 0, set$sigmaM) - set$sigmaM^2/2
    eH <- rnorm(years, 0, set$sigmaH) - set$sigmaH^2/2
    eMat <- rnorm(years, 0, set$sigmaMat) - set$sigmaMat^2/2
    eR0 <- rnorm(years, 0, set$sigmaR0) - set$sigmaR0^2/2
    ## some variables
    amax <- dat$amax + 1
    weight <- dat$weight
    R0 <- exp(dat$logR0)
    ## initiate matrices & vectors
    Nage <- CAA <- FAA <- matrix(nrow=amax,ncol=years+1)
    CW <- SP <- rep(NA, years)
    ## initialize starting values
    NAtmp <- NAtmp2 <- rep(R0, amax)
    Nage[,1] <- NAtmp
    ## loop
    for(y in 1:years){
        M <- dat$M * exp(eM[y])
        h <- dat$h * exp(eH[y])
        mat <- dat$mat * exp(eMat[y])
        Nage[,y] <- NAtmp
##        NAAmid <- Nage[,y] * exp(-M/2)
        FAA[,y] <- dat$sel * FM * exp(eF[y])
        CAA[,y] <- baranov(FAA[,y], M, Nage[,y]) ## FAA / Z * Nage[,y] * (1 - exp(-Z))
        CW[y] <- sum(CAA[,y] * dat$weight)
        SSBPR0 <- getSSBPR0(M, mat, fecun = 1, amax=amax)
        SSB <- sum(NAtmp * mat * dat$weight)
        R0 <- exp(dat$logR0) * exp(eR0[y])
        NAtmp2[1] <- recfunc(h = h, SSB = SSB,
                             SSBPR0 = SSBPR0, R0 = R0, method = dat$SR) * eR[y]
        Z <- M + FAA[,y]
        survival <- NAtmp * exp(-Z)
        NAtmp2[2:amax] <- survival[1:(amax-1)]
        NAtmp2[amax] <- NAtmp2[amax] + survival[amax]
        NAtmp <- NAtmp2
    }
    Bage <- apply(Nage, 2, function(x) x * dat$weight)
    SSBage <- apply(Nage, 2, function(x) x * dat$weight * dat$mat)
    ESBage <- apply(Nage, 2, function(x) x * dat$weight * dat$sel)

    ## surplus production (for reflevs)
    for(y in 1:(years-1)){
        SP[y] <- sum(Bage[,y+1]) - sum(Bage[,y]) + CW[y]
    }

    ## return
    return(list(CW=CW,
                TSB=apply(Bage[,1:years],2,sum,na.rm=TRUE),
                TSB1plus=apply(Bage[-1,1:years],2,sum,na.rm=TRUE),
                SP=SP,
                SSB=apply(SSBage[,1:years],2,sum,na.rm=TRUE),
                ESB=apply(ESBage[,1:years],2,sum,na.rm=TRUE)))
}


#' @name estRL
#' @export
estRL <- function(dat, set, fvec = seq(0,5,0.1)){

    ## Variables
    years <- set$refYears

    ## Remove variability
    set$sigmaF <- 0
    set$sigmaR <- 0
    set$sigmaM <- 0
    set$sigmaH <- 0
    set$sigmaMat <- 0
    set$sigmaR0 <- 0

    ## run over Fvec
    res <- sapply(fvec, function(x){
        with(simpop(x, dat, set),
             c(x, tail(TSB,1),
               head(tail(SP,2),1),
               tail(CW,1)))
    })

    ## refs
    msy <- res[4, which.max(res[3,])]
    fmsy <- fvec[which.max(res[3,])]
    bmsy <- res[2, which.max(res[3,])]
    binf <- res[2, which.max(res[2,])]

    ## B0
    unfished <- simpop(0, dat, set)
##    b0 <- tail(unfished$TSB1plus,1) ## CHECK: why?
    b0 <- tail(unfished$TSB,1)

    ## return
    return(list(refs = list(
                    Fmsy = fmsy,
                    Bmsy = bmsy,
                    MSY = msy,
                    Binf = binf,
                    B0 = b0),
        xvec = fvec,
        yvec = res[3,]
    ))
}
