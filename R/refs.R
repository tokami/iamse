## estimate reference levels for OM



#' @name simpopR
#' @export
simpopR <- function(FM, dat, set){

    ny <- set$refYears
    ns <- dat$nseasons
    nt <- ny * ns
    ## errors
    eF <- exp(rnorm(ny, 0, set$sigmaF) - set$sigmaF^2/2)
    eR <- genDevs(ny, set$sigmaR, set$rhoR)
    eM <- exp(rnorm(ny, 0, set$sigmaM) - set$sigmaM^2/2)
    eH <- exp(rnorm(ny, 0, set$sigmaH) - set$sigmaH^2/2)
    eMat <- exp(rnorm(ny, 0, set$sigmaMat) - set$sigmaMat^2/2)
    eR0 <- exp(rnorm(ny, 0, set$sigmaR0) - set$sigmaR0^2/2)
    eImp <- exp(rnorm(ny, 0, set$sigmaImp) - set$sigmaImp^2/2)
    ## some variables
    amax <- dat$amax + 1
    weight <- dat$weight
    weightF <- dat$weightF
    R0 <- dat$R0
    pzbm <- dat$pzbm
    sel <- dat$sel
    mat <- dat$mat
    h <- dat$h
    Ms <- dat$M / ns
    Fs <- FM * sel / ns
    initF <- dat$initF
    initN <- dat$initN
    ## containers
    NAA <- rep(0, amax)
    TSB <- TSB1plus <- SSB <- CW <- matrix(0, nrow=ny, ncol=ns)
    ## initialize starting values
    NAA[1] <- exp(initN[1]) * R0
    for(a in 2:amax)
        NAA[a] <- NAA[a-1] * exp(-Ms[a-1] + (initF*sel[a-1])/ns) * exp(initN[a])
    ## loop
    for(y in 1:ny){
        Fsy <- Fs * eF[y]
        Msy <- Ms * eM[y]
        Z <- Fsy + Msy
        maty <- mat * eMat[y]
        hy <- h * eH[y]
        R0y <- R0 * eR0[y]
        ## recruitment
        SSBtemp <- sum(NAA * weight * maty * exp(-pzbm * Z)) ## pre-recruitment mort
        SSBPR0 <- getSSBPR0(Msy, maty, 1, amax)
        rec <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSBtemp,
                       R0 = R0y, method = dat$SR)
        rec[rec<0] <- 1e-10
        NAA[1] <- rec * eR[y]
        ## seasons
        for(s in 1:ns){
            ## can't take more than what's there
            NAAmid <- NAA * exp(-Msy/2)
            Btemp <- sum(NAAmid * weight * sel * exp(-Msy/2))
            CAA <- baranov(Fsy, Msy, NAA)
            CW[y,s] <- sum(CAA * weight)
            ## maybe for refs too high F should overexploit the stock
            ## if(CW[y,s] > 0.99 * Btemp){
            ##     Fsy <- sel * min(set$maxF/ns,
            ##                            getFM(0.75 * Btemp, NAA = NAA[,y,s],
            ##                                  M = Msy, weight = weightF, sel = sel))
            ##     Z <- Msy + Fsy
            ##     CAA[,y,s] <- baranov(Fsy, Msy, NAA[,y,s])
            ##     CW[y,s] <- sum(weightF * CAA[,y,s])
            ## }
            ## TSB (use NAAmid?)
            TSB[y,s] <- sum(NAA * weight)
            TSB1plus[y,s] <- sum(NAA[-1] * weight[-1])
            ## SSB
            SSB[y,s] <- sum(NAA * weight * maty * exp(-pzbm * Z))
            NAA <- Ntemp <- NAA * exp(-Z)
            if(s == ns){
                NAA[amax] <- Ntemp[amax] + Ntemp[amax-1]
                for(a in 2:(amax-1)) NAA[a] <- Ntemp[a-1]
            }
        }
    }

    ## surplus production (for reflevs)
    TSBy <- apply(TSB, 1, mean)
    TSB1plusy <- apply(TSB1plus, 1, mean)
    SSBy <- apply(SSB, 1, mean)
    CWy <- apply(CW, 1, sum)
    SP <- rep(NA, ny)
    for(y in 1:(ny-1)){
        SP[y] <- TSBy[y+1] - TSBy[y] + CWy[y]
    }

    ## return
    return(list(CW=CWy,
                TSB=TSBy,
                TSB1plus=TSB1plusy,
                SP=SP,
                SSB=SSBy))
}


#' @name estRef
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @export
estRef <- function(dat, set, fvec = seq(0,5,0.1), ncores=parallel::detectCores()-1){

    ny <- dat$ny
    ns <- dat$nseasons
    nt <- ny * ns

    ## Remove variability
    set$sigmaF <- 0
    set$sigmaR <- 0
    set$rhoR <- 0
    set$sigmaM <- 0
    set$sigmaH <- 0
    set$sigmaMat <- 0
    set$sigmaR0 <- 0
    set$sigmaImp <- 0


    res0 <- parallel::mclapply(as.list(fvec),
                              function(x){
        with(simpopR(x, dat, set),
             c(x,
               tail(TSB,1),
               head(tail(SP,2),1),
               tail(CW,1)))
                              },
        mc.cores = ncores)

    ## refs
    res <- do.call(cbind, res0)
    msy <- res[4, which.max(res[3,])]
    fmsy <- res[1, which.max(res[3,])]
    bmsy <- res[2, which.max(res[3,])]
    binf <- res[2, which.max(res[2,])]

    ## B0
    unfished <- simpopR(0, dat, set)
    b0 <- tail(unfished$TSB1plus,1) ## CHECK: why?

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
