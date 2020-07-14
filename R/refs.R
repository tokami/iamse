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
    R0 <- exp(dat$logR0)
    pzbm <- dat$pzbm
    sel <- dat$sel
    mat <- dat$mat
    h <- dat$h
    Ms <- dat$M / ns
    Fs <- FM * sel / ns
    ## initiate matrices & vectors
    NAA <- CAA <- FAA <- array(0, c(amax, ny, ns),
                                dimnames = list(age = 0:(amax-1),
                                                years = 1:ny,
                                                seasons = 1:ns))
    ## containers
    TSB <- TSB1plus <- SSB <- ESB <- CW <- matrix(0, nrow=ny, ncol=ns)
    ## initialize starting values
    NAtmp <- NAtmp2 <- rep(R0, amax)
    NAA[,1,1] <- NAtmp
    ## loop
    for(y in 1:ny){
        Fsy <- Fs * eM[y]
        Msy <- Ms * eM[y]
        Z <- Fsy + Msy
        maty <- mat * eMat[y]
        hy <- h * eH[y]
        R0y <- R0 * eR0[y]
        ## recruitment
        SSBtemp <- sum(NAA[,y,1] * weight * maty * exp(-pzbm * Z)) ## pre-recruitment mort
        SSBPR0 <- getSSBPR0(Msy, maty, 1, amax)
        rec <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSBtemp,
                       R0 = R0y, method = dat$SR)
        rec[rec<0] <- 1e-10
        NAA[1,y,1] <- rec * eR[y]
        ## seasons
        for(s in 1:ns){
            ## can't take more than what's there
            Btemp <- sum(NAA[,y,s] * weight * sel * exp(-Msy/2))
            Ctemp <- sum(baranov(Fsy, Msy, NAA[,y,s]))
            if(Ctemp > 0.99 * Btemp){
                Ctemp <- 0.75 * Btemp
                Fsy <- sel * min(set$maxF/ns,
                                       getFM(Ctemp, NAA = NAA[,y,s],
                                             M = Msy, weight = weightF, sel = sel))
                Z <- Msy + Fsy
            }
            ## CAA
            CAA[,y,s] <- baranov(Fsy, Msy, NAA[,y,s])
            ## CW
            CW[y,s] <- sum(weightF * CAA[,y,s])
            ## TSB
            NAAmid <- NAA[,y,s] * exp(-Msy/2)
            TSB[y,s] <- sum(weight * NAAmid)
            TSB1plus[y,s] <- sum(weight[-1] * NAAmid[-1])
            ## SSB
            SSB[y,s] <- sum(NAAmid * weight * maty * exp(-pzbm * Z))
            ## ESB
            ESB[y,s] <- sum(NAAmid * weightF * sel)
            Ntemp <- NAA[,y,s] * exp(-Z)
            if(s < ns){
                NAA[,y,s+1] <- Ntemp
            }else if(y < ny){
                NAA[amax, y+1, 1] <- Ntemp[amax] + Ntemp[amax-1]
                for(a in 2:(amax-1)) NAA[a,y+1,1] <- Ntemp[a-1]
            }
        }
    }

    ## surplus production (for reflevs)
    TSBy <- apply(TSB, 1, sum)
    TSB1plusy <- apply(TSB1plus, 1, sum)
    SSBy <- apply(SSB, 1, sum)
    ESBy <- apply(ESB, 1, sum)
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
                SSB=SSBy,
                ESB=ESBy))
}


#' @name estRL
#' @export
estRL <- function(dat, set, fvec = seq(0,5,0.1)){

    ## Remove variability
    set$sigmaF <- 0
    set$sigmaR <- 0
    set$rhoR <- 0
    set$sigmaM <- 0
    set$sigmaH <- 0
    set$sigmaMat <- 0
    set$sigmaR0 <- 0
    set$sigmaImp <- 0

    ## run over Fvec
    res <- sapply(fvec, function(x){
        with(simpopR(x, dat, set),
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
