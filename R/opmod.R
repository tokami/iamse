
## input data:
## Number of years (ny)
## Maximum age of species (amax)
## Natural mortality (M)
## weight-at-age start (weightS)
## weight-at-age mid (weightH)
## maturity-at-age (mat)
## proportion of Z before maturation (pzbm)

#' @name initPop
#' @export
initPop <- function(specdat, set, reflev = TRUE){
    ## species data
    amax <- specdat$amax + 1  ## age 0
    M <- specdat$M
    weight <- specdat$weight
    weightF <- specdat$weightF
    mat <- specdat$mat
    pzbm <- specdat$pzbm
    initN <- specdat$initN
    recdev <- specdat$recdev
    SR <- specdat$SR
    logR0 <- specdat$logR0
    R0 <- exp(logR0)
    h <- specdat$h
    SSBPR0 <- specdat$SSBPR0
    rho <- specdat$rho
    Fvals <- specdat$Fvals
    sel <- specdat$sel
    initF <- specdat$initF
    ## indices
    ny <- specdat$ny
    ## errors
    eF <- set$eF
    eR <- set$eR
    eM <- set$eM
    eH <- set$eH
    eR0 <- set$eR0
    eMat <- set$eMat
    if(is.null(eF)) eF <- rnorm(ny, 0, set$sigmaF) - set$sigmaF^2/2
    if(is.null(eR)) {
        eR <- rnorm(ny, 0, set$sigmaR) - set$sigmaR^2/2
        if(set$autocor){  ## autocorrelated recruitment deviations
            eRar <- rep(NA,ny)
            eRar[1] <- eR[1]
            for (y in 2:ny){
                eRar[y] <- eRar[y-1] * rho + sqrt(1 - rho ^ 2) * eR[y]
            }
            eR <- eRar
        }
    }
    if(is.null(eM)) eM <- rnorm(ny, 0, set$sigmaM) - set$sigmaM^2/2
    if(is.null(eH)) eH <- rnorm(ny, 0, set$sigmaH) - set$sigmaH^2/2
    if(is.null(eR0)) eR0 <- rnorm(1, 0, set$sigmaR0) - set$sigmaR0^2/2
    if(is.null(eMat)) eMat <- rnorm(ny, 0, set$sigmaMat) - set$sigmaMat^2/2
    errs <- list(eF = eF,
                 eR = eR,
                 eM = eM,
                 eH = eH,
                 eR0 = eR0,
                 eMat = eMat)

    ## numbers-at-age
    NAA <- matrix(0, nrow=ny+1, ncol=amax)
    ## F-at-age
    FAA <- matrix(0, nrow=ny, ncol=amax)
    ## TSB
    TSB <- rep(0,ny)
    ## SSB
    SSB <- rep(0,ny)
    ## ESB
    ESB <- rep(0,ny)
    ## Z
    Z <- rep(0,amax)
    ## Predicted catch-at-age
    CAA <- matrix(0, nrow=ny, ncol=amax)
    ## Catch in weight
    CW <- rep(0,ny)
    ## TAC (for later)
    TACs <- rep(NA,ny)

    ## set up NAA
    NAA[1,1] <- R0 * exp(eR0) * exp(initN[1])
    for(a in 2:amax){
        NAA[1,a] <- NAA[1,a-1] * exp(-(M[a-1] + initF)) * exp(initN[a])
    }
    SSB0 <- sum(weight * mat * exp(-pzbm * Z) * NAA[1,])

    ## project forward
    for(y in 1:ny){
        ## FAA and Z
        FAA[y,] <- sel * Fvals[y] * exp(eF[y])
        M <- M * exp(eM[y])
        Z <- M + FAA[y,]
        ## CAA
        NAAmid <- NAA[y,] * exp(-M/2)
        CAA[y,] <- baranov(FAA[y,], M, NAAmid)
        ## CW
        CW[y] <- sum(weightF * CAA[y,])
        ## TSB
        TSB[y] <- sum(weight * NAA[y,])
        ## SSB
        maty <- mat * exp(eMat[y])
        SSB[y] <- sum(NAA[y,] * weight * maty * exp(-pzbm * Z))
        ## ESB
        ESB[y] <- sum(NAA[y,] * weightF * sel)
        ## remove Z
        Ntemp <- NAA[y,] * exp(-Z)
        ## update dynamics
        ##-----------------
        ## plus group
        NAA[y+1,amax] <- Ntemp[amax] + Ntemp[amax-1]
        ## other ages
        for(a in 2:(amax-1)) NAA[y+1,a] <- Ntemp[a-1]
        ## recruitment
        h <- h * exp(eH[y])
        SSBx <- c(SSB0,SSB)
        rec <- recfunc(h = h, SSBPR0 = SSBPR0, SSB = SSBx[y], R0 = R0, method = SR)
        rec[rec<0] <- 1e-10
        NAA[y+1,1] <- rec * exp(eR[y])
    }

    ## return
    out <- NULL
    out$NAA <- NAA
    out$TSB <- TSB
    out$SSB <- SSB
    out$ESB <- ESB
    out$CAA <- CAA
    out$FAA <- FAA
    out$FM <- apply(FAA, 1, function(x) mean(x / specdat$sel, na.rm=TRUE))
    out$CW <- CW
    out$TACs <- TACs
    out$errs <- errs
    return(out)
}



## advance population
#' @name advancePop
#' @export
advancePop <- function(specdat, hist, set, tacs){

    ##    if(nrow(hist$NAA) == 43) browser()

    ## parameters
    amax <- specdat$amax + 1  ## age 0
    pzbm <- specdat$pzbm
    SR <- specdat$SR
    logR0 <- specdat$logR0
    h <- specdat$h
    SSBPR0 <- specdat$SSBPR0
    rho <- specdat$rho
    ## parameters per age
    M <- specdat$M
    weight <- specdat$weight
    weightF <- specdat$weightF
    mat <- specdat$mat
    sel <- specdat$sel
    ## indices
    ny <- nrow(hist$FAA)
    y <- ny + 1
    yNAA <- y + 1  ## NAA one longer
    ysim <- y - specdat$ny
    ## errors
    eF <- set$eF[ysim]
    eR <- set$eR[ysim]
    eM <- set$eM[ysim]
    eH <- set$eH[ysim]
    eMat <- set$eMat[ysim]
    if(is.null(eR)) {
        eR <- rnorm(1, 0, set$sigmaR) - set$sigmaR^2/2
        if(set$autocor){  ## autocorrelated recruitment deviations
            eR <- set$eR[ysim-1] * rho + sqrt(1 - rho ^ 2) * eR
        }
    }
    if(is.null(eF)) eF <- rnorm(1, 0, set$sigmaF) - set$sigmaF^2/2
    if(is.null(eM)) eM <- rnorm(1, 0, set$sigmaM) - set$sigmaM^2/2
    if(is.null(eH)) eH <- rnorm(1, 0, set$sigmaH) - set$sigmaH^2/2
    if(is.null(eMat)) eMat <- rnorm(1, 0, set$sigmaMat) - set$sigmaMat^2/2
    if("errs" %in% names(hist)){
        errs <- list(eF = c(hist$errs$eF, eF),
                     eR = c(hist$errs$eR, eR),
                     eM = c(hist$errs$eM, eM),
                     eH = c(hist$errs$eH, eH),
                     eR0 = hist$errs$eR0,
                     eMat = c(hist$errs$eMat, eMat))
    }else{
        errs <- list(eF = eF,
                     eR = eR,
                     eM = eM,
                     eH = eH,
                     eMat = eMat)
    }

    ## numbers-at-age
    NAA <- hist$NAA
    ## F-at-age
    FAA <- hist$FAA
    ## TSB
    TSB <- hist$TSB
    ## SSB
    SSB <- hist$SSB
    ## ESB
    ESB <- hist$ESB
    ## Z
    Z <- hist$Z
    ## Predicted catch-at-age
    CAA <- hist$CAA
    ## Catch in weight
    CW <- hist$CW
    ## TAC
    TACs <- hist$TACs

    ## extend year dimension of matrices
    mattmp <- matrix(0, nrow=1, ncol=amax)
    NAA <- rbind(NAA, mattmp)
    FAA <- rbind(FAA, mattmp)
    TSB <- c(TSB,0)
    SSB <- c(SSB,0)
    ESB <- c(ESB,0)
    CAA <- rbind(CAA, mattmp)
    CW <- c(CW,0)
    TACs <- c(TACs,0)

    ## project forward
    ## Available mid year biomass
    M <- M * exp(eM)
    NAAmid <- NAA[yNAA-1,] * exp(-M/2)

    ## catch in weight = TAC
    if(tacs$id[nrow(tacs)] != "refFmsy" || is.na(tacs$id[nrow(tacs)])){
        TAC <- as.numeric(as.character(tacs$TAC[nrow(tacs)]))
        TACs[y] <- TAC
        FMtac <- min(set$maxF, getFM(TAC, NAA = NAAmid, M = M, weight = weightF, sel = sel))
        FAA[y,] <- sel * FMtac ## * exp(eF) ## that would be implementation error right?
        ## CAA[y,] <- baranov(FAA[y,], M, NAAmid) ## new: "To prevent
        ## a control from removing all exploitable biomass from the
        ## population in a year, we set the achieved catch to 75% of
        ## the midyear exploitable biomass in that year in cases where
        ## TAC exceeded the exploitable biomass."
        CAAtmp <- baranov(FAA[y,], M, NAAmid)
##        CAA[y,] <- CAAtmp
        CWreq <- sum(CAAtmp * weightF)
        CWpot <- sum(NAAmid * sel * weightF)
        ## print(paste0("Req:",CWreq))
        ## print(paste0("Pot:",CWpot))
        if(CWreq > CWpot){
            CAA[y,] <- 0.99 * NAAmid * sel
        }else{
            CAA[y,] <- CAAtmp
        }
    }else{
        FMtac <- hist$refs$Fmsy
        FAA[y,] <- sel * FMtac
        CAA[y,] <- baranov(FAA[y,], M, NAAmid)
        TAC <- sum(CAA[y,] * weightF, na.rm = TRUE)
        TACs[y] <- TAC
        tacs$TAC[nrow(tacs)] <- TAC
    }

    if(FALSE){
        ## CAA
        NAAtmp <- NAA[yNAA-1,] * weightF ## * sel
        CAAtmp <- TAC * (NAAtmp / sum(NAAtmp, na.rm=TRUE))## split TAC according to age distribution in pop
        ##    if(any(is.na(CAAtmp > NAAmid))) browser()
        CAAtmp[CAAtmp > NAAmid] <- NAAmid[CAAtmp > NAAmid]
        CAA[y,] <- CAAtmp / weightF

        ## FAA
        FAA[y,] <- CAA[y,] / (NAA[yNAA-1,] * exp(-M/2))
        FAA[y,] <- -log(1 - FAA[y,])
    }

    ## Z
    Z <- M + FAA[y,]

    ## NAA
    Ntmp <- NAA[yNAA-1,] * exp(-Z)
    ## plus group
    NAA[yNAA,amax] <- Ntmp[amax] + Ntmp[amax-1]
    for(a in 2:(amax-1)) NAA[yNAA,a] <- Ntmp[a-1]

    ## recruitment
    R0 <- exp(logR0)
    h <- h * exp(eH)
    rec <- recfunc(h = h, SSBPR0 = SSBPR0, SSB = SSB[y-1], R0 = R0, method = SR)
    rec[rec<0] <- 1e-10
    NAA[yNAA,1] <- rec * exp(eR)

    ## CW
    CW[y] <- sum(CAA[y,] * weightF, na.rm = TRUE)

    ## TSB
    TSB[y] <- sum(NAA[yNAA,] * weight, na.rm = TRUE)

    ## SSB
    maty <- mat * exp(eMat)
    SSB[y] <- sum(NAA[yNAA,] * weight * maty * exp(-pzbm * Z))

    ## ESB
    ESB[y] <- sum(NAA[yNAA,] * weightF * sel)

    ## return
    out <- NULL
    out$NAA <- NAA
    out$TSB <- TSB
    out$SSB <- SSB
    out$ESB <- ESB
    out$CAA <- CAA
    out$FAA <- FAA
    out$FM <- apply(FAA, 1, function(x) mean(x / specdat$sel, na.rm=TRUE))
    out$CW <- CW
    out$TACs <- TACs
    out$tacs <- tacs
    out$errs <- errs
    if("obs" %in% names(hist)) out$obs <- hist$obs
    if("refs" %in% names(hist)) out$refs <- hist$refs
    return(out)
}
