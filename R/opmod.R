
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
initPop <- function(specdat, set = NULL, refs = NULL, out.opt = 1, depl.quant = "B0"){
    ## indices
    ny <- specdat$ny
    ns <- specdat$nseason
    nt <- ny * ns
    ## species data
    amax <- specdat$amax + 1  ## age 0
    M <- specdat$M
    Ms <- M/ns
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
    ## errors
    if(!is.null(set)){
        eF <- set$eF
        eR <- set$eR
        eM <- set$eM
        eH <- set$eH
        eR0 <- set$eR0
        eMat <- set$eMat
        if(is.null(eF)) eF <- exp(rnorm(ny, 0, set$sigmaF) - set$sigmaF^2/2)
        if(is.null(eR)) {
            eR <- genDevs(ny, set$sigmaR, set$rhoR)
        }
        if(is.null(eM)) eM <- exp(rnorm(ny, 0, set$sigmaM) - set$sigmaM^2/2)
        if(is.null(eH)) eH <- exp(rnorm(ny, 0, set$sigmaH) - set$sigmaH^2/2)
        if(is.null(eR0)) eR0 <- exp(rnorm(1, 0, set$sigmaR0) - set$sigmaR0^2/2)
        if(is.null(eMat)) eMat <- exp(rnorm(ny, 0, set$sigmaMat) - set$sigmaMat^2/2)
    }else{
        eF <- exp(rnorm(ny, 0, 0))
        eR <- genDevs(ny, 0, 0)
        eM <- exp(rnorm(ny, 0, 0))
        eH <- exp(rnorm(ny, 0, 0))
        eR0 <- exp(rnorm(1, 0, 0))
        eMat <- exp(rnorm(ny, 0, 0))
    }
    errs <- list(eF = eF,
                 eR = eR,
                 eM = eM,
                 eH = eH,
                 eR0 = eR0,
                 eMat = eMat)
    ## numbers-at-age
    NAA <- array(0, c(amax, ny+1, ns),
                 dimnames = list(age = 0:(amax-1),
                                 years = 0:ny,
                                 seasons = 1:ns))
    ## F-at-age
    FAA <- array(0, c(amax, ny, ns),
                 dimnames = list(age = 0:(amax-1),
                                 years = 1:ny,
                                 seasons = 1:ns))
    ## Predicted catch-at-age
    CAA <- array(0, c(amax, ny, ns),
                 dimnames = list(age = 0:(amax-1),
                                 years = 1:ny,
                                 seasons = 1:ns))
    ## containers
    TSB <- SSB <- ESB <- CW <- matrix(0, nrow=ny, ncol=ns)
    ## TAC (for later)
    TACs <- rep(NA, ny)
    ## burnin period
    if(is.null(set)) burnin <- 20 else burnin <- set$burnin
    if(is.numeric(burnin) && burnin > 0){
        NAAbi <- array(0, c(amax, burnin, ns),
                       dimnames = list(age = 0:(amax-1),
                                       years = 1:burnin,
                                       seasons = 1:ns))
        CWbi <- matrix(0, nrow=ny, ncol=ns)
        NAAbi[1,1,1] <- R0 * exp(initN[1])
        Fbi <- (sel * Fvals[1]) / ns
        Zbi <- Ms + Fbi
        for(a in 2:amax) NAAbi[a,1,1] <- NAAbi[a-1,1,1] * exp(-(Zbi[a-1])) * exp(initN[a])
        for(y in 1:burnin){
            ## recruitment
            SSBtemp <- sum(NAAbi[,y,1] * weight * mat * exp(-pzbm * Zbi)) ## pre-recruitment mort
            rec <- recfunc(h = h, SSBPR0 = SSBPR0, SSB = SSBtemp,
                           R0 = R0, method = SR)
            rec[rec<0] <- 1e-10
            NAAbi[1,y,1] <- rec
            for(s in 1:ns){
                ## can't take more than what's there
                Btemp <- sum(NAAbi[,y,s] * weight * sel * exp(-M/2))
                CWbi[y,s] <- sum(baranov(Fbi, Ms, NAAbi[,y,s]))
                if(CWbi[y,s] > 0.99 * Btemp){
                    CWbi[y,s] <- 0.75 * Btemp
                    Fbi <- sel * min(set$maxF/ns,
                                     getFM(CWbi[y,s], NAA = NAAbi[,y,s],
                                           M = Ms, weight = weightF, sel = sel))
                    Zbi <- Ms + Fbi
                }
                Ntemp <- NAAbi[,y,s] * exp(-Zbi)
                if(s < ns){
                    NAAbi[,y,s+1] <- Ntemp
                }else if(y < burnin){
                    NAAbi[amax, y+1, 1] <- Ntemp[amax] + Ntemp[amax-1]
                    for(a in 2:(amax-1)) NAAbi[a,y+1,1] <- Ntemp[a-1]
                }
            }
        }
        NAA[,1,] <- NAAbi[,burnin,]
    }else{
        ## set up NAA
        NAA[1,1,1] <- R0 * exp(eR0) * exp(initN[1])
        for(a in 2:amax){
            NAA[a,1,1] <- NAA[a-1,1,1] * exp(-(Ms[a-1] + initF*sel[a-1])/ns) * exp(initN[a])
        } ## HERE: other seasons
    }
    ## main loop
    for(y in 1:ny){
        Msy <- Ms * eM[y]
        maty <- mat * eMat[y]
        hy <- h * eH[y]
        R0y <- R0 * eR0[y]
        ## recruitment
        SSB[y,1] <- sum(NAA[,y,1] * weight * maty * exp(-pzbm * (Msy + Fvals[y] * sel / ns))) ## pre-recruitment mort
        rec <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSB[y,1],
                       R0 = R0y, method = SR)
        rec[rec<0] <- 1e-10
        NAA[1,y,1] <- rec * eR[y]
        ## seasons
        for(s in 1:ns){
            ## FAA and Z
            FAA[,y,s] <- sel * Fvals[y] * eF[y] / ns
            Z <- Msy + FAA[,y,s]
            ## can't take more than what's there
            Btemp <- sum(NAA[,y,s] * weight * sel * exp(-Msy/2))
            Ctemp <- sum(baranov(FAA[,y,s], Msy, NAA[,y,s]))
            if(Ctemp > 0.99 * Btemp){
                Ctemp <- 0.75 * Btemp
                FAA[,y,s] <- sel * min(set$maxF/ns,
                                       getFM(Ctemp, NAA = NAA[,y,s],
                                             M = Msy, weight = weightF, sel = sel))
                Z <- Msy + FAA[,y,s]
            }
            ## CAA
            CAA[,y,s] <- baranov(FAA[,y,s], Msy, NAA[,y,s])
            ## CW
            CW[y,s] <- sum(weightF * CAA[,y,s])
            ## TSB
            NAAmid <- NAA[,y,s] * exp(-Msy/2)
            TSB[y,s] <- sum(weight * NAAmid)
            ## SSB
            SSB[y,s] <- sum(NAAmid * weight * maty * exp(-pzbm * Z))
            ## ESB
            ESB[y,s] <- sum(NAAmid * weightF * sel)
            Ntemp <- NAA[,y,s] * exp(-Z)
            if(s < ns){
                NAA[,y,s+1] <- Ntemp
            }else if(y < ny)){
                NAA[amax, y+1, 1] <- Ntemp[amax] + Ntemp[amax-1]
                for(a in 2:(amax-1)) NAA[a,y+1,1] <- Ntemp[a-1]
            }
        }
    }

    ## return
    out <- NULL
    if(out.opt == 1){
        out$NAA <- NAA
        out$TSB <- TSB
        out$SSB <- SSB
        out$ESB <- ESB
        out$CAA <- CAA
        out$FAA <- FAA
        out$FM <- apply(FAA, c(2,3), function(x) mean(x / specdat$sel, na.rm=TRUE))
        out$CW <- CW
        out$TACs <- TACs
        out$errs <- errs
    }else if(out.opt == 2){
        if(is.null(refs)){
            warning("The list with reference points is needed!")
        }else{
            out <- TSB[ny]/refs[[depl.quant]]
        }
    }
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
    eImp <- set$eImp[ysim]

    if(is.null(eR)) {
        eR <- genDevs(1, set$sigmaR, set$rhoR)
    }

    if(is.null(eF)) eF <- rnorm(1, 0, set$sigmaF) - set$sigmaF^2/2
    if(is.null(eM)) eM <- rnorm(1, 0, set$sigmaM) - set$sigmaM^2/2
    if(is.null(eH)) eH <- rnorm(1, 0, set$sigmaH) - set$sigmaH^2/2
    if(is.null(eMat)) eMat <- rnorm(1, 0, set$sigmaMat) - set$sigmaMat^2/2
    if(is.null(eImp)) eImp <- rnorm(1, 0, set$sigmaImp) - set$sigmaImp^2/2
    if("errs" %in% names(hist)){
        errs <- list(eF = c(hist$errs$eF, eF),
                     eR = c(hist$errs$eR, eR),
                     eM = c(hist$errs$eM, eM),
                     eH = c(hist$errs$eH, eH),
                     eR0 = hist$errs$eR0,
                     eMat = c(hist$errs$eMat, eMat),
                     eMat = c(hist$errs$eImp, eImp))
    }else{
        errs <- list(eF = eF,
                     eR = eR,
                     eM = eM,
                     eH = eH,
                     eMat = eMat,
                     eImp = eImp)
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
    ## NAAmid <- NAA[yNAA-1,] * exp(-M/2)

    ## catch in weight = TAC
    if(tacs$id[nrow(tacs)] != "refFmsy" || is.na(tacs$id[nrow(tacs)])){
        TAC <- as.numeric(as.character(tacs$TAC[nrow(tacs)]))
        TACs[y] <- TAC
        FMtac <- min(set$maxF, getFM(TAC, NAA = NAA[y,], M = M, weight = weightF, sel = sel))
        FAA[y,] <- sel * FMtac * exp(eImp)
        ## CAA[y,] <- baranov(FAA[y,], M, NAAmid) ## new: "To prevent
        ## a control from removing all exploitable biomass from the
        ## population in a year, we set the achieved catch to 75% of
        ## the midyear exploitable biomass in that year in cases where
        ## TAC exceeded the exploitable biomass."
        CAAtmp <- baranov(FAA[y,], M, NAA[y,])
##        CAA[y,] <- CAAtmp
        CWreq <- sum(CAAtmp * weightF)
        CWpot <- sum(NAA[y,] * sel * weightF)
        ## print(paste0("Req:",CWreq))
        ## print(paste0("Pot:",CWpot))
        if(CWreq > CWpot){
            CAA[y,] <- 0.99 * NAA[y,] * sel
        }else{
            CAA[y,] <- CAAtmp
        }
    }else{
        FMtac <- hist$refs$Fmsy
        FAA[y,] <- sel * FMtac * exp(eImp)
        CAA[y,] <- baranov(FAA[y,], M, NAA[y,])
        TAC <- sum(CAA[y,] * weightF, na.rm = TRUE)
        TACs[y] <- TAC
        tacs$TAC[nrow(tacs)] <- TAC
    }

    ## if(FALSE){
    ##     ## CAA
    ##     NAAtmp <- NAA[yNAA-1,] * weightF ## * sel
    ##     CAAtmp <- TAC * (NAAtmp / sum(NAAtmp, na.rm=TRUE))## split TAC according to age distribution in pop
    ##     ##    if(any(is.na(CAAtmp > NAAmid))) browser()
    ##     CAAtmp[CAAtmp > NAAmid] <- NAAmid[CAAtmp > NAAmid]
    ##     CAA[y,] <- CAAtmp / weightF

    ##     ## FAA
    ##     FAA[y,] <- CAA[y,] / (NAA[yNAA-1,] * exp(-M/2))
    ##     FAA[y,] <- -log(1 - FAA[y,])
    ## }

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
    NAA[yNAA,1] <- rec * eR

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
