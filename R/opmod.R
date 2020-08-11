
## input data:
## Number of years (ny)
## Maximum age of species (amax)
## Natural mortality (M)
## weight-at-age start (weightS)
## weight-at-age mid (weightH)
## maturity-at-age (mat)
## proportion of Z before maturation (pzbm)

#' @name initPop
#'
#' @export
initPop <- function(dat, set = NULL, out.opt = 1){

    ## indices
    if(is.null(set)) set <- checkSet()
    ny <- dat$ny
    ns <- dat$nseason
    nt <- ny * ns
    nsC <- set$catchSeasons
    nyhist <- set$nyhist
    idx <- (ny - nyhist + 1):ny
    surveyTimes <- set$surveyTimes
    nsurv <- length(surveyTimes)
    ## closest season
    seasonStart <- seq(0,1-1/ns,1/ns)
    idxS <- rep(0, nsurv)
    for(i in 1:nsurv){
        tmp <- seasonStart[seasonStart <= surveyTimes[i]]
        idxS[i] <- which.min((tmp - surveyTimes[i])^2)
    }

    ## species data
    amax <- dat$amax + 1  ## age 0
    Ms <- dat$Ms
    M <- dat$M
    Fs <- dat$Fs
    Fy <- dat$FM
    weights <- dat$weights
    weight <- dat$weight
    weightFs <- dat$weightFs
    weightF <- dat$weightF
    mats <- dat$mats
    mat <- dat$mat
    sels <- dat$sels
    sel <- dat$sel
    pzbm <- dat$pzbm
    R0 <- dat$R0
    h <- dat$h
    initN <- dat$initN
    q <- dat$q
    if(length(q) < nsurv) q <- rep(q, nsurv)

    ## errors
    if(!is.null(set)){
        eF <- set$eF
        eR <- set$eR
        eM <- set$eM
        eH <- set$eH
        eR0 <- set$eR0
        eMat <- set$eMat
        eImp <- set$eImp
        eC <- set$eC
        eI <- set$eI
        if(is.null(eF)) eF <- rlnorm(ny, muconv(1,set$sigmaF), sdconv(1,set$sigmaF))
        if(is.null(eR)) {
            eR <- genDevs(ny, set$sigmaR, set$rhoR)
        }
        if(is.null(eM)) eM <- rlnorm(ny, muconv(1,set$sigmaM), sdconv(1,set$sigmaM))
        if(is.null(eH)) eH <- rlnorm(ny, muconv(1,set$sigmaH), sdconv(1,set$sigmaH))
        if(is.null(eR0)) eR0 <- rlnorm(ny, muconv(1,set$sigmaR0), sdconv(1,set$sigmaR0))
        if(is.null(eMat)) eMat <- rlnorm(ny, muconv(1,set$sigmaMat), sdconv(1,set$sigmaMat))
        if(is.null(eImp)) eImp <- rlnorm(ny, muconv(1,set$sigmaImp), sdconv(1,set$sigmaImp))
        if(is.null(eC)) eC <- rlnorm(ny, muconv(1,set$CVC), sdconv(1,set$CVC))
        if(is.null(eI)){
            eI <- list()
            for(i in 1:nsurv){
                eI[[i]] <- rlnorm(ny, muconv(1,set$CVI), sdconv(1,set$CVI))
            }
        }
    }else{
        eF <- rep(1, ny)
        eR <- rep(1, ny)
        eM <- rep(1, ny)
        eH <- rep(1, ny)
        eR0 <- rep(1, ny)
        eMat <- rep(1, ny)
        eImp <- rep(1, ny)
        eC <- rep(1, ny)
        eI <- list()
        for(i in 1:nsurv){
            eI[[i]] <- rep(1, ny)
        }
    }
    errs <- list(eF = eF,
                 eR = eR,
                 eM = eM,
                 eH = eH,
                 eR0 = eR0,
                 eMat = eMat,
                 eImp = eImp,
                 eC = eC,
                 eI = eI)

    ## containers
    TSB <- TSB1plus <- ESB <- SSB <- CW <- FM <- matrix(0, nrow=ny, ncol=ns)
    TACs <- TSBfinal <- ESBfinal <- rec <- rep(NA, ny)
    obsI <- vector("list", nsurv)
    timeI <- vector("list", nsurv)
    ## burnin period
    if(is.null(set)) burnin <- 20 else burnin <- set$burnin
    if(is.numeric(burnin) && burnin > 0){
        NAAbi <- rep(NA, amax)
        NAAbi[1] <- R0 * exp(initN[1])
        for(a in 2:amax) NAAbi[a] <- NAAbi[a-1] * exp(-(M[a-1]+sel[a-1]*Fy[1])) * exp(initN[a])
        for(y in 1:burnin){
            ## recruitment
            Fbi <- sels * Fs[1]
            Zbi <- Ms + Fbi
            SSBtemp <- sum(NAAbi * weights[,1] * mats[,1] * exp(-pzbm * Zbi[,1])) ## pre-recruitment mort
            SSBPR0 <- getSSBPR(dat$M, dat$mat, dat$weight, 1, amax, dat$R0) ## annual M
            recbi <- recfunc(h = h, SSBPR0 = SSBPR0, SSB = SSBtemp,
                           R0 = R0, method = dat$SR, bp = dat$bp,
                           beta = dat$recBeta, gamma = dat$recGamma)
            recbi[recbi<0] <- 1e-10
            NAAbi[1] <- recbi
            for(s in 1:ns){
                ## can't take more than what's there
                Btemp <- sum(NAAbi * weights[,s] * sels[,s] * exp(-Ms[,s]/2))
                CWbi <- sum(baranov(Fbi[,s], Ms[,s], NAAbi) * weightFs[,s])
                if(CWbi > 0.99 * Btemp){
                    Fbi[,s] <- sels[,s] * getFM2(0.75 * Btemp, Btemp, 1/ns, Ms[,s], NAAbi,
                                                 weights[,s], weightFs[,s], sels[,s],
                                                 fmax = set$maxF/ns)
                    Zbi <- Ms + Fbi
                }
                ## ageing by season or year
                NAAbi <- Ntemp <- NAAbi * exp(-Zbi[,s])
                if(s == ns){
                    NAAbi[amax] <- Ntemp[amax] + Ntemp[amax-1]
                    for(a in 2:(amax-1)) NAAbi[a] <- Ntemp[a-1]
                }
            }
        }
        NAA <- NAAbi
    }else{
        ## set up NAA
        NAA <- rep(NA, amax)
        NAA[1] <- exp(initN[1]) * R0
        for(a in 2:amax)
            NAA[a] <- NAA[a-1] * exp(-(M[a-1] + Fy[1]*sel[a-1])) * exp(initN[a])
    }


    ## main loop
    for(y in 1:ny){
        ## Adding noise
        FM[y,] <- Fs[y] * eF[y]
        MAA <- Ms * eM[y]
        FAA <- FM[y,] * sels
        ZAA <- MAA + FAA
        maty <- mats * eMat[y]
        hy <- h * eH[y]
        R0y <- R0 * eR0[y]
        ## recruitment
        SSB[y,1] <- sum(NAA * weights[,1] * maty[,1] * exp(-pzbm * ZAA[,1])) ## pre-recruitment mort
        SSBPR0 <- getSSBPR(dat$M * eM[y], dat$mat * eMat[y], dat$weight, 1, amax, dat$R0) ## annual M
        rec[y] <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSB[y,1],
                       R0 = R0y, method = dat$SR, bp = dat$bp,
                       beta = dat$recBeta, gamma = dat$recGamma)
        rec[y]  <- ifelse(rec[y] < 0, 1e-10, rec[y])
        NAA[1] <- rec[y] * eR[y]
        ## seasons
        for(s in 1:ns){
            ## can't take more than what's there
            Btemp <- sum(NAA * weights[,s] * sels[,s] * exp(-MAA[,s]/2))
            CAA <- baranov(FAA[,s], MAA[,s], NAA)
            CW[y,s] <- sum(CAA * weightFs[,s])
            if(CW[y,s] > 0.99 * Btemp){
                FM[y,s] <- getFM2(0.75 * Btemp, Btemp, 1/ns, MAA[,s], NAA, weights[,s],
                                                 weightFs[,s], sels[,s], fmax = set$maxF/ns)
                FAA[,s] <- FM[y,s] * sels[,s]
                ZAA[,s] <- MAA[,s] + FAA[,s]
                CAA <- baranov(FAA[,s], MAA[,s], NAA)
                CW[y,s] <- sum(weightFs[,s] * CAA)
            }
            ## TSB
            TSB[y,s] <- sum(NAA * weights[,s])
            TSB1plus[y,s] <- sum(NAA[-1] * weights[-1,s])
            ## ESB
            ESB[y,s] <- sum(NAA * weights[,s] * sels[,s])
            ## SSB
            SSB[y,s] <- sum(NAA * weights[,s] * maty[,s] * exp(-pzbm * ZAA[,s]))

            ## index observations
            if(s %in% idxS){
                idxi <- which(idxS == s)
                for(i in 1:length(idxi)){
                    surveyTime <- surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                    NAAsurv <- exp(log(NAA) - ZAA[,s] * surveyTime)
                    ESBsurv <- sum(NAAsurv * weightFs[,s] * sels[,s])
                    obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * ESBsurv * eI[[idxi[i]]][y])
                    if(is.null(timeI[[idxi[i]]]))
                        timeIi <- 0 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                    timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + surveyTimes[idxi[i]])
                }
            }

            ## ageing by season
            NAA <- Ntemp <- NAA * exp(-ZAA[,s])
            if(s == ns){
                ## end of year biomass for risk P(B/Blim)
                TSBfinal[y] <- sum(NAA * weights[,ns])
                ESBfinal[y] <- sum(NAA * weights[,ns] * sels[,ns])
                ## Ageing by year
                NAA[amax] <- Ntemp[amax] + Ntemp[amax-1]
                for(a in 2:(amax-1)) NAA[a] <- Ntemp[a-1]
            }
        }
    }


    ## account for nyhist
    for(i in 1:nsurv){
        timeI[[i]] <- timeI[[i]][(ny-nyhist+1):ny]
        obsI[[i]] <- obsI[[i]][(ny-nyhist+1):ny]
    }

    ## catch observations
    if(ns > 1){
        if(nsC == 1){
            obsC <- apply(CW[idx,], 1, sum) * eC[idx]
            timeC <- idx
        }else if(nsC == ns){
            obsC <- as.numeric(t(CW[idx,] * eC[idx]))
            timeC <- rep(idx, each = ns) + rep(seasonStart, nyhist)
        }else{
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
    }else{
        if(nsC > 1) writeLines("Set dat$nseasons to > 1 for seasonal catches. Generating annual catches!")
        obsC <- CW[idx,] * eC[idx]
        timeC <- idx
    }

    inp <- list(obsC = obsC,
                timeC = timeC,
                obsI = obsI,
                timeI = timeI)

    ## return
    out <- NULL
    if(out.opt == 1){
        out$lastNAA <- NAA
        out$lastFAA <- FAA[,ns]
        out$TSBfinal <- TSBfinal
        out$ESBfinal <- ESBfinal
        out$rec <- rec
        out$TSB <- TSB
        out$TSB1plus <- TSB1plus
        out$ESB <- ESB
        out$SSB <- SSB
        out$FM <- FM
        out$CW <- CW
        out$TACs <- TACs
        out$errs <- errs
        out$inp <- inp
    }else if(out.opt == 2){
        refs <- dat$ref
        if(is.null(refs)){
            warning("The reference points are not part of dat! Use estRef to estimate them")
        }else{
            out <- TSBfinal[ny]/refs[[dat$depl.quant]]
        }
    }
    return(out)
}



## advance population
#' @name advancePop
#' @export
advancePop <- function(dat, hist, set, tacs){

    ## indices
    ny <- nrow(hist$TSB)
    y <- ny + 1
    ysim <- y - dat$ny
    ns <- dat$nseason
    nt <- ny * ns
    nsC <- set$catchSeasons

    ## survey
    nsurv <- length(set$surveyTimes)
    seasonStart <- seq(0,1-1/ns,1/ns)
    idxS <- rep(0, nsurv)
    for(i in 1:nsurv){
        tmp <- seasonStart[seasonStart < set$surveyTimes[i]]
        idxS[i] <- which.min((tmp - set$surveyTimes[i])^2)
    }

    ## parameters
    amax <- dat$amax + 1  ## age 0
    pzbm <- dat$pzbm
    R0 <- dat$R0
    h <- dat$h
    q <- dat$q
    if(length(q) < nsurv) q <- rep(q, nsurv)
    tacID <- tacs$id[nrow(tacs)]

    ## parameters per age

    Ms <- dat$Ms
    weight <- dat$weight
    weights <- dat$weights
    weightF <- dat$weightF
    weightFs <- dat$weightFs
    mat <- dat$mat
    mats <- dat$mats
    sel <- dat$sel
    sels <- dat$sels

    ## errors
    eF <- set$eF[ysim]
    eR <- set$eR[ysim]
    eM <- set$eM[ysim]
    eH <- set$eH[ysim]
    eR0 <- set$eR0[ysim]
    eMat <- set$eMat[ysim]
    eImp <- set$eImp[ysim]
    eC <- set$eC[ysim]
    eI <- list()
    for(i in 1:nsurv){
        eI[[i]] <- set$eI[[i]][ysim]
    }
    if(is.null(eR)) {
        eR <- genDevs(1, set$sigmaR, set$rhoR)
    }
    if(is.null(eF)) eF <- rlnorm(1, muconv(1,set$sigmaF), sdconv(1,set$sigmaF))
    if(is.null(eM)) eM <- rlnorm(1, muconv(1,set$sigmaM), sdconv(1,set$sigmaM))
    if(is.null(eR0)) eR0 <- rlnorm(1, muconv(1,set$sigmaR0), sdconv(1,set$sigmaR0))
    if(is.null(eH)) eH <- rlnorm(1, muconv(1,set$sigmaH), sdconv(1,set$sigmaH))
    if(is.null(eMat)) eMat <- rlnorm(1, muconv(1,set$sigmaMat), sdconv(1,set$sigmaMat))
    if(is.null(eImp)) eImp <- rlnorm(1, muconv(1,set$sigmaImp), sdconv(1,set$sigmaImp))
    if(is.null(eC)) eC <- rlnorm(1, muconv(1,set$CVC), sdconv(1,set$CVC))
    if(is.null(eI)){
        eI <- list()
        for(i in 1:nsurv){
            eI[[i]] <- rlnorm(1, muconv(1,set$CVI), sdconv(1,set$CVI))
        }
    }
    if("errs" %in% names(hist)){
        errs <- list(eF = c(hist$errs$eF, eF),
                     eR = c(hist$errs$eR, eR),
                     eM = c(hist$errs$eM, eM),
                     eH = c(hist$errs$eH, eH),
                     eR0 = c(hist$errs$eR0,eR0),
                     eMat = c(hist$errs$eMat, eMat),
                     eImp = c(hist$errs$eImp, eImp),
                     eC = c(hist$errs$eC, eC))
        errs$eI <- list()
        for(i in 1:nsurv){
            errs$eI[[i]] = c(hist$errs$eI[[i]], eI[[i]])
        }
    }else{
        errs <- list(eF = eF,
                     eR = eR,
                     eM = eM,
                     eH = eH,
                     eR0 = eR0,
                     eMat = eMat,
                     eImp = eImp,
                     eC = eC,
                     eI = eI)
    }

    ## Containers
    tmp <- matrix(0, 1, ns)
    TSB <- rbind(hist$TSB,tmp)
    SSB <- rbind(hist$SSB,tmp)
    ESB <- rbind(hist$ESB,tmp)
    CW  <- rbind(hist$CW,tmp)
    FM  <- rbind(hist$FM,tmp)
    NAA <- rep(0, amax)
    FAA <- ZAA <- MAA <- matrix(NA, amax, ns)
    TACs <- c(hist$TACs, NA)
    TSBfinal <- c(hist$TSBfinal, NA)
    ESBfinal <- c(hist$ESBfinal, NA)
    TACreal <- rep(NA, ns)
    ## for observations
    if(!is.null(hist$inp$obsC)){
        obsC <- hist$inp$obsC
    }
    if(!is.null(hist$inp$timeC)){
        timeC <- hist$inp$timeC
    }
    if(is.null(hist$inp$timeI)){
        timeI <- vector("list", nsurv)
    }else{
        timeI <- hist$inp$timeI
    }
    if(is.null(hist$inp$obsI)){
        obsI <- vector("list", nsurv)
    }else{
        obsI <- hist$inp$obsI
    }

    ## project forward
    MAA <- Ms * eM
    R0y <- R0 * eR0
    hy <- h * eH
    maty <- mats * eMat

    ## Survivors from previous season/year
    NAA <- hist$lastNAA
    Ztemp <- hist$lastFAA + dat$Ms[,1]  ## big assumptions that FAA in s=1 is equal to last FAA (pot s=4)

    ## recruitment
    SSBtemp <- sum(NAA * weights[,1] * maty[,1] * exp(-pzbm * Ztemp)) ## pre-recruitment mort
    SSBPR0 <- getSSBPR(dat$M * eM, dat$mat * eMat, dat$weight, 1, amax, dat$R0) ## annual M
    rec <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSBtemp, R0 = R0y,
                   method = dat$SR, bp = dat$bp,
                   beta = dat$recBeta, gamma = dat$recGamma)
    rec[rec<0] <- 1e-10
    NAA[1] <- rec * eR

    ## Define F/TACs
    if(tacID == "refFmsy"){
        ## Fmsy
        FMtac <- dat$ref$Fmsy / ns
    }else if(tacID == "noF"){
        ## noF
        FMtac <- 0
    }else{
        ## any other HCR
        TAC <- as.numeric(as.character(tacs$TAC[nrow(tacs)]))
        TACs[y] <- TAC
        TACreal <- rep(TAC/ns, ns)
    }

    ## seasons
    for(s in 1:ns){
        if(!(tacID %in% c("refFmsy","noF"))){
            FMtac <- min(set$maxF/ns,
                         getFM(TACreal[s], NAA = NAA, M = MAA[,s],
                               weight = weightFs[,s], sel = sels[,s]))
        }
        FM[y,s] <- FMtac * eImp
        FAA[,s] <- FM[y,s] * sels[,s]
        ZAA[,s] <- MAA[,s] + FAA[,s]
        ## can't take more than what's there
        Btemp <- sum(NAA * weights[,s] * sels[,s] * exp(-MAA[,s]/2))
        CAA <- baranov(FAA[,s], MAA[,s], NAA)
        CW[y,s] <- sum(CAA * weightFs[,s])
        if(CW[y,s] > 0.99 * Btemp){
            FM[y,s] <- getFM2(0.75 * Btemp, Btemp, 1/ns, MAA[,s], NAA, weights[,s],
                              weightFs[,s], sels[,s], fmax = set$maxF/ns)
            FAA[,s] <- FM[y,s] * sels[,s]
            ZAA[,s] <- MAA[,s] + FAA[,s]
            CAA <- baranov(FAA[,s], MAA[,s], NAA)
            CW[y,s] <- sum(CAA * weightFs[,s])
            if(s < ns){
                TACreal[s+1] <- TACreal[s+1] + TACreal[s] - CW[y,s]
            }else writeLines("Could not get full annual TAC.")
            TACreal[s] <- CW[y,s]
        }
        if(tacID == "refFmsy"){
            TACreal[s] <- sum(CAA * weightFs[,s], na.rm = TRUE)
            if(s == ns) TACs[y] <- tacs$TAC[nrow(tacs)] <- sum(TACreal)
        }
        ## TSB
        TSB[y,s] <- sum(NAA * weights[,s])
        ## ESB
        ESB[y,s] <- sum(NAA * weights[,s] * sels[,s])
        ## SSB
        SSB[y,s] <- sum(NAA * weights[,s] * maty[,s] * exp(-pzbm * ZAA[,s]))

        ## index observations
        if(s %in% idxS){
            idxi <- which(idxS == s)
            for(i in 1:length(idxi)){
                surveyTime <- set$surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                NAAsurv <- exp(log(NAA) - ZAA[,s] * surveyTime)
                ESBsurv <- sum(NAAsurv * dat$weightFs[,s] * dat$sels[,s])
                obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * ESBsurv * eI[[idxi[i]]])
                if(is.null(timeI[[idxi[i]]]))
                    timeIi <- ny-nyhist+1 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + set$surveyTimes[idxi[i]])
            }
        }

        ## ageing by season
        NAA <- Ntemp <- NAA * exp(-ZAA[,s])
        if(s == ns){
            ## end of year biomass for risk P(B/Blim)
            TSBfinal[y] <- sum(NAA * weights[,ns])
            ESBfinal[y] <- sum(NAA * weights[,ns] * sels[,ns])
            ## Ageing by year
            NAA[amax] <- Ntemp[amax] + Ntemp[amax-1]
            for(a in 2:(amax-1)) NAA[a] <- Ntemp[a-1]
        }
    }



    ## catch observations
    if(ns > 1){
        if(nsC == 1){
            obsC <- c(obsC, sum(CW[y,]) * eC)
            if(!is.null(timeC)) timeCi <- tail(timeC,1) else timeCi <- ny-nyhist+1
            timeC <- c(timeC, timeCi + 1)
        }else if(nsC == ns){
            obsC <- c(obsC, CW[y,] * eC)
            if(!is.null(timeC)) timeCi <- floor(tail(timeC,1)) else timeCi <- ny-nyhist+1
            timeC <- c(timeC, timeCi + 1 + rep(seasonStart, ny))
        }else{
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
    }else{
        if(nsC > 1) writeLines("Set dat$nseasons to > 1 for seasonal catches. Generating annual catches!")
        obsC <- c(obsC, sum(CW[y,]) * eC)
        if(!is.null(timeC)) timeCi <- tail(timeC,1) else timeCi <- ny-nyhist+1
        timeC <- c(timeC, timeCi + 1)
    }

    inp <- list(obsC = obsC,
                timeC = timeC,
                obsI = obsI,
                timeI = timeI)

    ## return
    out <- NULL
    out$lastNAA <- NAA
    out$lastFAA <- FAA[,ns]
    out$TSB <- TSB
    out$TSBfinal <- TSBfinal
    out$ESBfinal <- ESBfinal
    out$ESB <- ESB
    out$SSB <- SSB
    out$FM <- FM
    out$CW <- CW
    out$TACs <- TACs
    out$tacs <- tacs
    out$errs <- errs
    out$inp <- inp
    return(out)
}
