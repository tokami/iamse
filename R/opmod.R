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
    ns <- dat$ns
    nt <- ny * ns
    yvec <- dat$yvec
    svec <- dat$svec
    s1vec <- dat$s1vec
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
    Msels <- dat$Msels
    Msel <- dat$Msel
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
        eCmv <- set$eCmv
        eImv <- set$eImv
        if(is.null(eF)) eF <- genNoise(ny, set$noiseF[1], set$noiseF[2], bias.cor = set$noiseF[3])
        if(is.null(eR)) eR <- genNoise(ny, set$noiseR[1], set$noiseR[2], bias.cor = set$noiseR[3])
        if(is.null(eM)) eM <- genNoise(ny, set$noiseM[1], set$noiseM[2], bias.cor = set$noiseM[3])
        if(is.null(eH)) eH <- genNoise(ny, set$noiseH[1], set$noiseH[2], bias.cor = set$noiseH[3])
        if(is.null(eR0)) eR0 <- genNoise(ny, set$noiseR0[1], set$noiseR0[2], bias.cor = set$noiseR0[3])
        if(is.null(eMat)) eMat <- genNoise(ny, set$noiseMat[1], set$noiseMat[2], bias.cor = set$noiseMat[3])
        if(is.null(eImp)) eImp <- genNoise(ny, set$noiseImp[1], set$noiseImp[2], bias.cor = set$noiseImp[3])
        if(is.null(eC)) eC <- genNoise(ny, set$noiseC[1], set$noiseC[2], bias.cor = set$noiseC[3])
        if(is.null(eI)){
            eI <- list()
            for(i in 1:nsurv){
                eI[[i]] <- genNoise(ny, set$noiseI[1], set$noiseI[2], bias.cor = set$noiseI[3])
            }
        }
        if(is.null(eCmv)) eCmv <- genNoise(ny, set$noiseCmv[1], set$noiseCmv[2],
                                               bias.cor = set$noiseCmv[3],
                                               mv = TRUE, dat = dat)
        if(is.null(eImv)){
            eImv <- list()
            for(i in 1:nsurv){
                eImv[[i]] <- genNoise(ny, set$noiseImv[1], set$noiseImv[2],
                                    bias.cor = set$noiseImv[3],
                                    mv = TRUE, dat = dat)
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
        eCmv <- matrix(1, nrow=ny, ncol=amax)
        eImv <- list()
        for(i in 1:nsurv){
            eImv[[i]] <- matrix(1, nrow=ny, ncol=amax)
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
                 eI = eI,
                 eCmv = eCmv,
                 eImv = eImv)

    ## Flags
    mselFlag <- inherits(Msels, "list") && length(Msels) > 1

    ## containers
    TSB <- TSB1plus <- ESB <- SSB <- CW <- FM <- matrix(0, nrow=ny, ncol=ns)
    TACs <- TSBfinal <- SSBfinal <- ESBfinal <- rec <- rep(NA, ny)
    CAA <- vector("list", ny)
    for(i in 1:ny) CAA[[i]] <- matrix(NA, nrow = amax, ncol = ns)
    ## observations
    obsI <- timeI <- vector("list", nsurv)
    ## observations at age
    obsCA <- matrix(NA, nrow = length(idx), ncol = amax)
    obsIA <- vector("list", nsurv)
    for(i in 1:nsurv) obsIA[[i]] <- matrix(NA, nrow = ny, ncol = amax)

    ## burnin period
    if(is.null(set)) burnin <- 20 else burnin <- set$burnin
    if(is.numeric(burnin) && burnin > 0){
        NAAbi <- rep(NA, amax)
        NAAbi[1] <- R0 * exp(initN[1])
        for(a in 2:amax) NAAbi[a] <- NAAbi[a-1] * exp(-(Msel[[1]][a-1]+M[1]+sel[a-1]*Fy[1])) * exp(initN[a])
        for(y in 1:burnin){
            ## recruitment
            Fbi <- sels * Fs[1]
            Mbi <- t(t(Msels[[1]]) * Ms[1:ns])
            Zbi <- Mbi + Fbi
            SSBtemp <- sum(NAAbi * weights[,1] * mats[,1] * exp(-pzbm * Zbi[,1])) ## pre-recruitment mort
            SSBPR0 <- getSSBPR2(Mbi, dat$mats, dat$weights, fecun=1, amax, dat$R0,
                                ns = ns, recruitmentTiming = set$recruitmentTiming)
            recbi <- recfunc(h = h, SSBPR0 = SSBPR0, SSB = SSBtemp,
                           R0 = R0, method = dat$SR, bp = dat$bp,
                           beta = dat$recBeta, gamma = dat$recGamma)
            recbi[recbi<0] <- 1e-10
            NAAbi[1] <- recbi
            for(s in 1:ns){
                ## can't take more than what's there
                Btemp <- sum(NAAbi * weights[,s] * sels[,s] * exp(-(Mbi[,s])/2))
                CWbi <- sum(baranov(Fbi[,s], Mbi[,s], NAAbi) * weightFs[,s])
                if(CWbi > 0.99 * Btemp){
                    Fbi[,s] <- sels[,s] * getFM2(0.75 * Btemp, Btemp, 1/ns,
                                                 Mbi[,s], NAAbi,
                                                 weights[,s], weightFs[,s], sels[,s],
                                                 fmax = set$maxF/ns)
                    Zbi <- Mbi + Fbi
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
            NAA[a] <- NAA[a-1] * exp(-(M[1]*Msel[[1]][a-1] + Fy[1]*sel[a-1])) * exp(initN[a])
    }


    ## main loop
    for(y in 1:ny){
        ## Adding noise
        FM[y,] <- Fs[y] * eF[y]
        FAA <- FM[y,] * sels
        mselsInd <- ifelse(mselFlag, y, 1)
        MAA <- t(t(Msels[[mselsInd]]) * Ms[s1vec[y]:(s1vec[y]+ns-1)]) * eM[y]
        ZAA <- MAA + FAA
        maty <- mats * eMat[y]
        hy <- h * eH[y]
        R0y <- R0 * eR0[y]
        ## recruitment
        SSB[y,1] <- sum(NAA * weights[,1] * maty[,1] * exp(-pzbm * ZAA[,1])) ## pre-recruitment mort
        SSBPR0 <- getSSBPR2(MAA, maty, weights, fecun=1, amax, R0y, ns = ns,
                            recruitmentTiming = set$recruitmentTiming)
        rec[y] <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSB[y,1],
                       R0 = R0y, method = dat$SR, bp = dat$bp,
                       beta = dat$recBeta, gamma = dat$recGamma)
        rec[y]  <- ifelse(rec[y] < 0, 1e-10, rec[y])
        NAA[1] <- rec[y] * eR[y]
        ## seasons
        for(s in 1:ns){
            ## can't take more than what's there
            Btemp <- sum(NAA * weights[,s] * sels[,s] * exp(-MAA[,s]/2))
            CAA[[y]][,s] <- baranov(FAA[,s], MAA[,s], NAA)
            CW[y,s] <- sum(CAA[[y]][,s] * weightFs[,s])
            if(CW[y,s] > 0.99 * Btemp){
                FM[y,s] <- getFM2(0.75 * Btemp, Btemp, 1/ns, MAA[,s], NAA, weights[,s],
                                                 weightFs[,s], sels[,s], fmax = set$maxF/ns)
                FAA[,s] <- FM[y,s] * sels[,s]
                ZAA[,s] <- MAA[,s] + FAA[,s]
                CAA[[y]][,s] <- baranov(FAA[,s], MAA[,s], NAA)
                CW[y,s] <- sum(weightFs[,s] * CAA[[y]][,s])
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
                    ESBsurv <- NAAsurv * weightFs[,s] * sels[,s]
                    ## Total index (for spict)
                    obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * eI[[idxi[i]]][y])
                    if(is.null(timeI[[idxi[i]]]))
                        timeIi <- 0 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                    timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + surveyTimes[idxi[i]])
                    ## Index by age (sam, sms)
                    obsIA[[idxi[i]]][y,] <- q[idxi[i]] * NAAsurv * sels[,s] * eImv[[idxi[[i]]]][y,]
                }
            }

            ## ageing by season
            NAA <- Ntemp <- NAA * exp(-ZAA[,s])
            if(s == ns){
                ## end of year biomass for risk P(B/Blim)
                TSBfinal[y] <- sum(NAA * weights[,ns])
                ESBfinal[y] <- sum(NAA * weights[,ns] * sels[,ns])
                SSBfinal[y] <- sum(NAA * weights[,ns] * maty[,ns] * exp(-pzbm * ZAA[,ns]))
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
        obsIA[[i]] <- obsIA[[i]][(ny-nyhist+1):ny,]
        rownames(obsIA[[i]]) <- as.character((ny-nyhist+1):ny)
        colnames(obsIA[[i]]) <- as.character(0:dat$amax)
        attributes(obsIA[[i]]) <- c(attributes(obsIA[[i]]), list(time = set$surveyTimes))
    }

    ## catch observations
    ## ----------------
    if(ns > 1){
        ## seasonal OM
        ## ----------------
        if(nsC == 1){
            ## annual catches
            ## ----------------
            ## total catches (spict)
            obsC <- apply(CW[idx,], 1, sum) * eC[idx]
            timeC <- idx
            ## catch observations at age (SAM, SMS)
            for(y in 1:length(idx)) obsCA[y,] <- t(apply(CAA[[idx[y]]], 1, sum) * eCmv[y,])
        }else if(nsC == ns){
            ## seasonal catches (for each season in OM)
            ## ----------------
            obsC <- as.numeric(t(CW[idx,] * eC[idx]))
            timeC <- rep(idx, each = ns) + rep(seasonStart, nyhist)
            ## catch observations at age (seasonal CAA observations not yet implemented, what is target format?)
            for(y in 1:length(idx)) obsCA[y,] <- t(apply(CAA[[idx[y]]], 1, sum) * eCmv[y,])
        }else{
            ## seasonal catches not aligning with seasones in OM (not implemented)
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
    }else{
        ## annual OM
        if(nsC > 1) writeLines("Set dat$ns to > 1 for seasonal catches. Generating annual catches!")
        obsC <- CW[idx,] * eC[idx]
        timeC <- idx
        ## catch observations at age
        for(y in 1:length(idx)) obsCA[y,] <- t(CAA[[idx[y]]] * eCmv[y,])
    }
    rownames(obsCA) = as.character((ny-nyhist+1):ny)
    colnames(obsCA) <- as.character(0:dat$amax)

    ## other observations (required by SAM)
    ## ----------------

    ## natural mortality
    obsMAA <- as.matrix(dat$M[idx]) %*% t(as.matrix(dat$Msel[[1]]))
    rownames(obsMAA) = as.character((ny-nyhist+1):ny)
    colnames(obsMAA) = as.character(0:dat$amax)

    ## proportion mature
    propMature <- matrix(dat$mat, length(idx), amax, byrow = TRUE)
    rownames(propMature) = as.character((ny-nyhist+1):ny)
    colnames(propMature) = as.character(0:dat$amax)

    ## mean stock weight
    WAAs <- matrix(dat$weight, length(idx), amax, byrow = TRUE)
    rownames(WAAs) = as.character((ny-nyhist+1):ny)
    colnames(WAAs) = as.character(0:dat$amax)

    ## mean catch weight
    WAAc <- matrix(dat$weight, length(idx), amax, byrow = TRUE)
    rownames(WAAc) = as.character((ny-nyhist+1):ny)
    colnames(WAAc) = as.character(0:dat$amax)

    ## proportion females
    propFemale <- matrix(0.0, length(idx), amax)
    rownames(propFemale) <- as.character((ny-nyhist+1):ny)
    colnames(propFemale) <- as.character(0:dat$amax)

    ## landing fraction
    landFrac = matrix(1.0, length(idx), amax)
    rownames(landFrac) <- as.character((ny-nyhist+1):ny)
    colnames(landFrac) <- as.character(0:dat$amax)

    obs <- list(obsC = obsC,
                timeC = timeC,
                obsI = obsI,
                timeI = timeI,
                obsCA = obsCA,
                obsIA = obsIA,
                obsMAA = obsMAA,
                propMature = propMature,
                WAAs = WAAs,
                WAAc = WAAc,
                propFemale = propFemale,
                landFrac = landFrac)


    ## return
    out <- NULL
    if(out.opt == 1){
        out$lastNAA <- NAA
        out$lastFAA <- FAA[,ns]
        out$TSBfinal <- TSBfinal
        out$SSBfinal <- SSBfinal
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
        out$obs <- obs
    }else if(out.opt %in% c(2,3)){
        refs <- dat$ref
        if(is.null(refs)){
            warning("The reference points are not part of dat! Use estRef to estimate them")
        }else if(out.opt == 2){
            out <- TSBfinal[ny]/refs[[dat$depl.quant]]
        }else if(out.opt == 3){
            out <- SSBfinal[ny]/refs[[dat$depl.quant]]
        }
    }
    return(out)
}



## advance population
#' @name advancePop
#' @export
advancePop <- function(dat, hist, set, hcr, year){

    ## indices
    ny <- nrow(hist$TSB)
    y <- ny + 1
    ysim <- y - dat$ny
    ns <- dat$ns
    s1vec <- dat$s1vec
    nt <- ny * ns
    nsC <- set$catchSeasons
    nysim <- set$nysim
    assessYears <- seq(1, nysim, set$assessmentInterval)

    ## survey
    nsurv <- length(set$surveyTimes)
    seasonStart <- seq(0,1-1/ns,1/ns)
    idxS <- rep(0, nsurv)
    for(i in 1:nsurv){
        tmp <- seasonStart[seasonStart < set$surveyTimes[i]]
        idxS[i] <- which.min((tmp - set$surveyTimes[i])^2)
    }

    ## parameters
    tacs <- hist$tacs
    obs <- hist$obs
    refs <- dat$ref
    amax <- dat$amax + 1  ## age 0
    pzbm <- dat$pzbm
    R0 <- dat$R0
    h <- dat$h
    q <- dat$q
    if(length(q) < nsurv) q <- rep(q, nsurv)
    tacID <- hcr
    tacID2 <- unlist(strsplit(as.character(tacID), "-"))[1]
    Ms <- dat$Ms

    ## parameters per age
    weight <- dat$weight
    weights <- dat$weights
    weightF <- dat$weightF
    weightFs <- dat$weightFs
    mat <- dat$mat
    mats <- dat$mats
    sel <- dat$sel
    sels <- dat$sels
    Msel <- dat$Msel
    Msels <- dat$Msels

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
    eCmv <- set$eCmv[ysim,]
    eImv <- set$eImv[ysim,]
    if(is.null(eF)) eF <- genNoise(1, set$noiseF[1], set$noiseF[2], bias.cor = set$noiseF[3])
    if(is.null(eR)) eR <- genNoise(1, set$noiseR[1], set$noiseR[2], bias.cor = set$noiseR[3])
    if(is.null(eM)) eM <- genNoise(1, set$noiseM[1], set$noiseM[2], bias.cor = set$noiseM[3])
    if(is.null(eH)) eH <- genNoise(1, set$noiseH[1], set$noiseH[2], bias.cor = set$noiseH[3])
    if(is.null(eR0)) eR0 <- genNoise(1, set$noiseR0[1], set$noiseR0[2], bias.cor = set$noiseR0[3])
    if(is.null(eMat)) eMat <- genNoise(1, set$noiseMat[1], set$noiseMat[2], bias.cor = set$noiseMat[3])
    if(is.null(eImp)) eImp <- genNoise(1, set$noiseImp[1], set$noiseImp[2], bias.cor = set$noiseImp[3])
    if(is.null(eC)) eC <- genNoise(1, set$noiseC[1], set$noiseC[2], bias.cor = set$noiseC[3])
    if(is.null(eI)){
        eI <- list()
        for(i in 1:nsurv){
            eI[[i]] <- genNoise(1, set$noiseI[1], set$noiseI[2], bias.cor = set$noiseI[3])
        }
    }
    if(is.null(eCmv)) eCmv <- genNoise(1, set$noiseCmv[1], set$noiseCmv[2],
                                       bias.cor = set$noiseCmv[3],
                                       mv = TRUE, dat = dat)
    if(is.null(eImv)){
        eImv <- list()
        for(i in 1:nsurv){
            eImv[[i]] <- genNoise(1, set$noiseImv[1], set$noiseImv[2],
                                  bias.cor = set$noiseImv[3],
                                  mv = TRUE, dat = dat)
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
        errs$eCmv <- rbind(hist$errs$eCmv, eCmv)
        errs$eImv <- list()
        for(i in 1:nsurv){
            errs$eImv[[i]] = rbind(hist$errs$eImv[[i]], eImv[[i]])
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
                     eI = eI,
                     eCmv = eCmv,
                     eImv = eImv)
    }

    ## Flags
    mselFlag <- inherits(Msels, "list") && length(Msels) > 1

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
    rec <- c(hist$rec, NA)
    TSBfinal <- c(hist$TSBfinal, NA)
    SSBfinal <- c(hist$SSBfinal, NA)
    ESBfinal <- c(hist$ESBfinal, NA)

    ## project forward
    R0y <- R0 * eR0
    mselsInd <- ifelse(mselFlag, y, 1)
    MAA <- t(t(Msels[[mselsInd]]) * Ms[s1vec[y]:(s1vec[y]+ns-1)]) * eM
    hy <- h * eH
    maty <- mats * eMat

    ## TAC from previous year (if e.g. interim advice)
    ntac <- ns * set$assessmentInterval
    TACprev <- hist$TACprev ## TAC for full interval
    if(is.null(TACprev)){
        TACprev <- tail(as.numeric(t(hist$CW)),ntac)
    }
    ## Use only part before assessment from previous TAC in assessment year
    if(year %in% assessYears){
        TACreal <- rep(NA, ntac)
        if(set$assessmentTiming > 1){
            ind <- 1:(set$assessmentTiming-1)
            TACreal[ind] <- tail(TACprev, length(ind))
        }
    }else{
        ## Use full previous TAC if no assessment
        TACreal <- TACprev
    }


    NAA <- hist$lastNAA

    ## seasons
    for(s in 1:ns){

        ## Indices
        indtac <- (year - assessYears[max(which(assessYears <= year))]) * ns + s

        ## Recruitment
        if(s == set$recruitmentTiming){
            ## Survivors from previous season/year
            if(s == 1){
                Ztemp <- hist$lastFAA + MAA[,1]  ## i.e. FAA in s=1 is equal to last FAA (pot s=4)
            }else{
                Ztemp <- ZAA[,s-1]
            }
            NAAtemp <- NAA
            SSBtemp <- sum(NAAtemp * weights[,1] * maty[,1] * exp(-pzbm * Ztemp)) ## pre-recruitment mort
            SSBPR0 <- getSSBPR2(MAA, maty, weights, fecun=1, amax, R0y,
                                ns = ns, recruitmentTiming = set$recruitmentTiming)
            rec[y] <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSBtemp, R0 = R0y,
                              method = dat$SR, bp = dat$bp,
                              beta = dat$recBeta, gamma = dat$recGamma)
            rec[rec<0] <- 1e-10
            NAA[1] <- rec[y] * eR
        }

        ## Assessment
        if(year %in% assessYears && s == set$assessmentTiming){

            ## Estimate TAC
            if(hcr %in% c("refFmsy","consF")){
                ## Reference rule Fmsy
                tacs <- gettacs(tacs = tacs, id = hcr, TAC=NA)
            }else{
                ## True stock status
                TSBtmp <- sum(NAA * weights[,s])
                bbmsy <- TSBtmp/refs$Bmsy
                FMtmp <- as.numeric(t(FM))
                ffmsy <- sum(tail(FMtmp[1:((y*ns)-ns-(s-1))],4))/refs$Fmsy
                ## TAC
                tacs <- estTAC(obs = obs,
                               hcr = hcr,
                               tacs = tacs,
                               pars =
                                   list("ffmsy" = ffmsy,
                                        "bbmsy" = bbmsy))
                ## Split TAC to time steps
                TACs[y] <- as.numeric(as.character(tacs$TAC[nrow(tacs)])) * eImp
                TACreal <- rep(TACs[y] / ntac, ntac)
            }
        }

        ## Find F given TAC
        if(tacID2 == "refFmsy"){
            if(tacID == "refFmsy"){
                ## Fishing at Fmsy
                FMtac <- refs$Fmsy / ns
            }else{
                fraci <- as.numeric(unlist(strsplit(as.character(tacID), "-"))[2])
                ## Fishing at fraction of Fmsy
                FMtac <- (fraci * refs$Fmsy) / ns
            }
        }else if(tacID2 == "noF"){
            ## noF
            FMtac <- 0
        }else{
            ## any other HCR
            FMtac <- min(set$maxF/ns,
                         getFM(TACreal[indtac], NAA = NAA, M = MAA[,s],
                               weight = weightFs[,s], sel = sels[,s]))
        }

        ## Population dynamics
        FM[y,s] <- FMtac
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
            if(indtac < ntac){
                TACreal[indtac+1] <- TACreal[indtac+1] + TACreal[indtac] - CW[y,s]
            }else writeLines("Could not get full annual TAC.")
            TACreal[indtac] <- CW[y,s]
        }
        if(tacID2 == "refFmsy"){
            TACreal[indtac] <- sum(CAA * weightFs[,s], na.rm = TRUE)
            if(indtac == ntac) TACs[y] <- tacs$TAC[nrow(tacs)] <- sum(TACreal)
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
                ## survey observation: total catch in weight (spict)
                surveyTime <- set$surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                NAAsurv <- exp(log(NAA) - ZAA[,s] * surveyTime)
                ESBsurv <- NAAsurv * dat$weightFs[,s] * dat$sels[,s]
                obs$obsI[[idxi[i]]] <-
                    c(obs$obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * eI[[idxi[i]]])
                if(is.null(obs$timeI[[idxi[i]]])){
                    timeIi <- ny-nyhist+1
                }else timeIi <- floor(tail(obs$timeI[[idxi[i]]],1))
                obs$timeI[[idxi[i]]] <-
                    c(obs$timeI[[idxi[i]]], timeIi + 1 + set$surveyTimes[idxi[i]])
                ## survey observation: CAA
                obs$obsIA[[idxi[i]]] <- rbind(obs$obsIA[[idxi[i]]],
                                              q[idxi[i]] * NAAsurv * dat$sels[,s] * eImv[[idxi[[i]]]])
                rownames(obs$obsIA[[idxi[i]]])[nrow(obs$obsIA[[idxi[i]]])] <- as.character(y)  ## CHECK: instead of 1 (assuming one survey a year) this should allow for s surveys
                attributes(obs$obsIA[[i]]) <- c(attributes(obs$obsIA[[i]]), list(time = set$surveyTimes))
            }
        }

        ## ageing by season
        NAA <- Ntemp <- NAA * exp(-ZAA[,s])
        if(s == ns){
            ## end of year biomass for risk P(B/Blim)
            TSBfinal[y] <- sum(NAA * weights[,ns])
            ESBfinal[y] <- sum(NAA * weights[,ns] * sels[,ns])
            SSBfinal[y] <- sum(NAA * weights[,ns] * maty[,ns] * exp(-pzbm * ZAA[,ns]))
            ## Ageing by year
            NAA[amax] <- Ntemp[amax] + Ntemp[amax-1]
            for(a in 2:(amax-1)) NAA[a] <- Ntemp[a-1]
        }
    }

    ## catch observations
    if(ns > 1){
        if(nsC == 1){
            obs$obsC <- c(obs$obsC, sum(CW[y,]) * eC)
            if(!is.null(obs$timeC)) timeCi <- tail(obs$timeC,1) else timeCi <- ny-nyhist+1
            obs$timeC <- c(obs$timeC, timeCi + 1)
            ## catch observations at age (SAM, SMS)
            obs$obsCA <- rbind(obs$obsCA, t(apply(CAA, 1, sum) * eCmv))
        }else if(nsC == ns){
            obs$obsC <- c(obs$obsC, CW[y,] * eC)
            if(!is.null(obs$timeC))
                timeCi <- floor(tail(obs$timeC,1)) else timeCi <- ny-nyhist+1
            obs$timeC <- c(obs$timeC, timeCi + 1 + rep(seasonStart, ny))
            ## catch observations at age (SAM, SMS)
            obs$obsCA <- rbind(obs$obsCA, t(apply(CAA, 1, sum) * eCmv))
        }else{
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
    }else{
        if(nsC > 1) writeLines("Set dat$ns to > 1 for seasonal catches. Generating annual catches!")
        obs$obsC <- c(obs$obsC, sum(CW[y,]) * eC)
        if(!is.null(obs$timeC)) timeCi <- tail(obs$timeC,1) else timeCi <- ny-nyhist+1
        obs$timeC <- c(obs$timeC, timeCi + 1)
        ## catch observations at age (SAM, SMS)
        obs$obsCA <- rbind(obs$obsCA, t(CAA * eCmv))
    }
    rownames(obs$obsCA)[nrow(obs$obsCA)] <- as.character(y)

    ## natural mortalityp
    obs$obsMAA <- rbind(obs$obsMAA, y =
                    as.matrix(dat$M[y]) %*% t(as.matrix(dat$Msel[[1]])))
    rownames(obs$obsMAA)[nrow(obs$obsMAA)] = as.character(y)

    ## proportion mature
    obs$propMature <- rbind(obs$propMature, matrix(dat$mat, 1, amax, byrow = TRUE))
    rownames(obs$propMature)[nrow(obs$propMature)] = as.character(y)

    ## mean stock weight
    obs$WAAs <- rbind(obs$WAAs, matrix(dat$weight, 1, amax, byrow = TRUE))
    rownames(obs$WAAs)[nrow(obs$WAAs)] = as.character(y)

    ## mean catch weight
    obs$WAAc <- rbind(obs$WAAc, matrix(dat$weight, 1, amax, byrow = TRUE))
    rownames(obs$WAAc)[nrow(obs$WAAc)] = as.character(y)

    ## proportion females
    obs$propFemale <- rbind(obs$propFemale, matrix(0.0, 1, amax))
    rownames(obs$propFemale)[nrow(obs$propFemale)] <- as.character(y)

    ## landing fraction
    obs$landFrac = rbind(obs$landFrac, matrix(1.0, 1, amax))
    rownames(obs$landFrac)[nrow(obs$landFrac)] <- as.character(y)

    ## return
    out <- NULL
    out$lastNAA <- NAA
    out$lastFAA <- FAA[,ns]
    out$TSB <- TSB
    out$TSBfinal <- TSBfinal
    out$SSBfinal <- SSBfinal
    out$ESBfinal <- ESBfinal
    out$TACprev <- TACreal
    out$ESB <- ESB
    out$SSB <- SSB
    out$rec <- rec
    out$FM <- FM
    out$CW <- CW
    out$TACs <- TACs
    out$tacs <- tacs
    out$errs <- errs
    out$obs <- obs
    out$tacs <- tacs
    return(out)
}
