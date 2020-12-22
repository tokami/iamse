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
initPop <- function(dat, set = NULL, out.opt = 1, verbose = TRUE){
    ## indices
    if(is.null(set)) set <- checkSet()
    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    yvec <- dat$yvec
    svec <- dat$svec
    asvec <- dat$asvec
    savec <- dat$savec
    s1vec <- dat$s1vec
    s1avec <- dat$s1avec
    nsC <- set$catchSeasons
    nyhist <- set$nyhist
    if(nyhist > ny){
        nyhist <- ny
        if(verbose) writeLines("Period for observations ('set$nyhist') is larger than historical period ('dat$ny'). Setting nyhist == ny. ")
    }
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
    asmax <- amax * ns ## continuous max age
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
    spawning <- dat$spawning


    ## errors
    if(!is.null(set)){
        eF <- set$eF
        eR <- set$eR
        eM <- set$eM
        eH <- set$eH
        eR0 <- set$eR0
        eMat <- set$eMat
        eW <- set$eW
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
        if(is.null(eW)) eW <- genNoise(ny, set$noiseW[1], set$noiseW[2], bias.cor = set$noiseW[3])
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
        eW <- rep(1, ny)
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
                 eW = eW,
                 eImp = eImp,
                 eC = eC,
                 eI = eI,
                 eCmv = eCmv,
                 eImv = eImv)

    ## Flags
    mselFlag <- inherits(Msels, "list") && length(Msels) > 1

    ## containers
    TSB <- TSB1plus <- ESB <- SSB <- CW <- rec <- matrix(0, nrow=ny, ncol=ns)
    TACs <- TSBfinal <- SSBfinal <- ESBfinal <- rep(NA, ny)
    CAA <- vector("list", ny)
    for(i in 1:ny) CAA[[i]] <- matrix(NA, nrow = asmax, ncol = ns)
    ## observations
    obsI <- timeI <- vector("list", nsurv)
    ## observations at age
    obsCA <- matrix(NA, nrow = length(idx), ncol = amax)
    obsIA <- vector("list", nsurv)
    for(i in 1:nsurv) obsIA[[i]] <- matrix(NA, nrow = ny, ncol = amax)

    ## Initialise NAA
    hy <- h * eH[1]
    R0y <- R0 * eR0[1]
    maty <- as.numeric(t(mats)) * eMat[1]
    weighty <- as.numeric(t(weights)) * eW[1]
    sely <- as.numeric(t(sels))
    msely <- as.numeric(t(Msels[[1]]))
    NAAbi2 <- NAAbi <- matrix(0, asmax, ns)
    NAAbi[1,] <- R0y * spawning ## * exp(initN[1])
    ZAA <-  Ms[1:ns] * msely + Fs[1,] * sely
    ## each season
    for(as in 2:asmax)
        NAAbi[as,] <- NAAbi[as-1,] * exp(-ZAA[as-1])
    ## annual total mortality for plus group
    ## if(ns > 1){
    ##     for(s in 1:ns){
    ##         NAAtmp <- NAAbi[asmax,s]
    ##         s2 <- s
    ##         while(s2 < ns){
    ##             NAAtmp <- NAAtmp * exp(-ZAA[asmax-ns+s2])
    ##             s2 <- s2 + 1
    ##         }
    ##         NAAbi[asmax,s] <- NAAtmp
    ##     }
    ## }
    ## plusgroup <- NAAbi[asmax,] / (1 - exp(-sum(ZAA[(asmax-ns+1):asmax])))
    ## all seasons combined
    for(s in 1:ns){
        NAAbi2[seq(s,asmax,ns),s] <- NAAbi[seq(s,asmax,ns),s]
    }
    NAASbi <- rowSums(NAAbi2)
    NAASbi[1] <- 0
    ## plusgroup
## OLD:    NAASbi[asmax] <- sum(plusgroup)  ## REMOVE:
    plusgroup <- NAASbi[(asmax-ns+1):asmax] / (1 - exp(-sum(ZAA[(asmax-ns+1):asmax])))
    if(ns == 1){
        NAASbi[asmax] <- sum(plusgroup)
    }else{
        NAASbi[asmax] <- sum(plusgroup) - sum(NAASbi[(asmax-ns+1):asmax])
    }
    NAASbi2 <- NAASbi ## REMOVE:

    ## Burn-in period
    if(is.null(set)) burnin <- 5e2 else burnin <- set$burnin
    if(is.numeric(burnin) && burnin > 0){
        for(y in 1:burnin){
            Fbi <- Fs[1,] * sely
            Mbi <- Ms[1:ns] * msely
            Zbi <- Mbi + Fbi
            for(s in 1:ns){
                ## recruitment
                if(spawning[s] > 0){
                    SSBtmp <- sum(NAASbi * weighty * maty * exp(-pzbm * Zbi))
                    SSBPR0 <- getSSBPR(Zbi, maty, weighty, fecun=1, asmax, ns,
                                       spawning)
                    recbi <- spawning[s] * recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSBtmp,
                                                   R0 = R0y, method = dat$SR, bp = dat$bp,
                                                   beta = dat$recBeta, gamma = dat$recGamma)
                    recbi[recbi<0] <- 1e-10
                    NAASbi[1] <- recbi
                }
                CWbi <- sum(baranov(Fbi, Mbi, NAASbi) * as.numeric(t(weightFs)))
                ## can't take more than what's there
                Btmp <- sum(NAASbi * as.numeric(t(weights)) * as.numeric(t(sels)) * exp(-Mbi/2))
                if(is.na(CWbi) || is.na(Btmp)){
                    stop("Something went wrong in the burnin period. Catch or biomass is NA.")
                }
                if(CWbi > 0.99 * Btmp){
                    Fbi <- sely * min(getFM(0.75 * Btmp,
                                             NAA = NAASbi, MAA = Mbi,
                                             sels = sely,
                                             weights = weighty,
                                             seasons = s, ns = ns, y = y,
                                             h = hy, asmax = asmax, mats = maty,
                                             pzbm = pzbm, spawning = spawning,
                                             R0 = R0y, SR = dat$SR, bp = dat$pb, recBeta = dat$recBeta,
                                             recGamma = dat$recGamma, eR = eR,
                                             lastFM = 0.001),
                                      set$maxF/ns)
                    Zbi <- Mbi + Fbi
                }
                ## ageing by season or year
                NAASbi <- NAASbi * exp(-Zbi)
                NAASbi[asmax] <- NAASbi[asmax] + NAASbi[asmax-1]
                for(as in (asmax-1):2) NAASbi[as] <- NAASbi[as-1]
                NAASbi[1] <- 0
            }
        }
    }
    NAAS <- NAASbi

##    print(cbind(NAAbi,NAASbi2,NAASbi))

    ## main loop
    for(y in 1:ny){

        ## Adding noise
        sely <- as.numeric(t(sels))
        FAA <- Fs[y,] * eF[y] * sely
        mselsInd <- ifelse(mselFlag, y, 1)
        MAA <- Ms[yvec==y] * as.numeric(t(Msels[[mselsInd]])) * eM[y]
        ZAA <- MAA + FAA
        maty <- as.numeric(t(mats)) * eMat[y]
        weighty <- as.numeric(t(weights)) * eW[y]
        weightFy <- as.numeric(t(weightFs)) * eW[y]
        hy <- h * eH[y]
        R0y <- R0 * eR0[y]

        ## seasons
        for(s in 1:ns){
            ## recruitment
            if(spawning[s] > 0){
                ## Survivors from previous season/year
                SSB[y,s] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA)) ## pre-recruitment mort
                SSBPR0 <- getSSBPR(ZAA, maty, weighty, fecun=1, asmax,
                                   ns, spawning)
                rec[y,s] <-  spawning[s] * recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSB[y,s],
                                                   R0 = R0y, method = dat$SR, bp = dat$bp,
                                                   beta = dat$recBeta, gamma = dat$recGamma)
                rec[rec<0] <- 1e-10
                NAAS[1] <- rec[y,s] * eR[y]
            }

            ## catch
            CAA[[y]][,s] <- baranov(FAA, MAA, NAAS)
            CW[y,s] <- sum(CAA[[y]][,s] * weightFy)

            ## can't take more than what's there
            Btmp <- sum(NAAS * weighty * sely * exp(-MAA/2))
            if(CW[y,s] > 0.99 * Btmp){
                Fs[y,s] <- min(getFM(0.75 * Btmp,
                                             NAA = NAAS, MAA = MAA,
                                             sels = sely,
                                             weights = weighty,
                                             seasons = s, ns = ns, y = y,
                                             h = hy, asmax = asmax, mats = maty,
                                             pzbm = pzbm, spawning = spawning,
                                             R0 = R0y, SR = dat$SR, bp = dat$pb, recBeta = dat$recBeta,
                                             recGamma = dat$recGamma, eR = eR,
                                             lastFM = FM[y,s]),
                                      set$maxF/ns)
                FAA <- Fs[y,] * sely
                ZAA <- MAA + FAA
                CAA[[y]][,s] <- baranov(FAA, MAA, NAAS)
                CW[y,s] <- sum(weightFy * CAA[[y]][,s])
            }
            ## TSB
            TSB[y,s] <- sum(NAAS * weighty)
            TSB1plus[y,s] <- sum(NAAS[-1] * weighty[-1])
            ## ESB
            ESB[y,s] <- sum(NAAS * weighty * sely)
            ## SSB
            SSB[y,s] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))

            ## index observations
            if(s %in% idxS){
                idxi <- which(idxS == s)
                for(i in 1:length(idxi)){
                    surveyTime <- surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                    NAAsurv <- exp(log(NAAS) - ZAA * surveyTime)
                    ESBsurv <- NAAsurv * weightFy * sely
                    ## Total index (for spict)
                    obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * eI[[idxi[i]]][y])
                    if(is.null(timeI[[idxi[i]]]))
                        timeIi <- 0 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                    timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + surveyTimes[idxi[i]])
                    ## Index by age (sam, sms)
                    obsIA[[idxi[i]]][y,] <- q[idxi[i]] * as.numeric(by(NAAsurv*sely, asvec, sum)) *
                        eImv[[idxi[[i]]]][y,]
                }
            }

            ## Exponential decay
            NAAS <- NAAS * exp(-ZAA)
            if(s == ns){
                ## end of year biomasses
                TSBfinal[y] <- sum(NAAS * weighty)
                ESBfinal[y] <- sum(NAAS * weighty * sely)
                SSBfinal[y] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))
            }
            ## Continuous ageing
            NAAS[asmax] <- NAAS[asmax] + NAAS[asmax-1]
            for(as in (asmax-1):2) NAAS[as] <- NAAS[as-1]
            NAAS[1] <- 0
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
            for(y in 1:length(idx)) obsCA[y,] <- as.numeric(by(apply(CAA[[idx[y]]], 1, sum) * eCmv[y,],asvec,sum))
        }else if(nsC == ns){
            ## seasonal catches (for each season in OM)
            ## ----------------
            obsC <- as.numeric(t(CW[idx,] * eC[idx]))
            timeC <- rep(idx, each = ns) + rep(seasonStart, nyhist)
            ## catch observations at age (seasonal CAA observations not yet implemented, what is target format?)
            for(y in 1:length(idx)) obsCA[y,] <- as.numeric(by(apply(CAA[[idx[y]]], 1, sum) * eCmv[y,],asvec,sum))
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
        for(y in 1:length(idx)) obsCA[y,] <- as.numeric(by(CAA[[idx[y]]] * eCmv[y,],asvec,sum))
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
        out$lastNAA <- NAAS
        out$lastFAA <- FAA
        out$TSBfinal <- TSBfinal
        out$SSBfinal <- SSBfinal
        out$ESBfinal <- ESBfinal
        out$rec <- rec
        out$TSB <- TSB
        out$TSB1plus <- TSB1plus
        out$ESB <- ESB
        out$SSB <- SSB
        out$Fs <- Fs
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
    asvec <- dat$asvec
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
    asmax <- amax * ns
    pzbm <- dat$pzbm
    q <- dat$q
    if(length(q) < nsurv) q <- rep(q, nsurv)
    tacID <- hcr
    tacID2 <- unlist(strsplit(as.character(tacID), "_"))[1]
    Ms <- dat$Ms
    spawning <- dat$spawning
    assessmentTiming <- set$assessmentTiming
    if(ns == 1){
        assessSeasons <- assessmentTiming
    }else{
        assessSeasons <- rep(1:ns, 20)[assessmentTiming:(assessmentTiming+ns-1)]
    }

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
    ## TODO: put this in external wrapper function
    eF <- set$eF[ysim]
    eR <- set$eR[ysim]
    eM <- set$eM[ysim]
    eH <- set$eH[ysim]
    eR0 <- set$eR0[ysim]
    eMat <- set$eMat[ysim]
    eW <- set$eW[ysim]
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
    if(is.null(eW)) eW <- genNoise(1, set$noiseW[1], set$noiseW[2], bias.cor = set$noiseW[3])
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
                     eW = c(hist$errs$eW, eW),
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
                     eW = eW,
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
    Fs  <- rbind(hist$Fs,tmp)
    NAA <- rep(0, amax)
    FAA <- ZAA <- MAA <- CAA <- matrix(NA, asmax, ns)
    TACs <- c(hist$TACs, NA)
    rec <- rbind(hist$rec, tmp)
    TSBfinal <- c(hist$TSBfinal, NA)
    SSBfinal <- c(hist$SSBfinal, NA)
    ESBfinal <- c(hist$ESBfinal, NA)

    ## project forward
    R0y <- dat$R0 * eR0
    mselsInd <- ifelse(mselFlag, y, 1)
    MAA <- Ms[s1vec[y]:(s1vec[y]+ns-1)] * as.numeric(t(Msels[[mselsInd]])) * eM
    hy <- dat$h * eH
    maty <- as.numeric(t(mats)) * eMat
    sely <- as.numeric(t(sels))
    weighty <- as.numeric(t(weights)) * eW
    weightFy <- as.numeric(t(weightFs)) * eW

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


    NAAS <- hist$lastNAA

##    if(y == 22) browser()

    ## seasons
    for(s in 1:ns){

        ## Indices
        indtac <- (year - assessYears[max(which(assessYears <= year))]) * ns + s

        ## Recruitment
        if(spawning[s] > 0){
            ## Survivors from previous season/year
            if(s == 1){
                Ztmp <- hist$lastFAA + MAA  ## i.e. FAA in s=1 is equal to last FAA (pot s=4)
            }else{
                Ztmp <- ZAA
            }
            SSBtmp <- sum(NAAS * weighty * maty * exp(-pzbm * Ztmp)) ## pre-recruitment mort
            SSBPR0 <- getSSBPR(Ztmp, maty, weighty, fecun=1, asmax, ns, spawning)
            rec[y,s] <- spawning[s] * recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSBtmp, R0 = R0y,
                                              method = dat$SR, bp = dat$bp,
                                              beta = dat$recBeta, gamma = dat$recGamma)
            rec[rec<0] <- 1e-10
            NAAS[1] <- rec[y,s] * eR
        }

        ## Assessment
        if(year %in% assessYears && s == assessmentTiming){

            ## Estimate TAC
            if(tacID2 %in% c("refFmsy","consF")){
                ## Reference rule Fmsy
                tacs <- gettacs(tacs = tacs, id = hcr, TAC=NA)
            }else{
                ## True stock status
                TSBtmp <- sum(NAAS * weighty)
                bbmsy <- TSBtmp / refs$Bmsy
                FMtmp <- as.numeric(t(Fs))
                ffmsy <- sum(tail(FMtmp[1:((y*ns)-ns-(s-1))],ns)) / refs$Fmsy
                ## TAC
                tacs <- estTAC(obs = obs,
                               hcr = hcr,
                               tacs = tacs,
                               pars = list("ffmsy" = ffmsy,
                                           "bbmsy" = bbmsy))
                TACs[y] <- as.numeric(as.character(tacs$TAC[nrow(tacs)])) * eImp
                TACreal <- rep(TACs[y] / ntac, ntac)  ## CHECK: do not equally split tac among time steps!
            }

            ## Find F given TAC
            if(tacID2 == "refFmsy"){
                if(tacID == "refFmsy"){
                    ## Fishing at Fmsy
                    FMtac <- refs$Fmsy / ns
                }else{
                    fraci <- as.numeric(unlist(strsplit(as.character(tacID), "_"))[2])
                    ## Fishing at fraction of Fmsy
                    FMtac <- (fraci * refs$Fmsy) / ns
                }
            }else if(tacID2 == "noF"){
                ## noF
                FMtac <- 0
            }else if(tacID2 == "consF"){
                ## constant F
                fraci <- as.numeric(unlist(strsplit(as.character(tacID), "_"))[2])
                FMtac <- fraci / ns
            }else{
                ## any other HCR
                FMtac <- min(set$maxF/ns,
                             getFM(TACs[y],
                                    NAA = NAAS, MAA = MAA,
                                    sels = sely, weights = weighty,
                                    seasons = assessSeasons, ns = ns, y = y,
                                    h = hy, asmax = asmax, mats = maty,
                                    pzbm = pzbm, spawning = spawning,
                                    R0 = R0y, SR = dat$SR, bp = dat$pb, recBeta = dat$recBeta,
                                    recGamma = dat$recGamma, eR = eR,
                                    lastFM = Fs[y-1,s]
                                    )
                             )
            }
            Fs[y,] <- FMtac
        }

        ## Population dynamics
        FAA <- Fs[y,] * sely
        ZAA <- MAA + FAA
        CAA[,s] <- baranov(FAA, MAA, NAAS)
        CW[y,s] <- sum(CAA[,s] * weightFy)
        ## can't take more than what's there
        Btmp <- sum(NAAS * weighty * sely * exp(-MAA/2))
        if(CW[y,s] > 0.99 * Btmp){
            print("NOW")
            Fs[y,s] <- min(getFM(0.75 * Btmp,
                                         NAA = NAAS, MAA = MAA,
                                         sels = sely,
                                         weights = weighty,
                                         seasons = s, ns = ns, y = y,
                                         h = hy, asmax = asmax, mats = maty,
                                         pzbm = pzbm, spawning = spawning,
                                         R0 = R0y, SR = dat$SR, bp = dat$pb, recBeta = dat$recBeta,
                                         recGamma = dat$recGamma, eR = eR,
                                         lastFM = FM[y,s]),
                                  set$maxF/ns)
            FAA <- Fs[y,] * sely
            ZAA <- MAA + FAA
            CAA[,s] <- baranov(FAA, MAA, NAAS)
            CW[y,s] <- sum(CAA[,s] * weightFy)
            if(indtac < ntac){
                TACreal[indtac+1] <- TACreal[indtac+1] + TACreal[indtac] - CW[y,s]
            }else writeLines("Could not get full annual TAC.")
            TACreal[indtac] <- CW[y,s]
        }
        if(tacID2 == "refFmsy"){
            TACreal[indtac] <- sum(CAA[,s] * weightFs[,s], na.rm = TRUE)
            if(indtac == ntac) TACs[y] <- tacs$TAC[nrow(tacs)] <- sum(TACreal)
        }
        ## TSB
        TSB[y,s] <- sum(NAAS * weighty)
        ## ESB
        ESB[y,s] <- sum(NAAS * weighty * sely)
        ## SSB
        SSB[y,s] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))

        ## index observations
        if(s %in% idxS){
            idxi <- which(idxS == s)
            for(i in 1:length(idxi)){
                ## survey observation: total catch in weight (spict)
                surveyTime <- set$surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                NAAsurv <- exp(log(NAAS) - ZAA * surveyTime)
                ESBsurv <- NAAsurv * dat$weightFs[,s] * dat$sels[,s]
                obs$obsI[[idxi[i]]] <-
                    c(obs$obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * eI[[idxi[i]]])
                if(is.null(obs$timeI[[idxi[i]]])){
                    timeIi <- ny-nyhist+1
                }else timeIi <- floor(tail(obs$timeI[[idxi[i]]],1))
                obs$timeI[[idxi[i]]] <-
                    c(obs$timeI[[idxi[i]]], timeIi + 1 + set$surveyTimes[idxi[i]])
                ## survey observation at age
                obs$obsIA[[idxi[i]]] <- rbind(obs$obsIA[[idxi[i]]],
                                              q[idxi[i]] *
                                              as.numeric(by(NAAsurv * sely, asvec, sum)) * eImv[[idxi[[i]]]])
                rownames(obs$obsIA[[idxi[i]]])[nrow(obs$obsIA[[idxi[i]]])] <- as.character(y)  ## CHECK: instead of 1 (assuming one survey a year) this should allow for s surveys
                attributes(obs$obsIA[[i]]) <- c(attributes(obs$obsIA[[i]]), list(time = set$surveyTimes))
            }
        }

        ## Eponential decay
        NAAS <- NAAS * exp(-ZAA)
        if(s == ns){
            ## end of year biomasses
            TSBfinal[y] <- sum(NAAS * weighty)
            ESBfinal[y] <- sum(NAAS * weighty * sely)
            SSBfinal[y] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))
        }
        ## Continuous ageing
        NAAS[asmax] <- NAAS[asmax] + NAAS[asmax-1]
        for(as in (asmax-1):2) NAAS[as] <- NAAS[as-1]
        NAAS[1] <- 0
    }

    ## catch observations
    if(ns > 1){
        if(nsC == 1){
            obs$obsC <- c(obs$obsC, sum(CW[y,]) * eC)
            if(!is.null(obs$timeC)) timeCi <- tail(obs$timeC,1) else timeCi <- ny-nyhist+1
            obs$timeC <- c(obs$timeC, timeCi + 1)
            ## catch observations at age (SAM, SMS)
            obs$obsCA <- rbind(obs$obsCA, as.numeric(by(apply(CAA, 1, sum) * eCmv,asvec, sum)))
        }else if(nsC == ns){
            obs$obsC <- c(obs$obsC, CW[y,] * eC)
            if(!is.null(obs$timeC))
                timeCi <- floor(tail(obs$timeC,1)) else timeCi <- ny-nyhist+1
            obs$timeC <- c(obs$timeC, timeCi + 1 + rep(seasonStart, ny))
            ## catch observations at age (SAM, SMS)
            obs$obsCA <- rbind(obs$obsCA, as.numeric(by(apply(CAA, 1, sum) * eCmv,asvec, sum)))
        }else{
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
    }else{
        if(nsC > 1) writeLines("Set dat$ns to > 1 for seasonal catches. Generating annual catches!")
        obs$obsC <- c(obs$obsC, sum(CW[y,]) * eC)
        if(!is.null(obs$timeC)) timeCi <- tail(obs$timeC,1) else timeCi <- ny-nyhist+1
        obs$timeC <- c(obs$timeC, timeCi + 1)
        ## catch observations at age (SAM, SMS)
        obs$obsCA <- rbind(obs$obsCA, as.numeric(by(CAA * eCmv, asvec, sum)))
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
    out$lastNAA <- NAAS
    out$lastFAA <- FAA
    out$TSB <- TSB
    out$TSBfinal <- TSBfinal
    out$SSBfinal <- SSBfinal
    out$ESBfinal <- ESBfinal
    out$TACprev <- TACreal
    out$ESB <- ESB
    out$SSB <- SSB
    out$rec <- rec
    out$Fs <- Fs
    out$CW <- CW
    out$TACs <- TACs
    out$tacs <- tacs
    out$errs <- errs
    out$obs <- obs
    out$tacs <- tacs
    return(out)
}
