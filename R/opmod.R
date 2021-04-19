## input data:
## Number of years (ny)
## Maximum age of species (amax)
## Natural mortality (M)
## weight-at-age start (weightS)
## weight-at-age mid (weightH)
## maturity-at-age (mat)
## proportion of Z before maturation (pzbm)

#' @name initpop
#'
#' @export
initpop <- function(dat, set = NULL, out.opt = 1, verbose = TRUE){## indices
    if(is.null(set)) set <- check.set()
    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    yvec <- dat$yvec
    svec <- dat$svec
    s1vec <- dat$s1vec
    s1avec <- dat$s1avec
    as2a <- dat$as2a
    indage0 <- dat$indage0
    nsC <- dat$catchSeasons
    nyC <- dat$nyC
    if(is.numeric(nyC)){
        if(nyC > ny){
            nyC <- ny
            if(verbose) writeLines("Period for catch observations ('dat$nyC') is larger than historical period ('dat$ny'). Setting nyC == ny. ")
        }
        idxC <- (ny - nyC + 1):ny
    }
    nyI <- dat$nyI
    for(i in 1:length(nyI)){
        if(is.numeric(nyI) && any(nyI > ny) && verbose) writeLines("Period for index observations ('dat$nyI') is larger than historical period ('dat$ny'). Setting nyI == ny. ")
        if(is.numeric(nyI) && nyI[i] > ny){
            nyI[i] <- ny
        }
    }
    surveyTimes <- dat$surveyTimes
    nsurv <- length(surveyTimes)
    ## closest season
    seasonStart <- seq(0,1-1/ns,1/ns)
    if(is.numeric(nyI)){
        idxS <- rep(0, nsurv)
        for(i in 1:nsurv){
            tmp <- seasonStart[seasonStart <= surveyTimes[i]]
            idxS[i] <- which.min((tmp - surveyTimes[i])^2)
        }
    }
    ## effort
    nyE <- dat$nyE
    nsE <- dat$effortSeasons
    if(is.numeric(nyE)){
        if(nyE > ny){
            nyE <- ny
            if(verbose) writeLines("Period for effort observations ('dat$nyE') is larger than historical period ('dat$ny'). Setting nyE == ny. ")
        }
        idxE <- (ny - nyE + 1):ny
    }

    ## species data
    amax <- dat$amax + 1  ## age 0
    asmax <- amax * ns ## continuous max age
    M <- dat$M
    FM <- dat$FM
    weight <- dat$weight
    weightF <- dat$weightF
    mat <- dat$mat
    sel <- dat$sel
    Msel <- dat$Msel
    pzbm <- dat$pzbm
    R0 <- dat$R0
    h <- dat$h
    initN <- dat$initN
    q <- dat$q
    if(length(q) < nsurv) q <- rep(q, nsurv)
    qE <- dat$qE
    spawning <- dat$spawning


    ## errors
    if(!is.null(set)){
        eF <- set$eF
        eR <- set$eR
        eM <- set$eM
        eH <- set$eH
        eR0 <- set$eR0
        eMat <- set$eMat
        eSel <- set$eSel
        eW <- set$eW
        eImp <- set$eImp
        eC <- set$eC
        eI <- set$eI
        eCmv <- set$eCmv
        eImv <- set$eImv
        eE <- set$eE
        if(is.null(eF)) eF <- gen.noise(ny, set$noiseF[1], set$noiseF[2], bias.cor = set$noiseF[3])
        if(is.null(eR)) eR <- gen.noise(ny, set$noiseR[1], set$noiseR[2], bias.cor = set$noiseR[3])
        if(is.null(eM)) eM <- gen.noise(ny, set$noiseM[1], set$noiseM[2], bias.cor = set$noiseM[3])
        if(is.null(eH)) eH <- gen.noise(ny, set$noiseH[1], set$noiseH[2], bias.cor = set$noiseH[3])
        if(is.null(eR0)) eR0 <- gen.noise(ny, set$noiseR0[1], set$noiseR0[2], bias.cor = set$noiseR0[3])
        if(is.null(eMat)) eMat <- gen.noise(ny, set$noiseMat[1], set$noiseMat[2], bias.cor = set$noiseMat[3])
        if(is.null(eSel)) eSel <- gen.noise(ny, set$noiseSel[1], set$noiseSel[2], bias.cor = set$noiseSel[3])
        if(is.null(eW)) eW <- gen.noise(ny, set$noiseW[1], set$noiseW[2], bias.cor = set$noiseW[3])
        if(is.null(eImp)) eImp <- gen.noise(ny, set$noiseImp[1], set$noiseImp[2], bias.cor = set$noiseImp[3])
        if(is.null(eC)) eC <- gen.noise(ny, set$noiseC[1], set$noiseC[2], bias.cor = set$noiseC[3])
        if(is.null(eI)){
            eI <- list()
            for(i in 1:nsurv){
                eI[[i]] <- gen.noise(ny, set$noiseI[1], set$noiseI[2], bias.cor = set$noiseI[3])
            }
        }
        if(is.null(eCmv)) eCmv <- gen.noise(ny, set$noiseCmv[1], set$noiseCmv[2],
                                               bias.cor = set$noiseCmv[3],
                                               mv = TRUE, dat = dat)
        if(is.null(eImv)){
            eImv <- list()
            for(i in 1:nsurv){
                eImv[[i]] <- gen.noise(ny, set$noiseImv[1], set$noiseImv[2],
                                    bias.cor = set$noiseImv[3],
                                    mv = TRUE, dat = dat)
            }
        }
        if(is.null(eE)) eE <- gen.noise(ny, set$noiseE[1], set$noiseE[2], bias.cor = set$noiseE[3])
    }else{
        eF <- rep(1, ny)
        eR <- rep(1, ny)
        eM <- rep(1, ny)
        eH <- rep(1, ny)
        eR0 <- rep(1, ny)
        eMat <- rep(1, ny)
        eSel <- rep(1, ny)
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
        eE <- rep(1, ny)
    }
    errs <- list(eF = eF,
                 eR = eR,
                 eM = eM,
                 eH = eH,
                 eR0 = eR0,
                 eMat = eMat,
                 eSel = eSel,
                 eW = eW,
                 eImp = eImp,
                 eC = eC,
                 eI = eI,
                 eCmv = eCmv,
                 eImv = eImv,
                 eE = eE)

    ## Flags
    mselFlag <- inherits(Msel, "list") && length(Msel) > 1
    selFlag <- inherits(sel, "list") && length(sel) > 1

    ## containers
    TSB <- TSB1plus <- ESB <- SSB <- CW <- rec <- matrix(0, nrow=ny, ncol=ns)
    TACs <- TSBfinal <- SSBfinal <- ESBfinal <- rep(NA, ny)
    CAA <- vector("list", ny)
    for(i in 1:ny) CAA[[i]] <- matrix(NA, nrow = asmax, ncol = ns)
    ## observations
    if(all(!is.na(dat$surveyTimes)) && is.numeric(nyI)){
        obsI <- timeI <- vector("list", nsurv)
    }else{
        obsI <- timeI <- NULL
    }
    ## observations at age
    obsMAA <- obsCA <- matrix(0, nrow = ny, ncol = amax)
    obsCA <- matrix(0, nrow = length(idxC), ncol = amax)
    obsIA <- vector("list", nsurv)
    for(i in 1:nsurv) obsIA[[i]] <- matrix(0, nrow = ny, ncol = amax)

    ## Initialise NAA
    hy <- h * eH[1]
    R0y <- R0 * eR0[1]
    maty <- as.numeric(t(mat)) * eMat[1]
    weighty <- as.numeric(t(weight)) * eW[1]
    weightFy <- as.numeric(t(weightF)) * eW[1]
    sely <- as.numeric(t(sel[[1]])) * eSel[1]
    msely <- as.numeric(t(Msel[[1]]))
    NAAbi2 <- NAAbi <- matrix(0, asmax, ns)
    NAAbi[indage0,] <- R0y * spawning ## * exp(initN[1])
    Mbi <- M[1,] * eM[1] * msely
    Fbi <- FM[1,] * eF[1] * sely ## TODO: F,M different in various seasons?
    ZAA <-  Mbi + Fbi

    NAASbi <- initdistR(Mbi, Fbi, ns, asmax, indage0, spawning, R0y)

    ## Burn-in period
    if(is.null(set)) burnin <- 5e2 else burnin <- set$burnin
    if(is.numeric(burnin) && burnin > 0){
        for(y in 1:burnin){
            Fbi <- FM[1,] * eF[1] * sely
            Mbi <- M[1,] * eM[1] * msely
            Zbi <- Mbi + Fbi
            for(s in 1:ns){
                ## recruitment
                if(spawning[s] > 0){
                    SSBtmp <- sum(NAASbi * weighty * maty * exp(-pzbm * Zbi))
                    SSB0 <- get.ssb0(Mbi, maty, weighty, fecun = 1, asmax, ns,
                                     spawning, indage0 = indage0,
                                     R0 = R0y, season = s, FM = 0)
  ##                  print(paste0("SSB0: ",round(SSB0), " - SSBt: ",round(SSBtmp)))
                    recbi <- spawning[s] * recfunc(h = hy, SSBPR0 = SSB0/R0y, SSB = SSBtmp,
                                                   R0 = R0y, method = dat$SR, bp = dat$bp,
                                                   beta = dat$recBeta, gamma = dat$recGamma)
                    recbi[recbi<0] <- 1e-10
                    NAASbi[indage0] <- recbi
                }
                ## decay model
                NAASbi <- NAASbi * exp(-Zbi)
                ## ageing by season or year
                NAASbi[asmax] <- NAASbi[asmax] + NAASbi[asmax-1]
                for(as in (asmax-1):2) NAASbi[as] <- NAASbi[as-1]
                NAASbi[indage0] <- 0
            }
        }
    }
    NAAS <- NAASbi


    ## main loop
    for(y in 1:ny){

        ## Adding noise
        selInd <- ifelse(selFlag, y, 1)
        sely <- as.numeric(t(sel[[selInd]])) * eSel[y]
        FAA <- FM[y,] * eF[y] * sely
        mselInd <- ifelse(mselFlag, y, 1)
        MAA <- M[y,] * as.numeric(t(Msel[[mselInd]])) * eM[y]
        ZAA <- MAA + FAA
        maty <- as.numeric(t(mat)) * eMat[y]
        weighty <- as.numeric(t(weight)) * eW[y]
        weightFy <- as.numeric(t(weightF)) * eW[y]
        hy <- h * eH[y]
        R0y <- R0 * eR0[y]

        ## seasons
        for(s in 1:ns){
            ## recruitment
            if(spawning[s] > 0){
                ## Survivors from previous season/year
                SSB[y,s] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA)) ## pre-recruitment mort
##                print(SSB[y,s])
                SSB0 <- get.ssb0(MAA, maty, weighty, fecun=1, asmax,
                                 ns, spawning, indage0 = indage0, R0=R0y,
                                 season = s)
                rec[y,s] <-  spawning[s] * recfunc(h = hy, SSBPR0 = SSB0/R0y, SSB = SSB[y,s],
                                                   R0 = R0y, method = dat$SR, bp = dat$bp,
                                                   beta = dat$recBeta, gamma = dat$recGamma) * eR[y]
                rec[rec<0] <- 1e-10
                NAAS[indage0] <- rec[y,s]
            }

            ## catch
            CAA[[y]][,s] <- baranov(FAA, MAA, NAAS)
            CW[y,s] <- sum(CAA[[y]][,s] * weightFy)

            ## can't take more than what's there
            ## Btmp <- sum(NAAS * weighty * sely * exp(-MAA/2))
            ## if(CW[y,s] > 0.99 * Btmp){
            ##     FM[y,s] <- min(get.f(0.75 * Btmp,
            ##                                  NAA = NAAS, MAA = MAA,
            ##                                  sel = sely,
            ##                                  weight = weighty,
            ##                                  seasons = s, ns = ns, y = y,
            ##                                  h = hy, asmax = asmax, mat = maty,
            ##                                  pzbm = pzbm, spawning = spawning,
            ##                                  R0 = R0y, SR = dat$SR, bp = dat$pb, recBeta = dat$recBeta,
            ##                                  recGamma = dat$recGamma, eR = eR, indage0 = indage0,
            ##                                  lastFM = FM[y,s]),
            ##                           set$maxF/ns)
            ##     FAA <- FM[y,] * sely
            ##     ZAA <- MAA + FAA
            ##     CAA[[y]][,s] <- baranov(FAA, MAA, NAAS)
            ##     CW[y,s] <- sum(weightFy * CAA[[y]][,s])
            ## }
            ## TSB
            TSB[y,s] <- sum(NAAS * weighty)
            TSB1plus[y,s] <- sum(NAAS[-1] * weighty[-1])
            ## ESB
            ESB[y,s] <- sum(NAAS * weighty * sely)
            ## SSB
            SSB[y,s] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))

            ## index observations
            if(is.numeric(nyI)){
                if(s %in% idxS){
                    idxi <- which(idxS == s)
                    for(i in 1:length(idxi)){
                        surveyTime <- surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                        NAAsurv <- exp(log(NAAS) - ZAA * surveyTime)
                        ESBsurv <- NAAsurv * weightFy * as.numeric(t(dat$selI[[idxi[i]]]))
                        ## Total index (for spict)
                        obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * eI[[idxi[i]]][y])
                        if(is.null(timeI[[idxi[i]]]))
                            timeIi <- 0 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                        timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + surveyTimes[idxi[i]])
                        ## Index by age (sam, sms)
                        obsIA[[idxi[i]]][y,] <- q[idxi[i]] * as.numeric(by(NAAsurv *
                                                                           as.numeric(t(dat$selI[[idxi[i]]])),
                                                                           as2a, sum)) * eImv[[idxi[[i]]]][y,]
                    }
                }
            }
            ## observing annual M
            obsMAA[y,] <- obsMAA[y,] + MAA[seq(s,asmax,ns)]

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
            NAAS[indage0] <- 0
        }
    }

    ## account for nyhist (nyC, nyI)
    if(is.numeric(nyI) && all(!is.na(dat$surveyTimes))){
        for(i in 1:nsurv){
            timeI[[i]] <- timeI[[i]][(ny-nyI[i]+1):ny]
            obsI[[i]] <- obsI[[i]][(ny-nyI[i]+1):ny]
            obsIA[[i]] <- obsIA[[i]][(ny-nyI[i]+1):ny,]
            rownames(obsIA[[i]]) <- as.character((ny-nyI[i]+1):ny)
            colnames(obsIA[[i]]) <- as.character(0:dat$amax)
            attributes(obsIA[[i]]) <- c(attributes(obsIA[[i]]), list(time = dat$surveyTimes))
        }
    }
    obsMAA <- obsMAA[(ny-nyC+1):ny,]
    rownames(obsMAA) <- as.character((ny-nyC+1):ny)
    colnames(obsMAA) <- as.character(0:dat$amax)


    ## catch & effort observations
    ## ----------------
    if(ns > 1){
        ## Catch - seasonal OM
        ## ----------------
        if(nsC == 1){
            ## annual catches
            ## ----------------
            ## total catches (spict)
            obsC <- apply(CW[idxC,], 1, sum) * eC[idxC]
            timeC <- idxC
            ## catch observations at age (SAM, SMS)
            for(y in 1:length(idxC)) obsCA[y,] <- as.numeric(by(apply(CAA[[idxC[y]]], 1, sum) *
                                                                eCmv[y,],as2a,sum))
        }else if(nsC == ns){
            ## seasonal catches (for each season in OM)
            ## ----------------
            obsC <- as.numeric(t(CW[idxC,] * eC[idxC]))
            timeC <- rep(idxC, each = ns) + rep(seasonStart, nyC)
            ## catch observations at age (seasonal CAA observations not yet implemented, what is target format?)
            for(y in 1:length(idxC)) obsCA[y,] <- as.numeric(by(apply(CAA[[idxC[y]]], 1, sum) *
                                                                eCmv[y,],as2a,sum))
        }else{
            ## seasonal catches not aligning with seasones in OM (not implemented)
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
        ## Effort - seasonal OM
        ## ----------------
        if(is.numeric(nyE)){
            if(nsE == 1){
                ## annual effort
                ## ----------------
                obsE <- apply(FM[idxE,], 1, sum) * eF[idxE] * qE * eE[idxE]
                timeE <- idxE
            }else if(nsE == ns){
                ## seasonal effort (for each season in OM)
                ## ----------------
                obsE <- as.numeric(t(FM[idxE,] * eF[idxE] * qE * eE[idxE]))
                timeE <- rep(idxE, each = ns) + rep(seasonStart, nyE)
            }else{
                ## seasonal effort not aligning with seasones in OM (not implemented)
                stop("Effort observation seasons and operating model seasons do not match. Not yet implemented!")
            }
        }else{
            obsE <- NULL
            timeE <- NULL
        }
    }else{
        ## annual OM
        ## Catch
        if(nsC > 1) writeLines("Set dat$ns to > 1 for seasonal catches. Generating annual catches!")
        obsC <- CW[idxC,] * eC[idxC]
        timeC <- idxC
        ## catch observations at age
        for(y in 1:length(idxC)) obsCA[y,] <- as.numeric(by(CAA[[idxC[y]]] * eCmv[y,],as2a,sum))

        ## Effort
        if(is.numeric(nyE)){
            if(nsE > 1) writeLines("Set dat$ns to > 1 for seasonal effort. Generating annual catches!")
            obsE <- FM[idxE,] * eF[idxE] * qE * eE[idxE]
            timeE <- idxE
        }else{
            obsE <- NULL
            timeE <- NULL
        }
    }
    rownames(obsCA) = as.character((ny-nyC+1):ny)
    colnames(obsCA) <- as.character(0:dat$amax)

    ## other observations (required by SAM)
    ## ----------------
    ## natural mortality
    ## obsMAA <- as.matrix(dat$M[idxC]) %*% t(dat$Msel[[1]])
    ## rownames(obsMAA) = as.character((ny-nyC+1):ny)
    ## colnames(obsMAA) = as.character(0:dat$amax)

    ## proportion mature
    propMature <- matrix(rowSums(dat$mat), length(idxC), amax, byrow = TRUE)
    rownames(propMature) = as.character((ny-nyC+1):ny)
    colnames(propMature) = as.character(0:dat$amax)

    ## mean stock weight
    WAAs <- matrix(rowSums(dat$weight), length(idxC), amax, byrow = TRUE)
    rownames(WAAs) = as.character((ny-nyC+1):ny)
    colnames(WAAs) = as.character(0:dat$amax)

    ## mean catch weight
    WAAc <- matrix(rowSums(dat$weightF), length(idxC), amax, byrow = TRUE)
    rownames(WAAc) = as.character((ny-nyC+1):ny)
    colnames(WAAc) = as.character(0:dat$amax)

    ## proportion females
    propFemale <- matrix(0.0, length(idxC), amax)
    rownames(propFemale) <- as.character((ny-nyC+1):ny)
    colnames(propFemale) <- as.character(0:dat$amax)

    ## landing fraction
    landFrac = matrix(1.0, length(idxC), amax)
    rownames(landFrac) <- as.character((ny-nyC+1):ny)
    colnames(landFrac) <- as.character(0:dat$amax)

    obs <- list(obsC = obsC,
                timeC = timeC,
                obsI = obsI,
                timeI = timeI,
                obsE = obsE,
                timeE = timeE,
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
        out$FM <- FM * eF
        out$CW <- CW
        out$TACs <- TACs
        out$errs <- errs
        out$obs <- obs
    }else if(out.opt %in% c(2,3)){
        refs <- dat$ref
        if(is.null(refs)){
            warning("The reference points are not part of dat! Use est.ref.levels to estimate them")
        }else if(out.opt == 2){
            out <- TSBfinal[ny]/refs[[dat$depl.quant]][ny]
        }else if(out.opt == 3){
            out <- SSBfinal[ny]/refs[[dat$depl.quant]][ny]
        }
    }else if(out.opt == 4){
        out <- TSBfinal[ny]
    }
    return(out)
}



## advancepopulation
#' @name advancepop
#' @export
advancepop <- function(dat, hist, set, hcr, year, verbose = TRUE){

    ## indices
    ny <- nrow(hist$TSB)
    y <- ny + 1
    ysim <- y - dat$ny
    ns <- dat$ns
    s1vec <- dat$s1vec
    as2a <- dat$as2a
    indage0 <- dat$indage0
    nt <- ny * ns
    nyC <- dat$nyC
    nsC <- dat$catchSeasons
    nysim <- set$nysim
    assessYears <- seq(1, nysim, set$assessmentInterval)
    nyE <- dat$nyE
    nsE <- dat$effortSeasons

    ## survey
    nsurv <- length(dat$surveyTimes)
    seasonStart <- seq(0,1-1/ns,1/ns)
    idxS <- rep(0, nsurv)
    if(all(!is.na(dat$surveyTimes))){
        for(i in 1:nsurv){
            tmp <- seasonStart[seasonStart < dat$surveyTimes[i]]
            idxS[i] <- which.min((tmp - dat$surveyTimes[i])^2)
        }
    }

    ## parameters
    tacs <- hist$tacs
    obs <- hist$obs
    refs <- hist$ref
    amax <- dat$amax + 1  ## age 0
    asmax <- amax * ns
    pzbm <- dat$pzbm
    q <- dat$q
    if(length(q) < nsurv) q <- rep(q, nsurv)
    qE <- dat$qE
    tacID <- hcr
    tacID2 <- unlist(strsplit(as.character(tacID), "_"))[1]
    M <- dat$M
    spawning <- dat$spawning
    assessmentTiming <- set$assessmentTiming
    if(ns == 1){
        assessSeasons <- as.list(1)
    }else{
        if(length(assessmentTiming) == 1){
            assessSeasons <- list(rep(1:ns, 20)[assessmentTiming:(assessmentTiming+ns-1)])
        }else{
            if(verbose) writeLines("You specified several assessmentTimings. I hope you know what you are doing!")
            assessSeasons <- as.list(assessmentTiming)
        }

    }

    bbmsySD <- set$bbmsySD
    if(is.null(bbmsySD) || !is.numeric(bbmsySD)) bbmsySD <- 0
    bbmsyBias <- set$bbmsyBias
    if(is.null(bbmsyBias) || !is.numeric(bbmsyBias)) bbmsyBias <- 0
    ffmsySD <- set$ffmsySD
    if(is.null(ffmsySD) || !is.numeric(ffmsySD)) ffmsySD <- 0
    ffmsyBias <- set$ffmsyBias
    if(is.null(ffmsyBias) || !is.numeric(ffmsyBias)) ffmsyBias <- 0
    fmsyBias <- set$fmsyBias
    if(is.null(fmsyBias) || !is.numeric(fmsyBias)) fmsyBias <- 0
    tacSD <- set$tacSD
    if(is.null(tacSD) || !is.numeric(tacSD)) tacSD <- 0

    ## parameters per age
    weight <- dat$weight
    weightF <- dat$weightF
    mat <- dat$mat
    sel <- dat$sel
    Msel <- dat$Msel

    ## errors
    ## TODO: put this in external wrapper function
    eF <- set$eF[ysim]
    eR <- set$eR[ysim]
    eM <- set$eM[ysim]
    eH <- set$eH[ysim]
    eR0 <- set$eR0[ysim]
    eMat <- set$eMat[ysim]
    eSel <- set$eSel[ysim]
    eW <- set$eW[ysim]
    eImp <- set$eImp[ysim]
    eC <- set$eC[ysim]
    eI <- list()
    for(i in 1:nsurv){
        eI[[i]] <- set$eI[[i]][ysim]
    }
    eCmv <- set$eCmv[ysim,]
    eImv <- set$eImv[ysim,]
    eE <- set$eE[ysim]
    if(is.null(eF)) eF <- gen.noise(1, set$noiseF[1], set$noiseF[2], bias.cor = set$noiseF[3])
    if(is.null(eR)) eR <- gen.noise(1, set$noiseR[1], set$noiseR[2], bias.cor = set$noiseR[3])
    if(is.null(eM)) eM <- gen.noise(1, set$noiseM[1], set$noiseM[2], bias.cor = set$noiseM[3])
    if(is.null(eH)) eH <- gen.noise(1, set$noiseH[1], set$noiseH[2], bias.cor = set$noiseH[3])
    if(is.null(eR0)) eR0 <- gen.noise(1, set$noiseR0[1], set$noiseR0[2], bias.cor = set$noiseR0[3])
    if(is.null(eMat)) eMat <- gen.noise(1, set$noiseMat[1], set$noiseMat[2], bias.cor = set$noiseMat[3])
    if(is.null(eSel)) eSel <- gen.noise(1, set$noiseSel[1], set$noiseSel[2], bias.cor = set$noiseSel[3])
    if(is.null(eW)) eW <- gen.noise(1, set$noiseW[1], set$noiseW[2], bias.cor = set$noiseW[3])
    if(is.null(eImp)) eImp <- gen.noise(1, set$noiseImp[1], set$noiseImp[2], bias.cor = set$noiseImp[3])
    if(is.null(eC)) eC <- gen.noise(1, set$noiseC[1], set$noiseC[2], bias.cor = set$noiseC[3])
    if(is.null(eI)){
        eI <- list()
        for(i in 1:nsurv){
            eI[[i]] <- gen.noise(1, set$noiseI[1], set$noiseI[2], bias.cor = set$noiseI[3])
        }
    }
    if(is.null(eCmv)) eCmv <- gen.noise(1, set$noiseCmv[1], set$noiseCmv[2],
                                       bias.cor = set$noiseCmv[3],
                                       mv = TRUE, dat = dat)
    if(is.null(eImv)){
        eImv <- list()
        for(i in 1:nsurv){
            eImv[[i]] <- gen.noise(1, set$noiseImv[1], set$noiseImv[2],
                                  bias.cor = set$noiseImv[3],
                                  mv = TRUE, dat = dat)
        }
    }
    if(is.null(eE)) eE <- gen.noise(1, set$noiseE[1], set$noiseE[2], bias.cor = set$noiseE[3])
    if("errs" %in% names(hist)){
        errs <- list(eF = c(hist$errs$eF, eF),
                     eR = c(hist$errs$eR, eR),
                     eM = c(hist$errs$eM, eM),
                     eH = c(hist$errs$eH, eH),
                     eR0 = c(hist$errs$eR0,eR0),
                     eMat = c(hist$errs$eMat, eMat),
                     eSel = c(hist$errs$eSel, eSel),
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
        errs$eE <- c(hist$errs$eE, eE)
    }else{
        errs <- list(eF = eF,
                     eR = eR,
                     eM = eM,
                     eH = eH,
                     eR0 = eR0,
                     eMat = eMat,
                     eSel = eSel,
                     eW = eW,
                     eImp = eImp,
                     eC = eC,
                     eI = eI,
                     eCmv = eCmv,
                     eImv = eImv,
                     eE = eE)
    }

    ## Flags
    mselFlag <- inherits(Msel, "list") && length(Msel) > 1
    selFlag <- inherits(sel, "list") && length(sel) > 1

    ## Containers
    tmp <- matrix(0, 1, ns)
    TSB <- rbind(hist$TSB,tmp)
    SSB <- rbind(hist$SSB,tmp)
    ESB <- rbind(hist$ESB,tmp)
    CW  <- rbind(hist$CW,tmp)
    FM  <- rbind(hist$FM,tmp)
    NAA <- obsMAAtmp <- rep(0, amax)
    FAA <- ZAA <- MAA <- CAA <- matrix(NA, asmax, ns)
    TACs <- c(hist$TACs, NA)
    rec <- rbind(hist$rec, tmp)
    TSBfinal <- c(hist$TSBfinal, NA)
    SSBfinal <- c(hist$SSBfinal, NA)
    ESBfinal <- c(hist$ESBfinal, NA)

    ## project forward
    R0y <- dat$R0 * eR0
    mselInd <- ifelse(mselFlag, y, 1)
##    MAA <- M[s1vec[y]:(s1vec[y]+ns-1)] * as.numeric(t(Msel[[mselInd]])) * eM
    MAA <- M[y,] * as.numeric(t(Msel[[mselInd]])) * eM
    hy <- dat$h * eH
    maty <- as.numeric(t(mat)) * eMat
    selInd <- ifelse(selFlag, y, 1)
    sely <- as.numeric(t(sel[[selInd]])) * eSel
    weighty <- as.numeric(t(weight)) * eW
    weightFy <- as.numeric(t(weightF)) * eW

    ## TAC from previous year (if e.g. interim advice)
    ntac <- ns * set$assessmentInterval
    TACprev <- hist$TACprev ## TAC for full interval
    if(is.null(TACprev)){
        TACprev <- tail(as.numeric(t(hist$CW)),ntac)
    }
    ## Use only part before assessment from previous TAC in assessment year
    if(year %in% assessYears){
        TACreal <- rep(NA, ntac)
        if(min(set$assessmentTiming) > 1){
            ind <- 1:(min(set$assessmentTiming)-1)
            TACreal[ind] <- tail(TACprev, length(ind))
        }
    }else{
        ## Use full previous TAC if no assessment
        TACreal <- TACprev
    }

    NAAS <- hist$lastNAA


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
            SSB0 <- get.ssb0(MAA, maty, weighty, fecun=1, asmax, ns,
                             spawning, indage0 = indage0,
                             R0 = R0y, season = s)
            rec[y,s] <- spawning[s] * recfunc(h = hy, SSBPR0 = SSB0/R0y, SSB = SSBtmp, R0 = R0y,
                                              method = dat$SR, bp = dat$bp,
                                              beta = dat$recBeta, gamma = dat$recGamma) * eR
            rec[rec<0] <- 1e-10
            NAAS[indage0] <- rec[y,s]
        }

        ## Assessment
        if(year %in% assessYears && s %in% assessmentTiming){

            ## Estimate TAC
            if(tacID2 %in% c("refFmsy","consF")){
                ## Reference rule Fmsy
                tacs <- gettacs(tacs = tacs, id = hcr, TAC=NA)
            }else{
                ## True stock status
                TSBtmp <- sum(NAAS * weighty)
                ESBtmp <- sum(NAAS * weighty * sely)
                bbmsy <- TSBtmp / refs$Bmsy[y + s - 1]
                FMtmp <- as.numeric(t(FM))
                ffmsy <- sum(tail(FMtmp[1:((y*ns)-ns-(s-1))],ns)) / refs$Fmsy[y + s - 1]
                ## TAC
                tacs <- est.tac(obs = obs,
                                hcr = hcr,
                                tacs = tacs,
                                pars = list("ffmsy" = ffmsy,
                                            "bbmsy" = bbmsy,
                                            "f" = FMtmp,
                                            "b" = TSBtmp,
                                            "fmsy" = refs$Fmsy[y + s - 1],
                                            "bmsy" = refs$Bmsy[y + s - 1],
                                            "sel" = sely,
                                            "weight" = weightFy,
                                            "m" = MAA,
                                            "n" = NAAS,
                                            "bbmsySD" = bbmsySD,
                                            "bbmsyBias" = bbmsyBias,
                                            "ffmsySD" = ffmsySD,
                                            "ffmsyBias" = ffmsyBias,
                                            "fmsyBias" = fmsyBias,
                                            "tacSD" = tacSD,
                                            "ns" = ns
                                            ))

                TACs[y] <- as.numeric(as.character(tacs$TAC[nrow(tacs)])) * eImp
                TACreal <- rep(TACs[y] / ntac, ntac)
                ## CHECK: do not equally split tac among time steps! (only used for TACprev)
            }

            ## Find F given TAC
            if(tacID2 == "refFmsy"){
                if(tacID == "refFmsy"){
                    ## Fishing at Fmsy
                    FMtac <- refs$Fmsy[y + s - 1] / ns
                }else{
                    fraci <- as.numeric(unlist(strsplit(as.character(tacID), "_"))[2])
                    ## Fishing at fraction of Fmsy
                    FMtac <- (fraci * refs$Fmsy[y]) / ns
                    ## ## percentile
                    ## fraci <- as.numeric(unlist(strsplit(as.character(tacID), "_"))[2])
                    ## ## uncertainty
                    ## fraci2 <- as.numeric(unlist(strsplit(as.character(tacID), "_"))[3])
                    ## ## bias
                    ## fraci3 <- as.numeric(unlist(strsplit(as.character(tacID), "_"))[4])
                    ## ## Fishing at fraction of Fmsy
                    ## ## FMtac <- (fraci * refs$Fmsy[y]) / ns
                    ## FMtac <- exp(qnorm(fraci,log((refs$Fmsy[y] * fraci3)/ns), sd = fraci2))

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

                ## CHECK: shouldn't this give a vector FM of length ns?
                ## what if F cannot be assumed to be constant throughout the year? seasonal fishing effort?
                ## but which pattern to use? from last year? or in sel?

                FMtac <- min(set$maxF/ns,
                             get.f(TACs[y],
                                   NAA = NAAS, MAA = MAA,
                                   sel = sely, weight = weighty,
                                   seasons = assessSeasons[[s]], ns = ns, y = y,
                                   h = hy, asmax = asmax, mat = maty,
                                   pzbm = pzbm, spawning = spawning,
                                   R0 = R0y, SR = dat$SR, bp = dat$pb, recBeta = dat$recBeta,
                                   recGamma = dat$recGamma, eR = eR,
                                   indage0 = indage0,
                                   lastFM = 0.01, ## FM[y-1,s],  ## CHECK: had some issues with too high FM
                                   upper = log(set$maxF/ns)
                                   )
                             )
            }
            FM[y,] <- FMtac
        }

        ## Population dynamics
        FAA <- FM[y,] * sely
        ZAA <- MAA + FAA
        CAA[,s] <- baranov(FAA, MAA, NAAS)
        CW[y,s] <- sum(CAA[,s] * weightFy)
        ## can't take more than what's there
        Btmp <- sum(NAAS * weighty * sely * exp(-MAA/2))
        ## if(CW[y,s] > 0.99 * Btmp){
        ##     FM[y,s] <- min(get.f(0.75 * Btmp,
        ##                                  NAA = NAAS, MAA = MAA,
        ##                                  sel = sely,
        ##                                  weight = weighty,
        ##                                  seasons = s, ns = ns, y = y,
        ##                                  h = hy, asmax = asmax, mat = maty,
        ##                                  pzbm = pzbm, spawning = spawning,
        ##                                  R0 = R0y, SR = dat$SR, bp = dat$pb, recBeta = dat$recBeta,
        ##                                  recGamma = dat$recGamma, eR = eR,
        ##                                  indage0 = indage0,
        ##                                  lastFM = FM[y,s]),
        ##                           set$maxF/ns)
        ##     FAA <- FM[y,] * sely
        ##     ZAA <- MAA + FAA
        ##     CAA[,s] <- baranov(FAA, MAA, NAAS)
        ##     CW[y,s] <- sum(CAA[,s] * weightFy)
        ##     if(indtac < ntac){
        ##         TACreal[indtac+1] <- TACreal[indtac+1] + TACreal[indtac] - CW[y,s]
        ##     }else writeLines("Could not get full annual TAC.")
        ##     TACreal[indtac] <- CW[y,s]
        ## }
        if(tacID2 == "refFmsy"){
            TACreal[indtac] <- sum(CAA[,s] * weightF[,s], na.rm = TRUE)
            if(indtac == ntac) TACs[y] <- tacs$TAC[nrow(tacs)] <- sum(TACreal)
        }
        ## TSB
        TSB[y,s] <- sum(NAAS * weighty)
        ## ESB
        ESB[y,s] <- sum(NAAS * weighty * sely)
        ## SSB
        SSB[y,s] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))

        ## index observations
        if(all(!is.na(dat$surveyTimes))){
            if(s %in% idxS){
                idxi <- which(idxS == s)
                for(i in 1:length(idxi)){
                    ## survey observation: total catch in weight (spict)
                    surveyTime <- dat$surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                    NAAsurv <- exp(log(NAAS) - ZAA * surveyTime)
                    ESBsurv <- NAAsurv * weightFy * as.numeric(t(dat$selI[[idxi[i]]]))
                    obs$obsI[[idxi[i]]] <-
                        c(obs$obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * eI[[idxi[i]]])
                    if(is.null(obs$timeI[[idxi[i]]])){
                        timeIi <- ny-nyI[idxi[i]]+1
                    }else timeIi <- floor(tail(obs$timeI[[idxi[i]]],1))
                    obs$timeI[[idxi[i]]] <-
                        c(obs$timeI[[idxi[i]]], timeIi + 1 + dat$surveyTimes[idxi[i]])
                    ## survey observation at age
                    obs$obsIA[[idxi[i]]] <- rbind(obs$obsIA[[idxi[i]]],
                                                  q[idxi[i]] *
                                                  as.numeric(by(NAAsurv * as.numeric(t(dat$selI[[idxi[i]]])), as2a,
                                                                sum)) *
                                                  eImv[[idxi[[i]]]])
                    rownames(obs$obsIA[[idxi[i]]])[nrow(obs$obsIA[[idxi[i]]])] <- as.character(y)  ## CHECK: instead of 1 (assuming one survey a year) this should allow for s surveys
                }
            }
            for(i in 1:length(idxS)) attributes(obs$obsIA[[i]]) <-
                                         c(attributes(obs$obsIA[[i]]),
                                           list(time = dat$surveyTimes))
        }
        ## observing annual M
        obsMAAtmp <- obsMAAtmp + MAA[seq(s,asmax,ns)]


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
        NAAS[indage0] <- 0
    }

    ## catch & effort observations
    if(ns > 1){
        ## Catch
        if(nsC == 1){
            obs$obsC <- c(obs$obsC, sum(CW[y,]) * eC)
            if(!is.null(obs$timeC)) timeCi <- tail(obs$timeC,1) else timeCi <- ny-nyC+1
            obs$timeC <- c(obs$timeC, timeCi + 1)
            ## catch observations at age (SAM, SMS)
            obs$obsCA <- rbind(obs$obsCA, as.numeric(by(apply(CAA, 1, sum) * eCmv, as2a, sum)))
        }else if(nsC == ns){
            obs$obsC <- c(obs$obsC, CW[y,] * eC)
            if(!is.null(obs$timeC))
                timeCi <- floor(tail(obs$timeC,1)) else timeCi <- ny-nyC+1
            obs$timeC <- c(obs$timeC, timeCi + 1 + rep(seasonStart, ny))
            ## catch observations at age (SAM, SMS)
            obs$obsCA <- rbind(obs$obsCA, as.numeric(by(apply(CAA, 1, sum) * eCmv, as2a, sum)))
        }else{
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
        ## Effort
        if(is.numeric(nyE)){
            if(nsE == 1){
                obs$obsE <- c(obs$obsE, sum(FM[y,]) * qE * eE)
                if(!is.null(obs$timeE)) timeEi <- tail(obs$timeE,1) else timeEi <- ny-nyE+1
                obs$timeE <- c(obs$timeE, timeEi + 1)
            }else if(nsE == ns){
                obs$obsE <- c(obs$obsE, FM[y,] * qE * eE)
                if(!is.null(obs$timeE))
                    timeEi <- floor(tail(obs$timeE,1)) else timeEi <- ny-nyE+1
                obs$timeE <- c(obs$timeE, timeEi + 1 + rep(seasonStart, ny))
            }else{
                stop("Effort observation seasons and operating model seasons do not match. Not yet implemented!")
            }
        }
    }else{
        ## annual obs
        ## Catch
        if(nsC > 1) writeLines("Set dat$ns to > 1 for seasonal catches. Generating annual catches!")
        obs$obsC <- c(obs$obsC, sum(CW[y,]) * eC)
        if(!is.null(obs$timeC)) timeCi <- tail(obs$timeC,1) else timeCi <- ny-nyC+1
        obs$timeC <- c(obs$timeC, timeCi + 1)
        ## catch observations at age (SAM, SMS)
        obs$obsCA <- rbind(obs$obsCA, as.numeric(by(CAA * eCmv, as2a, sum)))
        ## Effort
        if(is.numeric(nyE)){
            if(nsE > 1) writeLines("Set dat$ns to > 1 for seasonal effort data. Generating annual effort observations!")
            obs$obsE <- c(obs$obsE, sum(FM[y,]) * qE * eE)
            if(!is.null(obs$timeE)) timeEi <- tail(obs$timeE,1) else timeEi <- ny-nyE+1
            obs$timeE <- c(obs$timeE, timeEi + 1)
        }
    }
    rownames(obs$obsCA)[nrow(obs$obsCA)] <- as.character(y)

    ## natural mortality
    obs$obsMAA <- rbind(obs$obsMAA, y = obsMAAtmp)
    rownames(obs$obsMAA)[nrow(obs$obsMAA)] = as.character(y)

    ## proportion mature
    obs$propMature <- rbind(obs$propMature, matrix(rowSums(dat$mat), 1, amax, byrow = TRUE))
    rownames(obs$propMature)[nrow(obs$propMature)] = as.character(y)

    ## mean stock weight
    obs$WAAs <- rbind(obs$WAAs, matrix(rowSums(dat$weight), 1, amax, byrow = TRUE))
    rownames(obs$WAAs)[nrow(obs$WAAs)] = as.character(y)

    ## mean catch weight
    obs$WAAc <- rbind(obs$WAAc, matrix(rowSums(dat$weightF), 1, amax, byrow = TRUE))
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
    out$FM <- FM
    out$CW <- CW
    out$TACs <- TACs
    out$tacs <- tacs
    out$errs <- errs
    out$obs <- obs
    out$tacs <- tacs
    out$refs <- refs
    return(out)
}
