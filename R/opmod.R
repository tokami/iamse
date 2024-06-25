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
initpop <- function(dat, set = NULL, out.opt = 1, verbose = TRUE,
                    plot = FALSE, dbg = 0){## indices
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
    if(!inherits(weight,"matrix")){
        weight <- matrix(weight, asmax, ny)
    }
    weightF <- dat$weightF
    if(!inherits(weightF,"matrix")){
        weightF <- matrix(weightF, asmax, ny)
    }
    mat <- dat$mat
    sel <- dat$sel
    Msel <- dat$Msel
    pzbm <- dat$pzbm
    R0 <- check.par(dat$R0, ny)
    h <- check.par(dat$h, ny)
    alpha <- check.par(dat$recAlpha, ny)
    beta <- check.par(dat$recBeta, ny)
    initN <- dat$initN
    q <- dat$q
    if(length(q) < nsurv) q <- rep(q, nsurv)
    qE <- dat$qE
    spawning <- dat$spawning

    ## errors
    errs <- list()
    errs$time <- get.errs(dat, set, 1:dat$ny)
    errs$rep <- get.errs(dat, set, 1, rep = TRUE)

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
    obsIL <- vector("list", nsurv)
    for(i in 1:nsurv) obsIL[[i]] <- matrix(0, nrow = ny, ncol = length(dat$mids))

    ## Initialise NAA
    ## TODO: use mean for all!
    hy <- h[1] * errs$time$eH[1] * errs$rep$eH
    R0y <- R0[1] * mean(errs$time$eR0 * errs$rep$eR0) ## * errs$time$eR0[1] * errs$rep$eR0
    maty <- as.numeric(t(mat)) * errs$time$eMat[1] * errs$rep$eMat
    weighty <- as.numeric(t(weight[,1])) * errs$time$eW[1] * errs$rep$eW
    weightFy <- as.numeric(t(weightF[,1])) * errs$time$eW[1] * errs$rep$eW
    sely <- as.numeric(t(sel[[1]])) * errs$time$eSel[1] * errs$rep$eSel
    msely <- as.numeric(t(Msel[[1]]))
    NAAbi2 <- NAAbi <- matrix(0, asmax, ns)
    NAAbi[indage0,] <- R0y * spawning ## * exp(initN[1])
    Mbi <- M[1,dat$as2s] * msely * mean(errs$time$eM * errs$rep$eM) ## BEFORE: * errs$time$eM[1] * errs$rep$eM
    Fbi <- FM[1,dat$as2s] * sely  * mean(errs$time$eF * errs$rep$eF) ## TODO: F,M different in various seasons?
    ZAA <-  Mbi + Fbi
    alphay <- alpha[1] * errs$time$eAlpha[1] * errs$rep$eAlpha
    betay <- beta[1] * errs$time$eBeta[1] * errs$rep$eBeta

    NAASbi <- initdistR(Mbi, Fbi, ns, asmax, indage0, spawning, R0y * mean(errs$time$eR * errs$rep$eR))

    ## Burn-in period
    if(is.null(set)) burnin <- 5e2 else burnin <- set$burnin
    if(is.numeric(burnin) && burnin > 0){
        for(y in 1:burnin){
            Fbi <- FM[1,] * sely * mean(errs$time$eF * errs$rep$eF)
            Mbi <- M[1,] * msely * mean(errs$time$eM * errs$rep$eM)
            Zbi <- Mbi + Fbi
            for(s in 1:ns){
                if(spawning[s] > 0){
                    SSBtmp <- sum(NAASbi * weighty * maty * exp(-pzbm * Zbi))
                    SSB0 <- get.ssb0(Mbi, maty, weighty, fecun = 1, asmax, ns,
                                     spawning, R0 = R0y * mean(errs$time$eR * errs$rep$eR),
                                     indage0 = indage0,
                                     season = s, FM = 0)
                    ##                  print(paste0("SSB0: ",round(SSB0), " - SSBt: ",round(SSBtmp)))
                    recbi <- spawning[s] * recfunc(h = hy, SSBPR0 = SSB0/(R0y * mean(errs$time$eR * errs$rep$eR)),
                                                   SSB = SSBtmp,
                                                   R0 = R0y, method = dat$SR,
                                                   bp = dat$bp,
                                                   beta = betay,
                                                   gamma = dat$recGamma,
                                                   alpha = alphay) * mean(errs$time$eR * errs$rep$eR)
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

    if(dbg > 0){
        Mbi <- M[1,dat$as2s] * msely * mean(errs$time$eM * errs$rep$eM) ## BEFORE: * errs$time$eM[1] * errs$rep$eM
        Fbi <- FM[1,dat$as2s] * sely  * mean(errs$time$eF * errs$rep$eF) ## TODO: F,M different in various seasons?
        ## initdist  only correct if F=0, so okay for get.ssb0
        print(data.frame(initR = initdistR(Mbi, Fbi, dat$ns, dat$asmax,
                                           dat$indage0, dat$spawning,
                                           dat$R0 * mean(errs$time$eR * errs$rep$eR)),
                         initC = initdist(Mbi, Fbi,
                                          dat$R0 * mean(errs$time$eR * errs$rep$eR),
                                          dat$spawning,  dat$indage0),
                         burnIn = NAASbi))
    }

    SSB0vec <- rep(NA,ny)
    SPRvec <- rep(NA,ny)

    ## main loop
    for(y in 1:ny){

        ## Adding noise
        selInd <- ifelse(selFlag, y, 1)
        sely <- as.numeric(t(sel[[selInd]])) * errs$time$eSel[y] * errs$rep$eSel
        FAA <- FM[y,dat$as2s] * errs$time$eF[y] * errs$rep$eF * sely
        mselInd <- ifelse(mselFlag, y, 1)
        MAA <- M[y,dat$as2s] * as.numeric(t(Msel[[mselInd]])) * errs$time$eM[y] * errs$rep$eM
        ZAA <- MAA + FAA
        maty <- as.numeric(t(mat)) * errs$time$eMat[y] * errs$rep$eMat
        weighty <- as.numeric(t(weight[,y])) * errs$time$eW[y] * errs$rep$eW
        weightFy <- as.numeric(t(weightF[,y])) * errs$time$eW[y] * errs$rep$eW
        hy <- h[y] * errs$time$eH[y] * errs$rep$eH
        R0y <- R0[y] * errs$time$eR0[y] * errs$rep$eR0
        alphay <- alpha[y] * errs$time$eAlpha[1] * errs$rep$eAlpha
        betay <- beta[y] * errs$time$eBeta[1] * errs$rep$eBeta

        if(plot){
            par(mfrow = n2mfrow(ns),
                mar = c(2,2,1,1), oma = c(4,4,3,1))
        }


        ## seasons
        for(s in 1:ns){

            if(s == 1 && set$recordLast == 0){
                ## beginning of year biomasses
                TSBfinal[y] <- sum(NAAS * weighty)
                ESBfinal[y] <- sum(NAAS * weighty * sely)
                SSBfinal[y] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))
            }

            if(plot) plot(as.vector(NAAS %*% dat$plba), ty= 'b',
                          xlab = "", ylab = "")

            ## index observations
            if(is.numeric(nyI) && set$surveyBeforeRec){
                if(s %in% idxS){
                    idxi <- which(idxS == s)
                    for(i in 1:length(idxi)){
                        surveyTime <- surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                        NAAsurv <- exp(log(NAAS) - ZAA * surveyTime)
                        ESBsurv <- NAAsurv * weightFy * as.numeric(t(dat$selI[[idxi[i]]]))
                        ## Total index (for spict)
                        obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * errs$time$eI[[idxi[i]]][y] * errs$rep$eI[[idxi[i]]])
                        if(is.null(timeI[[idxi[i]]]))
                            timeIi <- 0 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                        timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + surveyTimes[idxi[i]])
                        ## Index by age (sam, sms)
                        obsIA[[idxi[i]]][y,] <- q[idxi[i]] * as.numeric(by(NAAsurv *
                                                                           as.numeric(t(dat$selI[[idxi[i]]])),
                                                                           as2a, sum)) * errs$time$eImvA[[idxi[[i]]]][y,] * as.numeric(errs$rep$eImvA[[idxi[[i]]]])
                        ## Index by length
                        obsIL[[idxi[i]]][y,] <- q[idxi[i]] * (NAAsurv * as.numeric(t(dat$selI[[idxi[i]]])) * errs$time$eImvL[[idxi[[i]]]][y,] * as.numeric(errs$rep$eImvL[[idxi[[i]]]])) %*% dat$plba
                    }
                }
            }
            ## observing annual M
            obsMAA[y,] <- obsMAA[y,] + MAA[seq(s,asmax,ns)]

            ## recruitment
            if(spawning[s] > 0){
                ## Survivors from previous season/year
                SSB[y,s] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA)) ## pre-recruitment mort
##                print(SSB[y,s])
                SSB0 <- get.ssb0(MAA, maty, weighty, fecun=1, asmax,
                                 ns, spawning, indage0 = indage0,
                                 R0 = R0y * errs$time$eR[y] * errs$rep$eR,
                                 season = s, FM = 0)

                SSB0vec[y] <- SSB0
                SPRvec[y] <- SSB0/(R0y * errs$time$eR[y] * errs$rep$eR)

                rec[y,s] <-  spawning[s] * recfunc(h = hy,
                                                   SSBPR0 = SSB0/(R0y * errs$time$eR[y] * errs$rep$eR),
                                                   SSB = SSB[y,s],
                                                   R0 = R0y, method = dat$SR,
                                                   bp = dat$bp,
                                                   beta = betay,
                                                   gamma = dat$recGamma,
                                                   alpha = alphay) * errs$time$eR[y] * errs$rep$eR
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

            ##
            if(is.numeric(nyI) && !set$surveyBeforeRec){
                if(s %in% idxS){
                    idxi <- which(idxS == s)
                    for(i in 1:length(idxi)){
                        surveyTime <- surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                        NAAsurv <- exp(log(NAAS) - ZAA * surveyTime)
                        ESBsurv <- NAAsurv * weightFy * as.numeric(t(dat$selI[[idxi[i]]]))
                        ## Total index (for spict)
                        obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * errs$time$eI[[idxi[i]]][y] * errs$rep$eI[[idxi[i]]])
                        if(is.null(timeI[[idxi[i]]]))
                            timeIi <- 0 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                        timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + surveyTimes[idxi[i]])
                        ## Index by age (sam, sms)
                        obsIA[[idxi[i]]][y,] <- q[idxi[i]] * as.numeric(by(NAAsurv *
                                                                           as.numeric(t(dat$selI[[idxi[i]]])),
                                                                           as2a, sum)) * errs$time$eImvA[[idxi[[i]]]][y,] * as.numeric(errs$rep$eImvA[[idxi[[i]]]])
                        ## Index by length
                        obsIL[[idxi[i]]][y,] <- q[idxi[i]] * (NAAsurv * as.numeric(t(dat$selI[[idxi[i]]])) * errs$time$eImvL[[idxi[[i]]]][y,] * as.numeric(errs$rep$eImvL[[idxi[[i]]]])) %*% dat$plba
                    }
                }
            }

            ## Exponential decay
            NAAS <- NAAS * exp(-ZAA)
            if(s == ns && set$recordLast == 1){
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
            obsIL[[i]] <- obsIL[[i]][(ny-nyI[i]+1):ny,]
            rownames(obsIA[[i]]) <- as.character((ny-nyI[i]+1):ny)
            colnames(obsIA[[i]]) <- as.character(0:dat$amax)
            attributes(obsIA[[i]]) <- c(attributes(obsIA[[i]]), list(time = dat$surveyTimes))
            rownames(obsIL[[i]]) <- as.character((ny-nyI[i]+1):ny)
            colnames(obsIL[[i]]) <- as.character(dat$mids)
            attributes(obsIL[[i]]) <- c(attributes(obsIL[[i]]), list(time = dat$surveyTimes))
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
            obsC <- apply(CW[idxC,], 1, sum) *
                errs$time$eC[idxC] * errs$rep$eC
            timeC <- idxC
            ## catch observations at age (SAM, SMS)
            for(y in 1:length(idxC))
                obsCA[y,] <- as.numeric(by(apply(CAA[[idxC[y]]],
                                                 1, sum), as2a, sum)) *
                    errs$time$eCmv[y,] * as.numeric(errs$rep$eCmv)
        }else if(nsC == ns){
            ## seasonal catches (for each season in OM)
            ## ----------------
            obsC <- as.numeric(t(CW[idxC,] *
                                 errs$time$eC[idxC] * errs$rep$eC))
            timeC <- rep(idxC, each = ns) + rep(seasonStart, nyC)
            ## catch observations at age (seasonal CAA observations not yet implemented, what is target format?)
            for(y in 1:length(idxC))
                obsCA[y,] <- as.numeric(by(apply(CAA[[idxC[y]]],
                                                 1, sum),as2a,sum)) *
                    errs$time$eCmv[y,] * as.numeric(errs$rep$eCmv)
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
                obsE <- apply(FM[idxE,], 1, sum) * errs$time$eF[idxE] * errs$rep$eF * qE * errs$time$eE[idxE] * errs$rep$eE
                timeE <- idxE
            }else if(nsE == ns){
                ## seasonal effort (for each season in OM)
                ## ----------------
                obsE <- as.numeric(t(FM[idxE,] * errs$time$eF[idxE] * errs$rep$eF * qE * errs$time$eE[idxE] * errs$rep$eE))
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
        obsC <- CW[idxC,] * errs$time$eC[idxC] * errs$rep$eC
        timeC <- idxC
        ## catch observations at age
        for(y in 1:length(idxC)) obsCA[y,] <- as.numeric(by(CAA[[idxC[y]]], as2a, sum)) * errs$time$eCmv[y,] * as.numeric(errs$rep$eCmv)

        ## Effort
        if(is.numeric(nyE)){
            if(nsE > 1) writeLines("Set dat$ns to > 1 for seasonal effort. Generating annual catches!")
            obsE <- FM[idxE,] * errs$time$eF[idxE] * errs$rep$eF * qE * errs$time$eE[idxE] * errs$rep$eE
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
    propMature <- matrix(as.numeric(by(dat$mat,dat$as2a,mean)),
                         length(idxC), amax, byrow = TRUE) ## matrix(rowSums(dat$mat), length(idxC), amax, byrow = TRUE)
    rownames(propMature) = as.character((ny-nyC+1):ny)
    colnames(propMature) = as.character(0:dat$amax)

    ## mean stock weight
    tmp <- sapply(1:ny, function(x) as.numeric(by(weight[,x],dat$as2a,mean)))
    WAAs <- t(tmp)[idxC,]
        ## matrix(, length(idxC), amax, byrow = TRUE) ## matrix(rowSums(dat$weight), length(idxC), amax, byrow = TRUE)
    rownames(WAAs) = as.character((ny-nyC+1):ny)
    colnames(WAAs) = as.character(0:dat$amax)

    ## mean catch weight
    tmp <- sapply(1:ny, function(x) as.numeric(by(weightF[,x],dat$as2a,mean)))
    WAAc <- t(tmp)[idxC,]
##    WAAc <- matrix(as.numeric(by(dat$weightF,dat$as2a,mean)), length(idxC), amax, byrow = TRUE) ## matrix(rowSums(dat$weightF), length(idxC), amax, byrow = TRUE)
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
                obsIL = obsIL,
                obsMAA = obsMAA,
                propMature = propMature,
                WAAs = WAAs,
                WAAc = WAAc,
                propFemale = propFemale,
                landFrac = landFrac)

    atti <- attributes(obs$obsCA)
    for(j in 1:length(atti$dimnames)){
        maxi <- max(nchar(atti$dimnames[[j]]))
        attributes(obs$obsCA)$dimnames[[j]][ nchar(atti$dimnames[[j]]) < maxi] <- paste0(sapply(maxi - nchar(atti$dimnames[[j]])[nchar(atti$dimnames[[j]]) < maxi], function(x)
            rep("0",x)),
            atti$dimnames[[j]][ nchar(atti$dimnames[[j]]) < maxi])
    }
    for(i in 1:length(obs$obsIA)){
        atti <- attributes(obs$obsIA[[i]])
        for(j in 1:length(atti$dimnames)){
            maxi <- max(nchar(atti$dimnames[[j]]))
            attributes(obs$obsIA[[i]])$dimnames[[j]][ nchar(atti$dimnames[[j]]) < maxi] <- paste0(sapply(maxi - nchar(atti$dimnames[[j]])[nchar(atti$dimnames[[j]]) < maxi], function(x)
                rep("0",x)),
                atti$dimnames[[j]][ nchar(atti$dimnames[[j]]) < maxi])
        }
    }


    ## return
    out <- NULL
    if(out.opt == 1){
        out$lastNAA <- NAAS
        out$lastFAA <- FAA
        out$TSBfinal <- TSBfinal
        out$SSBfinal <- SSBfinal
        out$ESBfinal <- ESBfinal
        out$rec <- rec
        out$SPR <- SPRvec
        out$SSB0 <- SSB0vec
        out$TSB <- TSB
        out$TSB1plus <- TSB1plus
        out$ESB <- ESB
        out$SSB <- SSB
        out$FM <- FM * errs$time$eF * errs$rep$eF
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
    }else if(out.opt == 5){
        out <- CW[ny]
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
            tmp <- min(seasonStart[seasonStart <= dat$surveyTimes[i]])
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
            ## CHECK: not sure if that's correct! used as assessSeasons[[s]] later, but here always in position 1
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
    errs <- list()
    errs$time <- get.errs(dat, set, ysim, hist)
    errs$rep <- get.errs(dat, set, 1, hist, rep = TRUE)

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
    R0y <- dat$R0 * errs$time$eR0[y] * errs$rep$eR0
    mselInd <- ifelse(mselFlag, y, 1)
    ##    MAA <- M[s1vec[y]:(s1vec[y]+ns-1)] * as.numeric(t(Msel[[mselInd]])) * eM
    MAA <- M[y,] * as.numeric(t(Msel[[mselInd]])) * errs$time$eM[y] * errs$rep$eM
    hy <- dat$h * errs$time$eH[y]
    maty <- as.numeric(t(mat)) * errs$time$eMat[y] * errs$rep$eMat
    selInd <- ifelse(selFlag, y, 1)
    sely <- as.numeric(t(sel[[selInd]])) * errs$time$eSel[y] * errs$rep$eSel
    weighty <- as.numeric(t(weight)) * errs$time$eW[y] * errs$rep$eW
    weightFy <- as.numeric(t(weightF)) * errs$time$eW[y] * errs$rep$eW
    alphay <- dat$recAlpha * errs$time$eAlpha[y] * errs$rep$eAlpha
    betay <- dat$recBeta * errs$time$eBeta[y] * errs$rep$eBeta

    ## TAC from previous year (if e.g. interim advice)
    ntac <- ns
    TACprev <- hist$TACprev ## TAC for full interval
    if(is.null(TACprev)){
        TACprev <- tail(as.numeric(t(hist$CW)),ntac)
    }
    ## Use only part before assessment from previous TAC in assessment year
    if(year %in% assessYears){ ## TODO: here also account for intY?
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

        if(s == 1 && set$recordLast == 0){
            ## beginning of year biomasses
            TSBfinal[y] <- sum(NAAS * weighty)
            ESBfinal[y] <- sum(NAAS * weighty * sely)
            SSBfinal[y] <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))
        }

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
                             R0 = R0y * errs$time$eR[y] * errs$rep$eR, season = s)
            rec[y,s] <- spawning[s] * recfunc(h = hy, SSBPR0 = SSB0/(R0y * errs$time$eR[y] * errs$rep$eR),
                                              SSB = SSBtmp, R0 = R0y,
                                              method = dat$SR, bp = dat$bp,
                                              beta = betay,
                                              gamma = dat$recGamma,
                                              alpha = alphay) * errs$time$eR[y] * errs$rep$eR
            rec[rec<0] <- 1e-10
            NAAS[indage0] <- rec[y,s]
        }


        ## index observations before assessment (not implemented for obsIA yet)
        if(all(!is.na(dat$surveyTimes))){
            if(s %in% idxS){
                idxi <- which(idxS == s)
                for(i in 1:length(idxi)){
                    if(dat$surveyBeforeAssessment[idxi[i]]){
                        ## survey observation: total catch in weight (spict)
                        surveyTime <- dat$surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                        NAAsurv <- exp(log(NAAS) - ZAA * surveyTime)
                        ESBsurv <- NAAsurv * weightFy * as.numeric(t(dat$selI[[idxi[i]]]))
                        obs$obsI[[idxi[i]]] <-
                            c(obs$obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * errs$time$eI[[idxi[i]]][y] * errs$rep$eI[[idxi[i]]])
                        if(is.null(obs$timeI[[idxi[i]]])){
                            timeIi <- ny-nyI[idxi[i]]+1
                        }else timeIi <- floor(tail(obs$timeI[[idxi[i]]],1))
                        obs$timeI[[idxi[i]]] <-
                            c(obs$timeI[[idxi[i]]], timeIi + 1 + dat$surveyTimes[idxi[i]])
                    }
                }
            }
        }



        ## Assessment
        if(any(year > assessYears)){
            yearsAfterAssessment <- year - max(assessYears[assessYears < year])
        }else{
            yearsAfterAssessment <- 1
        }
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
                tacs <- est.tac(obs. = obs,
                                hcr. = hcr,
                                tacs. = tacs,
                                pars. = list("ffmsy" = ffmsy,
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

                if(set$assessmentIntYear > 0){
                    if(year > set$assessmentIntYear){
                        tmp <- tacs[nrow(tacs)-set$assessmentIntYear,
                                    paste0("TAC", yearsAfterAssessment)]
                        TACs[y] <- as.numeric(as.character(tmp)) * errs$time$eImp[y] * errs$rep$eImp
                        TACreal <- rep(TACs[y] / ntac, ntac)
                    }else{
                        TACs[y] <-  TACprev * errs$time$eImp[y] * errs$rep$eImp  ## * set$assessmentIntYear
                        TACreal <- rep(TACs[y] / ntac, ntac)
                    }

                }else{
                    tmp <- tacs[nrow(tacs), "TAC1"]
                    TACs[y] <- as.numeric(as.character(tmp)) * errs$time$eImp[y] * errs$rep$eImp
                    TACreal <- rep(TACs[y] / ntac, ntac)
                }

            }

            ## CHECK: do not equally split tac among time steps! (only used for TACprev)

        }else{

            if(s %in% assessmentTiming){

                tmp <- tacs[nrow(tacs), paste0("TAC", yearsAfterAssessment)]
                TACs[y] <- as.numeric(as.character(tmp)) * errs$time$eImp[y] * errs$rep$eImp
                TACreal <- rep(TACs[y] / ntac, ntac)

            }
        }
        ## TODO: not tested for ns > 1 & assessmentIntYear > 1

        if(year %in% assessYears && s %in% assessmentTiming){
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

                ## TODO: this does not consider assessmentInterval!! (e.g. if TAC corresponds to 2 years!)

                FMtac <- min(set$maxF/ns,
                             get.f(TACs[y],
                                   NAA = NAAS, MAA = MAA,
                                   sel = sely, weight = weighty,
                                   seasons = assessSeasons[[s]], ns = ns, y = y,
                                   h = hy, asmax = asmax, mat = maty,
                                   pzbm = pzbm, spawning = spawning,
                                   R0 = R0y, SR = dat$SR, bp = dat$pb,
                                   recBeta = dat$recBeta,
                                   recGamma = dat$recGamma, eR = errs$time$eR[y] * errs$rep$eR,
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
            TACreal[indtac] <- sum(CAA[,s] * weightF, na.rm = TRUE) ## weightF[,s]
            if(indtac == ntac) TACs[y] <- tacs$TAC1[nrow(tacs)] <- sum(TACreal)
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
                    if(!dat$surveyBeforeAssessment[idxi[i]]){
                        ## survey observation: total catch in weight (spict)
                        surveyTime <- dat$surveyTimes[idxi[i]] - seasonStart[idxS[idxi[i]]]
                        NAAsurv <- exp(log(NAAS) - ZAA * surveyTime)
                        ESBsurv <- NAAsurv * weightFy * as.numeric(t(dat$selI[[idxi[i]]]))
                        obs$obsI[[idxi[i]]] <-
                            c(obs$obsI[[idxi[i]]], q[idxi[i]] * sum(ESBsurv) * errs$time$eI[[idxi[i]]][y] * errs$rep$eI[[idxi[i]]])
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
                                                      errs$time$eImvA[[idxi[[i]]]][ysim,] * errs$rep$eImvA[[idxi[[i]]]])
                        rownames(obs$obsIA[[idxi[i]]])[nrow(obs$obsIA[[idxi[i]]])] <- as.character(y)  ## CHECK: instead of 1 (assuming one survey a year) this should allow for s surveys
                    }
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
        if(s == ns && set$recordLast == 1){
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
            obs$obsC <- c(obs$obsC, sum(CW[y,]) *
                                    errs$time$eC[y] * errs$rep$eC)
            if(!is.null(obs$timeC))
                timeCi <- tail(obs$timeC,1) else timeCi <- ny-nyC+1
            obs$timeC <- c(obs$timeC, timeCi + 1)
            ## catch observations at age (SAM, SMS)
            obs$obsCA <- rbind(obs$obsCA, as.numeric(by(apply(CAA, 1, sum), as2a, sum)) * errs$time$eCmv[ysim,] * as.numeric(errs$rep$eCmv))
        }else if(nsC == ns){
            obs$obsC <- c(obs$obsC, CW[y,] *
                                    errs$time$eC[y] * errs$rep$eC)
            if(!is.null(obs$timeC))
                timeCi <- floor(tail(obs$timeC,1)) else timeCi <- ny-nyC+1
            obs$timeC <- c(obs$timeC, timeCi + 1 + rep(seasonStart, ny))
            ## catch observations at age (SAM, SMS)
            obs$obsCA <- rbind(obs$obsCA, as.numeric(by(apply(CAA, 1, sum), as2a, sum)) * errs$time$eCmv[ysim,] * as.numeric(errs$rep$eCmv))
        }else{
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
        ## Effort
        if(is.numeric(nyE)){
            if(nsE == 1){
                obs$obsE <- c(obs$obsE, sum(FM[y,]) * qE * errs$time$eE[y] * errs$rep$eE)
                if(!is.null(obs$timeE)) timeEi <- tail(obs$timeE,1) else timeEi <- ny-nyE+1
                obs$timeE <- c(obs$timeE, timeEi + 1)
            }else if(nsE == ns){
                obs$obsE <- c(obs$obsE, FM[y,] * qE * errs$time$eE[y] * errs$rep$eE)
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
        obs$obsC <- c(obs$obsC, sum(CW[y,]) *
                                errs$time$eC[y] * errs$rep$eC)
        if(!is.null(obs$timeC)) timeCi <- tail(obs$timeC,1) else timeCi <- ny-nyC+1
        obs$timeC <- c(obs$timeC, timeCi + 1)
        ## catch observations at age (SAM, SMS)
        obs$obsCA <- rbind(obs$obsCA, as.numeric(by(CAA, as2a, sum)) * errs$time$eCmv[ysim,] * as.numeric(errs$rep$eCmv))
        ## Effort
        if(is.numeric(nyE)){
            if(nsE > 1) writeLines("Set dat$ns to > 1 for seasonal effort data. Generating annual effort observations!")
            obs$obsE <- c(obs$obsE, sum(FM[y,]) * qE * errs$time$eE[y] * errs$rep$eE)
            if(!is.null(obs$timeE)) timeEi <- tail(obs$timeE,1) else timeEi <- ny-nyE+1
            obs$timeE <- c(obs$timeE, timeEi + 1)
        }
    }
    rownames(obs$obsCA)[nrow(obs$obsCA)] <- as.character(y)

    ## natural mortality
    obs$obsMAA <- rbind(obs$obsMAA, y = obsMAAtmp)
    rownames(obs$obsMAA)[nrow(obs$obsMAA)] = as.character(y)

    ## proportion mature
    obs$propMature <- rbind(obs$propMature,
                            matrix(as.numeric(by(dat$mat,dat$as2a,mean)),
                                   1, amax, byrow = TRUE))
    ## obs$propMature <- rbind(obs$propMature, matrix(rowSums(dat$mat), 1, amax, byrow = TRUE))
    rownames(obs$propMature)[nrow(obs$propMature)] = as.character(y)

    ## mean stock weight
    tmp <- sapply(year,
                  function(x) as.numeric(by(dat$weight,dat$as2a,mean))) ## weight[,x]
    WAAs <- t(tmp)
    obs$WAAs <- rbind(obs$WAAs, matrix(WAAs, 1, amax, byrow = TRUE))
    ## obs$WAAs <- rbind(obs$WAAs, matrix(rowSums(dat$weight), 1, amax, byrow = TRUE))
    rownames(obs$WAAs)[nrow(obs$WAAs)] = as.character(y)

    ## mean catch weight
    tmp <- sapply(year,
                  function(x) as.numeric(by(dat$weightF,dat$as2a,mean))) ## weightF[,x]
    WAAc <- t(tmp)
    obs$WAAc <- rbind(obs$WAAc, matrix(WAAc, 1, amax, byrow = TRUE))
    ## obs$WAAc <- rbind(obs$WAAc, matrix(rowSums(dat$weightF), 1, amax, byrow = TRUE))
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
