
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
#' @importFrom abind abind
#'
#' @export
initPop <- function(specdat, set = NULL, refs = NULL, out.opt = 1, depl.quant = "B0"){

    ## indices
    ny <- specdat$ny
    ns <- specdat$nseason
    nt <- ny * ns
    nsC <- set$catchSeasons
    if(is.null(nsC)) nsC <- 1
    nyhist <- set$nyhist
    if(is.null(nyhist)) nyhist <- ny
    idx <- (ny - nyhist + 1):ny
    surveyTimes <- set$surveyTimes
    if(is.null(surveyTimes)) surveyTimes <- 0
    nsurv <- length(surveyTimes)
    ## closest season
    seasonStart <- seq(0,1-1/ns,1/ns)
    idxS <- rep(0, nsurv)
    for(i in 1:nsurv){
        tmp <- seasonStart[seasonStart <= surveyTimes[i]]
        idxS[i] <- which.min((tmp - surveyTimes[i])^2)
    }

    ## species data
    amax <- specdat$amax + 1  ## age 0
    M <- specdat$M
    Ms <- M / ns
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
        if(is.null(eF)) eF <- exp(rnorm(ny, 0, set$sigmaF) - set$sigmaF^2/2)
        if(is.null(eR)) {
            eR <- genDevs(ny, set$sigmaR, set$rhoR)
        }
        if(is.null(eM)) eM <- exp(rnorm(ny, 0, set$sigmaM) - set$sigmaM^2/2)
        if(is.null(eH)) eH <- exp(rnorm(ny, 0, set$sigmaH) - set$sigmaH^2/2)
        if(is.null(eR0)) eR0 <- exp(rnorm(ny, 0, set$sigmaR0) - set$sigmaR0^2/2)
        if(is.null(eMat)) eMat <- exp(rnorm(ny, 0, set$sigmaMat) - set$sigmaMat^2/2)
        if(is.null(eImp)) eImp <- exp(rnorm(ny, 0, set$sigmaImp) - set$sigmaImp^2/2)
        if(is.null(eC)) eC <- exp(rnorm(nyhist, 0, set$CVC) - set$CVC^2/2)
        if(is.null(eI)) eI <- exp(rnorm(nyhist, 0, set$CVI) - set$CVI^2/2)
    }else{
        eF <- exp(rnorm(ny, 0, 0))
        eR <- genDevs(ny, 0, 0)
        eM <- exp(rnorm(ny, 0, 0))
        eH <- exp(rnorm(ny, 0, 0))
        eR0 <- exp(rnorm(ny, 0, 0))
        eMat <- exp(rnorm(ny, 0, 0))
        eImp <- exp(rnorm(ny, 0, 0))
        eC <- exp(rnorm(nyhist, 0, 0))
        eI <- exp(rnorm(nyhist, 0, 0))
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
    TSB <- ESB <- SSB <- CW <- FM <- matrix(0, nrow=ny, ncol=ns)
    TACs <- rep(NA, ny)
    obsI <- vector("list", nsurv)
    timeI <- vector("list", nsurv)
    ## burnin period
    if(is.null(set)) burnin <- 20 else burnin <- set$burnin
    if(is.numeric(burnin) && burnin > 0){
        NAAbi <- rep(NA, amax)
        NAAbi[1] <- R0 * exp(initN[1])
        Fbi <- (sel * Fvals[1]) / ns
        Zbi <- Ms + Fbi
        for(a in 2:amax) NAAbi[a] <- NAAbi[a-1] * exp(-(Zbi[a-1])) * exp(initN[a])
        for(y in 1:burnin){
            ## recruitment
            SSBtemp <- sum(NAAbi * weight * mat * exp(-pzbm * Zbi)) ## pre-recruitment mort
            rec <- recfunc(h = h, SSBPR0 = SSBPR0, SSB = SSBtemp,
                           R0 = R0, method = SR)
            rec[rec<0] <- 1e-10
            NAAbi[1] <- rec
            for(s in 1:ns){
                ## can't take more than what's there
                Btemp <- sum(NAAbi * weight * sel * exp(-M/2))
                CWbi <- sum(baranov(Fbi, Ms, NAAbi) * weight)
                if(CWbi > 0.99 * Btemp){
                    Fbi <- sel * min(set$maxF/ns,
                                     getFM(0.75 * Btemp, NAA = NAAbi,
                                           M = Ms, weight = weightF, sel = sel))
                    Zbi <- Ms + Fbi
                }
                Ntemp <- NAAbi * exp(-Zbi)
                if(s < ns){
                    NAAbi <- Ntemp
                }else if(y < burnin){
                    NAAbi[amax] <- Ntemp[amax] + Ntemp[amax-1]
                    for(a in 2:(amax-1)) NAAbi[a] <- Ntemp[a-1]
                }
            }
        }
        NAA <- NAAbi
    }else{
        ## set up NAA
        NAA[1] <- exp(initN[1]) * R0  * eR0
        for(a in 2:amax)
            NAA[a] <- NAA[a-1] * exp(-Ms[a-1] + (initF*sel[a-1])/ns) * exp(initN[a])
    }

    ## main loop
    for(y in 1:ny){
        Msy <- Ms * eM[y]
        maty <- mat * eMat[y]
        hy <- h * eH[y]
        R0y <- R0 * eR0[y]
        ## recruitment
        SSB[y,1] <- sum(NAA * weight * maty * exp(-pzbm * (Msy + Fvals[y] * sel / ns))) ## pre-recruitment mort
        SSBPR0 <- getSSBPR0(Msy, maty, 1, amax)
        rec <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSB[y,1],
                       R0 = R0y, method = SR)
        rec[rec<0] <- 1e-10
        NAA[1] <- rec * eR[y]
        ## seasons
        for(s in 1:ns){
            ## FAA and Z
            FM[y,s] <- Fvals[y] * eF[y] / ns
            FAA <- sel * FM[y,s]
            Z <- Msy + FAA
            ## can't take more than what's there
            NAAmid <- NAA * exp(-Msy/2)
            Btemp <- sum(NAAmid * weight * sel * exp(-Msy/2))
            CAA <- baranov(FAA, Msy, NAA)
            CW[y,s] <- sum(CAA * weightF)
            if(CW[y,s] > 0.99 * Btemp){
                FM[y,s] <- min(set$maxF/ns,
                               getFM(0.75 * Btemp, NAA = NAA,
                                     M = Msy, weight = weightF, sel = sel))
                FAA <- sel * FM[y,s]
                Z <- Msy + FAA
                CAA <- baranov(FAA, Msy, NAA)
                CW[y,s] <- sum(weightF * CAA)
            }
            ## TSB
            TSB[y,s] <- sum(NAA * weight)
            ## ESB
            ESB[y,s] <- sum(NAA * weight * sel)
            ## SSB
            SSB[y,s] <- sum(NAA * weight * maty * exp(-pzbm * Z))

            ## index observations
            if(s %in% idxS){
                idxi <- which(idxS == s)
                for(i in 1:length(idxi)){
                    surveyTime <- set$surveyTimes[idxi[i]] - seasonStart[idxi[i]]
                    NAAsurv <- exp(log(NAA) - Z * surveyTime)
                    ESBsurv <- sum(NAAsurv * specdat$weightF * specdat$sel)
                    obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * ESBsurv * eI[y])
                    if(is.null(timeI[[idxi[i]]])) timeIi <- 0 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                    timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + set$surveyTimes[idxi[i]])
                }
            }

            ## ageing by season or year
            NAA <- Ntemp <- NAA * exp(-Z)
            if(s == ns){
                NAA[amax] <- Ntemp[amax] + Ntemp[amax-1]
                for(a in 2:(amax-1)) NAA[a] <- Ntemp[a-1]
            }
        }
    }


    ## catch observations
    if(ns > 1){
        if(nsC == 1){
            obsC <- apply(CW[idx,], 1, sum) * eC
            timeC <- 1:ny
        }else if(nsC == ns){
            obsC <- as.numeric(t(CW[idx,] * eC))
            timeC <- rep(1:ny, each = ns) + rep(seasonStart, ny)
        }else{
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
    }else{
        if(nsC > 1) writeLines("Set dat$nseasons to > 1 for seasonal catches. Generating annual catches!")
        obsC <- CW[idx,] * eC
        timeC <- 1:ny
    }

    ## return
    out <- NULL
    if(out.opt == 1){
        out$lastNAA <- NAA
        out$lastFAA <- FAA
        out$TSB <- TSB
        out$ESB <- ESB
        out$SSB <- SSB
        out$FM <- FM
        out$CW <- CW
        out$TACs <- TACs
        out$errs <- errs
        out$inp <- list()
        out$inp$obsC <- obsC
        out$inp$timeC <- timeC
        out$inp$obsI <- obsI
        out$inp$timeI <- timeI
    }else if(out.opt == 2){
        if(is.null(refs)){
            warning("The list with reference points is needed!")
        }else{
            out <- mean(TSB[ny,])/refs[[depl.quant]]
        }
    }
    return(out)
}



## advance population
#' @name advancePop
#' @export
advancePop <- function(specdat, hist, set, tacs){

    ## indices
    ny <- nrow(hist$TSB)
    y <- ny + 1
    ysim <- y - specdat$ny
    ns <- specdat$nseason
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
    amax <- specdat$amax + 1  ## age 0
    pzbm <- specdat$pzbm
    SR <- specdat$SR
    logR0 <- specdat$logR0
    h <- specdat$h
    SSBPR0 <- specdat$SSBPR0
    q <- dat$q
    if(length(q) < nsurv) q <- rep(q, nsurv)
    tacID <- tacs$id[nrow(tacs)]

    ## parameters per age
    Ms <- specdat$M / ns
    weight <- specdat$weight
    weightF <- specdat$weightF
    mat <- specdat$mat
    sel <- specdat$sel

    ## errors
    eF <- set$eF[ysim]
    eR <- set$eR[ysim]
    eM <- set$eM[ysim]
    eH <- set$eH[ysim]
    eR0 <- set$eR0[ysim]
    eMat <- set$eMat[ysim]
    eImp <- set$eImp[ysim]
    eC <- set$eC[ysim]
    eI <- set$eI[ysim]
    if(is.null(eR)) {
        eR <- genDevs(1, set$sigmaR, set$rhoR)
    }
    if(is.null(eF)) eF <- exp(rnorm(1, 0, set$sigmaF) - set$sigmaF^2/2)
    if(is.null(eM)) eM <- exp(rnorm(1, 0, set$sigmaM) - set$sigmaM^2/2)
    if(is.null(eR0)) eR0 <- exp(rnorm(1, 0, set$sigmaR0) - set$sigmaR0^2/2)
    if(is.null(eH)) eH <- exp(rnorm(1, 0, set$sigmaH) - set$sigmaH^2/2)
    if(is.null(eMat)) eMat <- exp(rnorm(1, 0, set$sigmaMat) - set$sigmaMat^2/2)
    if(is.null(eImp)) eImp <- exp(rnorm(1, 0, set$sigmaImp) - set$sigmaImp^2/2)
    if(is.null(eC)) eC <- exp(rnorm(1, 0, set$CVC) - set$CVC^2/2)
    if(is.null(eI)) eI <- exp(rnorm(1, 0, set$CVI) - set$CVI^2/2)
    if("errs" %in% names(hist)){
        errs <- list(eF = c(hist$errs$eF, eF),
                     eR = c(hist$errs$eR, eR),
                     eM = c(hist$errs$eM, eM),
                     eH = c(hist$errs$eH, eH),
                     eR0 = c(hist$errs$eR0,eR0),
                     eMat = c(hist$errs$eMat, eMat),
                     eImp = c(hist$errs$eImp, eImp),
                     eC = c(hist$errs$eC, eC),
                     eI = c(hist$errs$eI, eI))
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
    NAA <- rep(0, amax)
    TSB <- rbind(hist$TSB,tmp)
    SSB <- rbind(hist$SSB,tmp)
    ESB <- rbind(hist$ESB,tmp)
    CW  <- rbind(hist$CW,tmp)
    FM  <- rbind(hist$FM,tmp)
    TACs <- c(hist$TACs, NA)
    TACreal <- rep(NA, ns)
    obsI <- vector("list", nsurv)
    timeI <- vector("list", nsurv)

    ## project forward
    Msy <- Ms * eM
    R0y <- exp(logR0) * eR0
    hy <- h * eH
    maty <- mat * eMat

    ## Survivors from previous season/year
    Ztemp <- hist$lastFAA + specdat$M / ns
##    Ntemp <- hist$lastNAA ##* exp(-Ztemp)
    ## NAA[amax] <- Ntemp[amax] + Ntemp[amax-1]
    ## for(a in 2:(amax-1)) NAA[a] <- Ntemp[a-1]
    NAA <- hist$lastNAA ##* exp(-Ztemp)


    ## recruitment
    SSBtemp <- sum(NAA * weight * maty * exp(-pzbm * Ztemp)) ## pre-recruitment mort
    SSBPR0 <- getSSBPR0(Msy, maty, 1, amax)
    rec <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSBtemp, R0 = R0y, method = SR)
    rec[rec<0] <- 1e-10
    NAA[1] <- rec * eR

    ## Define F/TACs
    if(tacID == "refFmsy"){
        ## Fmsy
        FMtac <- hist$refs$Fmsy / ns
    }else if(tacID == "noF"){
        ## noF
        FMtac <- 0
    }else{
        ## any other HCR
        TAC <- as.numeric(as.character(tacs$TAC[nrow(tacs)]))
        TACs[y] <- TAC
        TACreal <- rep(TAC/ns, ns)
    }

    for(s in 1:ns){
        if(!(tacID %in% c("refFmsy","noF"))){
            FMtac <- min(set$maxF/ns,
                         getFM(TACreal[s], NAA = NAA, M = Msy,
                               weight = weightF, sel = sel))
        }
        FM[y,s] <- FMtac * eImp
        Fsy <- sel * FM[y,s]
        Z <- Msy + Fsy
        ## can't take more than what's there
        NAAmid <- NAA * exp(-Msy/2)
        Btemp <- sum(NAAmid * weight * sel * exp(-Msy/2))
        CAA <- baranov(Fsy, Msy, NAA)
        CW[y,s] <- sum(CAA * weightF)
        if(CW[y,s] > 0.99 * Btemp){
            FM[y,s] <- min(set$maxF/ns,
                             getFM(0.75 * Btemp, NAA = NAA,
                                   M = Msy, weight = weightF, sel = sel))
            Fsy <- FM[y,s] * sel
            Z <- Msy + Fsy
            CAA <- baranov(Fsy, Msy, NAA)
            CW[y,s] <- sum(CAA * weightF)
            if(s < ns){
                TACreal[s+1] <- TACreal[s+1] + TACreal[s] - CW[y,s]
            }else writeLines("Could not get full annual TAC.")
            TACreal[s] <- CW[y,s]
        }
        if(tacID == "refFmsy"){
            TACreal[s] <- sum(CAA * weightF, na.rm = TRUE)
            if(s == ns) TACs[y] <- tacs$TAC[nrow(tacs)] <- sum(TACreal)
        }
        ## TSB
        TSB[y,s] <- sum(NAA * weight)
        ## ESB
        ESB[y,s] <- sum(NAA * weight * sel)
        ## SSB
        SSB[y,s] <- sum(NAA * weight * maty * exp(-pzbm * Z))

        ## index observations
        if(s %in% idxS){
            idxi <- which(idxS == s)
            for(i in 1:length(idxi)){
                surveyTime <- set$surveyTimes[idxi[i]] - seasonStart[idxi[i]]
                NAAsurv <- exp(log(NAA) - Z * surveyTime)
                ESBsurv <- sum(NAAsurv * specdat$weightF * specdat$sel)
                obsI[[idxi[i]]] <- c(obsI[[idxi[i]]], q[idxi[i]] * ESBsurv * eI)
                if(is.null(timeI[[idxi[i]]])) timeIi <- 0 else timeIi <- floor(tail(timeI[[idxi[i]]],1))
                timeI[[idxi[i]]] <- c(timeI[[idxi[i]]], timeIi + 1 + set$surveyTimes[idxi[i]])
            }
        }

        ## ageing by season or year
        NAA <- Ntemp <- NAA * exp(-Z)
        if(s == ns){
            NAA[amax] <- Ntemp[amax] + Ntemp[amax-1]
            for(a in 2:(amax-1)) NAA[a] <- Ntemp[a-1]
        }
    }

    ## catch observations
    if(ns > 1){
        if(nsC == 1){
            obsC <- sum(CW[y,]) * eC
            timeC <- tail(hist$inp$timeC,1) + 1
        }else if(nsC == ns){
            obsC <- CW[y,] * eC
            timeC <- floor(tail(hist$inp$timeC,1)) + 1 + rep(seasonStart, ny)
        }else{
            stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
        }
    }else{
        if(nsC > 1) writeLines("Set dat$nseasons to > 1 for seasonal catches. Generating annual catches!")
        obsC <- sum(CW[y,]) * eC
        timeC <- tail(hist$inp$timeC,1) + 1
    }

    inp <- list(obsC = c(hist$obsC, obsC),
                timeC = c(hist$timeC, timeC),
                obsI = list(),
                timeI = list())
    for(i in 1:nsurv){
        inp$obsI[[i]] <- c(hist$obsI[[i]], obsI[[i]])
        inp$timeI[[i]] <- c(hist$timeI[[i]], timeI[[i]])
    }


    ## return
    out <- NULL
    out$lastNAA <- NAA
    out$lastFAA <- Fsy
    out$TSB <- TSB
    out$ESB <- ESB
    out$SSB <- SSB
    out$FM <- FM
    out$CW <- CW
    out$TACs <- TACs
    out$tacs <- tacs
    out$errs <- errs
    out$inp <- inp
    if("refs" %in% names(hist)) out$refs <- hist$refs
    return(out)
}
