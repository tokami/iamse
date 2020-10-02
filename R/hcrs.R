## Harvest control rules (HCRs)
##-----------------------------


#' @name gettacs
gettacs <- function(tacs=NULL, id="", TAC=NA, inp=NULL){
    if(!is.null(inp) && is.list(inp$obsI))
        nis <- length(inp$obsI) else nis <- 1
    tactmp <- data.frame(TAC=TAC, id=id,
                         hitSC=NA,
                         red=NA,
                         barID=FALSE,
                         sd=NA,
                         conv = FALSE,               ## model convergence
                         fmfmsy.est=NA,fmfmsy.sd=NA, ## bunch of spict estimates (point estimate + sd)
                         bpbmsy.est=NA,bpbmsy.sd=NA,
                         cp.est=NA,cp.sd=NA,
                         fmsy.est=NA,fmsy.sd=NA,
                         bmsy.est=NA,bmsy.sd=NA,
                         sdb.est=NA,sdb.sd=NA)
    sdi <- rep(c(NA,NA),nis)
    names(sdi) <- paste0(rep(c("sdi.est","sdi.sd"),nis),rep(1:nis,each=nis))
    tactmp <- data.frame(c(tactmp,
                           sdi,
                           sdf.est=NA,sdf.sd=NA,
                           sdc.est=NA,sdc.sd=NA,
                           bmbmsy.est=NA,bmbmsy.sd=NA,
                           n.est=NA,n.sd=NA,
                           K.est=NA,K.sd=NA,
                           m.est=NA,m.sd=NA,
                           indBref = NA,
                           bmID=NA, assessInt=NA,
                           medbpbref=NA, bpbref=NA))
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}




#' @name defHCRref
#'
#' @param consF either numeric indicating constant F level or "fmsy" for fishing at fmsy
#'
#' @export
#'
defHCRref <- function(consF = 0,
                      fracFmsy = NULL,
                      env = globalenv()
                      ){

    id <- NULL

    if(!is.null(fracFmsy)){
        id <- paste0("refFmsy-", as.character(fracFmsy))
        template  <- expression(paste0(
                '
structure(
    function(inp, tacs=NULL, pars=NULL){
        inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- gettacs(tacs, id="',id,'", TAC=NA, inp=inp)
        return(tacs)
    },
    class="hcr"
)
'))
    }else{
        if(!is.numeric(consF) && !(consF %in% c("fmsy","Fmsy","FMSY","refFmsy","refFMSY","reffmsy")))
            stop("'consF' has to be either the absolute F (numeric) or 'fmsy'.")

        if(is.numeric(consF)){
            ## TODO: don't think this is implemented, check mse.R
            id <- paste0("consF-",consF)
            template  <- expression(paste0(
                '
structure(
    function(inp, tacs=NULL, pars=NULL){
        inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- gettacs(tacs, id="',id,'", TAC=NA, inp=inp)
        return(tacs)
    },
    class="hcr"
)
'))
        }

        if(consF == "fmsy" || consF == "Fmsy" || consF == "FMSY"){
            id <- "refFmsy"
            template  <- expression(paste0(
                '
structure(
    function(inp, tacs=NULL, pars=NULL){
        inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- gettacs(tacs, id="',id,'", TAC=NA, inp=inp)
        return(tacs)
    },
    class="hcr"
)
'))
        }

        if(consF == 0){
            id <- "noF"
            template  <- expression(paste0(
                '
structure(
    function(inp, tacs=NULL, pars=NULL){
        inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- gettacs(tacs, id="',id,'", TAC=0, inp=inp)
        return(tacs)
    },
    class="hcr"
)
'))
        }
    }


    ## create HCR as functions
    templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}



#' @name defHCRconscat
#' @title Define harvest control rule
#'
#' @export
#'
defHCRconscat <- function(id = "conscat",
                          constantC = NULL,
                           clyears = 1,
                           red = NA,
                          redyears = 2,
                          assessmentInterval = 1,
                           env = globalenv()
                          ){

    if(is.null(constantC)) constantC = NA

    template  <- expression(paste0(
        '
structure(
    function(inp, tacs = NULL, pars=NULL){
        inp <- spict::check.inp(inp, verbose = FALSE)
        TAC <- ',constantC,'
        if(!is.numeric(TAC)){
            annualcatch <- spict:::annual(inp$timeC, inp$obsC/inp$dtc, type = "mean") ## CHECK: why not sum?
            TAC <- mean(tail(annualcatch$annvec, ',clyears,'))
        }

        ## bianunal reduction (usually 0.2)
        if(is.numeric(red)){
            if(is.null(tacs)){
                TAC <- TAC * (1-red)
                barID <- TRUE
            }else{
                idx1 <- ifelse(nrow(tacs) > (redyears-1), (nrow(tacs)-(redyears-2)), 1)
                idx <- idx1:nrow(tacs)
                if(all(as.logical(tacs$barID[idx]) == FALSE)){
                    TAC <- TAC * (1-red)
                    barID <- TRUE
                }else{
                    barID <- FALSE
                }
            }
        }else barID <- FALSE

        ## Account for non-annual assessments
        TAC <- TAC * ',assessmentInterval,'

        tacs <- gettacs(tacs, id = "',id,'", TAC = TAC, inp=inp)
        tacs$hitSC[nrow(tacs)] <- NA
        tacs$barID[nrow(tacs)] <- barID
        tacs$red[nrow(tacs)] <- red
        tacs$assessInt[nrow(tacs)] <- assessmentInterval
        return(tacs)
    },
class="hcr"
)
'))

    ## create HCR as functions
    templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}




#' @name defHCRindex
#' @title Define harvest control rule
#' @details This function allows to define harvest control rules (HCRs) which can be incorporated into a
#' management strategy evaluation framework (DLMtool package). HCRs are saved with a
#' generic name to the global environment and the names of the HCRs are returned if results of the
#' function are assigned to an object. HCR runs a SPiCT assessment using catch and
#' relative biomass index observations. Stock status estimates are used to set the TAC
#' for the next year. TAC can be based on the distribution of predicted catches (percentileC)
#' and/or the distribution of the Fmsy reference level (percentileFmsy).
#' Additionally, a cap can be applied to account for low biomass levels (below Bmsy).
#' Arguments of returned function are 'x' - the position in a data-limited mehods data object,
#' 'Data' - the data-limited methods data object (see DLMtool), and 'reps' - the number of
#' stochastic samples of the TAC recommendation (not used for this HCR).
#' One or several arguments of the function can be provided as vectors to generate several
#' HCRs at once (several vectors have to have same length).
#'
#' @export
#'
defHCRindex <- function(id = "r23",
                        x = 2,
                        y = 3,
                        stab = TRUE,
                        lower = 0.8,
                        upper = 1.2,
                        clType = "TAC",
                        clyears = 1,
                        red = NA,
                        redyears = 3,
                        redAlways = FALSE,
                        ffmsySD = 0,
                        bbtriggerSD = 0,
                        rightRef=1,
                        assessmentInterval = 1,
                        env = globalenv()
                        ){

    template  <- expression(paste0(
        '
structure(
    function(inp, tacs = NULL, pars = NULL){
        x <- ',x,'
        y <- ',y,'
        stab <- ',stab,'
        lower <- ',lower,'
        upper <- ',upper,'
        clyears <- ',clyears,'
        clType <- "',clType,'"
        red <- ',red,'
        redyears <- ',redyears,'

        ffmsy <- rnorm(1, pars$ffmsy, ',ffmsySD,')
        ## ffmsy <- runif(1, pars$ffmsy * ',ffmsySD,', pars$ffmsy)
        ffmsy[ffmsy < 0] <- 0
        bbtrigger <- rnorm(1, pars$bbmsy*0.5, ',bbtriggerSD,')
        ## bbtrigger <- runif(1, pars$bbmsy*0.5, pars$bbmsy*0.5 * ',bbtriggerSD,')
        bbtrigger[bbtrigger < 0] <- 0

        inp <- spict::check.inp(inp, verbose = FALSE)
        indBref <- inp$indBref[1]
        inds <- inp$obsI
        if(length(inds) > 1){
            ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
            indtab <- do.call(rbind, inds)
            ind <- apply(indtab, 2, mean)
        }else{
            ind <- unlist(inds)
        }
        ninds <- length(ind)
        inum <- ind[(ninds-(x-1)):ninds]
        iden <- ind[(ninds-(x+y-1)):(ninds-x)]
        r <- r0 <- mean(inum, na.rm = TRUE)/mean(iden, na.rm = TRUE)
        ## account for seasonal and annual catches
        ## cl <- sum(tail(inp$obsC, tail(1/inp$dtc,1))) ## CHECK: dtc required?
        if(clType == "observed"){
            cl <- mean(tail(inp$obsC, clyears))
        }else if(clType == "TAC"){
            if(is.null(tacs)){
                cl <- mean(tail(inp$obsC, 3))
            }else{
                cl <- tacs$TAC[nrow(tacs)]
            }
        }
        tac <- cl * r * 1 * 1 ## Clast * r * f * b
        ## uncertainty cap
        if(stab){
            cllo <- cl * lower
            clup <- cl * upper
            if(any(tac < cllo)) hitSC <- 1 else hitSC <- 0
            if(any(tac > clup)) hitSC <- 2 else hitSC <- 0
            tac[tac < cllo] <- cllo
            tac[tac > clup] <- clup
        }else hitSC <- 0
        ## PA buffer (e.g. 0.2 reduction of TAC) if B < Btrigger proxy or F > Fmsy
        if(is.numeric(red)){
            if(is.null(tacs)){
                ## apply in first year
                barID <- TRUE
            }else if(any(as.logical(tail(tacs$barID,(redyears-1))),na.rm=TRUE)){
                ## do not apply if applied during last x years (redyears)
                barID <- FALSE
            }else{
                if(',redAlways,'){
                    barID <- TRUE
                }else{
                    right <- ifelse(runif(1) <= ',rightRef,', TRUE, FALSE)
                    if((ffmsy > 1 || bbtrigger < 1) && right){
                        ## apply if any ref indicates overexploitation
                        barID <- TRUE
                    }else barID <- FALSE
                }
            }
        }else barID <- FALSE
        ## apply reduction
        if(barID){
            tac <- tac * (1-red)
        }
        ## Account for non-annual assessments
        tac <- tac * ',assessmentInterval,'

        tacs <- gettacs(tacs, id = "',id,'", TAC = tac, inp = inp)
        tacs$hitSC[nrow(tacs)] <- hitSC
        tacs$barID[nrow(tacs)] <- barID
        tacs$red[nrow(tacs)] <- red
        tacs$fmfmsy.est[nrow(tacs)] <- ffmsy
        tacs$bpbmsy.est[nrow(tacs)] <- bbtrigger
        tacs$fmfmsy.sd[nrow(tacs)] <- ffmsySD
        tacs$bpbmsy.sd[nrow(tacs)] <- bbtriggerSD
        tacs$n.est[nrow(tacs)] <- r0
        tacs$indBref[nrow(tacs)] <- indBref
        tacs$assessInt[nrow(tacs)] <- assessmentInterval
        return(tacs)
    },
    class="hcr"
)
'))

    ## create HCR as functions
    templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}





#' @name defHCRspict
#' @title Define harvest control rule
#' @details This function allows to define harvest control rules (HCRs) which can be incorporated into a
#' management strategy evaluation framework (DLMtool package). HCRs are saved with a
#' generic name to the global environment and the names of the HCRs are returned if results of the
#' function are assigned to an object. HCR runs a SPiCT assessment using catch and
#' relative biomass index observations. Stock status estimates are used to set the TAC
#' for the next year. TAC can be based on the distribution of predicted catches (percentileC)
#' and/or the distribution of the Fmsy reference level (percentileFmsy).
#' Additionally, a cap can be applied to account for low biomass levels (below Bmsy).
#' Arguments of returned function are 'x' - the position in a data-limited mehods data object,
#' 'Data' - the data-limited methods data object (see DLMtool), and 'reps' - the number of
#' stochastic samples of the TAC recommendation (not used for this HCR).
#' One or several arguments of the function can be provided as vectors to generate several
#' HCRs at once (several vectors have to have same length).
#'
#' @importFrom doBy which.minn
#'
#' @export
#'
defHCRspict <- function(id = "spict-msy",
                        fractiles = list(catch=0.5,
                                         ffmsy=0.5,
                                         bbmsy=0.5,
                                         bmsy = 0.5,
                                         fmsy = 0.5),
                        breakpointB = 0.0,
                        safeguard = list(limitB = 0, prob = 0.95),
                        dteuler = 1/4,
                        reportmode = 1,
                        stabilise = 0,
                        priorlogn = c(log(2),2,1),
                        priorlogalpha = c(log(1),2,1),
                        priorlogbeta = c(log(1),2,1),
                        fixn = FALSE,
                        bfac = NA,
                        bref = "current", ## lowest or "lowest5" or "average"
                        brefType = "target",
                        manstartdY = 0,
                        assessmentInterval = 1,
                        intC = NA,
                        nonconvHCR = "conscat",
                        clType = "TAC",
                        clyears = 1,
                        stab = FALSE,
                        lower = 0.8,
                        upper = 1.2,
                        bm = FALSE,  ## e.g. 5 => re-defining Bref at every benchmark
                        env = globalenv()
                        ){

    frc <- fractiles$catch
    if(is.null(frc)) frc <- 0.5
    frff <- fractiles$ffmsy
    if(is.null(frff)) frff <- 0.5
    frbb <- fractiles$bbmsy
    if(is.null(frbb)) frbb <- 0.5
    frb <- fractiles$bmsy
    if(is.null(frb)) frb <- 0.5
    frf <- fractiles$fmsy
    if(is.null(frf)) frf <- 0.5
    limitB <- safeguard$limitB
    if(is.null(limitB)) limitB <- 0
    prob <- safeguard$prob
    if(is.null(prob)) prob <- 0.95

    template  <- expression(paste0(
        '
structure(
    function(inp, tacs = NULL, pars=NULL){
        func <- get(nonconvHCR)
        inp$reportmode <- ',reportmode,'
        inp$dteuler <- ',dteuler,'
        inp$stabilise <- ',stabilise,'
        bfac <- ',bfac,'
        bm <- ',bm,'
        fixn <- "',fixn,'"
        prob <- ',prob,'
        ## Intermediate year
        manstart <- inp$timeC[length(inp$timeC)] + 1 + ',manstartdY,' ## assumes annual catches
        inp$maninterval <- c(manstart, manstart + ',assessmentInterval,')
        inp$maneval <- max(inp$maninterval)
        ## Check inp
        inp <- spict::check.inp(inp, verbose = FALSE)
        ## priors
        inp$priors$logn <- c(',priorlogn[1],',',priorlogn[2],',',priorlogn[3],')
        inp$priors$logalpha <- c(',priorlogalpha[1],',',priorlogalpha[2],',',priorlogalpha[3],')
        inp$priors$logbeta <- c(',priorlogbeta[1],',',priorlogbeta[2],',',priorlogbeta[3],')
        ## Catch for intermediate year
        if(is.na(',intC,')){
            intC2 <- NULL
        }else{
            if(is.null(tacs)){
                intC2 <- tail(inp$obsC,1) ## does not account for seasonal catches!
            }else{
                intC2 <- tacs$TAC[nrow(tacs)]
            }
        }
        if(fixn == "schaefer"){
            inp$phases$logn <- -1
            inp$ini$logn <- log(2)
        }else if(fixn == "thorsonMean"){
            inp$phases$logn <- -1
            inp$ini$logn <- log(1.478)
        }else if(fixn == "thorsonClupeids"){
            inp$phases$logn <- -1
            inp$ini$logn <- log(0.599)
        }
        if(is.list(inp$obsI)) nis <- length(inp$obsI)
        if(is.null(tacs)){
            indBref2 <- inp$indBref[1]
        }else{
            indBref2 <- tacs$indBref[nrow(tacs)]
        }
        medbpbref <- NA
        bpbref <- NA
        ## benchmark (assuming bm always in first year)
        if(is.null(tacs)){
            bmID <- TRUE
        }else{
            if(!is.numeric(bm) ||
               any(as.logical(tail(tacs$bmID,(bm-1))),na.rm=TRUE)){
                bmID <- FALSE
            }else bmID <- TRUE
        }
        rep <- try(fit.spict(inp), silent=TRUE)
        if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
            tacs <- func(inp, tacs=tacs, pars=pars)
            tacs$conv[nrow(tacs)] <- FALSE
            tacs$indBref[nrow(tacs)] <- indBref2
            tacs$bmID[nrow(tacs)] <- bmID
            tacs$assessInt[nrow(tacs)] <- assessmentInterval
        }else{
            fmfmsy <- try(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(!is.numeric(fmfmsy)) print(paste0("fmfmsy not numeric. fmfmsy: ", fmfmsy))
            if(all(is.numeric(fmfmsy))){
                fmfmsy <- round(fmfmsy,2)
            }else fmfmsy <- c(NA, NA)
            names(fmfmsy) <- c("fmfmsy.est","fmfmsy.sd")
            bpbmsy <- try(get.par("logBpBmsy",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(bpbmsy))){
                bpbmsy <- round(bpbmsy,2)
            }else bpbmsy <- c(NA,NA)
            names(bpbmsy) <- c("bpbmsy.est","bpbmsy.sd")
            cp <- try(get.par("logCp",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(cp))){
                cp <- round(cp,2)
            }else cp <- c(NA,NA)
            names(cp) <- c("cp.est","cp.sd")
            ##
            fmsy <- try(get.par("logFmsy",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(fmsy))){
                fmsy <- round(fmsy,2)
            }else fmsy <- c(NA,NA)
            names(fmsy) <- c("fmsy.est","fmsy.sd")
            bmsy <- try(get.par("logBmsy",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(bmsy))){
                bmsy <- round(bmsy,2)
            }else bmsy <- c(NA,NA)
            names(bmsy) <- c("bmsy.est","bmsy.sd")
            sdb <- try(get.par("logsdb",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(sdb))){
                sdb <- round(sdb,2)
            }else sdb <- c(NA,NA)
            names(sdb) <- c("sdb.est","sdb.sd")
            sdi <- try(get.par("logsdi",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(sdi))){
                sdi <- round(sdi,2)
                sdi <- as.numeric(t(sdi))
            }else{
                sdi <- rep(c(NA,NA),nis)
            }
            names(sdi) <- paste0(rep(c("sdi.est","sdi.sd"),nis),rep(1:nis,each=nis))
            sdf <- try(get.par("logsdf",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(sdf))){
                sdf <- round(sdf,2)
            }else sdf <- c(NA,NA)
            names(sdf) <- c("sdf.est","sdf.sd")
            sdc <- try(get.par("logsdc",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(sdc))){
                sdc <- round(sdc,2)
            }else sdc <- c(NA,NA)
            names(sdc) <- c("sdc.est","sdc.sd")
            bmbmsy <- try(get.par("logBmBmsy",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(bmbmsy))){
                bmbmsy <- round(bmbmsy,2)
            }else bmbmsy <- c(NA,NA)
            names(bmbmsy) <- c("bmbmsy.est","bmbmsy.sd")
            nest <- try(get.par("logn",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(nest))){
                nest <- round(nest,2)
            }else nest <- c(NA,NA)
            names(nest) <- c("n.est","n.sd")
            Kest <- try(get.par("logK",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(Kest))){
                Kest <- round(Kest,2)
            }else Kest <- c(NA,NA)
            names(Kest) <- c("K.est","K.sd")
            mest <- try(get.par("logm",rep, exp=TRUE)[,c(2,4)],silent=TRUE)
            if(all(is.numeric(mest))){
                mest <- round(mest,2)
            }else mest <- c(NA,NA)
            names(mest) <- c("m.est","m.sd")
            ##
            quantstmp <- c(fmfmsy, bpbmsy, cp, fmsy, bmsy, sdb, sdi, sdf, sdc, bmbmsy, nest, Kest, mest)
            if(!(fixn %in% c("schaefer","thorsonMean","thorsonClupeids")) &&
               reportmode %in% c(0,1) && any(is.na(quantstmp))){
               ## n.sd is NA if schaefer (n fixed), reportmode 3 does not report most quantities
                tacs <- func(inp, tacs=tacs, pars=pars)
                tacs$conv[nrow(tacs)] <- FALSE
                tacs$indBref[nrow(tacs)] <- indBref2
                tacs$bmID[nrow(tacs)] <- bmID
                tacs$assessInt[nrow(tacs)] <- assessmentInterval
                tacs$medbpbref[nrow(tacs)] <- medbpbref
                tacs$bpbref[nrow(tacs)] <- bpbref
            }else{
                ## resetting brefs at benchmark
                if(bmID){
                    logB <- rep$obj$report(rep$obj$env$last.par.best)$logB[inp$indest]
                    if(bref == "current"){
                        indBref <- inp$indlastobs
                    }else if(bref == "lowest"){
                        indBref <- which.min(logB)
                    }else if(bref == "lowest5"){
                        indBref <- doBy::which.minn(logB, 5)
                    }else if(bref == "average"){
                        indBref <- 1:length(logB)
                    }else stop(paste0("bref = ",bref, " not known! Either current, lowest, or lowest5."))
                }else{
                    indBref <- tail(tacs$indBref,1)
                }
                rep <- try(set.bref(rep, indBref = indBref),silent=TRUE)
                ## get TAC
                if(inherits(rep, "try-error")){
                    tacs <- func(inp, tacs=tacs, pars=pars)
                    tacs$conv[nrow(tacs)] <- FALSE
                    tacs$indBref[nrow(tacs)] <- indBref2
                    tacs$bmID[nrow(tacs)] <- bmID
                    tacs$assessInt[nrow(tacs)] <- assessmentInterval
                    tacs$medbpbref[nrow(tacs)] <- medbpbref
                    tacs$bpbref[nrow(tacs)] <- bpbref
                }else{
                    indBref2 <- rep$inp$indBref[1]
                    logBpBref <- get.par("logBpBref", rep, exp = FALSE)
                    medbpbref <- exp(logBpBref[,2])
                    bpbref <- exp(qnorm(1-prob, logBpBref[2], logBpBref[4]))
                    tac <- try(spict:::get.TAC(rep = rep,
                                               bfac = bfac,
                                               bref.type = "',brefType,'",
                                               fractiles = list(catch = ',frc,',
                                                                ffmsy = ',frff,',
                                                                bbmsy = ',frbb,',
                                                                bmsy  = ',frb,',
                                                                fmsy  = ',frf,'),
                                               breakpointB = ',breakpointB,',
                                               safeguardB = list(limitB = ',limitB,',prob = prob),
                                               intermediatePeriodCatch = intC2,
                                               verbose = FALSE),
                               silent = TRUE)
                    if(inherits(tac, "try-error") || !is.numeric(tac) || is.na(tac)){
                        tacs <- func(inp, tacs=tacs, pars=pars)
                        tacs$conv[nrow(tacs)] <- FALSE
                        tacs$indBref[nrow(tacs)] <- indBref2
                        tacs$bmID[nrow(tacs)] <- bmID
                        tacs$assessInt[nrow(tacs)] <- assessmentInterval
                        tacs$medbpbref[nrow(tacs)] <- medbpbref
                        tacs$bpbref[nrow(tacs)] <- bpbref
                    }else{
                        if(',stab,'){
                            if("',clType,'" == "observed"){
                                cl <- mean(tail(inp$obsC, ',clyears,'))
                            }else if("',clType,'" == "estimated"){
                                cl <- mean(tail(get.par("logCpred",rep, exp=TRUE)[,2], ',clyears,'))
                            }else if("',clType,'" == "TAC"){
                                if(is.null(tacs)){
                                    cl <- mean(tail(inp$obsC, 3))
                                }else{
                                    cl <- tacs$TAC[nrow(tacs)]
                                }
                            }
                            cllo <- cl * ',lower,'
                            clup <- cl * ',upper,'
                            if(any(tac < cllo)) hitSC <- 1 else hitSC <- 0
                            if(any(tac > clup)) hitSC <- 2 else hitSC <- 0
                            tac[tac < cllo] <- cllo
                            tac[tac > clup] <- clup
                        }else hitSC <- 0

                        tactmp <- data.frame(TAC=tac, id="',id,'", hitSC=hitSC,
                                             red=NA, barID=NA, sd=NA, conv = TRUE)
                        tactmp <- data.frame(c(tactmp, quantstmp,
                                               indBref = indBref2, bmID=bmID,
                                               assessInt = assessmentInterval,
                                               medbpbref = medbpbref, bpref = bpref))
                        if(is.null(tacs)){
                            tacs <- tactmp
                        }else{
                            tacs <- rbind(tacs, tactmp)
                        }
                    }
                }
            }
        }
        return(tacs)
    },
    class="hcr")
'))

    ## create HCR as functions
    templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}
