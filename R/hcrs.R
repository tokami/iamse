## Harvest control rules (HCRs)
##-----------------------------


#' @name gettacs
gettacs <- function(tacs=NULL, id="", TAC=NA, inp=NULL){
    if(!is.null(inp) && is.list(inp$obsI))
        nis <- length(inp$obsI) else nis <- 1
    tactmp <- data.frame(TAC=TAC, id=id, hitSC=NA, red=NA,
                         barID=FALSE, sd=NA, conv = FALSE,
                         fmfmsy.est=NA,fmfmsy.sd=NA,
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
                           indBpBx = NA))
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
    function(inp, tacs=NULL){
        inp <- spict::check.inp(inp)
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
    function(inp, tacs=NULL){
        inp <- spict::check.inp(inp)
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
    function(inp, tacs=NULL){
        inp <- spict::check.inp(inp)
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
    function(inp, tacs=NULL){
        inp <- spict::check.inp(inp)
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
                           env = globalenv()
                          ){

    if(is.null(constantC)) constantC = NA

    template  <- expression(paste0(
        '
structure(
    function(inp, tacs = NULL){
        inp <- spict::check.inp(inp)
        if(is.null(tacs)){
            indBpBx <- inp$indBpBx
        }else{
            indBpBx <- tacs$indBpBx[nrow(tacs)]
        }
        inp$indBpBx <- indBpBx

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

        tacs <- gettacs(tacs, id = "',id,'", TAC = TAC, inp=inp)
        tacs$hitSC <- NA
        tacs$barID <- barID
        tacs$red <- red
        tacs$indBpBx <- indBpBx
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
defHCRindex <- function(id = "r2_3",
                        x = 2,
                        y = 3,
                        stab = FALSE,
                        lower = 0.8,
                        upper = 1.2,
                        clyears = 1,
                        red = NA,
                        redyears = 2,
                        env = globalenv()
                        ){

    template  <- expression(paste0(
        'structure(
    function(inp, tacs = NULL){
        x <- ',x,'
        y <- ',y,'
        stab <- ',stab,'
        lower <- ',lower,'
        upper <- ',upper,'
        clyears <- ',clyears,'
        red <- ',red,'
        redyears <- ',redyears,'

        inp <- spict::check.inp(inp)
        if(is.null(tacs)){
            indBpBx <- inp$indBpBx
        }else{
            indBpBx <- tacs$indBpBx[nrow(tacs)]
        }
        inp$indBpBx <- indBpBx

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
        r <- mean(inum, na.rm = TRUE)/mean(iden, na.rm = TRUE)
        ## uncertainty cap
        if(stab){
            r[r<lower] <- lower
            r[r>upper] <- upper
            if(any(r < lower) || any(r > upper)) hitSC <- TRUE else hitSC <- FALSE
        }else hitSC <- FALSE
        ## account for seasonal and annual catches
        ## Cl <- sum(tail(inpin$obsC, tail(1/inpin$dtc,1))) ## CHECK: dtc required?
        Cl <- mean(tail(inp$obsC, clyears))
        TAC <- Cl * r * 1 * 1 ## Clast * r * f * b
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

        tacs <- gettacs(tacs, id = "',id,'", TAC = TAC, inp = inp)
        tacs$hitSC <- hitSC
        tacs$barID <- barID
        tacs$red <- red
        tacs$indBpBx <- indBpBx
        return(tacs)
    },
class="hcr"
)'))

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
                        schaefer = 0,
                        bfac = NA,
                        manstartdY = 1,
                        intC = NA,
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
    function(inp, tacs = NULL){
        inp$reportmode <- ',reportmode,'
        inp$dteuler <- ',dteuler,'
        inp$stabilise <- ',stabilise,'
        ## Intermediate year
        manstart <- inp$timeC[length(inp$timeC)] + ',manstartdY,'
        inp$maninterval <- c(manstart, manstart + 1)
        ## Check inp
        inp <- spict::check.inp(inp)
        ## priors
        inp$priors$logn <- c(',priorlogn[1],',',priorlogn[2],',',priorlogn[3],')
        inp$priors$logalpha <- c(',priorlogalpha[1],',',priorlogalpha[2],',',priorlogalpha[3],')
        inp$priors$logbeta <- c(',priorlogbeta[1],',',priorlogbeta[2],',',priorlogbeta[3],')
        ## Bfac rule
        if(is.null(tacs)){
            indBpBx <- inp$indBpBx
        }else{
            indBpBx <- tacs$indBpBx[nrow(tacs)]
        }
        inp$indBpBx <- indBpBx
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
        if(',schaefer,'){
            inp$phases$logn <- -1
            inp$ini$logn <- log(2)
        }
        if(is.list(inp$obsI)) nis <- length(inp$obsI)
        ## inp$optimiser.control <- list(iter.max = 1e3, eval.max = 1e3)
        rep <- try(fit.spict(inp), silent=TRUE)
        if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
            tacs <- conscat(inp, tacs=tacs)
            tacs$conv[nrow(tacs)] <- FALSE
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
            if(any(is.na(quantstmp))){
                tacs <- conscat(inp, tacs=tacs)
                tacs$conv[nrow(tacs)] <- FALSE
            }else{
                tac <- try(spict:::get.TAC(rep = rep,
                                           bfac = ',bfac,',
                                           fractiles = list(catch = ',frc,',
                                                            ffmsy = ',frff,',
                                                            bbmsy = ',frbb,',
                                                            bmsy  = ',frb,',
                                                            fmsy  = ',frf,'),
                                           breakpointB = ',breakpointB,',
                                           safeguardB = list(limitB = ',limitB,',prob = ',prob,'),
                                           intermediatePeriodCatch = intC2,
                                           verbose = FALSE),
                           silent = TRUE)
                if(inherits(tac, "try-error") || !is.numeric(tac)){
                    tacs <- conscat(inp, tacs=tacs)
                    tacs$conv[nrow(tacs)] <- FALSE
                }else{
                    tactmp <- data.frame(TAC=tac, id="',id,'", hitSC=NA,
                                         red=NA, barID=NA, sd=NA, conv = TRUE)
                    tactmp <- data.frame(c(tactmp, quantstmp,
                                           indBpBx = indBpBx))
                    if(is.null(tacs)){
                        tacs <- tactmp
                    }else{
                        tacs <- rbind(tacs, tactmp)
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
