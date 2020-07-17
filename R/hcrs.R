## Harvest control rules (HCRs)
##-----------------------------


#' @name gettacs
gettacs <- function(tacs=NULL, id="", TAC=NA){
    tactmp <- data.frame(TAC=TAC, id=id, hitSC=NA, red=NA,
                         barID=FALSE, sd=NA, conv = FALSE,
                         fmfmsy.est=NA,fmfmsy.sd=NA,
                         bpbmsy.est=NA,bpbmsy.sd=NA,
                         cp.est=NA,cp.sd=NA,
                         fmsy.est=NA,fmsy.sd=NA,
                         bmsy.est=NA,bmsy.sd=NA,
                         indBpBx = NA)
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}



#' @name conscat
#' @export
conscat <- structure(
    function(inp, tacs=NULL){
        inp <- check.inp(inp)
        if(is.null(tacs)){
            indBpBx <- inp$indBpBx
        }else{
            indBpBx <- tacs$indBpBx[nrow(tacs)]
        }
        tacs <- gettacs(tacs, id="cc",
                        TAC=sum(tail(inp$obsC, 1)))
        tacs$indBpBx[nrow(tacs)] <- indBpBx
        return(tacs)
    },
    class="hcr"
)




#' @name defHCRref
#'
#' @param consF either numeric indicating constant F level or "fmsy" for fishing at fmsy
#'
#' @export
#'
defHCRref <- function(consF = 0,
                      env = globalenv()
                      ){

    id <- NULL

    if(!is.numeric(consF) && !(consF %in% c("fmsy","Fmsy","FMSY","refFmsy","refFMSY","reffmsy")))
        stop("'consF' has to be either the absolute F (numeric) or 'fmsy'.")

    if(is.numeric(consF)){
        id <- paste0("consF",consF)
        template  <- expression(paste0(
            '
structure(
    function(inp, tacs=NULL){
        tacs <- gettacs(tacs, id="',id,'", TAC=NA)
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
        tacs <- gettacs(tacs, id="',id,'", TAC=NA)
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
        tacs <- gettacs(tacs, id="',id,'", TAC=0)
        return(tacs)
    },
    class="hcr"
)
'))
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
                           stab = FALSE,
                           lower = 0.8,
                           upper = 1.2,
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

        inp <- check.inp(inp)
        if(is.null(tacs)){
            indBpBx <- inp$indBpBx
        }else{
            indBpBx <- tacs$indBpBx[nrow(tacs)]
        }
        inp$indBpBx <- indBpBx

        TAC <- ',constantC,'
        if(is.na(TAC)) TAC <- sum(tail(inp$obsC, 1))
## TODO: account for seasonal catches sum(tail(inp$obsC, seasons)) ## where get seasons from?

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

        tacs <- gettacs(tacs, id = "',id,'", TAC = TAC)
        tacs$hitSC <- hitSC
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

        inp <- check.inp(inp)
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

        tacs <- gettacs(tacs, id = "',id,'", TAC = TAC)
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
        'structure(
        function(inp, tacs = NULL){
            inp$reportmode <- ',reportmode,'
            inp$dteuler <- ',dteuler,'
            inp$stabilise <- ',stabilise,'
            inp <- check.inp(inp)
            inp$priors$logn <- c(',priorlogn[1],',',priorlogn[2],',',priorlogn[3],')
            inp$priors$logalpha <- c(',priorlogalpha[1],',',priorlogalpha[2],',',priorlogalpha[3],')
            inp$priors$logbeta <- c(',priorlogbeta[1],',',priorlogbeta[2],',',priorlogbeta[3],')
            if(is.null(tacs)){
                indBpBx <- inp$indBpBx
            }else{
                indBpBx <- tacs$indBpBx[nrow(tacs)]
            }
            inp$indBpBx <- indBpBx
            if(',schaefer,'){
               inp$phases$logn <- -1
               inp$ini$logn <- log(2)
            }
            rep <- try(fit.spict(inp), silent=TRUE)
            if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                tacs <- conscat(inp, tacs=tacs)
            }else{
                fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4)],2)
                names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                bpbmsy <- round(get.par("logBpBmsy",rep, exp=TRUE)[,c(2,4)],2)
                names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4)],2)
                names(cp) <- paste0("cp-",names(cp))
                ##
                fmsy <- round(get.par("logFmsy",rep, exp=TRUE)[,c(2,4)],2)
                names(fmsy) <- paste0("fmsy-",names(fmsy))
                bmsy <- round(get.par("logBmsy",rep, exp=TRUE)[,c(2,4)],2)
                names(bmsy) <- paste0("bmsy-",names(bmsy))
                ##
                tac <- try(spict:::get.TAC(rep = rep,
                                           bfac = ',bfac,',
                                           fractiles = list(catch = ',frc,',
                                                            ffmsy = ',frff,',
                                                            bbmsy = ',frbb,',
                                                            bmsy  = ',frb,',
                                                            fmsy  = ',frf,'),
                                           breakpointB = ',breakpointB,',
                                           safeguardB = list(limitB = ',limitB,',prob = ',prob,'),
                                           verbose = FALSE),
                           silent = TRUE)
                if(inherits(tac, "try-error")){
                    tacs <- conscat(inp, tacs=tacs)
                }else{
                    tactmp <- data.frame(TAC=tac, id="',id,'", hitSC=NA,
                                         red=NA, barID=NA, sd=NA, conv = NA)
                    tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                           fmsy, bmsy, indBpBx = indBpBx))
                    if(is.null(tacs)){
                        tacs <- tactmp
                    }else{
                        tacs <- rbind(tacs, tactmp)
                    }
                }
            }
            return(tacs)
        },
        class="hcr")'))

    ## create HCR as functions
    templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}


















## 2/3
##-----------------------------
#' @name r23
#' @export
r23 <- function(inpin, stab=FALSE, lower=0.8, upper=1.2,
                red=NA, y_red = 2, tacs=NULL,
                targetC = FALSE){
    reps <- 1
    inds <- inpin$obsI
    if(length(inds) > 1){
        ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
        indtab <- do.call(rbind, inds)
        ind <- apply(indtab, 2, mean)
    }else{
        ind <- unlist(inds)
    }
    ninds <- length(ind)
    inum <- ind[(ninds-1):ninds]
    iden <- ind[(ninds-4):(ninds-2)]
    r23 <- mean(inum, na.rm = TRUE)/mean(iden, na.rm = TRUE)
    ## uncertainty cap
    if(stab){
        r23 <- spict:::stabilityClause(r23, lower, upper)
        if(any(r23 < lower) || any(r23 > upper)) hitSC <- TRUE else hitSC <- FALSE
    }else hitSC <- FALSE
    ## account for seasonal and annual catches
    ## Cl <- sum(tail(inpin$obsC, tail(1/inpin$dtc,1))) ## CHECK: dtc required?
    Cl <- sum(tail(inpin$obsC, 1))
    if(targetC) Cl <- mean(tail(inpin$obsC, 5))
    TACi <- Cl * r23 * 1 * 1 ## Clast * r * f * b
    ## bienniel reduction (usually 0.2)
    if(is.null(red)) red <- NA
    if(is.numeric(red)){
        if(is.null(tacs)){
            TACi <- TACi * (1-red)
            barID <- TRUE
        }else{
            idx1 <- ifelse(nrow(tacs) > (y_red-1), (nrow(tacs)-(y_red-2)), 1)
            idx <- idx1:nrow(tacs)
            if(all(as.logical(tacs$barID[idx]) == FALSE)){
                TACi <- TACi * (1-red)
                barID <- TRUE
            }else{
                barID <- FALSE
            }
        }
    }else barID <- FALSE
    TAC <- rep(TACi, reps)
    tactmp <- data.frame(TAC=TAC, id="23", hitSC=hitSC, red=red, barID=barID, sd=NA, conv = FALSE)
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}
##-----------------------------

## 12
##-----------------------------
#' @name r12
#' @export
r12 <- function(inpin, stab=FALSE, lower=0.8, upper=1.2,
                red=NA, y_red = 2, tacs=NULL,
                targetC = FALSE){
    reps <- 1
    inds <- inpin$obsI
    if(length(inds) > 1){
        ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
        indtab <- do.call(rbind, inds)
        ind <- apply(indtab, 2, mean)
    }else{
        ind <- unlist(inds)
    }
    ninds <- length(ind)
    inum <- ind[ninds]
    iden <- ind[(ninds-2):(ninds-1)]
    r23 <- mean(inum, na.rm = TRUE)/mean(iden, na.rm = TRUE)
    ## uncertainty cap
    if(stab){
        r23 <- spict:::stabilityClause(r23, lower, upper)
        if(any(r23 < lower) || any(r23 > upper)) hitSC <- TRUE else hitSC <- FALSE
    }else hitSC <- FALSE
    ## account for seasonal and annual catches
##    Cl <- sum(tail(inpin$obsC, tail(1/inpin$dtc,1)))
    Cl <- sum(tail(inpin$obsC, 1))
    if(targetC) Cl <- mean(tail(inpin$obsC, 5))
    TACi <- Cl * r23 * 1 * 1 ## Clast * r * f * b
    ## bienniel reduction (usually 0.2)
    if(is.null(red)) red <- NA
    if(is.numeric(red)){
        if(is.null(tacs)){
            TACi <- TACi * (1-red)
            barID <- TRUE
        }else{
            idx1 <- ifelse(nrow(tacs) > (y_red-1), (nrow(tacs)-(y_red-2)), 1)
            idx <- idx1:nrow(tacs)
            if(all(as.logical(tacs$barID[idx]) == FALSE)){
                TACi <- TACi * (1-red)
                barID <- TRUE
            }else{
                barID <- FALSE
            }
        }
    }else barID <- FALSE
    TAC <- rep(TACi, reps)
    tactmp <- data.frame(TAC=TAC, id="12", hitSC=hitSC, red=red, barID=barID, sd=NA, conv = FALSE)
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}
##-----------------------------

## 35
##-----------------------------
#' @name r35
#' @export
r35 <- function(inpin, stab=FALSE, lower=0.8, upper=1.2,
                red=NA, y_red = 2, tacs=NULL,
                targetC = FALSE){
    reps <- 1
    inds <- inpin$obsI
    if(length(inds) > 1){
        ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
        indtab <- do.call(rbind, inds)
        ind <- apply(indtab, 2, mean)
    }else{
        ind <- unlist(inds)
    }
    ninds <- length(ind)
    inum <- ind[(ninds-2):ninds]
    iden <- ind[(ninds-7):(ninds-3)]
    r23 <- mean(inum, na.rm = TRUE)/mean(iden, na.rm = TRUE)
    ## uncertainty cap
    if(stab){
        r23 <- spict:::stabilityClause(r23, lower, upper)
        if(any(r23 < lower) || any(r23 > upper)) hitSC <- TRUE else hitSC <- FALSE
    }else hitSC <- FALSE
    ## account for seasonal and annual catches
##    Cl <- sum(tail(inpin$obsC, tail(1/inpin$dtc,1)))
    Cl <- sum(tail(inpin$obsC, 1))
    if(targetC) Cl <- mean(tail(inpin$obsC, 5))
    TACi <- Cl * r23 * 1 * 1 ## Clast * r * f * b
    ## bienniel reduction (usually 0.2)
    if(is.null(red)) red <- NA
    if(is.numeric(red)){
        if(is.null(tacs)){
            TACi <- TACi * (1-red)
            barID <- TRUE
        }else{
            idx1 <- ifelse(nrow(tacs) > (y_red-1), (nrow(tacs)-(y_red-2)), 1)
            idx <- idx1:nrow(tacs)
            if(all(as.logical(tacs$barID[idx]) == FALSE)){
                TACi <- TACi * (1-red)
                barID <- TRUE
            }else{
                barID <- FALSE
            }
        }
    }else barID <- FALSE
    TAC <- rep(TACi, reps)
    tactmp <- data.frame(TAC=TAC, id="12", hitSC=hitSC, red=red, barID=barID, sd=NA, conv = FALSE)
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}
##-----------------------------


## derivative method
##-----------------------------
#' @name derimeth
#' @export
derimeth <- function(inpin, stab=FALSE, lower=0.8, upper=1.2,
                red=NA, y_red = 2, tacs=NULL, lambda = 0.4, cfac = 0.8){
    reps <- 1
    inds <- inpin$obsI
    if(length(inds) > 1){
        ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
        indtab <- do.call(rbind, inds)
        ind <- apply(indtab, 2, mean)
    }else{
        ind <- unlist(inds)
    }
    ninds <- length(ind)
    indi <- ind[(ninds-4):ninds]
    x <- 1:5
    mod <- summary(lm(indi ~ x))
    slope <- mod$coef[2]
    ## uncertainty cap
    if(stab){
        slope <- spict:::stabilityClause(slope, lower, upper)
        if(any(slope < lower) || any(slope > upper)) hitSC <- TRUE else hitSC <- FALSE
    }else hitSC <- FALSE
    ## account for seasonal and annual catches
    if(is.null(tacs)){
        taclast <- cfac * mean(tail(inpin$obsC, 5))
    }else{
        taclast <- tacs$TAC[nrow(tacs)]
    }
    TACi <- taclast * (1 + lambda*slope)
    ## bienniel reduction (usually 0.2)
    if(is.null(red)) red <- NA
    if(is.numeric(red)){
        if(is.null(tacs)){
            TACi <- TACi * (1-red)
            barID <- TRUE
        }else{
            idx1 <- ifelse(nrow(tacs) > (y_red-1), (nrow(tacs)-(y_red-2)), 1)
            idx <- idx1:nrow(tacs)
            if(all(as.logical(tacs$barID[idx]) == FALSE)){
                TACi <- TACi * (1-red)
                barID <- TRUE
            }else{
                barID <- FALSE
            }
        }
    }else barID <- FALSE
    TAC <- rep(TACi, reps)
    tactmp <- data.frame(TAC=TAC, id="dm", hitSC=hitSC, red=red, barID=barID, sd=NA, conv = FALSE)
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}
##-----------------------------


## DCV
##-----------------------------
## variance for 2/3 rule
#' @name var23
#' @export
var23 <- function(index, sd){  ## sdI list with associated sd for each index element (needed for var23)
    n <- length(index)
    idxnum <- c(n-1,n)
    idxden <- c(n-4,n-3,n-2)
    munum <- mean(index[idxnum])
    muden <- mean(index[idxden])
    sdnum <- sd[idxnum]
    sdden <- sd[idxden]
    (1/muden^2) * ((1/2)^2 *sum(sdnum^2)) + (munum^2)/(muden^4) * ((1/3)^2 *sum(sdden^2))
}
## r23 or spict-dl
#' @name spictDCV
#' @export
spictDCV <- function(repin, stab=FALSE, lower=0.8, upper=1.2, red=NA,
                     tacs=NULL,verbose=FALSE){
    inpin <- repin$inp
    reps <- 1
    inds <- inpin$obsI
    if(length(inds) > 1){
        ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
        indtab <- do.call(rbind, inds)
        ind <- apply(indtab, 2, mean)
        sdtab <- do.call(rbind, inpin$obsIsd)
        sdI <- apply(sdtab, 2, mean)
    }else{
        ind <- unlist(inds)
        sdI <- unlist(repin$inp$obsIsd)
    }
    ## compare sds
    sd23 <- sqrt(var23(ind, sd = sdI))
    sdspict <- get.par("logsdi",repin,exp=TRUE)[,2]
    if(verbose) print(ifelse(sdspict <= sd23, "spictsdl", "23"))
    if(sdspict <= sd23){
        tactmp <- get.TAC(repin, hcr="dl", bfrac=1, prob = 0.5,
                          stab=stab, lower=lower, upper=upper)
        tactmp$sd <- sdspict
    }else{
        ninds <- length(ind)
        inum <- ind[(ninds-1):ninds]
        iden <- ind[(ninds-4):(ninds-2)]
        r23 <- mean(inum, na.rm = TRUE)/mean(iden, na.rm = TRUE)
        ## uncertainty cap
        if(stab){
            r23 <- spict:::stabilityClause(r23, lower, upper)
            if(any(r23 < lower) || any(r23 > upper)) hitSC <- TRUE else hitSC <- FALSE
        }else hitSC <- FALSE
        ## account for seasonal and annual catches
        Cl <- sum(tail(inpin$obsC, tail(1/inpin$dtc,1)))
        TACi <- Cl * r23 * 1 * 1  ## Clast * r * f * b
        TAC <- rep(TACi, reps)
        tactmp <- data.frame(TAC=TAC, id="r23", hitSC=hitSC, red=red, barID=barID, sd=sd23, conv=FALSE)
    }
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}
##-----------------------------




## WCV
##-----------------------------
## r23 and spict-dl weighted
#' @name spictWCV
#' @export
spictWCV <- function(repin, stab=FALSE, lower=0.8, upper=1.2, red=NA,
                     tacs=NULL, verbose=FALSE){
    inpin <- repin$inp
    reps <- 1
    inds <- inpin$obsI
    if(length(inds) > 1){
        ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
        indtab <- do.call(rbind, inds)
        ind <- apply(indtab, 2, mean)
        sdtab <- do.call(rbind, inpin$obsIsd)
        sdI <- apply(sdtab, 2, mean)
    }else{
        ind <- unlist(inds)
        sdI <- unlist(inpin$obsIsd)
    }
    ## est spict tac
    tacspict <- get.TAC(repin, hcr="dl", bfrac=1, prob = 0.5,
                        stab=stab, lower=lower, upper=upper)$TAC
    ## est 23 tac
    ninds <- length(ind)
    inum <- ind[(ninds-1):ninds]
    iden <- ind[(ninds-4):(ninds-2)]
    r23 <- mean(inum, na.rm = TRUE)/mean(iden, na.rm = TRUE)
    ## uncertainty cap
    if(stab){
        r23 <- spict:::stabilityClause(r23, lower, upper)
        if(any(r23 < lower) || any(r23 > upper)) hitSC <- TRUE else hitSC <- FALSE
    }else hitSC <- FALSE
    ## account for seasonal and annual catches
    Cl <- sum(tail(inpin$obsC, tail(1/inpin$dtc,1)))
    tac23 <- Cl * r23 * 1 * 1  ## Clast * r * f * b
    ## est variances
    sd23 <- sqrt(var23(ind, sd = sdI))
    sdspict <- get.par("logsdi",repin,exp=TRUE)[,2]
    if(verbose) print(ifelse(sdpspict <= sd23, "spictsdl", "23"))
    ## weigh according to variances
    omega <- sdspict / (sd23 + sdspict)
    tac <- omega*tacspict + (1-omega)*tac23
    return(list(TAC=tac, id="", hitSC=hitSC, red = red, barID = NA, sd="", conv = FALSE))
    ## output
    tactmp <- data.frame(TAC=tac, id="wcv", hitSC=hitSC,
                         red=red, barID=NA, sd=NA, conv=TRUE)
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}
##-----------------------------
