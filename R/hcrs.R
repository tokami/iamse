## management procedures
##-----------------------------


## constant catch
##-----------------------------
#' @name conscat
#' @export
conscat <- function(inpin, tacs=NULL){
    inpin <- check.inp(inpin)
    if(is.null(tacs)){
        indBpBx <- inpin$indBpBx
    }else{
        indBpBx <- tacs$indBpBx[nrow(tacs)]
    }
    TAC <- sum(tail(inpin$obsC, 1))
    tactmp <- data.frame(TAC=TAC, id="cc", hitSC=FALSE, red=FALSE,
                         barID=FALSE, sd=NA, conv = FALSE,
                         fmfmsy.est=NA,fmfmsy.sd=NA,
                         bpbmsy.est=NA,bpbmsy.sd=NA,
                         cp.est=NA,cp.sd=NA,
                         fmsy.est=NA,fmsy.sd=NA,
                         bmsy.est=NA,bmsy.sd=NA,
                         indBpBx=indBpBx)
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}
##-----------------------------


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
