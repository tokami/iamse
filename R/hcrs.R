## Harvest control rules (HCRs)
##-----------------------------

#' @name est.tac
#' @export
est.tac <- function(obs, hcr, tacs=NULL, pars=NULL){
    func <- get(hcr)
    tacs <- func(obs, tacs, pars)
    return(tacs)
}


#' @name gettacs
gettacs <- function(tacs=NULL, id="", TAC=NA, obs=NULL){
    if(!is.null(obs) && is.list(obs$obsI) && length(obs$obsI) != 0)
        nis <- length(obs$obsI) else nis <- 1
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
                           medbpbref=NA, bpbref=NA,
                           rwF=NA))
    if(is.null(tacs)){
        tacs <- tactmp
    }else{
        tacs <- rbind(tacs, tactmp)
    }
    return(tacs)
}




#' @name def.hcr.ref
#'
#' @param consF either numeric indicating constant F level or "fmsy" for fishing at fmsy
#'
#' @export
#'
def.hcr.ref <- function(consF = 0,
                        fracFmsy = NULL,
                        env = globalenv()
                        ){

    id <- NULL

    if(!is.null(fracFmsy)){
        id <- paste0("refFmsy_", as.character(fracFmsy))
        template  <- expression(paste0(
            '
structure(
    function(obs, tacs=NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }
        ## inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- gettacs(tacs, id="',id,'", TAC=NA, obs=inp)
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
            id <- paste0("consF_",consF)
            template  <- expression(paste0(
                '
structure(
    function(obs, tacs=NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }
        ## inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- gettacs(tacs, id="',id,'", TAC=NA, obs=inp)
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
    function(obs, tacs=NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }
        ## inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- gettacs(tacs, id="',id,'", TAC=NA, obs=inp)
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
    function(obs, tacs=NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }
        ## inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- gettacs(tacs, id="',id,'", TAC=0, obs=inp)
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



#' @name def.hcr.conscat
#' @title Define harvest control rule
#'
#' @export
#'
def.hcr.conscat <- function(id = "conscat",
                            constantC = NULL,
                            clyears = 1,
                            red = NA,
                            redyears = 2,
                            redAlways = FALSE,
                            assessmentInterval = 1,
                            ffmsySD = 0,
                            bbtriggerSD = 0,
                            rightRef=1,
                            env = globalenv()
                            ){

    if(is.null(constantC)) constantC = NA

    template  <- expression(paste0(
        '
structure(
    function(obs, tacs = NULL, pars=NULL){
        red <- ',red,'
        redyears <- ',redyears,'
        assessInt <- ',assessmentInterval,'

        ffmsy <- rnorm(1, pars$ffmsy, ',ffmsySD,')
        ## ffmsy <- runif(1, pars$ffmsy * ',ffmsySD,', pars$ffmsy)
        ffmsy[ffmsy < 0] <- 0
        bbtrigger <- rnorm(1, pars$bbmsy*2, ',bbtriggerSD,')
        ## bbtrigger <- runif(1, pars$bbmsy*2, pars$bbmsy*2 * ',bbtriggerSD,')
        bbtrigger[bbtrigger < 0] <- 0

        ## obs <- spict::check.inp(obs, verbose = FALSE)
        obs$obsC <- as.vector(na.omit(obs$obsC))
        obs$obsI <- lapply(obs$obsI,function(x) as.vector(na.omit(x)))
        obs$dtc <- min(diff(obs$timeC))
        ## if(is.null(tacs)){
        ##     indBref <- obs$indBref
        ## }else{
        ##     indBref <- as.numeric(as.character(unlist(strsplit(as.character(tacs$indBref[nrow(tacs)]), "-"))))
        ## }
        ## indBref2 <- paste(indBref, collapse="-")
        tac <- ',constantC,'
        if(!is.numeric(tac)){
            annualcatch <- spict:::annual(obs$timeC, obs$obsC/obs$dtc, type = "mean") ## CHECK: why not sum? TODO: remove spict dependency
            tac <- mean(tail(annualcatch$annvec, ',clyears,'))
            ## Account for non-annual assessments
            tac <- tac * assessInt
        }

        ## PA buffer (e.g. 0.2 reduction of TAC) if B < Btrigger proxy or F > Fmsy
        if(is.numeric(red)){
            if(is.null(tacs)){
                if(',redAlways,'){
                    barID <- TRUE
                }else{
                    if(ffmsy > 1 || bbtrigger < 1){
                        barID <- TRUE
                    }else barID <- FALSE
                    right <- ifelse(runif(1) <= ',rightRef,', barID, !barID)
                }
            }else if(any(as.logical(tail(tacs$barID,ceiling(redyears/assessInt-1))),na.rm=TRUE)){
                ## do not apply if applied during last x years (redyears)
                barID <- FALSE
            }else{
                if(',redAlways,'){
                    barID <- TRUE
                }else{
                    if(ffmsy > 1 || bbtrigger < 1){
                        barID <- TRUE
                    }else barID <- FALSE
                    right <- ifelse(runif(1) <= ',rightRef,', barID, !barID)
                }
            }
        }else barID <- FALSE
        ## apply reduction
        if(barID){
            tac <- tac * (1-red)
        }
        tacs <- gettacs(tacs, id = "',id,'", TAC = tac, obs=obs)
        tacs$hitSC[nrow(tacs)] <- NA
        tacs$barID[nrow(tacs)] <- barID
        tacs$red[nrow(tacs)] <- red
##        tacs$indBref[nrow(tacs)] <- indBref2
        tacs$assessInt[nrow(tacs)] <- assessInt
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




#' @name def.hcr.index
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
def.hcr.index <- function(id = "r23",
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
                          env = globalenv(),
                          dbg = FALSE
                          ){

    template  <- expression(paste0(
        '
structure(
    function(obs, tacs = NULL, pars = NULL){
        x <- ',x,'
        y <- ',y,'
        stab <- ',stab,'
        lower <- ',lower,'
        upper <- ',upper,'
        clyears <- ',clyears,'
        clType <- "',clType,'"
        red <- ',red,'
        redyears <- ',redyears,'
        assessInt <- ',assessmentInterval,'

        ffmsy <- rnorm(1, pars$ffmsy, ',ffmsySD,')
        ## ffmsy <- runif(1, pars$ffmsy * ',ffmsySD,', pars$ffmsy)
        ffmsy[ffmsy < 0] <- 0
        bbtrigger <- rnorm(1, pars$bbmsy*2, ',bbtriggerSD,')
        ## bbtrigger <- runif(1, pars$bbmsy*2, pars$bbmsy*2 * ',bbtriggerSD,')
        bbtrigger[bbtrigger < 0] <- 0

        ## obs <- spict::check.inp(obs, verbose = FALSE)
        obs$obsC <- as.vector(na.omit(obs$obsC))
        obs$obsI <- lapply(obs$obsI,function(x) as.vector(na.omit(x)))
        obs$dtc <- min(diff(na.omit(obs$timeC)))

        ## if(is.null(tacs)){
        ##     indBref <- obs$indBref
        ## }else{
        ##     indBref <- as.numeric(as.character(unlist(strsplit(as.character(tacs$indBref[nrow(tacs)]), "-"))))
        ## }
        ## indBref2 <- paste(indBref, collapse="-")
        ## benchmark (only benchmark in spict)
        bmID <- FALSE
        inds <- obs$obsI
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

if(',dbg,') print(paste0("r: ", round(r,2)))

        ## account for seasonal and annual catches
        ## cl <- sum(tail(obs$obsC, tail(1/obs$dtc,1))) ## CHECK: dtc required?
        if(clType == "observed"){
            cl <- mean(tail(obs$obsC, clyears))
            ## Account for non-annual assessments
            cl <- cl * assessInt
        }else if(clType == "TAC"){
            if(is.null(tacs)){
                cl <- mean(tail(obs$obsC, 3))
                ## Account for non-annual assessments
                cl <- cl * assessInt
            }else{
                cl <- tacs$TAC[nrow(tacs)]
            }
        }
        tac <- cl * r * 1 * 1 ## Clast * r * f * b
if(',dbg,') print(paste0("tac: ", round(tac,2)))
        ## uncertainty cap
        if(stab){
            cllo <- cl * lower
            clup <- cl * upper
if(',dbg,') print(paste0("cl: ", round(cl,2)))
            if(any(tac < cllo)) hitSC <- 1 else hitSC <- 0
            if(any(tac > clup)) hitSC <- 2 else hitSC <- 0
            tac[tac < cllo] <- cllo
            tac[tac > clup] <- clup
        }else hitSC <- 0
        ## PA buffer (e.g. 0.2 reduction of TAC) if B < Btrigger proxy or F > Fmsy
        if(is.numeric(red)){
            if(is.null(tacs)){
                if(',redAlways,'){
                    barID <- TRUE
                }else{
                    if(ffmsy > 1 || bbtrigger < 1){
                        barID <- TRUE
                    }else barID <- FALSE
                    barID <- ifelse(runif(1) <= ',rightRef,', barID, !barID)
                }
            }else if(any(as.logical(tail(tacs$barID,ceiling(redyears/assessInt-1))),na.rm=TRUE)){
                ## do not apply if applied during last x years (redyears)
                barID <- FALSE
            }else{
                if(',redAlways,'){
                    barID <- TRUE
                }else{
                    if(ffmsy > 1 || bbtrigger < 1){
                        barID <- TRUE
                    }else barID <- FALSE
                    barID <- ifelse(runif(1) <= ',rightRef,', barID, !barID)
                }
            }
        }else barID <- FALSE
        ## apply reduction
        if(barID){
            tac <- tac * (1-red)
        }
if(',dbg,') print(paste0("tac: ", round(tac,2)))
        tacs <- gettacs(tacs, id = "',id,'", TAC = tac, obs = obs)
        tacs$hitSC[nrow(tacs)] <- hitSC
        tacs$barID[nrow(tacs)] <- barID
        tacs$red[nrow(tacs)] <- red
        tacs$fmfmsy.est[nrow(tacs)] <- ffmsy
        tacs$bpbmsy.est[nrow(tacs)] <- bbtrigger
        tacs$fmfmsy.sd[nrow(tacs)] <- ffmsySD
        tacs$bpbmsy.sd[nrow(tacs)] <- bbtriggerSD
        tacs$n.est[nrow(tacs)] <- r0
        ## tacs$indBref[nrow(tacs)] <- indBref2
        tacs$bmID[nrow(tacs)] <- bmID
        tacs$assessInt[nrow(tacs)] <- assessInt
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




#' @name def.hcr.hr
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
def.hcr.hr <- function(id = "prop_8.75",
                       prop = 0.0875,
                       env = globalenv()
                       ){

    template  <- expression(paste0(
        '
structure(
    function(obs, tacs = NULL, pars = NULL){
        prop <- ',prop,'

        ## obs <- spict::check.inp(obs, verbose = FALSE)
        obs$obsC <- as.vector(na.omit(obs$obsC))
        obs$obsI <- lapply(obs$obsI,function(x) as.vector(na.omit(x)))
        obs$dtc <- min(diff(na.omit(obs$timeC)))

        ## if(is.null(tacs)){
        ##     indBref <- obs$indBref
        ## }else{
        ##     indBref <- as.numeric(as.character(unlist(strsplit(as.character(tacs$indBref[nrow(tacs)]), "-"))))
        ## }
        ## indBref2 <- paste(indBref, collapse="-")
        ## benchmark (only benchmark in spict)
        bmID <- FALSE
        inds <- obs$obsI
        if(length(inds) > 1){
            ## WHAT TO DO IF SEVERAL INDICES AVAILABLE? ## for now: mean
            indtab <- do.call(rbind, inds)
            ind <- apply(indtab, 2, mean)
        }else{
            ind <- unlist(inds)
        }
        ninds <- length(ind)

        tac <- ind[ninds] * prop

        tacs <- gettacs(tacs, id = "',id,'", TAC = tac, obs = obs)
        ## tacs$indBref[nrow(tacs)] <- indBref2
        tacs$bmID[nrow(tacs)] <- bmID
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





#' @name def.hcr.spict
#' @title Define harvest control rule
#'
#' @param id Name/ID of HCR. Default: "spict-msy"
#' @param fractiles Fractiles. List
#' @param breakpointB breakpointb
#' @param safeguard safeguard. List
#' @param dteuler Time step of Forward Euler scheme.
#' @param reportmode Reportmode
#' @param stabilise stabilise option of spict. Default: FALSE.
#' @param priorlogn prior
#' @param priorlogalpha prior
#' @param priorlogbeta prior
#' @param priorlogbkfrac prior
#' @param fixn fixn
#' @param bfac bfac
#' @param bref bref
#' @param brefType brefType
#' @param nyBref nyBref
#' @param btar btar
#' @param probtar probtar
#' @param brule Brule
#' @param red Reduction
#' @param redyears Reduction years
#' @param redAlways Reduction always? Default: FALSE.
#' @param rai Raising factor. Default: 0.2.
#' @param manstartdY Management start year. Default: 0.
#' @param assessmentInterval Assessment interval. Default: 1.
#' @param intC Interval C. Default: NA.
#' @param nonconvHCR HCR if SPiCT does not converge. Default: "conscat" (constant catch).
#' @param clType Catch type for uncertainty cap. Default: "TAC".
#' @param clyears Number of years for uncertainty cap. Default: 1.
#' @param stab Uncertainty cap. Default: FALSE.
#' @param lower Upper bound of uncertainty cap. Default: 0.8.
#' @param upper Upper bound of uncertainty cap. Default: 1.2.
#' @param bm Benchmark. Default: FALSE,  ## e.g. 5 => re-defining Bref at every benchmark
#' @param env Environment. Default: globalenv()
#'
#'
#' @importFrom doBy which.minn
#'
#' @export
#'
def.hcr.spict <- function(id = "spict-msy",
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
                          priorlogsdf = c(5,2,0),
                          priorlogsdc = c(log(0.2),2,0),
                          priorlogalpha = c(log(1),2,1),
                          priorlogbeta = c(log(1),2,1),
                          priorlogbkfrac = c(log(0.5),2,0),
                          fixn = NA,
                          bfac = NA,
                          bref = "current", ## lowest or "highest" or "average" or "last"
                          brefType = "target",
                          nyBref = 5,
                          btar = "bmsy",
                          probtar = 0.4,
                          brule = 0,
                          red = NA,
                          redyears = 3,
                          redAlways = FALSE,
                          rai = 0.2,
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
    breakpointB1 <- breakpointB[1]
    breakpointB2 <- breakpointB[2]

    template  <- expression(paste0(' structure(function(obs, tacs = NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }

        func <- get("',nonconvHCR,'")
        inp$reportmode <- ',reportmode,'
        inp$dteuler <- ',dteuler,'
        inp$stabilise <- ',stabilise,'
        bfac <- ',bfac,'
        bm <- ',bm,'
        fixn <- ',fixn,'
        prob <- ',prob,'
        lower <- ',lower,'
        upper <- ',upper,'
        probtar <- ',probtar,'
        red <- ',red,'
        redyears <- ',redyears,'
        rai <- ',rai,'
        brule <- ',brule,'
        assessInt <- ',assessmentInterval,'
        clyears <- ',clyears,'
        clType <- "',clType,'"
        nyBref <- ',nyBref,'
        ## Intermediate year
        manstart <- inp$timeC[length(inp$timeC)] + 1 + ',manstartdY,' ## assumes annual catches
        inp$maninterval <- c(manstart, manstart + assessInt)
        inp$maneval <- max(inp$maninterval)
        ## Check inp
        inp <- spict::check.inp(inp, verbose = FALSE)
        ## priors
        inp$priors$logn <- c(',priorlogn[1],',',priorlogn[2],',',priorlogn[3],')
        inp$priors$logsdf <- c(',priorlogsdf[1],',',priorlogsdf[2],',',
        priorlogsdf[3],')
        inp$priors$logsdc <- c(',priorlogsdc[1],',',priorlogsdc[2],',',
        priorlogsdc[3],')
        inp$priors$logalpha <- c(',priorlogalpha[1],',',priorlogalpha[2],',',
        priorlogalpha[3],')
        inp$priors$logbeta <- c(',priorlogbeta[1],',',priorlogbeta[2],',',
        priorlogbeta[3],')
        inp$priors$logbkfrac <- c(',priorlogbkfrac[1],',',priorlogbkfrac[2],',',
        priorlogbkfrac[3],')
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
        if(is.numeric(fixn) && !is.na(fixn)){
            inp$phases$logn <- -1
            inp$ini$logn <- log(fixn)
        }
        if(is.list(inp$obsI) && length(inp$obsI) != 0) nis <- length(inp$obsI) else nis <- 1
        if(is.null(tacs)){
            indBref <- inp$indBref
        }else{
            indBref <- as.numeric(as.character(unlist(strsplit(as.character(tacs$indBref[nrow(tacs)]), "-"))))
        }
        indBref2 <- paste(indBref, collapse="-")
        medbpbref <- NA
        bpbref <- NA
        ## benchmark (assuming bm always in first year)
        if(is.null(tacs)){
            bmID <- TRUE
        }else{
            if(!is.numeric(bm) ||
               any(as.logical(tail(tacs$bmID,ceiling(bm/assessInt-1))),na.rm=TRUE)){
                bmID <- FALSE
            }else bmID <- TRUE
        }

        ## ## if TAC/catch very small -> adjust random walk
        ## if(!is.null(tacs) && any(tacs$TAC  < 10)){
        ##     inp$priors$logbeta <- c(0,0,0)
        ##     inp$priors$logsdf <- c(1,0.5,1)
        ##     inp$priors$logsdc <- c(log(0.2),0.2,1)
        ##     rwF <- TRUE
        ## }else{
        ##     rwF <- FALSE
        ## }

        ## fit spict
        rwF  <- FALSE
##  inp$reportmode <- 0
        fit <- try(fit.spict(inp), silent=TRUE)
##  try(plot(fit),silent=TRUE)
        ## non-convergence
        if(class(fit) == "try-error" || fit$opt$convergence != 0 || any(is.infinite(fit$sd))){
            tacs <- func(inp, tacs=tacs, pars=pars)
            tacs$conv[nrow(tacs)] <- FALSE
            tacs$indBref[nrow(tacs)] <- indBref2
            tacs$bmID[nrow(tacs)] <- bmID
            tacs$assessInt[nrow(tacs)] <- assessInt
            tacs$rwF[nrow(tacs)] <- rwF
            return(tacs)
        }
        ## last years catch/tac
        if(clType == "observed"){
            cl <- mean(tail(inp$obsC, clyears))
        }else if(clType == "estimated"){
            cl <- mean(tail(get.par("logCpred",fit, exp=TRUE)[,2], clyears))
        }else if(clType == "TAC"){
            if(is.null(tacs)){
                cl <- mean(tail(inp$obsC, clyears))
                cl <- cl * assessInt
            }else{
                cl <- tacs$TAC[nrow(tacs)]
            }
        }
        fmfmsy <- try(get.par("logFmFmsynotS",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(!is.numeric(fmfmsy)) print(paste0("fmfmsy not numeric. fmfmsy: ", fmfmsy))
        if(all(is.numeric(fmfmsy))){
            fmfmsy <- round(fmfmsy,2)
        }else fmfmsy <- c(NA, NA)
        names(fmfmsy) <- c("fmfmsy.est","fmfmsy.sd")
        bpbmsy <- try(get.par("logBpBmsy",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(bpbmsy))){
            bpbmsy <- round(bpbmsy,2)
        }else bpbmsy <- c(NA,NA)
        names(bpbmsy) <- c("bpbmsy.est","bpbmsy.sd")
        cp <- try(get.par("logCp",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(cp))){
            cp <- round(cp,2)
        }else cp <- c(NA,NA)
        names(cp) <- c("cp.est","cp.sd")
        ##
        fmsy <- try(get.par("logFmsy",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(fmsy))){
            fmsy <- round(fmsy,2)
        }else fmsy <- c(NA,NA)
        names(fmsy) <- c("fmsy.est","fmsy.sd")
        bmsy <- try(get.par("logBmsy",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(bmsy))){
            bmsy <- round(bmsy,2)
        }else bmsy <- c(NA,NA)
        names(bmsy) <- c("bmsy.est","bmsy.sd")
        sdb <- try(get.par("logsdb",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(sdb))){
            sdb <- round(sdb,2)
        }else sdb <- c(NA,NA)
        names(sdb) <- c("sdb.est","sdb.sd")
        sdi <- try(get.par("logsdi",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(sdi))){
            sdi <- round(sdi,2)
            sdi <- as.numeric(t(sdi))
        }else{
            sdi <- rep(c(NA,NA),nis)
        }
        names(sdi) <- paste0(rep(c("sdi.est","sdi.sd"),nis),rep(1:nis,each=nis))
        sdf <- try(get.par("logsdf",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(sdf))){
            sdf <- round(sdf,2)
        }else sdf <- c(NA,NA)
        names(sdf) <- c("sdf.est","sdf.sd")
        sdc <- try(get.par("logsdc",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(sdc))){
            sdc <- round(sdc,2)
        }else sdc <- c(NA,NA)
        names(sdc) <- c("sdc.est","sdc.sd")
        bmbmsy <- try(get.par("logBmBmsy",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(bmbmsy))){
            bmbmsy <- round(bmbmsy,2)
        }else bmbmsy <- c(NA,NA)
        names(bmbmsy) <- c("bmbmsy.est","bmbmsy.sd")
        nest <- try(get.par("logn",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(nest))){
            nest <- round(nest,2)
        }else nest <- c(NA,NA)
        names(nest) <- c("n.est","n.sd")
        Kest <- try(get.par("logK",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(Kest))){
            Kest <- round(Kest,2)
        }else Kest <- c(NA,NA)
        names(Kest) <- c("K.est","K.sd")
        mest <- try(get.par("logm",fit, exp=TRUE)[,c(2,4)],silent=TRUE)
        if(all(is.numeric(mest))){
            mest <- round(mest,2)
        }else mest <- c(NA,NA)
        names(mest) <- c("m.est","m.sd")
        ##
        quantstmp <- c(fmfmsy, bpbmsy, cp, fmsy, bmsy, sdb, sdi, sdf, sdc, bmbmsy, nest, Kest, mest)

        if(is.na(fixn) &&
           reportmode %in% c(0,1) && any(is.na(quantstmp)[-which(names(quantstmp) == "sdi.sd1")])){
            ## n.sd is NA if schaefer (n fixed), reportmode 3 does not report most quantities
            tacs <- func(inp, tacs=tacs, pars=pars)
            tacs$conv[nrow(tacs)] <- FALSE
            tacs$indBref[nrow(tacs)] <- indBref2
            tacs$bmID[nrow(tacs)] <- bmID
            tacs$assessInt[nrow(tacs)] <- assessInt
            tacs$medbpbref[nrow(tacs)] <- medbpbref
            tacs$bpbref[nrow(tacs)] <- bpbref
            tacs$rwF[nrow(tacs)] <- rwF
            return(tacs)
        }

        ## resetting brefs at benchmark
        if(bmID){
            logB <- fit$obj$report(fit$obj$env$last.par.best)$logB[inp$indest]
            logB[1:(2/inp$dteuler)] <- NA    ## hack: remove first year, because first B est often outlier
            logB[is.infinite(logB)] <- NA    ## hack: remove first year, because first B est often outlier
            if(bref == "current"){
                indBref <- inp$indlastobs
            }else if(bref == "lowest"){
                timeX <- inp$time[inp$indest]
                ann <- spict:::annual(timeX,logB)
                ann2 <- ann$anntime[sort(doBy::which.minn(ann$annvec, nyBref))]
                indBref <- which(floor(timeX) %in% ann2)
            }else if(bref == "highest"){
                timeX <- inp$time[inp$indest]
                ann <- spict:::annual(timeX,logB)
                ann2 <- ann$anntime[sort(doBy::which.maxn(ann$annvec, nyBref))]
                indBref <- which(floor(timeX) %in% ann2)
            }else if(bref == "average"){
                indBref <- (2/inp$dteuler):length(logB)
            }else if(bref == "last"){
                indBref <- (length(logB)-(nyBref * 1/inp$dteuler)):length(logB)
            }else stop(paste0("bref = ",bref, " not known! Either current, lowest, highest, average, or last."))
        }else{
            indBref <- as.numeric(as.character(unlist(strsplit(as.character(tacs$indBref[nrow(tacs)]), "-"))))
        }
        if(any(is.na(indBref)) || length(indBref) < 1 || any(is.infinite(indBref))){
            tacs <- func(inp, tacs=tacs, pars=pars)
            tacs$conv[nrow(tacs)] <- FALSE
            tacs$indBref[nrow(tacs)] <- indBref2
            tacs$bmID[nrow(tacs)] <- bmID
            tacs$assessInt[nrow(tacs)] <- assessInt
            tacs$medbpbref[nrow(tacs)] <- medbpbref
            tacs$bpbref[nrow(tacs)] <- bpbref
            tacs$rwF[nrow(tacs)] <- rwF
            return(tacs)
        }else{
            fit <- try(set.bref(fit, indBref = indBref),silent=TRUE)
        }

        ## get TAC
        if(inherits(fit, "try-error")){
            tacs <- func(inp, tacs=tacs, pars=pars)
            tacs$conv[nrow(tacs)] <- FALSE
            tacs$indBref[nrow(tacs)] <- indBref2
            tacs$bmID[nrow(tacs)] <- bmID
            tacs$assessInt[nrow(tacs)] <- assessInt
            tacs$medbpbref[nrow(tacs)] <- medbpbref
            tacs$bpbref[nrow(tacs)] <- bpbref
            tacs$rwF[nrow(tacs)] <- rwF
            return(tacs)
        }

        ## Indicators
        ## -----------------------
        if("',btar,'" == "bmsy"){
            logBpBtar <- get.par("logBpBmsy", fit, exp = FALSE)
        }else if("',btar,'" == "btrigger"){
            logBpBtar <- get.par("logBpBtrigger", fit, exp = FALSE)
        }else stop("btar not known. Use either bmsy or btrigger")
        logFmFtar <- get.par("logFmFmsynotS", fit, exp = FALSE)
        bindi <- exp(qnorm(probtar, logBpBtar[2], logBpBtar[4]))
        findi <- exp(qnorm(1-probtar, logFmFtar[2], logFmFtar[4]))
        indBref2 <- paste(fit$inp$indBref, collapse="-")
        logBpBref <- get.par("logBpBref", fit, exp = FALSE)
        medbpbref <- as.numeric(exp(logBpBref[,2]))
        bpbref <- exp(qnorm(1-prob, logBpBref[2], logBpBref[4]))
        barID <- FALSE
        if(brule == 0){
            ## standard rule
            tac <- try(spict:::get.TAC(fit,
                                       ## bfac = bfac,                    ## not in pubFF yet
                                       ## bref.type = "',brefType,'",     ## not in pubFF yet
                                       fractiles = list(catch = ',frc,',
                                                        ffmsy = ',frff,',
                                                        bbmsy = ',frbb,',
                                                        bmsy  = ',frb,',
                                                        fmsy  = ',frf,'),
                                       breakpointB = c(',breakpointB1,',',
breakpointB2,'),
                                       safeguardB = list(limitB = ',limitB,',prob = prob),
                                       intermediatePeriodCatch = intC2,
                                       verbose = FALSE),
                       silent = TRUE)
        }else if(brule == 1){
            ## standard bref rule + pa buffer
            if(!is.numeric(bindi) || is.na(bindi) || !is.numeric(findi) || is.na(findi) ||
               !is.numeric(bpbref) || is.na(bpbref)){
                tacs <- func(inp, tacs=tacs, pars=pars)
                tacs$conv[nrow(tacs)] <- FALSE
                tacs$indBref[nrow(tacs)] <- indBref2
                tacs$bmID[nrow(tacs)] <- bmID
                tacs$assessInt[nrow(tacs)] <- assessInt
                tacs$medbpbref[nrow(tacs)] <- medbpbref
                tacs$bpbref[nrow(tacs)] <- bpbref
                tacs$rwF[nrow(tacs)] <- rwF
                return(tacs)
            }
            tac <- try(spict:::get.TAC(fit,
                                       bfac = bfac,
                                       bref.type = "',brefType,'",
                                       fractiles = list(catch = ',frc,',
                                                        ffmsy = ',frff,',
                                                        bbmsy = ',frbb,',
                                                        bmsy  = ',frb,',
                                                        fmsy  = ',frf,'),
                                       breakpointB = c(',breakpointB1,',',
breakpointB2,'),
                                       safeguardB = list(limitB = ',limitB,',prob = prob),
                                       intermediatePeriodCatch = intC2,
                                       verbose = FALSE),
                       silent = TRUE)

            ## PA buffer (e.g. 0.2 reduction of TAC) if B < Btrigger or F > Fmsy
            if(is.numeric(red)){
                if(any(as.logical(tail(tacs$barID,ceiling(redyears/assessInt-1))),na.rm=TRUE)){
                    ## do not apply if applied during last x years (redyears)
                    barID <- FALSE
                }else{
                    if(',redAlways,'){
                        barID <- TRUE
                    }else{
                        if((bindi - 1) < -1e-3 || (1 - findi) < -1e-3){
                            ## apply if any ref indicates overexploitation
                            barID <- TRUE
                        }else barID <- FALSE
                    }
                }
            }else barID <- FALSE

        }else if(brule == 2){
            ## decision tree using spict reference levels qualitatively
            if(!is.numeric(bindi) || is.na(bindi) || !is.numeric(findi) || is.na(findi) ||
               !is.numeric(bpbref) || is.na(bpbref)){
                tacs <- func(inp, tacs=tacs, pars=pars)
                tacs$conv[nrow(tacs)] <- FALSE
                tacs$indBref[nrow(tacs)] <- indBref2
                tacs$bmID[nrow(tacs)] <- bmID
                tacs$assessInt[nrow(tacs)] <- assessInt
                tacs$medbpbref[nrow(tacs)] <- medbpbref
                tacs$bpbref[nrow(tacs)] <- bpbref
                tacs$rwF[nrow(tacs)] <- rwF
                return(tacs)
            }
            ## 4 stock status categories
            ## -------------------------
            if((bpbref - bfac) < -1e-3){
                ## Overfished
                ## -----------------------
                ## -> find F that meets pi
                tac <- try(spict:::get.TAC(fit,
                                           bfac = bfac,
                                           bref.type = "',brefType,'",
                                           fractiles = list(catch = ',frc,',
                                                            ffmsy = ',frff,',
                                                            bbmsy = ',frbb,',
                                                            bmsy  = ',frb,',
                                                            fmsy  = ',frf,'),
                                           breakpointB = c(',breakpointB1,',',
breakpointB2,'),
                                           safeguardB = list(limitB = ',limitB,',prob = prob),
                                           intermediatePeriodCatch = intC2,
                                           verbose = FALSE),
                           silent = TRUE)
            }else if((bindi - 1) < -1e-3 || (1 - findi) < -1e-3){
                ## Indication of overfishing
                ## -----------------------
                ## -> reduce F/TAC by red%
                tac <- try(spict:::get.TAC(fit,
                                           ffac = 1 * (1-red),
                                           intermediatePeriodCatch = intC2,
                                           verbose = FALSE),
                           silent = TRUE)
            }else{
                ## No indication of overfishing
                ## -----------------------
                ## -> raise F/TAC by rai%
                tac <- try(spict:::get.TAC(fit,
                                           ffac = 1 * (1+rai),
                                           intermediatePeriodCatch = intC2,
                                           verbose = FALSE),
                           silent = TRUE)
            }
            ## if((bindi - 1) < -1e-3 || (1 - findi) < -1e-3){
            ##     ## Indication of overfishing
            ##     ## -----------------------
            ##     ## -> reduce TAC by red%
            ##     tac <- ',red,' * cl
            ## }else{
            ##     ## No indication of overfishing
            ##     ## -----------------------
            ##     ## -> raise TAC by rai%
            ##     tac <- ',rai,' * cl
            ## }
        }

        if(inherits(tac, "try-error") || !is.numeric(tac) || is.na(tac)){
            tacs <- func(inp, tacs=tacs, pars=pars)
            tacs$conv[nrow(tacs)] <- FALSE
            tacs$indBref[nrow(tacs)] <- indBref2
            tacs$bmID[nrow(tacs)] <- bmID
            tacs$assessInt[nrow(tacs)] <- assessInt
            tacs$medbpbref[nrow(tacs)] <- medbpbref
            tacs$bpbref[nrow(tacs)] <- bpbref
            tacs$rwF[nrow(tacs)] <- rwF
            return(tacs)
        }

        if(',stab,'){
            cllo <- cl * lower
            clup <- cl * upper
            if(any(tac < cllo)) hitSC <- 1 else hitSC <- 0
            if(any(tac > clup)) hitSC <- 2 else hitSC <- 0
            tac[tac < cllo] <- cllo
            tac[tac > clup] <- clup
        }else hitSC <- 0

        ## apply reduction
        if(barID){
            tac <- tac * (1-red)
        }
        tactmp <- data.frame(TAC=tac, id="',id,'", hitSC=hitSC,
                             red=red, barID=barID, sd=NA, conv = TRUE)
        tactmp <- data.frame(c(tactmp, quantstmp,
                               indBref = indBref2, bmID=bmID,
                               assessInt = assessInt,
                               medbpbref = medbpbref, bpbref = bpbref,
                               rwF = rwF))
        if(is.null(tacs)){
            tacs <- tactmp
        }else{
            tacs <- rbind(tacs, tactmp)
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




#' @name def.hcr.sam
#' @title Define harvest control rule for SAM
#'
#' @param id Name/ID of HCR. Default: "sam"
#' @param nonconvHCR HCR if SAM does not converge. Default: "conscat" (constant catch).
#' @param silent silent
#' @param verbose verbose
#' @param env Environment. Default: globalenv()
#'
#'
#' @export
#'
def.hcr.sam <- function(id = "sam",
                        nonconvHCR = "conscat",
                        silent = TRUE,
                        verbose = FALSE,
                        env = globalenv()
                        ){

    template  <- expression(paste0('
structure(function(obs, tacs = NULL, pars=NULL){
    silent <- ',silent,'
    verbose <- ',verbose,'
    func <- get("',nonconvHCR,'")

    ## setup SAM data
    dat <- stockassessment::setup.sam.data(surveys = obs$obsIA,
                                           residual.fleet = obs$obsCA,
                                           prop.mature = obs$propMature,
                                           stock.mean.weight = obs$WAAs,
                                           catch.mean.weight = obs$WAAc,
                                           dis.mean.weight = obs$WAAc,
                                           land.mean.weight = obs$WAAc,
                                           prop.f = obs$propFemale,
                                           prop.m = obs$propFemale,
                                           natural.mortality = obs$obsMAA,
                                           land.frac = obs$landFrac)

    ## configurations
    conf <- stockassessment::defcon(dat)

    ## starting values
    par <- stockassessment::defpar(dat, conf)

    ## fit SAM (to make faster re-code sam.fit)
    fit <- try(stockassessment::sam.fit(dat, conf, par, ignore.parm.uncertainty = TRUE,
                                        rel.tol = 1e-6,
                                        silent=silent),
               silent=TRUE)

    if(class(fit) == "try-error"){
        if(verbose) cat("Error in model fitting.\n")
        tacs <- func(obs, tacs=tacs, pars=pars)
        tacs$conv[nrow(tacs)] <- FALSE
        return(tacs)
    }else{
        ## reference levels (using simEQ)
        ## estimating fval
        ## fval = Fmsy based on simEQ

        ## sam forecast given Fval
        ## fore = forecast(fit, nextssb = c(NA,ssbref), fval = c(1,NA))  ## doesnt work
        fore <- try(stockassessment::forecast(fit, fval = c(1,1)), silent = TRUE)

        if(class(fore) == "try-error"){
            if(verbose) cat("Error in model fitting.\n")
            tacs <- func(obs, tacs=tacs, pars=pars)
            tacs$conv[nrow(tacs)] <- FALSE
            return(tacs)

        }else{

            ## predicted catches
            shorttab <- attr(fore,"shorttab")
            tac <- shorttab[4,2]

            ## plots
            ## ssbplot(fit)
            ## catchplot(fit)
            ## dataplot(fit)
            ## fitplot(fit)

            ## write output object
            tacs <- gettacs(tacs, id="',id,'", TAC=tac, obs=obs)
            tacs$conv[nrow(tacs)] <- TRUE
            return(tacs)
        }
    }
},
class="hcr")
'))

    ## create HCR as functions
    templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}



#' @name def.hcr.pseudo
#' @title Define harvest control rule with pseudo assessment
#'
#' @param id Name/ID of HCR. Default: "pseudo-msy"
#' @param fractiles Fractiles. List
#' @param breakpointB breakpointb
#' @param safeguard safeguard. List
#' @param dteuler Time step of Forward Euler scheme.
#' @param reportmode Reportmode
#' @param stabilise stabilise option of spict. Default: FALSE.
#' @param priorlogn prior
#' @param priorlogalpha prior
#' @param priorlogbeta prior
#' @param priorlogbkfrac prior
#' @param fixn fixn
#' @param bfac bfac
#' @param bref bref
#' @param brefType brefType
#' @param nyBref nyBref
#' @param btar btar
#' @param probtar probtar
#' @param brule Brule
#' @param red Reduction
#' @param redyears Reduction years
#' @param redAlways Reduction always? Default: FALSE.
#' @param rai Raising factor. Default: 0.2.
#' @param manstartdY Management start year. Default: 0.
#' @param assessmentInterval Assessment interval. Default: 1.
#' @param intC Interval C. Default: NA.
#' @param nonconvHCR HCR if SPiCT does not converge. Default: "conscat" (constant catch).
#' @param clType Catch type for uncertainty cap. Default: "TAC".
#' @param clyears Number of years for uncertainty cap. Default: 1.
#' @param stab Uncertainty cap. Default: FALSE.
#' @param lower Upper bound of uncertainty cap. Default: 0.8.
#' @param upper Upper bound of uncertainty cap. Default: 1.2.
#' @param bm Benchmark. Default: FALSE,  ## e.g. 5 => re-defining Bref at every benchmark
#' @param env Environment. Default: globalenv()
#'
#'
#' @importFrom doBy which.minn
#'
#' @export
#'
def.hcr.pseudo <- function(id = "pseudo-msy",
                          fractiles = list(catch=0.5,
                                           ffmsy=0.5,
                                           bbmsy=0.5,
                                           bmsy = 0.5,
                                           fmsy = 0.5),
                          breakpointB = 0.0,
                          clType = "TAC",
                          clyears = 1,
                          stab = FALSE,
                          lower = 0.8,
                          upper = 1.2,
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
    breakpointB <- sort(breakpointB)
    if(length(breakpointB) > 1){
        blim <- breakpointB[1]
        btrigger <- breakpointB[2]
    }else{
        blim <- 0
        btrigger <- breakpointB[1]
    }
    ## Knife-edge hockey-stick
    flagKE <- ifelse(blim == btrigger, TRUE, FALSE)

    template  <- expression(paste0(
        '
structure(
    function(obs, tacs = NULL, pars=NULL){


        inp <- obs[c("obsC","timeC")]
        bbmsy <- pars[["bbmsy"]]
        bbmsySD <- pars[["bbmsySD"]]
        bbmsyBias <- pars[["bbmsyBias"]]
        ffmsy <- pars[["ffmsy"]]
        ffmsySD <- pars[["ffmsySD"]]
        ffmsyBias <- pars[["ffmsyBias"]]
        tacSD <- pars[["tacSD"]]
        fmsy <- pars[["fmsy"]]
        fmsyBias <- pars[["fmsyBias"]]

        lower <- ',lower,'
        upper <- ',upper,'
        clyears <- ',clyears,'
        clType <- "',clType,'"

        ## pseudo-assessment-hcr
        ffmsyi <- exp(qnorm(1 - frf, log(ffmsy + ffmsy * ffmsyBias), ffmsySD))
        ffmsy5 <- exp(qnorm(0.5, log(ffmsy + ffmsy * ffmsyBias), ffmsySD))
        fred <- ffmsy5 / ffmsyi
        if(!flagKE){
        if(btrigger > 0){
           hsSlope <- 1/(btrigger-blim)
           hsIntercept <- - hsSlope * blim
           bbmsyi <- hsSlope * exp(qnorm(frb, log(bbmsy + bbmsy * bbmsyBias),
                                         bbmsySD)) + hsIntercept
           fred <- fred * min(1, max(0,bbmsyi))
        }
}else{
        if(btrigger > 0){
           bbmsyi <- 1/blim * exp(qnorm(frb, log(bbmsy + bbmsy * bbmsyBias),
                                         bbmsySD))
           fred <- fred * ifelse(bbmsyi < 1, 0, 1)
        }
}
        targetF <- (fred + 1e-8) * (fmsy + fmsy * fmsyBias) / pars[["ns"]]

        faa <- targetF * pars[["sel"]]
        caa <- baranov(faa, pars[["m"]], pars[["n"]])
        tac <- sum(caa * pars[["weight"]])
        tac <- exp(qnorm(frc, log(tac), tacSD))

        ## last years catch/tac
        if(clType == "observed"){
            cl <- mean(tail(inp$obsC, clyears))
        }else if(clType == "TAC"){
            if(is.null(tacs)){
                cl <- mean(tail(inp$obsC, clyears))
                cl <- cl ## * assessInt
            }else{
                cl <- tacs$TAC[nrow(tacs)]
            }
        }

        if(',stab,'){
            cllo <- cl * lower
            clup <- cl * upper
            if(any(tac < cllo)) hitSC <- 1 else hitSC <- 0
            if(any(tac > clup)) hitSC <- 2 else hitSC <- 0
            tac[tac < cllo] <- cllo
            tac[tac > clup] <- clup
        }else hitSC <- 0

        tacs <- gettacs(tacs, id="',id,'", TAC=tac, obs=obs)
        tacs$hitSC[nrow(tacs)] <- hitSC
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
