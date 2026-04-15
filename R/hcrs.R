
#' est.tac
#'
#' @param obs. observation
#' @param hcr.fun HCR function
#' @param tacs. tac daa frame
#' @param pars. parameters
#'
#' @export
est.tac <- function(obs., hcr.fun, tacs.=NULL, pars.=NULL) {

    ## func <- get(hcr.)
    ## res <- func(obs., tacs., pars.)
    ## return(res)
    hcr.fun(obs., tacs., pars.)

}



#' gettacs
#'
#' @param tacs. info
#' @param id. info
#' @param TAC. info
#' @param obs. info
#'
gettacs <- function(tacs.=NULL, id.="", TAC.=NA, obs.=NULL) {

    if(!is.null(obs.) && is.list(obs.$obsI) && length(obs.$obsI) != 0)
        nis <- length(obs.$obsI) else nis <- 1
    names(TAC.) <- paste0("TAC", seq_len(length(TAC.)))
    tactmp <- data.frame(t(TAC.), id=id.,
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
    if(is.null(tacs.)){
        tacs. <- tactmp
    }else{
        tacs. <- rbind(tacs., tactmp)
    }
    return(tacs.)
}




#' def.hcr.ref
#'
#' @param consF either numeric indicating constant F level or "fmsy" for fishing at fmsy
#' @param fracFmsy fraction of fmsy
#' @param set. settings
#' @param env environment
#'
#' @export
#'
def.hcr.ref <- function(consF = 0,
                        fracFmsy = NULL,
                        set. = check.set(),
                        env = globalenv()
                        ){

    id <- NULL

    if(!is.null(fracFmsy)){
        id <- paste0("refFmsy_", as.character(fracFmsy))
        template  <- expression(paste0(
            'function(obs, tacs=NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }
        ## inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- iamse:::gettacs(tacs.=tacs, id.="',id,'",
                        TAC.=rep(0, ',set.$assessmentInterval,'), obs.=inp)
        return(tacs)
    }'))
    }else{
        if(!is.numeric(consF) && !(consF %in% c("fmsy","Fmsy","FMSY","refFmsy","refFMSY","reffmsy")))
            stop("'consF' has to be either the absolute F (numeric) or 'fmsy'.")

        if(is.numeric(consF)){
            ## TODO: don't think this is implemented, check mse.R
            id <- paste0("consF_",consF)
            template  <- expression(paste0(
                'function(obs, tacs=NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }
        ## inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- iamse:::gettacs(tacs.=tacs, id.="',id,'",
                        TAC.=rep(0, ',set.$assessmentInterval,'), obs.=inp)
        return(tacs)
    }'))
        }

        if(consF == "fmsy" || consF == "Fmsy" || consF == "FMSY"){
            id <- "refFmsy"
            template  <- expression(paste0(
                'function(obs, tacs=NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }
        ## inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- iamse:::gettacs(tacs.=tacs, id.="',id,'",
                        TAC.=rep(NA, ',set.$assessmentInterval,'), obs.=inp)
        return(tacs)
    }'))
        }

        if(consF == 0){
            id <- "noF"
            template  <- expression(paste0(
                'function(obs, tacs=NULL, pars=NULL){
        if(is.null(obs$timeE)){
             inp <- obs[c("obsC","timeC","obsI","timeI")]
        }else if(is.null(obs$timeI)){
             inp <- obs[c("obsC","timeC","obsE","timeE")]
        }else{
             inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
        }
        ## inp <- spict::check.inp(inp, verbose = FALSE)
        tacs <- iamse:::gettacs(tacs.=tacs, id.="',id,'",
                        TAC.=rep(0, ',set.$assessmentInterval,'), obs.=inp)
        return(tacs)
    }'))
        }
    }


    ## create HCR as functions
    ## templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    ## assign(value=templati, x=id, envir=env)

    templati <- eval(parse(text = eval(template)))
    class(templati) <- c(class(templati), "hcr")
    attributes(templati)$id <- id
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}



#' Define a constant-catch harvest control rule
#'
#' `def.hcr.conscat()` defines a simple harvest control rule (HCR) that
#' recommends a constant catch (TAC) over time. The baseline TAC can be
#' supplied directly via `constantC` or estimated from recent catches, and
#' optional reduction rules can be used to step down the TAC when performance
#' is poor (e.g. low biomass or repeated assessment failure).
#'
#' The function is usually called for its side effects: it registers a new HCR
#' with identifier `id` in the `iamse` settings object `set.` and, by default,
#' also in the supplied environment `env`. The HCR can then be included in
#' `set$hcr` and evaluated with [run.mse()].
#'
#' @param id Character string giving the identifier for this HCR (for example
#'   `"conscat"`). This label is used when selecting HCRs in `set$hcr` and in
#'   summaries and plots.
#' @param constantC Optional numeric giving the constant catch level (e.g. in
#'   tonnes). If `NULL` (default), the constant catch is derived internally
#'   from recent catches, typically using the last `clyears` years.
#' @param clyears Integer specifying the number of historical years used to
#'   calculate the baseline catch level when `constantC` is `NULL`.
#' @param red Optional numeric reduction factor applied to the constant catch
#'   (e.g. `0.8` corresponds to a 20% reduction). If `NA` (default), no
#'   automatic reduction is applied and the TAC remains constant.
#' @param redyears Integer giving the number of years over which the reduction
#'   criteria are evaluated. The exact interpretation depends on the internal
#'   implementation (e.g. how many years of poor performance trigger a
#'   reduction).
#' @param redAlways Logical; if `TRUE`, the reduction rule is evaluated in all
#'   years where the criteria are met. If `FALSE`, the reduction may be applied
#'   only once or under more restrictive conditions, depending on the internal
#'   implementation.
#' @param assessmentInterval Positive integer giving the interval (in years)
#'   between assessments and TAC updates. The default `1` corresponds to an
#'   annual advice cycle.
#' @param ffmsySD Non-negative numeric giving the standard deviation of the
#'   (log) assessment error in \eqn{F/F_{MSY}} used when simulating the assessed
#'   status for this HCR.
#' @param bbtriggerSD Non-negative numeric giving the standard deviation of the
#'   (log) assessment error in biomass-related indicators (e.g. \eqn{B/B_{trigger}})
#'   used for this HCR.
#' @param rightRef Integer index specifying which reference-point set should be
#'   used when evaluating this HCR (e.g. an index into `dat$ref`).
#' @param set. An `iamse` settings list as returned by [check.set()]. This
#'   object is updated to include the newly defined HCR.
#' @param env Environment in which the HCR function is created and stored.
#'   Defaults to [globalenv()]. In most cases the default should be used.
#'
#' @return Invisibly returns the updated `set.` object. The function is called
#'   for its side effects (defining and registering the HCR).
#'
#' @export
def.hcr.conscat <- function(id = "conscat",
                            constantC = NULL,
                            clyears = 1,
                            red = NA,
                            redyears = 2,
                            redAlways = FALSE,
                            assessmentInterval = 1,
                            ffmsySD = 0,
                            bbtriggerSD = 0,
                            rightRef = 1,
                            set. = check.set(),
                            env = globalenv()
                            ){

    if (!requireNamespace("spict", quietly = TRUE)) {
        stop("The 'spict' package is required for this function. ",
             "Please install it from GitHub: DTUAqua/spict/spict.", call. = FALSE)
    }

    if(is.null(constantC)) constantC = NA

    template  <- expression(paste0(
        'function(obs, tacs = NULL, pars=NULL){
        red <- ',red,'
        redyears <- ',redyears,'
        assessInt <- ',set.$assessmentInterval,'

        ffmsy <- rnorm(1, pars$ffmsy, ',ffmsySD,')
        ## ffmsy <- runif(1, pars$ffmsy * ',ffmsySD,', pars$ffmsy)
        ffmsy[ffmsy < 0] <- 0
        bbtrigger <- rnorm(1, pars$bbmsy*2, ',bbtriggerSD,')
        ## bbtrigger <- runif(1, pars$bbmsy*2, pars$bbmsy*2 * ',bbtriggerSD,')
        bbtrigger[bbtrigger < 0] <- 0

        obs <- spict::check.inp(obs, verbose = FALSE)
        if(is.null(tacs)){
            indBref <- obs$indBref
        }else{
            indBref <- as.numeric(as.character(unlist(strsplit(as.character(tacs$indBref[nrow(tacs)]), "-"))))
        }
        indBref2 <- paste(indBref, collapse="-")


        tac <- ',constantC,'
        if(!is.numeric(tac)){
            annualcatch <- spict:::annual(obs$timeC, obs$obsC/obs$dtc, type = "mean") ## CHECK: why not sum?
            tac <- mean(tail(annualcatch$annvec, ',clyears,'))
            ## accounted for by rep(tac, assessInt) below
            ## ## Account for non-annual assessments
            ## tac <- tac * assessInt
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

        tac <- rep(tac, ',set.$assessmentInterval,')
        tacs <- iamse:::gettacs(tacs.=tacs, id. = "',id,'", TAC. = tac, obs.=obs)
        tacs$hitSC[nrow(tacs)] <- NA
        tacs$barID[nrow(tacs)] <- barID
        tacs$red[nrow(tacs)] <- red
        tacs$indBref[nrow(tacs)] <- indBref2
        tacs$assessInt[nrow(tacs)] <- assessInt
        return(tacs)
    }'))

    ## create HCR as functions
    ## templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    ## assign(value=templati, x=id, envir=env)

    templati <- eval(parse(text = eval(template)))
    class(templati) <- c(class(templati), "hcr")
    attributes(templati)$id <- id
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}




#' Define an index-based harvest control rule
#'
#' `def.hcr.index()` defines one or more harvest control rules (HCRs) that use
#' changes in a relative biomass index to set catch advice (TAC) or fishing
#' mortality. The function constructs HCR functions and assigns them to the
#' specified environment, typically the global environment, so they can be
#' used within an MSE framework.
#'
#' The generated HCR functions usually take arguments of the form
#' `function(x, Data, reps, ...)`, where `Data` contains the catch and index
#' time series and `reps` is the number of stochastic TAC draws (if used).
#' Internally, the rule uses recent values of a biomass index, applies
#' stabilisation bounds on year-to-year changes, and can optionally reduce
#' catch when biomass is low.
#'
#' Several arguments (e.g. `id`, `x`, `y`) can be provided as vectors of equal
#' length, in which case multiple HCRs are created in a single call.
#'
#' @param id Character vector giving the identifier(s) for the HCR(s)
#'   (e.g. `"r23"`). These names are used when referring to the HCR in the
#'   MSE setup and outputs.
#' @param x Numeric scalar or vector used to parameterise the index-based
#'   rule. Together with `y`, it controls how changes in the biomass index
#'   translate into changes in TAC (e.g. slope, threshold, or tuning factor).
#' @param y Numeric scalar or vector used alongside `x` to tune the response
#'   of the HCR to the biomass index. If `x` and `y` are vectors, they must
#'   have the same length.
#' @param stab Logical; if `TRUE` (default), apply stabilisation bounds to
#'   limit year-to-year changes in the TAC or fishing mortality.
#' @param lower Numeric value giving the lower stabilisation bound, typically
#'   interpreted as the minimum allowed ratio \eqn{\text{TAC}_{y+1} /
#'   \text{TAC}_y}. The default is `0.8` (maximum 20\% decrease).
#' @param upper Numeric value giving the upper stabilisation bound, typically
#'   interpreted as the maximum allowed ratio \eqn{\text{TAC}_{y+1} /
#'   \text{TAC}_y}. The default is `1.2` (maximum 20\% increase).
#' @param clType Character string specifying the type of control variable
#'   returned by the HCR. Common options include `"TAC"` for catch limits or
#'   `"F"` for fishing mortality targets. The default is `"TAC"`.
#' @param clyears Integer specifying the number of recent years used to
#'   calculate the reference catch level when deriving TAC advice (e.g. mean
#'   catch over the last `clyears` years).
#' @param red Optional numeric reduction factor applied to the catch when
#'   biomass is low (e.g. `0.8` for a 20\% reduction). If `NA` (default), no
#'   automatic biomass-based reduction is applied.
#' @param redyears Integer giving the number of years over which the biomass
#'   condition for reduction is evaluated (e.g. how many years below a
#'   threshold before the reduction is triggered).
#' @param redAlways Logical; if `TRUE`, the reduction rule is evaluated and
#'   applied whenever the criteria are met. If `FALSE`, the reduction may only
#'   be applied under more restrictive conditions, depending on the internal
#'   implementation.
#' @param ffmsySD Non-negative numeric specifying the standard deviation of
#'   the (log) assessment error in \eqn{F/F_{MSY}} associated with this HCR,
#'   when simulating perceived stock status.
#' @param bbtriggerSD Non-negative numeric specifying the standard deviation
#'   of the (log) assessment error in biomass-related indicators (e.g.
#'   \eqn{B/B_{trigger}}) associated with this HCR.
#' @param rightRef Integer index specifying which set of reference points
#'   should be used when evaluating this HCR (e.g. an index into a reference
#'   table).
#' @param assessmentInterval Positive integer giving the interval (in years)
#'   between full assessments and updates of the advice. The default `1`
#'   corresponds to annual updates.
#' @param env Environment in which the generated HCR function(s) are created
#'   and stored. Defaults to [globalenv()].
#'
#' @return A character vector with the name(s) of the HCR function(s) created.
#'   The main purpose of the function is its side effect of defining these HCRs
#'   in `env`.
#'
#' @examples
#' \dontrun{
#'   ## Define a single index-based HCR with default tuning
#'   hcr_names <- def.hcr.index(id = "index_HCR")
#'   hcr_names
#'
#'   ## Define several index-based HCRs at once with different tuning
#'   hcr_names2 <- def.hcr.index(
#'     id = c("index_low", "index_high"),
#'     x  = c(0.8, 1.0),
#'     y  = c(1.2, 1.5)
#'   )
#'   hcr_names2
#'
#'   ## The created functions are available in the chosen environment
#'   ls(pattern = "index_", envir = globalenv())
#' }
#'
#' @export
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
                          env = globalenv()
                          ){

    if (!requireNamespace("spict", quietly = TRUE)) {
        stop("The 'spict' package is required for this function. ",
             "Please install it from GitHub: DTUAqua/spict/spict.", call. = FALSE)
    }

    template  <- expression(paste0(
        'function(obs, tacs = NULL, pars = NULL){
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

        obs <- spict::check.inp(obs, verbose = FALSE)
        if(is.null(tacs)){
            indBref <- obs$indBref
        }else{
            indBref <- as.numeric(as.character(unlist(strsplit(as.character(tacs$indBref[nrow(tacs)]), "-"))))
        }
        indBref2 <- paste(indBref, collapse="-")
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
                cl <- tacs[nrow(tacs),paste0("TAC",assessInt)]
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
        tacs <- iamse:::gettacs(tacs.=tacs, id. = "',id,'", TAC. = tac, obs. = obs)
        tacs$hitSC[nrow(tacs)] <- hitSC
        tacs$barID[nrow(tacs)] <- barID
        tacs$red[nrow(tacs)] <- red
        tacs$fmfmsy.est[nrow(tacs)] <- ffmsy
        tacs$bpbmsy.est[nrow(tacs)] <- bbtrigger
        tacs$fmfmsy.sd[nrow(tacs)] <- ffmsySD
        tacs$bpbmsy.sd[nrow(tacs)] <- bbtriggerSD
        tacs$n.est[nrow(tacs)] <- r0
        tacs$indBref[nrow(tacs)] <- indBref2
        tacs$bmID[nrow(tacs)] <- bmID
        tacs$assessInt[nrow(tacs)] <- assessInt
        return(tacs)
    }'))

    ## create HCR as functions
    ## templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    ## assign(value=templati, x=id, envir=env)

    templati <- eval(parse(text = eval(template)))
    class(templati) <- c(class(templati), "hcr")
    attributes(templati)$id <- id
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}





#' Define a SPiCT-based harvest control rule
#'
#' `def.hcr.spict()` defines one or more harvest control rules (HCRs) that use
#' the SPiCT model (stochastic surplus production model in continuous time)
#' to generate catch advice. For each call, a function is created that:
#'
#' 1. Fits or updates a SPiCT assessment using catch and index data,
#' 2. Extracts stock status and predicted catch distributions, and
#' 3. Translates these into TAC recommendations according to specified
#'    fractiles, biomass breakpoints, safeguards, and stabilisation rules.
#'
#' The generated HCR function(s) are assigned to the environment `env`
#' (typically the global environment) and can then be used within a management
#' strategy evaluation (MSE) framework.
#'
#' @param id Character string giving the name/ID of the HCR. Default is
#'   `"spict-msy"`. If vectors are supplied for some arguments, multiple HCRs
#'   can be created in one call, each with its own `id`.
#' @param fractiles List of fractiles used to derive advice from SPiCT output.
#'   Typical elements include quantiles for predicted catch and reference
#'   points, e.g. `list(catch = 0.5, ffmsy = 0.5, bbmsy = 0.5, bmsy = 0.5,
#'   fmsy = 0.5)`. Values are probabilities in \[0, 1].
#' @param breakpointB Numeric biomass breakpoint (on the scale of
#'   \eqn{B/B_{MSY}} or related indicator) that can modify the rule below a
#'   given biomass level (e.g. more conservative advice when biomass is low).
#' @param evalBreakpointB Numeric or integer specifying how the breakpoint
#'   should be evaluated (e.g. over which period or using which biomass
#'   indicator). The exact interpretation depends on the internal
#'   implementation.
#' @param safeguard List describing an additional biomass safeguard, typically
#'   of the form `list(limitB = <threshold>, prob = <probability>)`. For
#'   example, the rule may reduce or cap TAC if the probability that biomass is
#'   below `limitB` exceeds `1 - prob`.
#' @param dteuler Time step of the forward Euler scheme used by SPiCT to solve
#'   the underlying continuous-time model. Default is `1/4` (quarterly).
#' @param reportmode Integer flag passed to SPiCT controlling the level of
#'   detail in reporting and output.
#' @param stabilise Numeric or logical argument passed to SPiCT's `stabilise`
#'   option, controlling how strongly the SPiCT fit is stabilised (e.g. in the
#'   presence of limited data). Default is `0` (no additional stabilisation).
#'
#' @param priorlogn Numeric vector giving the prior for the initial depletion
#'   (log scale), typically of the form `c(mean, sd, p)` where `p` is an
#'   indicator controlling whether the prior is active.
#' @param priorlogsdf Numeric vector giving the prior for the process-noise
#'   standard deviation on fishing mortality (log scale).
#' @param priorlogsdc Numeric vector giving the prior for the process-noise
#'   standard deviation on catches (log scale).
#' @param priorlogalpha Numeric vector giving the prior for the production
#'   function exponent \eqn{\alpha} (log scale).
#' @param priorlogbeta Numeric vector giving the prior for the production
#'   function exponent \eqn{\beta} (log scale).
#' @param priorlogbkfrac Numeric vector giving the prior for the fraction of
#'   carrying capacity at which biomass starts (log scale).
#' @param priorlogr Numeric vector giving the prior for the intrinsic growth
#'   rate \eqn{r} (log scale).
#'
#' @param fixn Optional value indicating whether the depletion parameter
#'   (`n`) should be fixed in the SPiCT fit. Default `NA` means it is
#'   estimated.
#' @param bfac Optional numeric factor used to scale biomass-based reference
#'   points when computing advice. Default `NA` means no additional scaling.
#' @param bref Character string specifying how the biomass reference level is
#'   defined, e.g. `"current"`, `"lowest"`, `"highest"`, `"average"` or
#'   `"last"`. This reference is used together with `brefType`, `nyBref`,
#'   and `btar` to construct biomass-based triggers.
#' @param brefType Character string specifying the type of biomass reference,
#'   typically `"target"` or `"limit"`. This affects how conservative the rule
#'   becomes when biomass is near or below the reference.
#' @param nyBref Integer giving the number of years used to define the biomass
#'   reference (e.g. average biomass over the last `nyBref` years).
#' @param btar Character string specifying the target biomass level used in
#'   the rule, e.g. `"bmsy"` for biomass at MSY.
#' @param probtar Numeric probability (in \[0, 1]) indicating how conservative
#'   the target rule should be (e.g. requiring that the probability of being
#'   above `btar` exceeds `probtar`).
#' @param brule Numeric code or scaling factor controlling how the biomass
#'   rule is applied below the target (e.g. linear reduction, more aggressive
#'   reductions, or no reduction).
#'
#' @param red Optional numeric reduction factor applied to the advised catch
#'   when additional criteria are met (e.g. low biomass, repeated assessment
#'   failure). For example, `0.8` corresponds to a 20\% reduction. If `NA`,
#'   no such reduction is applied.
#' @param redyears Integer indicating over how many years the reduction
#'   criteria are evaluated (e.g. number of years below a threshold before the
#'   reduction is triggered).
#' @param redAlways Logical; if `TRUE`, the reduction rule is evaluated and
#'   applied every time the criteria are met. If `FALSE`, the reduction may be
#'   applied more sparingly, depending on the internal implementation.
#' @param rai Numeric "raising factor" that controls how quickly the rule
#'   responds when moving back towards higher catches (e.g. after a reduction).
#'   Default is `0.2`.
#'
#' @param manstartdY Numeric management start year offset relative to the
#'   assessment time series (e.g. `0` for the last year, negative values for
#'   earlier years). Default is `0`.
#' @param assessmentInterval Positive integer giving the interval (in years)
#'   between full SPiCT assessments and advice updates. Default is `1`, i.e.
#'   annual assessments.
#' @param cpvec Logical; if `TRUE`, use a vector of catch profiles or control
#'   points when computing advice. Default is `FALSE`.
#' @param halfAssessInt Logical; if `TRUE`, the assessment is effectively
#'   shifted by half an assessment interval (e.g. mid-year updates). Default
#'   is `FALSE`.
#' @param intC Numeric or `NA`. If non-`NA`, this can be used to specify a
#'   predefined catch level in an initial period or over an interval. Default
#'   is `NA`.
#'
#' @param nonconvHCR Character string giving the name of an alternative HCR to
#'   be used if the SPiCT assessment does not converge. Default is `"conscat"`,
#'   which typically corresponds to a constant-catch rule defined by
#'   [def.hcr.conscat()].
#'
#' @param clType Character string specifying the type of control variable
#'   returned by the HCR, e.g. `"TAC"` for catch advice. Used in combination
#'   with the stabilisation cap.
#' @param clyears Integer giving the number of years used to compute the
#'   reference level for the stabilisation cap (e.g. mean catch over the last
#'   `clyears` years). Default is `1`.
#' @param stab Logical; if `TRUE`, apply an uncertainty or stabilisation cap
#'   that limits year-to-year changes in the advice (e.g. TAC). Default is
#'   `FALSE`.
#' @param lower Numeric lower bound for the stabilisation cap, typically
#'   interpreted as the minimum allowed ratio \eqn{\text{Advice}_{y+1} /
#'   \text{Advice}_y}. Default is `0.8` (maximum 20\% decrease).
#' @param upper Numeric upper bound for the stabilisation cap, typically
#'   interpreted as the maximum allowed ratio \eqn{\text{Advice}_{y+1} /
#'   \text{Advice}_y}. Default is `1.2` (maximum 20\% increase).
#'
#' @param bm Logical or numeric. If `FALSE` (default), no explicit benchmarks
#'   are applied. If a positive integer is supplied (e.g. `5`), the biomass
#'   reference `Bref` may be redefined every `bm` years (benchmarking).
#' @param env Environment in which the generated HCR function(s) are created
#'   and stored. Defaults to [globalenv()].
#' @param scenario Optional label (character or `NA`) used to identify the
#'   scenario associated with this HCR in subsequent summaries and plots.
#'
#' @return A character vector containing the name(s) of the HCR function(s)
#'   created. The main purpose of the function is its side effect of defining
#'   SPiCT-based HCRs in `env`.
#'
#' @details
#' This function assumes that the **spict** package is available and that the
#' underlying catch and index data can be passed to SPiCT in the expected
#' format. Many arguments (e.g. priors and Euler time step) are passed directly
#' to SPiCT's input list and control how the assessment behaves. For details on
#' SPiCT itself, see the **spict** documentation.
#'
#' @importFrom doBy which.minn
#'
#' @export
def.hcr.spict <- function(id = "spict-msy",
                          fractiles = list(catch=0.5,
                                           ffmsy=0.5,
                                           bbmsy=0.5,
                                           bmsy = 0.5,
                                           fmsy = 0.5),
                          breakpointB = 0.0,
                          evalBreakpointB = 0,
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
                          priorlogr = c(log(0.4),2,0),
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
                          cpvec = FALSE,
                          halfAssessInt = FALSE,
                          intC = NA,
                          nonconvHCR = "conscat",
                          clType = "TAC",
                          clyears = 1,
                          stab = FALSE,
                          lower = 0.8,
                          upper = 1.2,
                          bm = FALSE,  ## e.g. 5 => re-defining Bref at every benchmark
                          env = globalenv(),
                          scenario = NA
                          ) {

    if (!requireNamespace("spict", quietly = TRUE)) {
        stop("The 'spict' package is required for this function. ",
             "Please install it from GitHub: DTUAqua/spict/spict.", call. = FALSE)
    }

    if(is.null(scenario)) stop("Set scenario to NA!")

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

  template  <- expression(paste0('function(obs, tacs = NULL, pars=NULL){

if(is.null(obs$timeE)){

inp <- obs[c("obsC","timeC","obsI","timeI")]

}else if(is.null(obs$timeI)){

inp <- obs[c("obsC","timeC","obsE","timeE")]

}else{

inp <- obs[c("obsC","timeC","obsI","timeI","obsE","timeE")]
}

        func2 <- get("',nonconvHCR,'")
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
        if(',halfAssessInt,'){
            inp$maninterval <- c(manstart, manstart + 1)
        }else{
            inp$maninterval <- c(manstart, manstart + assessInt)
        }
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
        inp$priors$logr <- c(',priorlogr[1],',',priorlogr[2],',',priorlogr[3],')
        ## Catch for intermediate year
        if(is.na(',intC,')){
            intC2 <- NULL
        }else{
            if(is.null(tacs)){
                intC2 <- tail(inp$obsC,',manstartdY,') ## does not account for seasonal catches!
            }else{
                intC2 <- tacs[nrow(tacs),paste0("TAC",assessInt)] * ',
        manstartdY,'
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
## save.fits <- TRUE ## HERE:
## if(save.fits){
## inp$reportmode <- 0
## }
## if(',scenario,' == 1){
## if(inp$nobsC > 50) inp$stdevfacC[51] <- 0.1
## if(inp$nobsC > 51) inp$stdevfacC[52] <- 0.1
## }
## if(',scenario,' == 2){
## inp$priors$logsdf <- c(2,0.5,1)
## inp$priors$logsdc <- c(log(0.15),0.5,1)
## }
        fit <- try(fit.spict(inp), silent=TRUE)
## if(save.fits){
## save(fit, file = paste0("fit_",inp$nobsC,"_",',scenario,',".RData"))
## }
##  try(plot(fit),silent=TRUE)
        ## non-convergence
        if(class(fit) == "try-error" || fit$opt$convergence != 0 || any(is.infinite(fit$sd))){
            tacs <- func2(inp, tacs, pars)
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
                cl <- tacs[nrow(tacs),paste0("TAC",assessInt)]
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
            tacs <- func2(inp, tacs, pars)
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
            tacs <- func2(inp, tacs, pars)
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
            tacs <- func2(inp, tacs, pars)
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
                                       evalBreakpointB = ',evalBreakpointB,',
                                       cpvec = ',cpvec,',
                                       safeguardB = list(limitB = ',limitB,',prob = prob),
                                       intermediatePeriodCatch = intC2,
                                       verbose = FALSE),
                       silent = TRUE)
        }else if(brule == 1){
            ## standard bref rule + pa buffer
            if(!is.numeric(bindi) || is.na(bindi) || !is.numeric(findi) || is.na(findi) ||
               !is.numeric(bpbref) || is.na(bpbref)){
                tacs <- func2(inp, tacs=tacs, pars=pars)
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
                                       evalBreakpointB = ',evalBreakpointB,',
                                       cpvec = ',cpvec,',
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
                tacs <- func2(inp, tacs=tacs, pars=pars)
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
                                           evalBreakpointB = ',evalBreakpointB,',
                                           cpvec = ',cpvec,',
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

        if(any(inherits(tac, "try-error")) || any(!is.numeric(tac)) || any(is.na(tac))){
            tacs <- func2(inp, tacs=tacs, pars=pars)
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

        if(!',cpvec,' && assessInt > 1){
             tac <- rep(tac / assessInt, assessInt)
        }
        if(',halfAssessInt,'){
            tac <- rep(tac, assessInt)
        }

        names(tac) <- paste0("TAC", seq(length(tac)))
        tactmp <- data.frame(t(tac), id="',id,'", hitSC=hitSC,
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
    }'))

    ## create HCR as functions
    ## templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    ## assign(value=templati, x=id, envir=env)

    templati <- eval(parse(text = eval(template)))
    class(templati) <- c(class(templati), "hcr")
    attributes(templati)$id <- id
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}




#' def.hcr.sam
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

    template  <- expression(paste0('function(obs, tacs = NULL, pars=NULL){
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
            tacs <- iamse:::gettacs(tacs.=tacs, id.="',id,'", TAC.=tac, obs.=obs)
            tacs$conv[nrow(tacs)] <- TRUE
            return(tacs)
        }
    }
}'))

    ## create HCR as functions
    ## templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    ## assign(value=templati, x=id, envir=env)

    templati <- eval(parse(text = eval(template)))
    class(templati) <- c(class(templati), "hcr")
    attributes(templati)$id <- id
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}



#' def.hcr.pseudo
#'
#' @title Define harvest control rule with pseudo assessment
#'
#' @param id Name/ID of HCR. Default: "pseudo-msy"
#' @param fractiles Fractiles. List
#' @param breakpointB breakpointb
#' @param clType Catch type for uncertainty cap. Default: "TAC".
#' @param clyears Number of years for uncertainty cap. Default: 1.
#' @param stab Uncertainty cap. Default: FALSE.
#' @param lower Upper bound of uncertainty cap. Default: 0.8.
#' @param upper Upper bound of uncertainty cap. Default: 1.2.
#' @param env Environment. Default: globalenv()
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
        'function(obs, tacs = NULL, pars=NULL){


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
                cl <- tacs[nrow(tacs),paste0("TAC",assessInt)]
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

        tacs <- iamse:::gettacs(tacs.=tacs, id.="',id,'", TAC.=tac, obs.=obs)
        tacs$hitSC[nrow(tacs)] <- hitSC
        return(tacs)
    }'))

    ## create HCR as functions
    ## templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
    ## assign(value=templati, x=id, envir=env)

    templati <- eval(parse(text = eval(template)))
    class(templati) <- c(class(templati), "hcr")
    attributes(templati)$id <- id
    assign(value=templati, x=id, envir=env)

    ## allow for assigning names
    invisible(id)
}


#' def.hcr.ss3
#' @title Define harvest control rule for SS3
#'
#' @param id Name/ID of HCR. Default: "ss3"
#' @param nonconvHCR HCR if SAM does not converge. Default: "conscat" (constant catch).
#' @param silent silent
#' @param verbose verbose
#' @param env Environment. Default: globalenv()
#'
#'
#' @export
#'
def.hcr.ss3 <- function(id = "ss3",
                        nonconvHCR = "conscat",
                        silent = TRUE,
                        verbose = FALSE,
                        env = globalenv()
){
  
  
  template  <- expression(paste0('function(obs, tacs = NULL, pars=NULL){
    silent <- ',silent,'
    verbose <- ',verbose,'
    func <- get("',nonconvHCR,'")
    
    browser()

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
            tacs <- iamse:::gettacs(tacs.=tacs, id.="',id,'", TAC.=tac, obs.=obs)
            tacs$conv[nrow(tacs)] <- TRUE
            return(tacs)
        }
    }
}'))
  
  ## create HCR as functions
  ## templati <- eval(parse(text=paste(parse(text = eval(template)),collapse=" ")))
  ## assign(value=templati, x=id, envir=env)
  
  templati <- eval(parse(text = eval(template)))
  class(templati) <- c(class(templati), "hcr")
  attributes(templati)$id <- id
  assign(value=templati, x=id, envir=env)
  
  ## allow for assigning names
  invisible(id)
}