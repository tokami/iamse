#' @title check.dat
#' @description Checks data list and fills missing slots
#' @param dat - data list with all parameters
#' @param verbose - print informative messages
#' @return Updated data list
#' @export
check.dat <- function(dat = NULL, verbose = TRUE){

    if(is.null(dat)) dat <- list()

    ## Do not change processes (sel, mat, weight, weightF) when re-running check.dat()
    if(is.null(dat$fixProcs)) dat$fixProcs <- FALSE


    ## Number of historic years
    ##------------------
    if(!any(names(dat) == "ny")){
        dat$ny <- 35
    }
    ny <- dat$ny

    ## Number of seasons
    ##------------------
    if(!any(names(dat) == "ns")){
        dat$ns <- 1
    }
    ns <- dat$ns
    nt <- ny * ns

    ## ages
    ##------------------
    if(!any(names(dat) == "amax")) dat$amax <- 10 ## stop("Maximum age missing: amax")
    if(dat$amax%%1 != 0){
        if(verbose) writeLines(paste0("Maximum age (dat$amax = ",dat$amax,") is not natural number. This is not possible, setting amax equal to ",ceiling(dat$amax)))
        dat$amax <- ceiling(dat$amax)
    }

    amax <- dat$amax + 1
    asmax <- amax * ns
    dat$asmax <- asmax
    ages <- 0:(amax-1)
    ## OLD: ages <- t(t(matrix(rep(ages,ns),ncol=ns,nrow=amax)) + seq(0, 1-1/ns, 1/ns))
    ds05 <- 1/ns/2
    ages <- t(t(matrix(rep(ages,each = ns),ncol=ns,nrow=amax, byrow=TRUE)) + seq(ds05, 1-ds05, length.out=ns))
    dat$ages <- ages

    dat$yvec <- rep(1:ny, each = ns)
    dat$svec <- rep(1:ns, ny)
    dat$s1vec <- seq(1, nt, ns)
    dat$s1avec <- seq(1, asmax, ns)
    dat$as2a <- rep(1:amax, each = ns)
    dat$as2s <- rep(1:ns, amax)


    ## Parameterisation on life history parameters
    ##------------------
    if(!any(names(dat) == "Linf")){
        dat$Linf <- NULL
    }else Linf <- dat$Linf
    if(!any(names(dat) == "K")){
        dat$K <- NULL
    }else K <- dat$K
    if(!any(names(dat) == "t0")){
        dat$t0 <- NULL
    }else t0 <- dat$t0
    ## Parameters of the seasonalised VBGF
    if(!any(names(dat) == "C")){
        dat$C <- NULL
    }
    if(!any(names(dat) == "ts")){
        dat$ts <- NULL
    }
    if(is.null(dat$Linf) && any(!is.null(dat$K), !is.null(dat$t0))) dat$Linf <- 100 ## stop("dat$Linf not provided.")
    if(is.null(dat$K) && any(!is.null(dat$Linf), !is.null(dat$t0))) dat$K <- 0.1 ## stop("dat$K not provided.")
    if(is.null(dat$t0) && any(!is.null(dat$Linf), !is.null(dat$K))) dat$t0 <- -0.1 ## stop("dat$t0 not provided.")
    if(!any(is.null(dat$t0), is.null(dat$Linf), is.null(dat$K))){
        if(!is.null(dat$C) && !is.null(dat$ts) && dat$C > 0){
            LA <- dat$Linf * (1 - exp(-(dat$K * (ages - dat$t0) + (((dat$C * dat$K)/(2 * pi)) *
                                                                   sin(2 * pi * (ages - dat$ts))) -
                                        (((dat$C * dat$K)/(2 * pi)) * sin(2 * pi * (dat$t0 - dat$ts))))))
        }else{
            LA <- dat$Linf * (1 - exp(-dat$K * (ages - dat$t0)))
        }
    }else{
        LA <- NULL
    }
    dat$LA <- LA
    if(!any(names(dat) == "binwidth")){
        if(is.null(LA)){
            dat$binwidth <- NULL
        }else dat$binwidth <- 1
    }
    binwidth <- dat$binwidth
    if(!is.null(LA)){
        mids <- seq((binwidth/2), dat$Linf * 1.5, by = binwidth)
        dat$mids <- mids
        highs <- mids + binwidth/2
        lows <- mids - binwidth/2
        lbprobs <- function(mnl, sdl){
            return(pnorm(highs, mnl, sdl) - pnorm(lows, mnl, sdl))
        }
        vlprobs <- Vectorize(lbprobs, vectorize.args=c("mnl","sdl"))
        if(!any(names(dat) == "CVlen")){
            dat$CVlen <- 0.1
        }
        LA2 <- array(LA, dim = c(amax,1,ns))  ## NEW: LA2 <- array(t(LA), dim = c(asmax,1,ns)) ## this plba can translate age+season to length and back at every season, I think this might require changing all process matrices to smooth_age by season (not as currently: age by season)
        plba <- apply(apply(LA2, c(1,3), function(x) vlprobs(x, x * dat$CVlen)), c(1,3), t)
        for(i in 1:dim(plba)[3]){
            tmp <- rowSums(plba[,,i])
            plba[,,i] <- plba[,,i] / tmp
            plba[tmp == 0,,i] <- 0
        }
        ## plba <- t(vlprobs(LA, LA * dat$CVlen))
        ## plba <- plba / rowSums(plba)

        ## NEW:
        LA2 <- array(t(LA), dim = c(asmax,1,ns))
        ## this plba can translate age+season to length and back at every season, I think this might require changing all process matrices to smooth_age by season (not as currently: age by season)
        plba2 <- apply(apply(LA2, c(1,3), function(x) vlprobs(x, x * dat$CVlen)), c(1,3), t)
        for(i in 1:dim(plba)[3]){
            tmp <- rowSums(plba2[,,i])
            plba2[,,i] <- plba2[,,i] / tmp
            plba2[tmp == 0,,i] <- 0
        }

    }else{
        plba <- mids <- NULL
        plba2 <- NULL
    }
    dat$plba <- plba
    dat$mids <- mids
    dat$plba2 <- plba2


    ## weight at age
    ##------------------
    if(!any(names(dat) == "lwa")){
        lwa <- NULL
    }else lwa <- dat$lwa
    if(!any(names(dat) == "lwb")){
        lwb <- NULL
    }else lwb <- dat$lwb
    if(is.null(lwa) && !is.null(lwb)) stop("Length-weight coefficient missing: dat$lwa")
    if(is.null(lwb) && !is.null(lwa)) stop("Length-weight exponent missing: dat$lwb")
    if(is.null(LA) && !is.null(lwa) && verbose) writeLines("Parameters of the length-weight relationship not used because growth parameters for the age-length key not provided.")
    if(any(is.null(LA),is.null(lwb),is.null(lwa))){
        if(!any(names(dat) == "weight")) dat$weight <- matrix(NA, nrow = amax, ncol = ns)
    }else dat$weight <- lwa * LA ^ lwb



    ## weight at age (fishery)
    ##------------------
    if(!any(names(dat) == "lwaF")){
        lwaF <- NULL
    }else lwaF <- dat$lwaF
    if(!any(names(dat) == "lwbF")){
        lwbF <- NULL
    }else lwbF <- dat$lwbF
    if(is.null(lwaF) && !is.null(lwbF)) stop("Length-weight coefficient for catch weight-at-age missing: dat$lwaF")
    if(is.null(lwbF) && !is.null(lwaF)) stop("Length-weight exponent for catch weight-at-age missing: dat$lwbF")
    if(is.null(lwaF) && !is.null(lwa)){
        lwaF <- lwa
        if(verbose) writeLines("Length-weight coefficient for catch weight-at-age missing: dat$lwaF. Setting dat$lwaF equal to stock length-weight coefficient: dat$lwa")
    }
    if(is.null(lwbF) && !is.null(lwb)){
        lwbF <- lwb
        if(verbose) writeLines("Length-weight exponent for catch weight-at-age missing: dat$lwbF. Setting dat$lwbF equal to stock length-weight exponent: dat$lwb")
    }
    dat$lwaF <- lwaF
    dat$lwbF <- lwbF
    if(any(is.null(LA),is.null(lwbF),is.null(lwaF))){
        if(!any(names(dat) == "weightF")) dat$weightF <- matrix(NA, nrow = amax, ncol = ns)
    }else dat$weightF <- lwaF * LA ^ lwbF


    ## maturity
    ##------------------
    if(!any(names(dat) == "Lm50")){
        Lm50 <- NULL
    }else Lm50 <- dat$Lm50
    if(!any(names(dat) == "Lm95")){
        Lm95 <- NULL
    }else Lm95 <- dat$Lm95
    if(is.null(Lm50) && !is.null(Lm95)) stop("Length at 95% maturity missing: dat$Lm50")
    if(is.null(Lm95) && !is.null(Lm50)) stop("Length at 50% maturity missing: dat$Lm95")
    if(any(is.null(LA),is.null(Lm50),is.null(Lm95))){
        if(!any(names(dat) == "mat")) dat$mat <- matrix(NA, nrow = amax, ncol = ns)
    }else if(!dat$fixProcs) dat$mat <- getMat(Lm50, Lm95, mids, plba)


    ## selectivity
    ##------------------
    if(!any(names(dat) == "selType")) dat$selType <- "logistic"
    if(!any(names(dat) == "selAge")) dat$selAge <- FALSE
    if(!any(names(dat) == "LFS")) dat$LFS <- NULL
    if(!any(names(dat) == "sl")) dat$sl <- NULL
    if(!any(names(dat) == "sr")) dat$sr <- NULL
    if(dat$selType == "dnormal" && any(is.null(c(dat$LFS,dat$sl,dat$sr))))
        stop("Double-normal selectivity specified, but not all required parameters provided. Please provide 'LFS', 'sl', and 'sr'.")
    if(!any(names(dat) == "Ls50")){
        Ls50 <- NULL
    }else Ls50 <- dat$Ls50
    if(!any(names(dat) == "Ls95")){
        Ls95 <- NULL
    }else Ls95 <- dat$Ls95
    if(dat$selType == "logistic" && is.null(Ls50) && !is.null(Ls95))
        stop("Length at 50% selectivity missing: dat$Ls50")
    if(dat$selType == "logistic" && is.null(Ls95) && !is.null(Ls50))
        stop("Length at 95% selectivity missing: dat$Ls95")
    if(dat$selType == "logistic" && is.null(Ls50) && !is.null(Lm50)){
        Ls50 <- Lm50
        if(verbose) writeLines("Length at 50% selectivity missing: dat$Ls50. Setting dat$Ls50 equal to maturity parameter: dat$Lm50.")
    }
    if(dat$selType == "logistic" && is.null(Ls95) && !is.null(Lm95)){
        Ls95 <- Lm95
        if(verbose) writeLines("Length at 95% selectivity missing: dat$Ls95. Setting dat$Ls95 equal to maturity parameter: dat$Lm95.")
    }
    dat$Ls50 <- Ls50
    dat$Ls95 <- Ls95
    ## get selectivity
    if(dat$selType == "logistic" && any(is.null(LA),is.null(Ls50),is.null(Ls95))){
        if(!any(names(dat) == "sel")) dat$sel <- list(matrix(NA, nrow = amax, ncol = ns))
    }else if(dat$selType == "dnormal" &&
             any(is.null(LA),is.null(dat$LFS),is.null(dat$sl),is.null(dat$sr))){
        if(!any(names(dat) == "sel")) dat$sel <- list(matrix(NA, nrow = amax, ncol = ns))
    }else if(!dat$fixProcs) dat$sel <- getSel(dat, mids, plba, dat[["selType"]], dat[["selAge"]], ages)

    if(!inherits(dat$sel, "list") && dim(dat$sel) == c(amax,ns)){
        dat$sel <- list(dat$sel/max(dat$sel))
    }else if(inherits(dat$sel, "list") && length(dat$sel) != ny && length(dat$sel) != 1) stop("Selectivity at age ('dat$sel') has incorrect dimensions. dat$sel has to be a matrix with dimensions: maximum age + 1 (age 0) x the number of seasons (rows x columns)!")


    ## fecundity
    ##------------------
    if(!any(names(dat) == "fecun")) dat$fecun <- 1


    ## natural mortality over time
    ##------------------
    if(is.null(LA)){
        if(!any(names(dat) == "M")){
            dat$M <- matrix(NA, nrow = ny, ncol = ns)
        }else if(!inherits(dat$M, "matrix") && length(dat$M) == 1){
            dat$M <- matrix(dat$M/ns, nrow=ny, ncol=ns)
        }else if(!inherits(dat$M, "matrix") && length(dat$M) == ny){
            dat$M <- matrix(rep(dat$M/ns,each=ns), nrow=ny, ncol=ns)
        }else if(any(dim(dat$M) != c(ny,ns))){
            writeLines("dat$M does not have the correct dimensions. dat$M has to be a matrix with dimensions: number of years x the number of seasons (rows x columns). Extending M matrix.")
            dat$M <- rbind(dat$M, matrix(rep(dat$M[dim(dat$M)[1]],ny-dim(dat$M)[1]),
                                         nrow=ny-dim(dat$M)[1],ncol=ns))
            }
    }else{
        if(!any(names(dat) == "M")){
            if(verbose) writeLines("No natural mortality provided. Setting time-invariant M corresponding to annual M based Gislason's formula and growth parameters.")
            dat$M <- matrix(getM(Linf, K, mids) / ns, nrow = ny, ncol = ns)
        }else if(!inherits(dat$M, "matrix") && length(dat$M) == 1 && !dat$fixProcs){
            dat$M <- matrix(rep(dat$M/ns,ny), nrow = ny, ncol = ns, byrow=TRUE)
        }else if(nrow(dat$M) == 1 && !dat$fixProcs){
            dat$M <- matrix(rep(dat$M/ns,ny), nrow = ny, ncol = ns, byrow=TRUE)
        }else if(any(dim(dat$M) != c(ny,ns))){
            stop(paste0("Intra-annual natural mortality ('dat$M') has dimensions ",
                        paste0(dim(dat$M),collapse=" x "),
                        ". It should either be a single numeric or correspond to a matrix with dimensions: dat$ny x dat$ns."))
        }
    }


    ## natural mortality at length
    ##------------------
    if(is.null(LA)){
        if(!any(names(dat) == "Msel")) dat$Msel <- list(matrix(NA, nrow = amax, ncol = ns))
    }else{
        if(!any(names(dat) == "Msel")){
            if(verbose) writeLines("No natural mortality at age provided. Setting M-at-age based on the Gislason's (2010) empirical formula.")
            dat$Msel <- getMsel(Linf, K, mids, plba, below10 = dat$below10)
        }else if(!inherits(dat$Msel, "list") && dim(dat$Msel) == c(amax,ns) && !dat$fixProcs){
            dat$Msel <- list(dat$Msel/max(dat$Msel))
        }else if(inherits(dat$Msel, "list") && length(dat$Msel) != ny && length(dat$Msel) != 1) stop("Natural mortality at age (dat$Msel) has incorrect dimensions. dat$Msel has to be a matrix with dimensions: maximum age + 1 (age 0) x the number of seasons (rows x columns)!")
    }


    ## Historic fishing mortality
    ##------------------
    ## timeseries <- c(seq(1985,2015,5),2019)
    ## otter <- c(2400, 2500, 1900, 1600, 1000, 700, 500, 500)
    ## beam <- c(200, 1200, 1500, 1300, 1200, 1300, 900, 900)
    ## eff <- beam  + otter
    timeseries <- c(seq(1970,2015,5), 2019)
    eff <- c(0.1, 0.2, 0.35, 0.55, 0.7, 0.8, 0.9, 0.97, 1.0, 1.0, 0.9)
    effrel <- eff/max(eff)
    mod <- stats::smooth.spline(x=timeseries, y=effrel)

    ## Intra-annual (seasonl) FM (with sum = 1)
    if(!any(names(dat) == "iaFM")){
        dat$iaFM <- rep(1/ns, ns)
    }else if(length(dat$iaFM) != ns){
        if(verbose) writeLines(paste0(
                        "The intra-annual FM (dat$iaFM) does not have the correct length (dat$ns). ",
                        "Please revise dat$iaFM. Removing the intra-annual FM pattern for now."
                    ))
        dat$iaFM <- rep(1/ns, ns)
    }else if(all(dat$iaFM == 0)){
        dat$iaFM <- rep(1/ns, ns)
    }else if(sum(dat$iaFM) > 1 + 1e-8 || sum(dat$iaFM) < 1 - 1e-8){
        if(verbose) writeLines(paste0(
                        "The intra-annual FM (dat$iaFM) does not sum up to 1. Please revise dat$iaFM. ",
                        "Overwriting dat$iaFM with the correct values that sum up to 1."
                    ))
        dat$iaFM <- dat$iaFM / sum(dat$iaFM)
    }

    ## Combined FM pattern (historic years vs seasons)
    if(!any(names(dat) == "FM")){
        tmp <- stats::predict(mod, x = seq(1970, 2019, length.out = ny))$y
        dat$FM <- as.matrix(tmp) %*% t(as.matrix(dat$iaFM))
    }else if((!inherits(dat$FM, "matrix") && length(dat$FM) == ny) ||
             (inherits(dat$FM, "matrix") && any(dim(as.matrix(dat$FM)) == 1))){
        if(verbose) writeLines(paste0(
                        "Assuming that the provided fishing mortality (dat$FM) corresponds to the ",
                        "annual fishing mortality and, thus, extending it to the FM matrix by year ",
                        " and season using the intra-annual FM pattern (dat$iaFM)."
                    ))
        if(!inherits(dat$FM, "matrix")){
            dat$FM <- as.matrix(dat$FM) %*% t(as.matrix(dat$iaFM))
        }else if(ncol(dat$FM) == 1){
            dat$FM <- dat$FM %*% t(as.matrix(dat$iaFM))
        }else if(ncol(dat$FM) == 1){
            dat$FM <- t(dat$FM) %*% t(as.matrix(dat$iaFM))
        }
    }else if((!inherits(dat$FM, "matrix") && length(dat$FM) != ny) ||
             (inherits(dat$FM, "matrix") && ncol(dat$FM) != ns && nrow(dat$FM) != ny)){
        if(verbose) writeLines(paste0(
                        "The fishing mortality (dat$FM) does not have the correct length/dimensions. ",
                        "dat$FM has to be a matrix with dimensions: number of years x ",
                        "the number of seasons (rows x columns). Overwriting dat$FM with default values."
                    ))
        tmp <- predict(mod, x = seq(1970, 2019, length.out = ny))$y
        dat$FM <- as.matrix(tmp) %*% t(as.matrix(dat$iaFM))
    }


    ## Depletion level final year
    ##------------------
    if(!"depl" %in% names(dat)){
        dat$depl <- 0.5
    }

    ## Depletion relative to
    ##------------------
    if(!"depl.quant" %in% names(dat)){
        dat$depl.quant <- "B0"
    }

    ## Depletion level with probability
    ##------------------
    if(!"depl.quant" %in% names(dat)){
        dat$depl.prob <- 0.5
    }

    ## initial pop size
    ##------------------
    if(!"initN" %in% names(dat)){
        dat$initN <- rep(0, amax)
    }

    ## Fishing
    ##------------------
    if(!"initF" %in% names(dat)){
        dat$initF <- 0.2
    }

    ## Recruitment
    ##------------------
    ## Stock-recruitment model
    if(!"SR" %in% names(dat)){
        dat$SR <- "bevholt"
    }
    ## Steepness
    if(!"h" %in% names(dat)){
        dat$h <- 0.74              ## mean of Thorson 2018
    }
    ## Number of recruits
    if(!"R0" %in% names(dat)){
        dat$R0 <- 1e6
    }
    ## breakpoint (bp)
    if(!"bp" %in% names(dat)){
        dat$bp <- NA
    }
    if(!"recBeta" %in% names(dat)){
        dat$recBeta <- NA
    }
    ## measure of the radius of curvature near the breakpoint (as gamma -> 0 = sharp bent like hockey-stick)
    if(!"recGamma" %in% names(dat)){
        dat$recGamma <- NA
    }
    ## virgin biomass (vb)
    if(!"vb" %in% names(dat)){
        dat$vb <- NA
    }


    ## Spawning
    ##------------------
    if(is.null(dat$spawning)){
        if(ns == 1){
            dat$spawning <- 1
        }else{
            dat$spawning <- c(1, rep(0, ns-1))
        }
    }else if(length(dat$spawning) != ns){
        writeLines("The length of 'dat$spawning' is not equal to number of seasons (dat$ns). Assuming spawning in first season.")
        if(ns == 1){
            dat$spawning <- 1
        }else{
            dat$spawning <- c(1, rep(0, ns-1))
        }
    }
    dat$spawning <- dat$spawning / sum(dat$spawning)
    if(!"age0" %in% names(dat)){
        dat$age0 <- 0
    }
    ##    dat$indage0 <- which(ages == ages[which.min(abs(ages - 1/ns/2 - dat$age0))], arr.ind = TRUE)
    dat$indage0 <- which.min(abs(t(ages) - 1/ns/2 - dat$age0))

    ## pzbm
    ##------------------
    if(!"pzbm" %in% names(dat)){
        dat$pzbm <- 0
    }


    ## Survey information
    ##------------------
    ## North Sea - IBITS: the majority of countries have only carried
    ## out a survey twice a year; a first quarter survey (January-February) and a
    ## third quarter survey (August-September)
    if(is.null(dat$surveyTimes)) dat$surveyTimes <- c(1/12, 7/12)
    nsurv <- length(dat$surveyTimes)

    if(is.null(dat$surveyBeforeAssessment)) dat$surveyBeforeAssessment <- c(FALSE, FALSE)
    if(length(dat$surveyBeforeAssessment) != nsurv){
        writeLines(paste0("Length of dat$surveyBeforeAssessment does not match with dat$surveyTimes! Overwriting with FALSE."))
        dat$surveyBeforeAssessment <- rep(FALSE, nsurv)
    }

    ## survey selectivities and timings
    if(is.null(dat$selI)){
        dat$selI <- list()
        for(i in 1:nsurv){
            dat$selI[[i]] <- dat$sel[[1]]
        }
    }

    ## catchability
    if(!any(names(dat) == "q")) dat$q <- 0.005

    ## effort coefficient
    if(!any(names(dat) == "qE")) dat$qE <- 5e2

    ## years of data available to assessment
    if(is.null(dat$nyC[1])) dat$nyC <- ny
    if(is.null(dat$nyI[1])) dat$nyI <- rep(ny, nsurv)
    if(is.null(dat$nyE[1])) dat$nyE <- NA
    ## generate Length frequency distributions (LFD)?
    if(is.null(dat$obsLFD[1])) dat$obsLFD <- FALSE

    ## most recent years of data available to assessment ## TODO: RENAME:
    ## 0: last year available to assessment
    ## 1: last year not available
    ## 2: last 2 years not available
    ## etc.
    if(is.null(dat$lastyC[1])) dat$lastyC <- 0
    if(is.null(dat$lastyI[1])) dat$lastyI <- rep(0, nsurv)
    if(is.null(dat$lastyE[1])) dat$lastyE <- 0


    ## Seasonal catch observations
    if(is.null(dat$catchSeasons)) dat$catchSeasons <- 1

    ## Seasonal effort observations
    if(is.null(dat$effortSeasons)) dat$effortSeasons <- 1

    ## Timing of catch observation interval (default: calendar year => 1)
    if(is.null(dat$catchObsTiming)) dat$catchObsTiming <- 1

    ## return
    return(dat)
}
