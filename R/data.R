#' @name checkDat
#' @description Checks data list and fills missing slots
#' @param dat - data list with all parameters
#' @param verbose - print informative messages
#' @return Updated data list
#' @export
checkDat <- function(dat, verbose = TRUE){

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
    if(!any(names(dat) == "amax")) stop("Maximum age missing: amax")
    amax <- (dat$amax)
    asmax <- (amax+1) * ns
    ages <- 0:amax
    ## OLD: ages <- t(t(matrix(rep(ages,ns),ncol=ns,nrow=amax+1)) + seq(0, 1-1/ns, 1/ns))
    ds05 <- 1/ns/2
    ages <- t(t(matrix(rep(ages,each = ns),ncol=ns,nrow=amax+1, byrow=TRUE)) + seq(ds05, 1-ds05, length.out=ns))
    dat$ages <- ages

    dat$yvec <- rep(1:ny, each = ns)
    dat$svec <- rep(1:ns, ny)
    dat$savec <- rep(1:ns, amax+1)
    dat$asvec <- rep(1:(amax+1), each = ns)
    dat$s1vec <- seq(1, nt, ns)
    dat$s1avec <- seq(1, asmax, ns)

    ## growth at age
    ##------------------
    if(!any(names(dat) == "a0")){
        dat$a0 <- 0
    }
    if(!any(names(dat) == "linf")) stop("Asymptotic length missing: linf")
    if(!any(names(dat) == "k")) stop("Growth coefficient missing: k")
    LA <- dat$linf * (1 - exp(-dat$k * (ages - dat$a0)))
    dat$LA <- LA
    if(!any(names(dat) == "binwidth")){
        dat$binwidth <- 1
    }
    binwidth <- dat$binwidth
    mids <- seq((binwidth/2), dat$linf * 1.2, by = binwidth)
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
    LA <- array(LA, dim = c(amax+1,1,ns))
    plba <- apply(apply(LA, c(1,3), function(x) vlprobs(x, x * dat$CVlen)), c(1,3), t)
    for(i in 1:dim(plba)[3]){
        tmp <- rowSums(plba[,,i])
        plba[,,i] <- plba[,,i] / tmp
        plba[tmp == 0,,i] <- 0
    }
    ## plba <- t(vlprobs(LA, LA * dat$CVlen))
    ## plba <- plba / rowSums(plba)
    dat$plba <- plba


    ## weight at age
    ##------------------
    if(!any(names(dat) == "lwa")) stop("Length-weight coefficient missing: lwa")
    if(!any(names(dat) == "lwb")) stop("Length-weight exponent missing: lwb")
    dat$weights <- dat$lwa * dat$LA ^ dat$lwb
    dat$weight <- rowMeans(dat$weights)


    ## weight at age (fishery)
    ##------------------
    if(!any(names(dat) == "lwaF")){
        dat$lwaF <- dat$lwa
    }
    if(!any(names(dat) == "lwbF")){
        dat$lwbF <- dat$lwb
    }
    dat$weightFs <- dat$lwaF * dat$LA ^ dat$lwbF
    dat$weightF <- rowMeans(dat$weightFs)



    ## maturity
    ##------------------
    if(!any(names(dat) == "Lm50")) stop("Length at 50% maturity missing: Lm50")
    if(!any(names(dat) == "Lm95")) stop("Length at 95% maturity missing: Lm95")
    if(!is.na(dat$Lm50) && !is.na(dat$Lm95)){
        dat$mats <- getMat(dat$Lm50, dat$Lm95, dat$mids, dat$plba)
        ##        dat$mat <- getMat(dat$Lm50, dat$Lm95, dat$mids, dat$plba[,,1, drop=FALSE]) ## or rowMeans?
        dat$mat <- rowMeans(dat$mats)
    }


    ## selectivity
    ##------------------
    if(!any(names(dat) == "L50")) dat$L50 <- dat$Lm50
    if(!any(names(dat) == "L95")) dat$L95 <- dat$Lm95
    if(!is.na(dat$L50) && !is.na(dat$L95)){
        dat$sels <- getSel(dat$L50, dat$L95, dat$mids, dat$plba)
##        dat$sel <- getSel(dat$L50, dat$L95, dat$mids, dat$plba[,,1, drop=FALSE]) ## or rowMeans?
        dat$sel <- rowMeans(dat$sels)
    }


    ## fecundity
    ##------------------
    if(!any(names(dat) == "fecun")) dat$fecun <- 1


    ## natural mortality over time
    ##------------------
    ## per year
    if(!any(names(dat) == "M")){
        dat$M <- rep(getM(dat$linf, dat$k, dat$mids), dat$ny)
        if(verbose) writeLines("No natural mortality (M) provided. Setting time-invariant M.")
    }else if(length(dat$M) == 1){
        dat$M <- rep(dat$M, dat$ny)
    }else if(length(dat$M) < dat$ny){
        stop(paste0("Natural mortality ('dat$M') has length ",length(dat$M),". It should either be a single numeric or correspond at least to the number of historical years ('dat$ny')."))
    }
    ## per season
    if(!any(names(dat) == "Ms")){
        dat$Ms <- rep(dat$M / ns, each = ns)
        if(verbose) writeLines("No in-year natural mortality (Ms) provided. Setting time-invariant Ms corresponding to annual M (adjusted by in-year time steps).")
    }else if(length(dat$Ms) == 1){
        dat$Ms <- rep(dat$Ms, dat$ny * ns)
    }else if(length(dat$Ms) < dat$ny * ns){
        stop(paste0("In-year natural mortality ('dat$Ms') has length ",length(dat$Ms),". It should either be a single numeric or correspond at least to the number of historical time steps ('dat$ny' * 'dat$ns')."))
    }


    ## natural mortality at length
    ##------------------
    if(!any(names(dat) == "Msel")){
        if(verbose) writeLines("No natural mortality at age provided. Setting M-at-age based on the Gislason's (2010) empirical formula.")
        dat$Msels <- getMsel(dat$linf, dat$k, dat$mids, dat$plba)
        dat$Msel <- lapply(dat$Msels, rowMeans)
    }else if(!inherits(dat$Msel, "list") && length(dat$Msel) == amax+1){
        dat$Msel <- list(dat$Msel/max(dat$Msel))
        dims <- dim(dat$plba)
        dat$Msels <- list(matrix(dat$Msel, ncol=dims[3], nrow=dims[1]))
    }else if(inherits(dat$Msel, "list") && (length(dat$Msel) == dat$ny || length(dat$Msel) == 1)){
    }else stop("Natural mortality at age ('dat$Msel') has incorrect length. Length has to be equal to maximum age + 1 (age 0)!")
    ## dat$Mtot <- matrix(NA, nrow=amax+1, ncol=dat$ny)
    ## for(y in 1:dat$ny){
    ##     dat$Mtot[,y] <- rowSums(dat$Ms[y] * dat$Msels)
    ## }


    ## historic fishing mortality
    ##------------------
    ## timeseries <- c(seq(1985,2015,5),2019)
    ## otter <- c(2400, 2500, 1900, 1600, 1000, 700, 500, 500)
    ## beam <- c(200, 1200, 1500, 1300, 1200, 1300, 900, 900)
    ## eff <- beam  + otter
    timeseries <- c(seq(1970,2015,5),2019)
    eff <- c(0.1, 0.2, 0.35, 0.55, 0.7, 0.8, 0.9, 0.97, 1.0, 1.0, 0.9)
    effrel <- eff/max(eff)
    mod <- smooth.spline(x=timeseries, y=effrel)
    if(!"FM" %in% names(dat)){
        dat$FM <- predict(mod, x = seq(1970, 2019, length.out = dat$ny))$y
    }else if(length(dat$FM) < dat$ny){
        warning("Length of FM not equal to ny. Overwriting FM.")
        dat$FM <- predict(mod, x = seq(1970, 2019, length.out = dat$ny))$y
    }
    ## account for seasons
    dat$Fs <- rep(dat$FM / ns, each = ns)

    ## Depletion level final year
    ##------------------
    if(!"depl" %in% names(dat)){
        dat$depl <- 0.5
    }

    ## Depletion relative to:
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
        dat$initN <- rep(0, amax + 1)
    }

    ## Fishing
    ##------------------
    if(!"initF" %in% names(dat)){
        dat$initF <- 0.2
    }

    ## recruitment
    ##------------------
    if(!"SR" %in% names(dat)){
        dat$SR <- "bevholt"
    }
    if(!"h" %in% names(dat)){
        dat$h <- 0.74              ## mean of Thorson 2018
    }
    if(!"R0" %in% names(dat)){
        dat$R0 <- 1e6
    }
    if(!"bp" %in% names(dat)){
        dat$bp <- 0
    }
    if(!"recBeta" %in% names(dat)){
        dat$recBeta <- 0
    }
    if(!"recGamma" %in% names(dat)){
        dat$recGamma <- 0  ## measure of the radius of curvature near the breakpoint (as gamma -> 0 = sharp bent like hockey-stick)
    }
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


    ## pzbm
    ##------------------
    if(!"pzbm" %in% names(dat)){
        dat$pzbm <- 0
    }


    ## return
    return(dat)
}
