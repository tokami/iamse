#' @name checkDat
#' @description Checks data list and fills missing slots
#' @param dat - data list with all parameters
#' @return Updated data list
#' @export
checkDat <- function(dat){

    ## Number of historic years
    ##------------------
    if(!any(names(dat) == "ny")){
        dat$ny <- 35
    }

    ## Number of seasons
    ##------------------
    if(!any(names(dat) == "nseasons")){
        dat$nseasons <- 1
    }

    ## ages
    ##------------------
    amax <- (dat$amax)
    ages <- 0:amax
    dat$ages <- ages

    ## growth at age
    ##------------------
    if(!"sel" %in% names(dat) | !"mat" %in% names(dat)){
        if("a0" %in% names(dat)){
            a0 <- dat$a0
        }else{
            a0 <- 0
        }
        LA <- dat$linf * (1 - exp(-dat$k * (ages - a0)))
        dat$LA <- LA
        binwidth <- dat$binwidth
        mids <- seq((binwidth/2), dat$linf * 1.2, by = binwidth)
        dat$mids <- mids
        highs <- mids + binwidth/2
        lows <- mids - binwidth/2
        lbprobs <- function(mnl, sdl){
            return(pnorm(highs, mnl, sdl) - pnorm(lows, mnl, sdl))
        }
        vlprobs <- Vectorize(lbprobs, vectorize.args=c("mnl","sdl"))
        plba <- t(vlprobs(LA, LA * dat$CVlen))
        plba <- plba / rowSums(plba)
        dat$plba <- plba
    }

    ## weight at age
    ##------------------
    if(!"weight" %in% names(dat)){
        dat$weight <- dat$lwa * dat$LA ^ dat$lwb
    }

    ## weight at age (fishery)
    ##------------------
    if(!"lwaF" %in% names(dat)){
        dat$lwaF <- dat$lwa
    }
    if(!"lwbF" %in% names(dat)){
        dat$lwbF <- dat$lwb
    }
    if(!"weightF" %in% names(dat)){
        dat$weightF <- dat$lwaF * dat$LA ^ dat$lwbF
    }

    ## selectivity
    ##------------------
    if(!"sel" %in% names(dat)){
        dat$sel <- getSel(dat$L50, dat$L95, dat$mids, dat$plba)
    }

    ## maturity
    ##------------------
    if(!"mat" %in% names(dat)){
        dat$mat <- getMat(dat$Lm50, dat$Lm95, dat$mids, dat$plba)
    }

    ## fecundity
    ##------------------
    if(!"fecun" %in% names(dat)){
        dat$fecun <- 1
    }

    ## natural mortality
    ##------------------
    if(length(dat$M) == 1){
        dat$M <- rep(dat$M, amax + 1)
    }

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
    if(!"Fvals" %in% names(dat)){
        dat$Fvals <- predict(mod, x = seq(1970, 2019, length.out = dat$ny))$y
    }else if(length(dat$Fvals) != dat$ny){
        warning("Length of dat$Fvals not equal to dat$ny. Overwriting dat$Fvals.")
        dat$Fvals <- predict(mod, x = seq(1970, 2019, length.out = dat$ny))$y
    }

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

    ## spwaning potential ratio
    ##------------------
    if(!"SSBPR0" %in% names(dat)){
        dat$SSBPR0 <- getSSBPR0(M = dat$M, mat = dat$mat, fecun = dat$fecun, amax = amax + 1)
    }

    ## recruitment
    ##------------------
    if(!"SR" %in% names(dat)){
        dat$SR <- "bevholt"
    }
    if(!"h" %in% names(dat)){
        dat$h <- 0.71              ## median over all species Myers 1999
    }
    if(!"R0" %in% names(dat)){
        dat$R0 <- 1e6
    }

    ## pzbm
    ##------------------
    if(!"pzbm" %in% names(dat)){
        dat$pzbm <- 0
    }


    ## return
    return(dat)
}
