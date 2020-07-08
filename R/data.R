#' @name checkDat
#' @description Checks data list and fills missing slots
#' @param dat - data list with all parameters
#' @return Updated data list
#' @export
checkDat <- function(dat){

    ## ages
    ##------------------
    amax <- (dat$amax+1)
    ages <- 1:amax
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
        dat$M <- rep(dat$M, amax)
    }

    ## historic fishing mortality
    ##------------------
    if(!"Fvals" %in% names(dat)){
        dat$Fvals <- seq(0.01, 1.5, length.out = dat$ny)
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
        dat$initN <- rep(0, amax)
    }

    ## Fishing
    ##------------------
    if(!"initF" %in% names(dat)){
        dat$initF <- 0.2
    }

    ## spwaning potential ratio
    ##------------------
    if(!"SSBPR0" %in% names(dat)){
        dat$SSBPR0 <- getSSBPR0(M = dat$M, mat = dat$mat, fecun = dat$fecun, amax = amax)
    }

    ## recruitment
    ##------------------
    if(!"SR" %in% names(dat)){
        dat$SR <- "bevholt"
    }
    if(!"h" %in% names(dat)){
        dat$h <- 0.71              ## median over all species Myers 1999
    }

    ## return
    return(dat)
}
