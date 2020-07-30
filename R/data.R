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
    ns <- dat$nseasons

    ## ages
    ##------------------
    if(!any(names(dat) == "amax")) stop("Maximum age missing: amax")
    amax <- (dat$amax)
    ages <- 0:amax
    ages <- t(t(matrix(rep(ages,ns),ncol=ns,nrow=amax+1)) + seq(0, 1-1/ns, 1/ns))
    dat$ages <- ages


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
    dat$weight <- dat$lwa * dat$LA[,1] ^ dat$lwb ## or rowMeans

    ## weight at age (fishery)
    ##------------------
    if(!any(names(dat) == "lwaF")){
        dat$lwaF <- dat$lwa
    }
    if(!any(names(dat) == "lwbF")){
        dat$lwbF <- dat$lwb
    }
    dat$weightFs <- dat$lwaF * dat$LA ^ dat$lwbF
    dat$weightF <- dat$lwaF * dat$LA[,1] ^ dat$lwbF ## or rowMeans



    ## maturity
    ##------------------
    if(!any(names(dat) == "Lm50")) stop("Length at 50% maturity missing: Lm50")
    if(!any(names(dat) == "Lm95")) stop("Length at 95% maturity missing: Lm95")
    dat$mats <- getMat(dat$Lm50, dat$Lm95, dat$mids, dat$plba)
    dat$mat <- getMat(dat$Lm50, dat$Lm95, dat$mids, dat$plba[,,1, drop=FALSE]) ## or rowMeans?


    ## selectivity
    ##------------------
    if(!any(names(dat) == "L50")) dat$L50 <- dat$Lm50
    if(!any(names(dat) == "L95")) dat$L95 <- dat$Lm95
    dat$sels <- getSel(dat$L50, dat$L95, dat$mids, dat$plba)
    dat$sel <- getSel(dat$L50, dat$L95, dat$mids, dat$plba[,,1, drop=FALSE]) ## or rowMeans?


    ## fecundity
    ##------------------
    if(!any(names(dat) == "fecun")) dat$fecun <- 1


    ## natural mortality
    ##------------------
    if(!any(names(dat) == "M") || dim(dat$Ms)[2] != tail(dim(dat$plba),1)){
        writeLines("M not defined. Using Gislason's empirical formula for M at length.")
        M <- exp(0.55 - 1.61 * log(dat$LA) + 1.44 * log(dat$linf) + log(dat$k))
        for(i in 1:dim(M)[2]){
            M[dat$LA[,i] < 10,i] <- exp(0.55 - 1.61 * log(10) + 1.44 * log(dat$linf) + log(dat$k))
        }
        ## account for seasons
        dat$Ms <- M / ns
        dat$M <- M[,1]
    }else if(length(dat$M) == 1){
        dat$M <- rep(dat$M, amax + 1)
        dat$Ms <- dat$M / ns
    }else if(length(dat$M) != amax+1) stop("Natural mortality has incorrect length. Length has to be equal to maximum age + 1 (age 0)!")



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
    }else if(length(dat$FM) != dat$ny){
        warning("Length of FM not equal to ny. Overwriting FM.")
        dat$FM <- predict(mod, x = seq(1970, 2019, length.out = dat$ny))$y
    }
    ## account for seasons
    dat$Fs <- dat$FM / ns

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

    ## pzbm
    ##------------------
    if(!"pzbm" %in% names(dat)){
        dat$pzbm <- 0
    }


    ## return
    return(dat)
}
