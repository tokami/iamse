
#' @name genDevs
#' @export
genDevs <- function(n, sd, rho=0){

    rnoise <- rnorm(n, 0, sd) - sd^2/2

    res <- numeric(n)
    res[1] <- rnoise[1]
    if(n > 1){
        for(i in 2:n) res[i] <- rho * res[i-1] + (1 - rho) * rnoise[i] ## sqrt(1-rho^2) * rnoise[i]
    }

    res <- exp(res)
    res <- res/mean(res)

    return(res)

}

#' @name estDepl
#' @export
estDepl <- function(dat, fmax = 10, verbose = TRUE){

    if(!any(names(dat) == "ref")) stop("Reference points are missing in dat. Use estRef to estimate reference points.")

    depl <- dat$depl
    depl.quant <- dat$depl.quant

    frel <- dat$FM/max(dat$FM)

    fn <- function(logfabs, frel, depl, dat, opt=1){
        datx <- dat
        fpat <- frel * exp(logfabs)
        datx$FM <- fpat
        datx$Fs <- fpat / datx$nseasons
        dreal <- initPop(datx, out.opt = 2, depl.quant = depl.quant)
        if(opt==1) return((depl - dreal)^2)
        if(opt==2) return(dreal)
    }

    opt <- optimize(fn, c(log(0.0001),log(fmax)), frel = frel, depl = depl, dat = dat)
    fabs <- exp(opt$minimum)
    fvals <- frel * fabs
    dreal <- round(fn(log(fabs), frel, depl, dat, opt=2),3)

    if(verbose){
        print(paste0("Required depletion relative to ", depl.quant,
                    ": ", depl, " -- Realised: ", dreal, " with abs F equal to ",round(fabs,3)))
    }

    dat$FM <- fvals
    dat$Fs <- fvals / dat$nseasons
    dat$depl <- dreal

    return(dat)
}


#' @name estProd
#' @export
estProd <- function(dat, set= NULL,
                    ny = 100,
                    fmax = 100,
                    tsSplit = 3,
                    plot = TRUE){

    dat$ny <- ny
    ns <- dat$nseasons
    nt <- ny * ns

    ## noise
    if(is.null(set)) set <- checkSet()
    set$sigmaF <- 0
    set$sigmaR <- 0
    set$rhoR <- 0
    set$sigmaR0 <- 0
    set$sigmaH <- 0
    set$sigmaM <- 0
    set$sigmaMat <- 0
    set$sigmaImp <- 0

    ##
    len1 <- len3 <- floor(ny/tsSplit)
    len2 <- ny - len1 - len3

    ## increasing effort
    dat$FM <- c(rep(0,len1),
                   seq(0, fmax, length.out = len2),
                rep(fmax,len3))
    dat$Fs <- dat$FM / ns
    pop1 <- initPop(dat, set)
    tsb1 <- pop1$TSB[,1]
    cw1 <- apply(pop1$CW,1,sum)
    prod1 <- rep(NA, ny)
    for(i in 1:(ny-1)){
        prod1[i] <- tsb1[i+1] - tsb1[i] + cw1[i]
    }

    ## decreasing effort
    dat$FM <- c(rep(fmax, len1),
                   seq(fmax, 0, length.out = len2),
                   rep(0, len3))
    pop2 <- initPop(dat, set)
    dat$Fs <- dat$FM / ns
    tsb2 <- pop2$TSB[,1]
    cw2 <- apply(pop2$CW,1,sum)
    prod2 <- rep(NA, ny)
    for(i in 1:(ny-1)){
        prod2[i] <- tsb2[i+1] - tsb2[i] + cw2[i]
    }

    if(plot){

        plot(tsb1, prod1, ty='n',
             xlim = range(0,tsb1,tsb2),
             ylim = range(0,prod1,prod2,na.rm=TRUE))
        lines(tsb1, prod1, ty='b',
              col = "dodgerblue2")
        lines(tsb2, prod2, ty='b',
              col = "darkgreen")
        legend("topright",
               title = "Effort",
               legend = c("increasing","decreasing"),
               lty=1, col = c("dodgerblue2","darkgreen"))

        if(FALSE){
        plot(tsb1/refs$B0, prod1/refs$MSY, ty='n',
             xlim = c(0,1.05), ylim = c(0,1.5))
        lines(tsb1/refs$B0, prod1/refs$MSY, ty='b',
              col = "dodgerblue2")
        lines(tsb2/refs$B0, prod2/refs$MSY, ty='b',
              col = "darkgreen")
        legend("topright",
               title = "Effort",
               legend = c("increasing","decreasing"),
               lty=1, col = c("dodgerblue2","darkgreen"))
        }

        ## abs plot
        ## plot(tsb1/refs$B0, prod1, ty='n',
        ##      xlim = c(0,1), ylim = range(prod1,prod2,na.rm=TRUE))
        ## lines(tsb1/refs$B0, prod1, ty='b',
        ##       col = "dodgerblue2")
        ## lines(tsb2/refs$B0, prod2, ty='b',
        ##       col = "darkgreen")
    }

    ## CHECK: different production curves as a results of different age/length composition of stock (not at equilibrium age composition at given F, because F changes to quickly. If F changes small -> two curves are the same!

    res <- list(
        incr = data.frame(tsb = tsb1,
                          cw = cw1,
                          prod = prod1),
        decr = data.frame(tsb = tsb2,
                          cw = cw2,
                          prod = prod2)
    )

    return(res)

}




#' @name checkSet
#' @export
checkSet <- function(set = NULL){

    ## empty list
    if(is.null(set)) set <- list()

    ## lognormal SDs of noise
    if(is.null(set$sigmaR)) set$sigmaR <- 0
    if(is.null(set$sigmaF)) set$sigmaF <- 0
    if(is.null(set$sigmaM)) set$sigmaM <- 0
    if(is.null(set$sigmaH)) set$sigmaH <- 0
    if(is.null(set$sigmaR0)) set$sigmaR0 <- 0
    if(is.null(set$sigmaMat)) set$sigmaMat <- 0
    if(is.null(set$sigmaImp)) set$sigmaImp <- 0

    ## auto-correlated recruitment devs
    if(is.null(set$rhoR)) set$rhoR <- 0

    ## maximum F for baranov solution for F given TAC
    if(is.null(set$maxF)) set$maxF <- 5

    ## for estimation of ref levels
    if(is.null(set$refN)) set$refN <- 1e6
    if(is.null(set$refYears)) set$refYears <- 300

    ## number of years available to assessment method
    if(is.null(set$nyhist)) set$nyhist <- 35
    if(is.null(set$nysim)) set$nysim <- 35
    if(is.null(set$nrep)) set$nrep <- 50

    ## observation noise
    if(is.null(set$CVC)) set$CVC <- 0
    if(is.null(set$CVI)) set$CVI <- 0

    ## HCR
    if(is.null(set$hcr)) set$hcr <- c(defHCRref(),defHCRref(consF = "fmsy"))
    ## define constant catch rule
    defHCRconscat()
    if(is.null(set$stab)) set$stab <- FALSE

    ## Index timing:
    ## North Sea - IBITS: the majority of countries have only carried
    ## out a survey twice a year; a first quarter survey (January-February) and a
    ## third quarter survey (August-September)
    if(is.null(set$surveyTimes)) set$surveyTimes <- c(1/12,7/12)

    ## Seasonal catch observations
    if(is.null(set$catchSeasons)) set$catchSeasons <- 1

    ## burn in period
    if(is.null(set$burnin)) set$burnin <- 20

    ## return
    return(set)
}


#' @name fpat
#' @export
fpat <- function(fmax, fscen = 1){
    fscen <- as.character(fscen)
    switch(fscen,
           "1" = {  ## flat
               fs <- c(rep(0,20),                               ## no fishing - burn-in
                       seq(0, fmax, length.out = 10),           ## increasing effort
                       rep(fmax, 20))                           ## flat
           },
           "2" = {  ## decreasing
               fs <- c(rep(0,20),                               ## no fishing - burn-in
                       seq(0, fmax, length.out = 10),           ## increasing effort
                       rep(fmax, 10),                           ## flat
                       seq(fmax, 0.6 * fmax, length.out = 10))  ## decreasing effort
           },
           "3" = {  ## increasing
               fs <- c(rep(0,20),                               ## no fishing - burn-in
                       seq(0, fmax, length.out = 10),           ## increasing effort
                       rep(fmax, 10),                           ## flat
                       seq(fmax, 1.4 * fmax, length.out = 10))  ## increasing effort
           })
    return(fs)
}



#' @name baranov
#' @export
baranov <- function(FAA,M,NAA){
    Z <- FAA + M
    return(FAA/Z * NAA * (1 - exp(-Z)))
}


#' @name getFAA
#' @export
getFM <- function(TAC, NAA, M, weight, sel){
    tacEst <- function(FM, NAA, M, sel, weight, TAC){
        (TAC - sum(baranov(exp(FM) * sel, M, NAA) * weight))^2
    }
    opt <- optimise(tacEst, c(-10,10), NAA = NAA, M = M, TAC = TAC,
                    weight = weight, sel = sel)
    return(exp(opt$minimum))
}



#' @name getSel
#' @description Function to estimate selectivity ogive
#' @param L50 - length at 50% selectivity
#' @param L95 - length at 95% selectivity
#' @param mids - midlengths
#' @param plba - probability of being in mids given age
getSel <- function(L50, L95, mids, plba){
    selL <- (1 /(1 + exp(-log(19)*(mids - L50)/(L95 - L50))))
    dims <- dim(plba)
    selA <- matrix(NA, ncol = dims[3], nrow = dims[1])
    for(i in 1:dim(plba)[3]){
        selA[,i] <- apply(t(plba[,,i]) * selL, 2, sum)
    }
##    selA <- apply(t(plba) * selL, 2, sum)
##    selA[1] <- 1e-9 # it should be zero for age 0
    return(selA)
}


#' @name getMat
#' @description Function to estimate maturity at age
#' @param Lm50 - length at 50% maturity
#' @param Lm95 - length at 95% maturity
#' @param mids - midlengths
#' @param plba - probability of being in mids given age
getMat <- function(Lm50, Lm95, mids, plba){
    ## maturity at length
    matL <- (1 /(1 + exp(-log(19)*(mids - Lm50)/(Lm95 - Lm50))))
    ## maturity at age
    dims <- dim(plba)
    matA <- matrix(NA, ncol = dims[3], nrow = dims[1])
    for(i in 1:dim(plba)[3]){
        matA[,i] <- apply(t(plba[,,i]) * matL, 2, sum)
    }
##    matA <- apply(t(plba) * matL, 2, sum)
##    matA <- c(1e-9,matA[-1])
    return(matA)
}


#' @name getSSBPR0
#' @description Function to calculate spawners per recruit
#' @param M - natural mortality
#' @param mat - maturity ogive
#' @param fecun - fecundity matrix
#' @param amax - number of age classes
#' @return spawning biomass per recruit
#' @export
getSSBPR0 <- function(M, mat, fecun=1, amax){
    N <- rep(NA, amax)
    N[1] <- 1
    for(age in 2:amax)
        N[age] <- N[age-1] * exp(-M[age-1])
    N[amax] <- N[amax] / (1 - exp(-M[amax-1]))
    SBPR <- sum(N * mat * fecun)
    return(SBPR)
}


#' @name recfunc
#' @description Function to calculate recruitment (Beverton - Holt)
#' @param h - steepness
#' @param R0 - recruitment in unfished population
#' @param SSBPR0 - spawning biomass produced by one recrut in its lifetime
#' @param SSB - spawning biomass
#' @export
recfunc <- function(h, SSBPR0, SSB,  R0 = 1000, method = "bevholt"){
    if(method == "bevholt"){
        alpha <- SSBPR0 * ((1-h)/(4*h))
        beta <- (5*h-1)/(4*h*R0)
        rec <- SSB / (alpha + beta * SSB)
    }else if(method == "ricker"){
        beta <- log(5 * h) / (0.8 * R0)
        alpha <- exp(beta * R0)/SSBPR0
        rec <- alpha * SSB * exp(-beta * SSB)
    }else if(method == "average"){
        rec <- R0
    }else print("method not known!")
    return (rec)
}



#' @name estTAC
#' @export
estTAC <- function(inp, hcr, tacs=NULL){
    func <- get(hcr)
    tacs <- func(inp, tacs)
    return(tacs)
}








## REMOVE:

#' @name estTACOLD
#' @export
estTACOLD <- function(inp, hcr, hist=NULL, stab=FALSE, tacs=NULL, tcv=NA){

    switch(hcr,

           ## Bfacs
           ## --------------------------------------
           "spict-bfac" = {
               inp$reportmode <- 3
               inp$dteuler <- 1/4
               inp$stabilise <- 0
               inp$priors$logn <- c(0,0,0)
               inp$priors$logalpha <- c(0,0,0)
               inp$priors$logbeta <- c(0,0,0)
               inp <- check.inp(inp)
               if(is.null(tacs)){
                   indBpBx <- inp$indBpBx
               }else{
                   indBpBx <- tacs$indBpBx[nrow(tacs)]
               }
               inp$indBpBx <- indBpBx
               inp$phases$logn <- -1
               inp$ini$logn <- log(2)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
               }else{
                   tac <- try(spict:::get.TAC(rep = rep,
                                              bfac = 1,
                                              safeguardB = list(prob = 0.5),
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5,
                                                               bmsy = 0.5, fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-bfac", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- as.data.frame(c(tactmp,
                                                 fmfmsy.est=NA,fmfmsy.sd=NA,
                                                 bpbmsy.est=NA,bpbmsy.sd=NA,
                                                 cp.est=NA,cp.sd=NA,
                                                 fmsy.est=NA,fmsy.sd=NA,
                                                 bmsy.est=NA,bmsy.sd=NA,
                                                 indBpBx = indBpBx))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-bfac60" = {
               inp$reportmode <- 3
               inp$dteuler <- 1/4
               inp$stabilise <- 0
               inp$priors$logn <- c(0,0,0)
               inp$priors$logalpha <- c(0,0,0)
               inp$priors$logbeta <- c(0,0,0)
               inp <- check.inp(inp)
               if(is.null(tacs)){
                   indBpBx <- inp$indBpBx
               }else{
                   indBpBx <- tacs$indBpBx[nrow(tacs)]
               }
               inp$indBpBx <- indBpBx
               inp$phases$logn <- -1
               inp$ini$logn <- log(2)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
               }else{
                   tac <- try(spict:::get.TAC(rep = rep,
                                              bfac = 1,
                                              safeguardB = list(prob = 0.6),
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5,
                                                               bmsy = 0.5, fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-bfac60", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- as.data.frame(c(tactmp,
                                                 fmfmsy.est=NA,fmfmsy.sd=NA,
                                                 bpbmsy.est=NA,bpbmsy.sd=NA,
                                                 cp.est=NA,cp.sd=NA,
                                                 fmsy.est=NA,fmsy.sd=NA,
                                                 bmsy.est=NA,bmsy.sd=NA,
                                                 indBpBx = indBpBx))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-bfac70" = {
               inp$reportmode <- 3
               inp$dteuler <- 1/4
               inp$stabilise <- 0
               inp$priors$logn <- c(0,0,0)
               inp$priors$logalpha <- c(0,0,0)
               inp$priors$logbeta <- c(0,0,0)
               inp <- check.inp(inp)
               if(is.null(tacs)){
                   indBpBx <- inp$indBpBx
               }else{
                   indBpBx <- tacs$indBpBx[nrow(tacs)]
               }
               inp$indBpBx <- indBpBx
               inp$phases$logn <- -1
               inp$ini$logn <- log(2)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
               }else{
                   tac <- try(spict:::get.TAC(rep = rep,
                                              bfac = 1,
                                              safeguardB = list(prob = 0.7),
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5,
                                                               bmsy = 0.5, fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-bfac70", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- as.data.frame(c(tactmp,
                                                 fmfmsy.est=NA,fmfmsy.sd=NA,
                                                 bpbmsy.est=NA,bpbmsy.sd=NA,
                                                 cp.est=NA,cp.sd=NA,
                                                 fmsy.est=NA,fmsy.sd=NA,
                                                 bmsy.est=NA,bmsy.sd=NA,
                                                 indBpBx = indBpBx))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-bfac80" = {
               inp$reportmode <- 3
               inp$dteuler <- 1/4
               inp$stabilise <- 0
               inp$priors$logn <- c(0,0,0)
               inp$priors$logalpha <- c(0,0,0)
               inp$priors$logbeta <- c(0,0,0)
               inp <- check.inp(inp)
               if(is.null(tacs)){
                   indBpBx <- inp$indBpBx
               }else{
                   indBpBx <- tacs$indBpBx[nrow(tacs)]
               }
               inp$indBpBx <- indBpBx
               inp$phases$logn <- -1
               inp$ini$logn <- log(2)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
               }else{
                   tac <- try(spict:::get.TAC(rep = rep,
                                              bfac = 1,
                                              safeguardB = list(prob = 0.8),
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5,
                                                               bmsy = 0.5, fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-bfac80", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- as.data.frame(c(tactmp,
                                                 fmfmsy.est=NA,fmfmsy.sd=NA,
                                                 bpbmsy.est=NA,bpbmsy.sd=NA,
                                                 cp.est=NA,cp.sd=NA,
                                                 fmsy.est=NA,fmsy.sd=NA,
                                                 bmsy.est=NA,bmsy.sd=NA,
                                                 indBpBx = indBpBx))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-bfac90" = {
               inp$reportmode <- 3
               inp$dteuler <- 1/4
               inp$stabilise <- 0
               inp$priors$logn <- c(0,0,0)
               inp$priors$logalpha <- c(0,0,0)
               inp$priors$logbeta <- c(0,0,0)
               inp <- check.inp(inp)
               if(is.null(tacs)){
                   indBpBx <- inp$indBpBx
               }else{
                   indBpBx <- tacs$indBpBx[nrow(tacs)]
               }
               inp$indBpBx <- indBpBx
               inp$phases$logn <- -1
               inp$ini$logn <- log(2)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
               }else{
                   tac <- try(spict:::get.TAC(rep = rep,
                                              bfac = 1,
                                              safeguardB = list(prob = 0.9),
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5,
                                                               bmsy = 0.5, fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-bfac90", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- as.data.frame(c(tactmp,
                                                 fmfmsy.est=NA,fmfmsy.sd=NA,
                                                 bpbmsy.est=NA,bpbmsy.sd=NA,
                                                 cp.est=NA,cp.sd=NA,
                                                 fmsy.est=NA,fmsy.sd=NA,
                                                 bmsy.est=NA,bmsy.sd=NA,
                                                 indBpBx = indBpBx))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-bfac99" = {
               inp$reportmode <- 3
               inp$dteuler <- 1/4
               inp$stabilise <- 0
               inp$priors$logn <- c(0,0,0)
               inp$priors$logalpha <- c(0,0,0)
               inp$priors$logbeta <- c(0,0,0)
               inp <- check.inp(inp)
               if(is.null(tacs)){
                   indBpBx <- inp$indBpBx
               }else{
                   indBpBx <- tacs$indBpBx[nrow(tacs)]
               }
               inp$indBpBx <- indBpBx
               inp$phases$logn <- -1
               inp$ini$logn <- log(2)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
               }else{
                   tac <- try(spict:::get.TAC(rep = rep,
                                              bfac = 1,
                                              safeguardB = list(prob = 0.99),
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5,
                                                               bmsy = 0.5, fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-bfac99", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- as.data.frame(c(tactmp,
                                                 fmfmsy.est=NA,fmfmsy.sd=NA,
                                                 bpbmsy.est=NA,bpbmsy.sd=NA,
                                                 cp.est=NA,cp.sd=NA,
                                                 fmsy.est=NA,fmsy.sd=NA,
                                                 bmsy.est=NA,bmsy.sd=NA,
                                                 indBpBx = indBpBx))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },



           ## spict MSY
           "spict-msy" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
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
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5,
                                                               ffmsy=0.5,
                                                               bbmsy=0.5,
                                                               bmsy = 0.5,
                                                               fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE),
                              silent = TRUE)
                   if(inherits(tac, "try-error")){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                            fmsy, bmsy, indBpBx = NA))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-msy-fmsy30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
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
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5,
                                                               ffmsy=0.5,
                                                               bbmsy=0.5,
                                                               bmsy = 0.5,
                                                               fmsy = 0.3),
                                              breakpointB = 0, verbose = FALSE),
                              silent = TRUE)
                   if(inherits(tac, "try-error")){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                            fmsy, bmsy, indBpBx = NA))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-msy-bmsy30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
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
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5,
                                                               ffmsy=0.5,
                                                               bbmsy=0.5,
                                                               bmsy = 0.3,
                                                               fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE),
                              silent = TRUE)
                   if(inherits(tac, "try-error")){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                            fmsy, bmsy, indBpBx = NA))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-hs50" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
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
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5,
                                                               ffmsy=0.5,
                                                               bbmsy=0.5,
                                                               bmsy = 0.5,
                                                               fmsy = 0.5),
                                              breakpointB = 0.5, verbose = FALSE),
                              silent = TRUE)
                   if(inherits(tac, "try-error")){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                            fmsy, bmsy, indBpBx = NA))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-msy-a30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
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
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.3,
                                                               ffmsy=0.3,
                                                               bbmsy=0.3,
                                                               bmsy = 0.5,
                                                               fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE),
                              silent = TRUE)
                   if(inherits(tac, "try-error")){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                            fmsy, bmsy, indBpBx = NA))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-msy-b30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
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
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5,
                                                               ffmsy=0.5,
                                                               bbmsy=0.3,
                                                               bmsy = 0.5,
                                                               fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE),
                              silent = TRUE)
                   if(inherits(tac, "try-error")){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                            fmsy, bmsy, indBpBx = NA))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-msy-f30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
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
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5,
                                                               ffmsy=0.3,
                                                               bbmsy=0.5,
                                                               bmsy = 0.5,
                                                               fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE),
                              silent = TRUE)
                   if(inherits(tac, "try-error")){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                            fmsy, bmsy, indBpBx = NA))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "spict-msy-c30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
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
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.3,
                                                               ffmsy=0.5,
                                                               bbmsy=0.5,
                                                               bmsy = 0.5,
                                                               fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE),
                              silent = TRUE)
                   if(inherits(tac, "try-error")){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- data.frame(c(tactmp, fmfmsy, bpbmsy, cp,
                                            fmsy, bmsy, indBpBx = NA))
                       if(is.null(tacs)){
                           tacs <- tactmp
                       }else{
                           tacs <- rbind(tacs, tactmp)
                       }
                   }
               }
               return(tacs)
           },
           "2/3" = {
##               inp <- check.inp(inp)
               tacs <- r23(inp, stab=stab, red = NA, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "2/3_2red02" = {
               inp <- check.inp(inp)
               tacs <- r23(inp, stab=stab, red = 0.2, y_red = 2, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "2/3_3red02" = {
               inp <- check.inp(inp)
               tacs <- r23(inp, stab=stab, red = 0.2, y_red = 3, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "2/3_3red02_stab" = {
               inp <- check.inp(inp)
               tacs <- r23(inp, stab=TRUE, red = 0.2, y_red = 3, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "2/3_stab" = {
               inp <- check.inp(inp)
               tacs <- r23(inp, stab=TRUE, red = NA, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "1/2" = {
##               inp <- check.inp(inp)
               tacs <- r12(inp, stab=FALSE, red = NA, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="12", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "1/2_stab" = {
##               inp <- check.inp(inp)
               tacs <- r12(inp, stab=TRUE, red = NA, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="12", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "3/5" = {
##               inp <- check.inp(inp)
               tacs <- r35(inp, stab=FALSE, red = NA, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="35", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "3/5_stab" = {
               inp <- check.inp(inp)
               tacs <- r35(inp, stab=TRUE, red = NA, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="35", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "2/3_4red02" = {
               inp <- check.inp(inp)
               tacs <- r23(inp, stab=stab, red = 0.2, y_red = 4, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "2/3_tarC" = {
               inp <- check.inp(inp)
               tacs <- r23(inp, stab=stab, red = NA, tacs=tacs, targetC = TRUE)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "dm" = {
               inp <- check.inp(inp)
               tacs <- derimeth(inp, stab=stab, red = NA, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "dm_tarC" = {
               inp <- check.inp(inp)
               tacs <- derimeth(inp, stab=stab, red = NA, tacs=tacs)##, targetC = TRUE)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="23", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(is.null(tacs)){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "conscat" = {
               inp <- check.inp(inp)
               tacs <- conscat(inp, tacs=tacs)
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="cc", hitSC=NA, red=NA,
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
               }
               return(tacs)
           },
           "refFmsy" = {
               tactmp <- data.frame(TAC=NA, id="refFmsy", hitSC=NA,
                                  red=NA, barID=NA, sd=NA, conv=NA)
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
               }
               return(tacs)
           },
           "noF" = {
               tactmp <- data.frame(TAC=0, id="noF", hitSC=NA,
                                  red=NA, barID=NA, sd=NA, conv=NA)
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
               }
               return(tacs)
           },
           stop("HCR not known!"))
}
