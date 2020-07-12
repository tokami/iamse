#' @name stocklist
#' @title Fisheries data included in Polacheck et al. (1993).
#' @details Fisheries data for south Atlantic albacore, northern Namibian hake, and New Zealand rock lobster.
#' @docType data
#' @keywords datasets
#' @usage data(stocklist)
#' @format Data are lists containing data
NULL

#' @name genDevs
#' @export
genDevs <- function(n, sd, rho=0){


    rnoise <- rnorm(n, 0, sd) - sd^2/2

    res <- numeric(n)
    res[1] <- rnoise[1]
    for(i in 2:n) res[i] <- rho * res[i-1] + (1 - rho) * rnoise[i] ## sqrt(1-rho^2) * rnoise[i]

    res <- exp(res)
    res <- res/mean(res)

    return(res)

}

#' @name estDepl
#' @export
estDepl <- function(dat, refs, verbose = TRUE){

    depl <- dat$depl
    depl.quant <- dat$depl.quant

    frel <- dat$Fvals/max(dat$Fvals)

    fn <- function(fabs, frel, depl, opt=1){
        fpat <- frel * fabs
        dat$Fvals <- fpat
        dreal <- initPop(dat, refs = ref$refs, out.opt = 2, depl.quant = depl.quant)
        if(opt==1) return((depl - dreal)^2)
        if(opt==2) return(dreal)
    }

    opt <- optimize(fn, c(0,20), frel = frel, depl = depl)
    fabs <- opt$minimum
    fvals <- frel * fabs

    dreal <- round(fn(fabs, frel, depl, opt=2),3)

    if(verbose){
        print(paste0("Required depletion relative to ", depl.quant,
                    ": ", depl, " -- Realised: ", dreal, " with abs F equal to ",round(fabs,3)))
    }

    dat$Fvals <- fvals
    dat$depl <- dreal

    return(dat)
}


#' @name estProd
#' @export
estProd <- function(dat, set, refs, ny = 100, plot = TRUE){

    dat$ny <- ny
    set$sigmaF <- 0
    set$sigmaR <- 0
    set$sigmaR0 <- 0
    set$sigmaH <- 0
    set$sigmaM <- 0
    set$sigmaMat <- 0

    ## increasing effort
    dat$Fvals <- c(rep(0,ny/4),
                   seq(0, 100*refs$Fmsy, length.out = ny/2),
                   rep(100*refs$Fmsy,ny/4))
    pop1 <- initPop(dat, set)
    tsb1 <- pop1$TSB
    cw1 <- pop1$CW
    prod1 <- rep(NA, ny)
    for(i in 1:(ny-1)){
        prod1[i] <- tsb1[i+1] - tsb1[i] + cw1[i]
    }

    ## decreasing effort
    dat$Fvals <- c(rep(100*refs$Fmsy, ny/4),
                   seq(100*refs$Fmsy, 0, length.out = ny/2),
                   rep(0, ny/4))
    pop2 <- initPop(dat, set)
    tsb2 <- pop2$TSB
    cw2 <- pop2$CW
    prod2 <- rep(NA, ny)
    for(i in 1:(ny-1)){
        prod2[i] <- tsb2[i+1] - tsb2[i] + cw2[i]
    }


    if(plot){

        plot(tsb1/refs$B0, prod1/refs$MSY, ty='n',
             xlim = c(0,1.05), ylim = c(0,1.5))
        lines(tsb1/refs$B0, prod1/refs$MSY, ty='b',
              col = "dodgerblue2")
        lines(tsb2/refs$B0, prod2/refs$MSY, ty='b',
              col = "darkgreen")

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
    if(is.null(set$refN)) set$refN <- 1e4
    if(is.null(set$refYears)) set$refYears <- 150

    ## number of years available to assessment method
    if(is.null(set$nyhist)) set$nyhist <- 35
    if(is.null(set$nysim)) set$nysim <- 35
    if(is.null(set$nrep)) set$nrep <- 50

    ## observation noise
    if(is.null(set$CVC)) set$CVC <- 0
    if(is.null(set$CVI)) set$CVI <- 0

    ## HCR
    if(is.null(set$hcr)) set$hcr <- c("refFmsy","noF")
    if(is.null(set$stab)) set$stab <- FALSE

    ## Index timing:
    ## North Sea - IBITS: the majority of countries have only carried
    ## out a survey twice a year; a first quarter survey (January-February) and a
    ## third quarter survey (August-September)
    if(is.null(set$surveyTimes)) set$surveyTimes <- c(1/12,7/12)

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
    selA <- apply(t(plba) * selL, 2, sum)
    selA[1] <- 1e-9 # it should be zero for age 0
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
    matA <- apply(t(plba) * matL, 2, sum)
    matA <- c(1e-9,matA[-1])
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
    }else print("method not known!")
    return (rec)
}


#' @name estTAC
#' @export
estTAC <- function(inp, hcr, hist=NULL, stab=FALSE, tacs=NULL, tcv=NA){

    switch(hcr,
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
                                            red=NA, barID=NA, sd=NA, conv = NA,
                                            fmfmsy=fmfmsy, bpbmsy=bpbmsy, cp=cp,
                                            fmsy=fmsy, bmsy=bmsy, indBpBx = NA)
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
                   if(class(tac) == "try-error"){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA,
                                            fmfmsy=fmfmsy, bpbmsy=bpbmsy, cp=cp,
                                            fmsy=fmsy, bmsy=bmsy,
                                            indBpBx = NA)
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
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5,
                                                               bmsy = 0.3, fmsy = 0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA,
                                            fmfmsy=fmfmsy, bpbmsy=bpbmsy, cp=cp,
                                            fmsy=fmsy, bmsy=bmsy,
                                            indBpBx = NA)
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
                   tactmp <- conscat(inp, tacs=tacs)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBpBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5),
                                              breakpointB = 0.5, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-hs50", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                       tactmp <- as.data.frame(c(tactmp, fmfmsy, bpbmsy, cp))
                   }
               }
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
               }
               rm(rep)
               return(tacs)
           },
           "spict-msy-a30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tactmp <- conscat(inp, tacs=tacs)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBpBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.3, ffmsy=0.3, bbmsy=0.3),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA,
                                            fmfmsy=fmfmsy, bpbmsy=bpbmsy, cp=cp,
                                            fmsy=fmsy, bmsy=bmsy,
                                            indBpBx = NA)
                   }
               }
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
               }
               return(tacs)
           },
           "spict-msy-b30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tactmp <- conscat(inp, tacs=tacs)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBpBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.3),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA,
                                            fmfmsy=fmfmsy, bpbmsy=bpbmsy, cp=cp,
                                            fmsy=fmsy, bmsy=bmsy,
                                            indBpBx = NA)
                   }
               }
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
               }
               return(tacs)
           },
           "spict-msy-f30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tactmp <- conscat(inp, tacs=tacs)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBpBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5, ffmsy=0.3, bbmsy=0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA,
                                            fmfmsy=fmfmsy, bpbmsy=bpbmsy, cp=cp,
                                            fmsy=fmsy, bmsy=bmsy,
                                            indBpBx = NA)
                   }
               }
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
               }
               return(tacs)
           },
           "spict-msy-c30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tactmp <- conscat(inp, tacs=tacs)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBpBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.3, ffmsy=0.5, bbmsy=0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA,
                                            fmfmsy=fmfmsy, bpbmsy=bpbmsy, cp=cp,
                                            fmsy=fmsy, bmsy=bmsy,
                                            indBpBx = NA)
                   }
               }
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
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
