#' @name stocklist
#' @title Fisheries data included in Polacheck et al. (1993).
#' @details Fisheries data for south Atlantic albacore, northern Namibian hake, and New Zealand rock lobster.
#' @docType data
#' @keywords datasets
#' @usage data(stocklist)
#' @format Data are lists containing data
NULL



#' @name checkSet
#' @export
checkSet <- function(set = NULL){
    if(is.null(set)) set <- list()

    ## lognormal SDs of noise
    if(is.null(set$sigmaR)) set$sigmaR <- 0
    if(is.null(set$sigmaF)) set$sigmaF <- 0
    if(is.null(set$sigmaM)) set$sigmaM <- 0
    if(is.null(set$sigmaH)) set$sigmaH <- 0
    if(is.null(set$sigmaR0)) set$sigmaR0 <- 0
    if(is.null(set$sigmaMat)) set$sigmaMat <- 0

    ## auto-correlated recruitment devs (not yet implemented)
    if(is.null(set$autocor)) set$autocor <- FALSE

    ## maximum F for baranov solution for F given TAC
    if(is.null(set$maxF)) set$maxF <- 5

    ## for estimation of ref levels
    if(is.null(set$refN)) set$refN <- 1e4
    if(is.null(set$refYears)) set$refYears <- 150

    if(is.null(set$nyhist)) set$nyhist <- 50
    if(is.null(set$nysim)) set$nysim <- 50
    if(is.null(set$nrep)) set$nrep <- 50

    ## observation noise
    if(is.null(set$CVC)) set$CVC <- 0
    if(is.null(set$CVI)) set$CVI <- 0

    ## HCR
    if(is.null(set$hcr)) set$hcr <- c("refFmsy","noF")
    if(is.null(set$stab)) set$stab <- FALSE

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
           "pbbfred" = {
               inp <- check.inp(inp)
               inp$phases$logn <- -1
               inp$ini$logn <- log(2)
               inp$MSEmode <- 3
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
               }else{
                   tacs <- spictPBBFred(rep, bfrac=1, prob = 0.5, stab=stab, tacs=tacs)
               }

               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbbfred", hitSC=NA,
                                        red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
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
                   tactmp <- conscat(inp, tacs=tacs)
                   fmfmsy <- c(NA,NA,NA)
                   bpbmsy <- c(NA,NA,NA)
                   cp <- c(NA,NA,NA)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBmBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                   }

               }
               tactmp <- as.data.frame(c(tactmp, fmfmsy, bpbmsy, cp))
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
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
                   fmfmsy <- c(NA,NA,NA)
                   bpbmsy <- c(NA,NA,NA)
                   cp <- c(NA,NA,NA)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBmBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
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
                   }

               }
               tactmp <- as.data.frame(c(tactmp, fmfmsy, bpbmsy, cp))
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
               }
               return(tacs)
           },
           "spict-msy-a30" = {
               inp$reportmode <- 1
               inp$dteuler <- 1/4
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 || any(is.infinite(rep$sd))){
                   tactmp <- conscat(inp, tacs=tacs)
                   fmfmsy <- c(NA,NA,NA)
                   bpbmsy <- c(NA,NA,NA)
                   cp <- c(NA,NA,NA)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBmBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.3, ffmsy=0.3, bbmsy=0.3),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy-a30", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                   }

               }
               tactmp <- as.data.frame(c(tactmp, fmfmsy, bpbmsy, cp))
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
                   fmfmsy <- c(NA,NA,NA)
                   bpbmsy <- c(NA,NA,NA)
                   cp <- c(NA,NA,NA)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBmBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5, ffmsy=0.5, bbmsy=0.3),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy-b30", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                   }

               }
               tactmp <- as.data.frame(c(tactmp, fmfmsy, bpbmsy, cp))
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
                   fmfmsy <- c(NA,NA,NA)
                   bpbmsy <- c(NA,NA,NA)
                   cp <- c(NA,NA,NA)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBmBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.5, ffmsy=0.3, bbmsy=0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy-f30", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                   }

               }
               tactmp <- as.data.frame(c(tactmp, fmfmsy, bpbmsy, cp))
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
                   fmfmsy <- c(NA,NA,NA)
                   bpbmsy <- c(NA,NA,NA)
                   cp <- c(NA,NA,NA)
               }else{
                   fmfmsy <- round(get.par("logFmFmsynotS",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(fmfmsy) <- paste0("fmfmsy-",names(fmfmsy))
                   bpbmsy <- round(get.par("logBmBmsy",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(bpbmsy) <- paste0("bpbmsy-",names(bpbmsy))
                   cp <- round(get.par("logCp",rep, exp=TRUE)[,c(2,4,5)],2)
                   names(cp) <- paste0("cp-",names(cp))
                   tac <- try(spict:::get.TAC(rep=rep,
                                              fractiles = list(catch=0.3, ffmsy=0.5, bbmsy=0.5),
                                              breakpointB = 0, verbose = FALSE), silent = TRUE)
                   if(class(tac) == "try-error"){
                       tactmp <- conscat(inp, tacs=tacs)
                   }else{
                       tactmp <- data.frame(TAC=tac, id="spict-msy-c30", hitSC=NA,
                                            red=NA, barID=NA, sd=NA, conv = NA)
                   }

               }
               tactmp <- as.data.frame(c(tactmp, fmfmsy, bpbmsy, cp))
               if(is.null(tacs)){
                   tacs <- tactmp
               }else{
                   tacs <- rbind(tacs, tactmp)
               }
               return(tacs)
           },
           "pbb065-msy025" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.65, stab=stab, tacs=tacs, tcv=0.25)
                   }
               }else{
                   tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.65, stab=stab, tacs=tacs, tcv=0.25)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb_msy", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb065-msy05" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.65, stab=stab, tacs=tacs, tcv=0.5)
                   }
               }else{
                   tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.65, stab=stab, tacs=tacs, tcv=0.5)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb_msy", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb065-msy10" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.65, stab=stab, tacs=tacs, tcv=1)
                   }
               }else{
                   tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.65, stab=stab, tacs=tacs, tcv=1)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb_msy", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb075-msy025" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.75, stab=stab, tacs=tacs, tcv=0.25)
                   }
               }else{
                   tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.75, stab=stab, tacs=tacs, tcv=0.25)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb_msy", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb075-msy05" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.75, stab=stab, tacs=tacs, tcv=0.5)
                   }
               }else{
                   tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.75, stab=stab, tacs=tacs, tcv=0.5)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb_msy", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb075-msy10" = {
##               browser()
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               ## plot(rep)
               ## (flfmsy <- get.par("logFlFmsy", rep, exp=TRUE)[5])
               ## (blbmsy <- get.par("logBlBmsy", rep, exp=TRUE)[5])
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                       tacs <- conscat(inp, tacs=tacs)
                   }else{
                       tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.75, stab=stab, tacs=tacs, tcv=1)
                   }
               }else{
                   tacs <- spictPBB_MSY(rep, bfrac=1, prob = 0.75, stab=stab, tacs=tacs, tcv=1)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb_msy", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb05" = {
               inp <- check.inp(inp)
##               if(!is.null(tacs) && nrow(tacs) == 4) browser()
##               inp$priors$logn <- c(0,0,0)
##               inp$priors$logalpha <- c(0,0,0)
##               inp$priors$logbeta <- c(0,0,0)
               rep <- try(fit.spict(inp), silent=TRUE)
               plot(rep)

               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
                   }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.5, stab=stab, tacs=tacs)
                   }
               }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.5, stab=stab, tacs=tacs)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb05", hitSC=NA,
                                        red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb05S" = {
##               if(!is.null(tacs) && nrow(tacs) == 2) browser()
               inp <- check.inp(inp)
               inp$priors$logn <- c(log(2),0.2,1)
               rep <- try(fit.spict(inp), silent=TRUE)
               btmp <- get.par("logB",rep,exp=TRUE)[,2]
               print(tail(btmp,17))
               print(tail(btmp,1) > tail(btmp,16)[1])
               plot(rep)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
                   }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.5, stab=stab, tacs=tacs)
                   }
               }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.5, stab=stab, tacs=tacs)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb05", hitSC=NA,
                                        red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb055" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
                   }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.55, stab=stab, tacs=tacs)
                   }
               }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.55, stab=stab, tacs=tacs)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb05", hitSC=NA,
                                        red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb06" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
##               plot(rep)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
                   }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.6, stab=stab, tacs=tacs)
                   }
               }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.6, stab=stab, tacs=tacs)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb06", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb07" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
                   }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.7, stab=stab, tacs=tacs)
                   }
               }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.7, stab=stab, tacs=tacs)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb07", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }
               }
               return(tacs)
           },
           "pbb08" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
                   }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.8, stab=stab, tacs=tacs)
                   }
               }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.8, stab=stab, tacs=tacs)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb08", hitSC=NA,
                                        red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
                   }

               }
               return(tacs)
           },
           "pbb09" = {
               inp <- check.inp(inp)
               rep <- try(fit.spict(inp), silent=TRUE)
               if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                  any(is.infinite(rep$sd))){
                   inp$priors$logn <- c(log(2),1e-4,1)
                   rep <- try(fit.spict(inp), silent=TRUE)
                   if(class(rep) == "try-error" || rep$opt$convergence != 0 ||
                      any(is.infinite(rep$sd))){
                   tacs <- conscat(inp, tacs=tacs)
                   }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.9, stab=stab, tacs=tacs)
                   }
               }else{
                   tacs <- spictPBB(rep, bfrac=1, prob = 0.9, stab=stab, tacs=tacs)
               }
               if(is.na(tacs$TAC[nrow(tacs)])){
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="pbb09", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
                       tacs <- tactmp
                   }else{
                       tacs <- rbind(tacs, tactmp)
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   if(nrow(tacs) == 1){
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
                   tactmp <- data.frame(TAC=inp$obsC[length(inp$obsC)], id="cc", hitSC=NA,
                                      red=NA, barID=NA, sd=NA, conv=NA)
                   if(nrow(tacs) == 1){
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
