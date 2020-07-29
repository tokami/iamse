

#' @name genConvs
#' @description Get converged simulates from a resMSE object
#' @export
getConvs <- function(mse, convyears = "all", convhcrs = "all", verbose = FALSE){

    nhcr <- length(mse)
    nrep <- length(mse[[1]])
    nquant <- length(mse[[1]][[1]])
    nysim <- nrow(mse[[1]][[1]]$tacs)
    dims <- dim(mse[[1]][[1]]$CW)
    ny <- dims[1] - nysim
    ns <- dims[2]

    if(convyears[1] == "all") convyears <- 1:ny
    if(convhcrs[1] == "all") convhcrs <- 1:nhcr


    indlist <- vector("list",nhcr)
    for(hcr in 1:nhcr){
        tmp <- do.call(rbind,lapply(mse[[hcr]], function(x) x[["tacs"]]$conv))
        indlist[[hcr]] <- apply(tmp[,convyears], 1, all)
    }

    ## across hcrs only
    indlist2 <- do.call(cbind,indlist[convhcrs])


    ## ## across hcrs and scens
    ## tmp <- lapply(indlist, function(x) do.call(cbind,x[convhcrs]))
    ## tmp2 <- do.call(cbind, tmp)
    ## print(5*length(which(apply(tmp2, 1, all))))

    inds <- which(apply(indlist2,1,all))
    if(verbose)
        writeLines(paste0("Converged reps: ",length(inds), " of ",nrep,
                          " reps = ",round(length(inds)/nrep*100),"%"))
    res <- vector("list",nhcr)
    for(hcr in 1:nhcr){
        res[[hcr]] <- mse[[hcr]][inds]
        names(res[[hcr]]) <- inds
    }

    ## return
    return(res)

}


#' @name genDevs
#' @export
genDevs <- function(n, sd, rho=0){

    rnoise <- rnorm(n, 0, sd) - sd^2/2

    res <- numeric(n)
    res[1] <- rnoise[1]
    if(n > 1){
        for(i in 2:n) res[i] <- rho * res[i-1] + sqrt(1 - rho^2) * rnoise[i]
    }

    res <- exp(res)
    res <- res/mean(res)

    return(res)

}

#' @name estDepl
#' @export
estDepl <- function(dat, set=NULL, fmax = 10, nrep = 100, verbose = TRUE){

    if(!any(names(dat) == "ref")) stop("Reference points are missing in dat. Use estRef to estimate reference points.")

    if(is.null(set)){
        set <- checkSet()
        nrep <- 1
    }

    depl <- dat$depl
    depl.quant <- dat$depl.quant
    depl.prob <- dat$depl.prob

    frel <- dat$FM/max(dat$FM)

    fn <- function(logfabs, frel, depl, depl.prob, nrep, dat, set, opt=1){
        datx <- dat
        fpat <- frel * exp(logfabs)
        datx$FM <- fpat
        datx$Fs <- fpat / datx$nseasons
        dreal <- sapply(1:nrep, function(x) initPop(datx, set, out.opt = 2))
        drealQ <- quantile(dreal, probs = depl.prob)
        if(opt==1) return((depl - drealQ)^2)
        if(opt==2) return(drealQ)
    }

    opt <- optimize(fn, c(log(0.0001),log(fmax)), frel = frel, depl = depl, depl.prob = depl.prob,
                    nrep = nrep, dat = dat, set=set)
    fabs <- exp(opt$minimum)
    fvals <- frel * fabs
    dreal <- round(fn(log(fabs), frel, depl, depl.prob, nrep, dat, set, opt=2),3)

    if(verbose){
        print(paste0("Target depletion level as ",depl.prob * 100, "% quantile of ", depl, " ", depl.quant,
                     " -- Realised depletion level at ", dreal, " ", depl.quant,
                     " with absolute F equal to ",round(fabs,3)))
    }

    dat$FM <- fvals
    dat$Fs <- fvals / dat$nseasons

    return(dat)
}


estDeplOLD <- function(dat, set=NULL, fmax = 10, verbose = TRUE){

    if(!any(names(dat) == "ref")) stop("Reference points are missing in dat. Use estRef to estimate reference points.")

    browser()

    if(is.null(set)) set <- checkSet()

    depl <- dat$depl
    depl.quant <- dat$depl.quant

    frel <- dat$FM/max(dat$FM)

    fn <- function(logfabs, frel, depl, dat, set, opt=1){
        datx <- dat
        fpat <- frel * exp(logfabs)
        datx$FM <- fpat
        datx$Fs <- fpat / datx$nseasons
        dreal <- initPop(datx, set, out.opt = 2, depl.quant = depl.quant)
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
                    fmax = 10,
                    tsSplit = 8,
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
    if(is.null(set$refN)) set$refN <- 1e4
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

    ## seed
    if(is.null(set$seed)) set$seed <- NA

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


#' @name getSSBPR
#' @description Function to calculate spawners per recruit
#' @param Z - total mortality
#' @param mat - maturity ogive
#' @param fecun - fecundity matrix
#' @param amax - number of age classes
#' @return spawning biomass per recruit
#' @export
getSSBPR <- function(M, mat, weight, fecun=1, amax, R0 = 1, FM = NULL){
    N <- rep(NA, amax)
    N[1] <- R0
    M <- cumsum(M)
    if(!is.null(FM)) Z = M + FM else Z = M
    N[2:(amax-1)] <- R0 * exp(-Z[1:(amax-2)])
    N[amax] <- R0 * exp(-Z[amax-1]) / (1-exp(-Z[amax]))
    SBPR <- sum(N * mat * weight * fecun)
    return(SBPR)
}


#' @name recfunc
#' @description Function to calculate recruitment (Beverton - Holt)
#' @param h - steepness
#' @param R0 - recruitment in unfished population
#' @param SSBPR0 - spawning biomass produced by one recrut in its lifetime
#' @param SSB - spawning biomass
#' @param bp - breakpoint for hockey-stick SR
#' @param method - SR type
#'
#' @export
recfunc <- function(h, SSBPR0, SSB,  R0 = 1e6, method = "bevholt", bp = 0,
                    beta = 0, gamma = 0){

    if(method == "bevholt"){
        ## alpha <- SSBPR0 * ((1-h)/(4*h))
        ## beta <- (5*h-1)/(4*h*R0)
        ## rec <- SSB / (alpha + beta * SSB)

        rec <- (4 * h * R0 * SSB / (SSBPR0 * (1-h) + SSB * (5*h-1)))

    }else if(method == "ricker"){
        beta <- log(5 * h) / (0.8 * R0)
        alpha <- exp(beta * R0)/SSBPR0
        rec <- alpha * SSB * exp(-beta * SSB)
    }else if(method == "average"){
        rec <- R0
    }else if(method == "hockey-stick"){
        rec <- ifelse(SSB > bp, R0, SSB * R0/bp)
    }else if(method == "bent-hyperbola"){  ## Watts-Bacon bent hyperbola
        rec <- beta * (SSB + sqrt(bp^2 + (gamma^2)/4) - sqrt((SSB-bp)^2 + (gamma^2)/4))
    }else print("Stock-recruitment method not known! Implemented methods: 'bevholt', 'ricker', 'average', and 'hockey-stick'.")

    return (rec)
}



#' @name estTAC
#' @export
estTAC <- function(inp, hcr, tacs=NULL){
    func <- get(hcr)
    tacs <- func(inp, tacs)
    return(tacs)
}
