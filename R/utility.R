

#' @name genConvs
#' @description Get converged simulates from a resMSE object
#' @export
getConvs <- function(mse, convyears = "all", convhcrs = "all", out = 0, verbose = FALSE){

    nhcr <- length(mse)
    nrep <- length(mse[[1]])
    nquant <- length(mse[[1]][[1]])
    nysim <- nrow(mse[[1]][[1]]$tacs)
    dims <- dim(mse[[1]][[1]]$CW)
    ny <- dims[1] - nysim
    ns <- dims[2]

    if(!is.na(convyears[1]) && convyears[1] == "all")
        convyears <- 1:ny
    if(!is.na(convhcrs[1]) && convhcrs[1] == "all")
        convhcrs <- 1:nhcr

    if(is.numeric(convyears[1]) && is.numeric(convhcrs[1])){
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
    }else{
        hcr <- 3
        res <- vector("list",nhcr)
        for(hcr in 1:nhcr){
            id <- unique(mse[[hcr]][[1]]$tacs$id)[1]
            id2 <- unlist(strsplit(as.character(id), "-"))[1]
            if(!(id2 %in% c("noF","refFmsy"))){
                tmp <- do.call(rbind,lapply(mse[[hcr]], function(x) x[["tacs"]]$conv))
                inds <- which(apply(tmp[,convyears], 1, all))
                as.numeric(inds)
            }else inds <- 1:nrep
            if(verbose)
                writeLines(paste0("Converged reps: ",length(inds), " of ",nrep,
                                  " reps = ",round(length(inds)/nrep*100),"%"))
            res[[hcr]] <- mse[[hcr]][inds]
            names(res[[hcr]]) <- inds
        }
    }


    ## return
    if(out == 0){
        return(res)
    }else if(out == 1){
        return(sapply(res, length))
    }

}


#' @name sdconv
#' @export
sdconv <- function(mu, sd) (log(1 + ((sd^2)/(mu^2))))^0.5


#' @name muconv
#' @export
muconv <- function(mu, sd) log(mu) - 0.5 * log(1 + ((sd^2)/(mu^2)))


#' @name genNoise
#' @export
genNoise <- function(n, sd, rho=0, bias.cor = 0){

    if(bias.cor == 0){
        rnoise <- rnorm(n, 0, sd)
    }else if(bias.cor == 1){
        rnoise <- rnorm(n, 0, sd) - sd^2/2
    }else if(bias.cor == 2){
        rnoise <- log(rlnorm(n, muconv(1,sd), sdconv(1,sd)))
    }

    if(rho > 0){
        res <- numeric(n)
        res[1] <- rnoise[1]
        if(n > 1){
            for(i in 2:n) res[i] <- rho * res[i-1] + sqrt(1 - rho^2) * rnoise[i]
        }

        res <- exp(res)
        res <- res/mean(res)

    }else{
        res <- exp(rnoise)
    }

    return(res)
}


#' @name estDepl
#' @export
estDepl <- function(dat, set=NULL, fmax = 10, nrep = 100, verbose = TRUE, method = "percentile"){

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
        if(method == "mean"){
            drealQ <- mean(dreal)
        }else if(method == "median"){
            drealQ <- quantile(dreal, probs = 0.5)
        }else if(method == "percentile"){
            drealQ <- quantile(dreal, probs = depl.prob)
        }
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
    set$noiseF <- c(0,0,0)
    set$noiseR <- c(0,0,0)
    set$noiseR0 <- c(0,0,0)
    set$noiseH <- c(0,0,0)
    set$noiseM <- c(0,0,0)
    set$noiseMat <- c(0,0,0)
    set$noiseImp <- c(0,0,0)

    ##
    len1 <- len3 <- floor(ny/tsSplit)
    len2 <- ny - len1 - len3

    ## increasing effort
    dat$FM <- c(rep(0,len1),
                   seq(0, fmax, length.out = len2),
                rep(fmax,len3))
    dat$Fs <- dat$FM / ns
    pop1 <- initPop(dat, set)
    tsb1 <- pop1$TSBfinal
    cw1 <- apply(pop1$CW,1,sum)
    prod1 <- rep(NA, ny)
    for(i in 2:ny){
        prod1[i] <- tsb1[i] - tsb1[i-1] + cw1[i]
    }

    ## decreasing effort
    dat$FM <- c(rep(fmax, len1),
                   seq(fmax, 0, length.out = len2),
                   rep(0, len3))
    pop2 <- initPop(dat, set)
    dat$Fs <- dat$FM / ns
    tsb2 <- pop2$TSBfinal
    cw2 <- apply(pop2$CW,1,sum)
    prod2 <- rep(NA, ny)
    for(i in 2:ny){
        prod2[i] <- tsb2[i] - tsb2[i-1] + cw2[i]
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

    ## noise vectors (SD, rho, bias correction)
    if(is.null(set$noiseR)){
        set$noiseR <- c(0,0,0)
    }else if(length(set$noiseR) != 3){
        writeLines("'set$noiseR' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseF)){
        set$noiseF <- c(0,0,0)
    }else if(length(set$noiseF) != 3){
        writeLines("'set$noiseF' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseM)){
        set$noiseM <- c(0,0,0)
    }else if(length(set$noiseM) != 3){
        writeLines("'set$noiseM' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseH)){
        set$noiseH <- c(0,0,0)
    }else if(length(set$noiseH) != 3){
        writeLines("'set$noiseH' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseR0)){
        set$noiseR0 <- c(0,0,0)
    }else if(length(set$noiseR0) != 3){
        writeLines("'set$noiseR0' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseMat)){
        set$noiseMat <- c(0,0,0)
    }else if(length(set$noiseMat) != 3){
        writeLines("'set$noiseMat' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseImp)){
        set$noiseImp <- c(0,0,0)
    }else if(length(set$noiseImp) != 3){
        writeLines("'set$noiseImp' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseC)){
        set$noiseC <- c(0,0,0)
    }else if(length(set$noiseC) != 3){
        writeLines("'set$noiseC' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseI)){
        set$noiseI <- c(0,0,0)
    }else if(length(set$noiseI) != 3){
        writeLines("'set$noiseI' needs to be a vector with 3 values corresponding to: sd, rho, bias.corr (see genNoise). Setting to c(0,0,0)!")
    }

    ## maximum F for baranov solution for F given TAC
    if(is.null(set$maxF)) set$maxF <- 5

    ## for estimation of ref levels
    if(is.null(set$refN)) set$refN <- 1e3
    if(is.null(set$refYears)){
        set$refYears <- 300
    }else if(set$refYears < 100){
        writeLines("set$refYears has to be at least 100 years, as the median surplus production over the last 50 years are used for reference estimation. Setting to 100.")
        set$refYears <- 100
    }
    if(is.null(set$refYearsMSY)){
        set$refYearsMSY <- 100
    }else if(set$refYearsMSY >= set$refYears){
        writeLines("set$refYearsMSY cannot be longer than 'set$refYears'. Setting to half of 'set$refYears'.")
        set$refYearsMSY <- floor(set$refYears / 2)
    }
    if(is.null(set$refMethod)){
        set$refMethod <- "mean"
    }else if(set$refYears < 100){
        writeLines("'set$refMethod' has to be mean or median. Setting to 'mean'.")
        set$refMethod <- "mean"
    }

    ## number of years available to assessment method
    if(is.null(set$nyhist)) set$nyhist <- 35
    if(is.null(set$nysim)) set$nysim <- 35
    if(is.null(set$nrep)) set$nrep <- 50

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


#' @name getFM
#' @export
getFM <- function(TAC, NAA, M, weight, sel){
    tacEst <- function(FM, NAA, M, sel, weight, TAC){
        (TAC - sum(baranov(exp(FM) * sel, M, NAA*exp(-M/2)) * weight))^2
    }
    opt <- optimise(tacEst, c(-10,10), NAA = NAA, M = M, TAC = TAC,
                    weight = weight, sel = sel)
    return(exp(opt$minimum))
}


#' @name getFM2
#' @description Hybrid method (Methot and Wetzel)
#' @export
getFM2 <- function(TAC, TSB, ds, M, NAA, weight, weightF, sel, fmax = 1){
    fout <- 0
    if(TAC > 0){
        tmp <- TAC/(TSB + 0.1*TAC)
        j <- (1 + exp(30 * (tmp - 0.95)))^-1
        tmp2 <- j * tmp + 0.95 * (1 - j)
        fout <- -log(1 - tmp2) / ds
        for(i in 1:4){
            Fs <- fout * sel
            Z <- M + Fs
            Alpha <- (1 - exp(-Z))
            Ctmp <- sum((Fs/Z) * (NAA * weightF * sel) * Alpha)

            Zadj <- TAC/(Ctmp + 0.0001)

            Zprime <- M + Zadj * (Z - M)
            Alpha <- (1 - exp(-Zprime))/(Zprime)

            tmp <- sum(NAA * weight * sel * Alpha)
            Ftmp <-  TAC/(tmp + 0.0001)
            j2 <- 1/(1 + exp(30 * (Ftmp - 0.95 * fmax)))

            fout <- j2 * Ftmp + (1 - j2)
        }
    }
    return(fout)
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
