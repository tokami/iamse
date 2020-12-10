

#' @name genConvs
#' @description Get converged simulates from a resMSE object
#' @export
getConvs <- function(mse, convyears = "all", convhcrs = "all", out = 0, verbose = FALSE){

    nhcr <- length(mse)
    hcrs <- names(mse)
    nrep <- length(mse[[1]])
    nquant <- length(mse[[1]][[1]])
    nysim <- nrow(mse[[1]][[1]]$tacs)
    dims <- dim(mse[[1]][[1]]$CW)
    ny <- dims[1] - nysim
    ns <- dims[2]

    if(!is.na(convyears[1]) && convyears[1] == "all")
        convyears <- 1:nysim
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
            if(hcrs[hcr] != "refFmsy"){
            res[[hcr]] <- mse[[hcr]][inds]
            names(res[[hcr]]) <- inds
            }else{
                res[[hcr]] <- mse[[hcr]]
            names(res[[hcr]]) <- 1:nrep
            }
        }
    }else if(!is.na(convyears[1])){
        res <- vector("list",nhcr)
        for(hcr in 1:nhcr){
            tmpid <- unlist(lapply(strsplit(as.character(mse[[hcr]][[1]]$tacs$id),"-"), "[[", 1))
            if(any(tmpid %in% c("Bref","Bref2"))){
                id <- mse[[hcr]][[1]]$tacs$id[which(tmpid %in% c("Bref","Bref2"))[1]]
            }else{
                if(length(unique(mse[[hcr]][[1]]$tacs$id)) > 1){ ## TODO: find better solution to this
                    id <- unique(mse[[hcr]][[1]]$tacs$id)[which(!unique(mse[[hcr]][[1]]$tacs$id) %in% c("conscat","r23","r12","r35"))]
                }else{
                    id <- unique(mse[[hcr]][[1]]$tacs$id)
                }
            }
            if(id == 1) browser()
            print(id)
            id2 <- unlist(strsplit(as.character(id), "-"))[1]
            if(!(id2 %in% c("noF","refFmsy","r11","r12","r23","r35"))){
                tmp <- do.call(rbind,lapply(mse[[hcr]], function(x) x[["tacs"]]$conv))
                inds <- which(apply(tmp[,convyears], 1, all))
            }else inds <- 1:nrep
            if(verbose)
                writeLines(paste0("Converged reps: ",length(inds), " of ",nrep,
                                  " reps = ",round(length(inds)/nrep*100),"%"))
            res[[hcr]] <- mse[[hcr]][inds]
            names(res[[hcr]]) <- inds
        }
    }else{
        res <- vector("list",nhcr)
        for(hcr in 1:nhcr){
            res[[hcr]] <- mse[[hcr]]
            names(res[[hcr]]) <- 1:nrep
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
genNoise <- function(n, sd, rho=0, bias.cor = 0, mv=FALSE, dat=NULL){

    if(mv){
        ## multivariate noise
        stopifnot(!is.null(dat))
        amax <- dat$amax + 1
        Sigma <- matrix(NA, amax, amax)
        for(i in 1:amax) for(j in 1:amax) Sigma[i,j] = rho^abs(i - j) * sd^2
        res <- MASS::mvrnorm(n, rep(0,ncol(Sigma)), Sigma)
        if(bias.cor == 1){
            for(i in 1:amax) res[,i] <- res[,i] - Sigma[i,i] / 2
        }else if(bias.cor != 0) stop("bias.cor has to be 0 or 1 for multivariate noise. Please check set$noiseCmult and set$noiseImult.")
        res <- exp(res)
    }else{
        if(bias.cor == 0){
            rnoise <- rnorm(n, 0, sd)
        }else if(bias.cor == 1){
            rnoise <- rnorm(n, 0, sd) - sd^2/2
        }else if(bias.cor == 2){
            rnoise <- log(rlnorm(n, muconv(1,sd), sdconv(1,sd)))
        }else stop("bias.cor has to be 0, 1, or 2. Please check the different set$noise* settings.")

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

    if(depl.quant %in% c("Bmsy","Blim")){
        outopt <- 2
    }else if(depl.quant %in% c("SSBmsy","SSBlim")){
        outopt <- 3
    }else stop("depl.quant not implemented. Please use Bmsy, Blim,SSBmsy or SSBlim. Or implement others.")

    frel <- dat$FM/max(dat$FM)

    fn <- function(logfabs, frel, depl, depl.prob, nrep, dat, set, outopt, optFn=1){
        datx <- dat
        fpat <- frel * exp(logfabs)
        datx$FM <- fpat
        datx$Fs <- fpat / datx$ns
        dreal <- sapply(1:nrep, function(x) initPop(datx, set, out.opt = outopt))
        if(method == "mean"){
            drealQ <- mean(dreal)
        }else if(method == "median"){
            drealQ <- quantile(dreal, probs = 0.5)
        }else if(method == "percentile"){
            drealQ <- quantile(dreal, probs = depl.prob)
        }
        if(optFn==1) return((depl - drealQ)^2)
        if(optFn==2) return(drealQ)
    }

    opt <- optimize(fn, c(log(0.0001),log(fmax)), frel = frel, depl = depl, depl.prob = depl.prob,
                    nrep = nrep, dat = dat, set=set, outopt = outopt, optFn = 1)
    fabs <- exp(opt$minimum)
    fvals <- frel * fabs
    dreal <- round(fn(log(fabs), frel, depl, depl.prob, nrep, dat, set, outopt=outopt, optFn = 2),3)

    if(verbose){
        print(paste0("Target depletion level as ",depl.prob * 100, "% quantile of ", depl, " ", depl.quant,
                     " -- Realised depletion level at ", dreal, " ", depl.quant,
                     " with absolute F equal to ",round(fabs,3)))
    }

    dat$FM <- fvals
    dat$Fs <- fvals / dat$ns

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
    ns <- dat$ns
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
    ## CHECK: how to estimate productivity with time variant M?
    dat$M <- mean(dat$M)
    dat$Ms <- mean(dat$Ms)
    dat <- checkDat(dat)
    pop1 <- initPop(dat, set)
    tsb1 <- pop1$TSBfinal
    esb1 <- pop1$ESBfinal
    cw1 <- apply(pop1$CW,1,sum)
    prod1 <- rep(NA, ny)
    if(set$spType == 0){
        for(i in 2:ny){
            prod1[i] <- tsb1[i] - tsb1[i-1] + cw1[i]
        }
    }else if(set$spType == 1){
        for(i in 2:ny){
            prod1[i] <- esb1[i] - esb1[i-1] + cw1[i]
        }
    }

    ## est blim as fraction of B corresponding to 0.5 MSY (ICES WKBUT 2013, Cadrin 1999)
    msy1 <- max(prod1, na.rm=TRUE)
    Blim1 <- tsb1[which.min(abs(prod1 - msy1/2))]

    ## decreasing effort
    dat$FM <- c(rep(fmax, len1),
                   seq(fmax, 0, length.out = len2),
                rep(0, len3))
    dat$Fs <- dat$FM / ns
    dat$M <- mean(dat$M)
    dat$Ms <- mean(dat$Ms)
    dat <- checkDat(dat)
    pop2 <- initPop(dat, set)
    dat$Fs <- dat$FM / ns
    tsb2 <- pop2$TSBfinal
    esb2 <- pop2$ESBfinal
    cw2 <- apply(pop2$CW,1,sum)
    prod2 <- rep(NA, ny)
    if(set$spType == 0){
        for(i in 2:ny){
            prod2[i] <- tsb2[i] - tsb2[i-1] + cw2[i]
        }
    }else if(set$spType == 1){
        for(i in 2:ny){
            prod2[i] <- esb2[i] - esb2[i-1] + cw2[i]
        }
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
        Blim = Blim1,
        incr = data.frame(tsb = tsb1,
                          esb = esb1,
                          cw = cw1,
                          prod = prod1),
        decr = data.frame(tsb = tsb2,
                          esb = esb2,
                          cw = cw2,
                          prod = prod2)
    )

    return(res)

}



#' @name estProd
#' @export
estProdStoch <- function(dat, set= NULL,
                         fmax = 10,
                         nf = 1e3,
                         prob = c(0.1,0.9),
                         ncores = parallel::detectCores()-1,
                         plot = TRUE){

    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    ## noise
    if(is.null(set)) set <- checkSet()
    nyref <- set$refYears
    nrep <- set$refN
    nyrefmsy <- set$refYearsMSY

    set$noiseR <- c(dat$sigmaR, dat$rhoR, 1)
    dist <- NULL
    if(!(set$refMethod %in% c("mean","median"))){
        stop("'set$refMethod' not known! Has to be 'mean' or 'median'!")
    }

    ## errors (have to be re-used for estimation of Bmsy)
    errs <- vector("list", nrep)
    for(i in 1:nrep){
        errs[[i]] <- vector("list", 7)
        errs[[i]]$eF <- genNoise(nyref, set$noiseF[1], set$noiseF[2], set$noiseF[3])
        errs[[i]]$eR <- genNoise(nyref, set$noiseR[1], set$noiseR[2], set$noiseR[3])
        errs[[i]]$eM <- genNoise(nyref, set$noiseM[1], set$noiseM[2], set$noiseM[3])
        errs[[i]]$eH <- genNoise(nyref, set$noiseH[1], set$noiseH[2], set$noiseH[3])
        errs[[i]]$eR0 <- genNoise(nyref, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
        errs[[i]]$eMat <- genNoise(nyref, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
        errs[[i]]$eImp <- genNoise(nyref, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
    }

    ##
    fms <- seq(0, fmax, length.out = nf)
    resList <- vector("list", nf)
    for(fx in 1:nf){
        tmp0 <- parallel::mclapply(as.list(1:nrep), function(x){
            setx <- c(set, errs[[x]])
            pop <- simpop(log(fms[fx]), dat, setx, out=0)
            tsb <- tail(pop$TSB,1)
            esb <- tail(pop$ESB,1)
            ssb <- tail(pop$SSB,1)
            cw <- tail(pop$CW,1)
            sp <- tail(pop$SP,1)
            return(c(TSB = tsb, SSB = ssb, ESB = esb, CW = cw, SP = sp))
        }, mc.cores = ncores)
        tmp <- do.call(rbind, tmp0)
        resList[[fx]] <- cbind(f = rep(fms[fx],nrep), tmp)
    }

    ## est blim as fraction of B corresponding to 0.5 MSY (ICES WKBUT 2013, Cadrin 1999)
    bs <- do.call(rbind, lapply(resList, function(x) x[,2]))
    sps <- do.call(rbind, lapply(resList, function(x) x[,6]))
    blims <- rep(NA, nrep)
    for(i in 1:nrep){
        msy <- max(sps[,i], na.rm=TRUE)
        blims[i] <- bs[which.min(abs(sps[,i] - msy/2)),i]
    }

    means <- as.data.frame(do.call(rbind,lapply(resList,
                                                function(x) apply(x,2, mean, na.rm=TRUE))))
    meds <- as.data.frame(do.call(rbind,lapply(resList,
                                               function(x) apply(x,2, median, na.rm=TRUE))))
    lo <- as.data.frame(do.call(rbind,lapply(resList,
                                             function(x) apply(x,2, quantile, prob=min(prob),
                                                                        na.rm=TRUE))))
    up <- as.data.frame(do.call(rbind,lapply(resList,
                                             function(x) apply(x,2, quantile, prob=max(prob),
                                                          na.rm=TRUE))))

    res <- list(meds = meds,
                means = means,
                lo = lo,
                up = up,
                blims = blims)
    return(res)

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
    tacEst <- function(logFM, NAA, M, sel, weight, TAC){
##        (TAC - sum(baranov(exp(logFM) * sel, M, NAA*exp(-M/2)) * weight))^2
        (TAC - sum(baranov(exp(logFM) * sel, M, NAA) * weight))^2
    }
    opt <- optimise(tacEst, c(-10,10), NAA = NAA, M = M, TAC = TAC,
                    weight = weight, sel = sel)
    return(exp(opt$minimum))
}


#' @name getFM3
#' @details get FM accounting for seasons
#' @export
getFM3 <- function(TAC, NAA, MAA, weights, sels, ns, spawning,
                   lastFM = 0.1){
    tacEst <- function(logFM, NAA, MAA, sels, weights, TAC, ns){
        Ctmp <- 0
        NAAtmp <- NAA
        for(s in 1:ns){
            ## recruitment
            if(spawning[s] > 0){
                ## Survivors from previous season/year
                FAA <- exp(logFM) * sels[,s]
                Ztemp <- ZAA[,s-1]

                ## HERE: incorporate spawning in this funtion
                ## TODO: allow to have any TAC intervals here (starting in middle of year being 2 years long)
                ##

                SSB[y,1] <- sum(NAA * weights[,1] * maty[,1] * exp(-pzbm * Ztemp)) ## pre-recruitment mort
                SSBPR0 <- getSSBPR2(MAA, maty, weights, fecun=1, amax, R0y, ns = ns,
                                    season = s)
                rec[y] <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSB[y,1],
                                  R0 = R0y, method = dat$SR, bp = dat$bp,
                                  beta = dat$recBeta, gamma = dat$recGamma)
                rec[y] <- ifelse(rec[y] < 0, 1e-10, rec[y])
                NAA[1] <- NAA[1] + spawning[s] * rec[y] * eR[y]
            }
            Ctmp <- Ctmp + sum(baranov(FAA, MAA[,s], NAAtmp) * weights[,s])
            NAAtmp <- NAAtmp * exp(-(MAA[,s] + FAA))
        }
        (TAC - Ctmp)^2
    }
    opt <- nlminb(log(lastFM), tacEst, lower = -10, upper = 10, NAA = NAA, MAA = MAA, TAC = TAC,
                  weights = weights, sels = sels, ns = ns, control = list(rel.tol = 1e-15))
    return(exp(opt$par))
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


#' @name getM
#' @description Function to estimate selectivity of natural mortality
#' @param linf - Linf of vBGF
#' @param k - K of vBGF
#' @param mids - midlengths
getM <- function(linf, k, mids, a = 0.55, b = 1.61, c = 1.44){
    n <- max(c(length(a),length(b),length(c)))
    maxM <- rep(NA, n)
    for(i in 1:n){
        selL <- exp(a[i] - b[i] * log(mids) + c[i] * log(linf) + log(k))
        selL[mids < 10] <- exp(a[i] - b[i] * log(10) + c[i] * log(linf) + log(k))
        selL <- round(selL, 3)
        maxM[i] <- max(selL)
    }
    return(maxM)
}


#' @name getMsel
#' @description Function to estimate selectivity of natural mortality
#' @param linf - Linf of vBGF
#' @param k - K of vBGF
#' @param mids - midlengths
#' @param plba - probability of being in mids given age
getMsel <- function(linf, k, mids, plba, a = 0.55, b = 1.61, c = 1.44){
    n <- max(c(length(a),length(b),length(c)))
    sels <- vector("list", n)
    for(i in 1:n){
        selL <- exp(a[i] - b[i] * log(mids) + c[i] * log(linf) + log(k))
        selL[mids < 10] <- exp(a[i] - b[i] * log(10) + c[i] * log(linf) + log(k))
        dims <- dim(plba)
        selA <- matrix(NA, ncol = dims[3], nrow = dims[1])
        for(j in 1:dim(plba)[3]){
            selA[,j] <- apply(t(plba[,,j]) * selL, 2, sum)
        }
        maxM <- max(selA)
        sels[[i]] <- selA/maxM
    }
    return(sels)
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

#' @name getSSBPR2
#' @description Function to calculate spawners per recruit
#' @param Z - total mortality
#' @param mat - maturity ogive
#' @param fecun - fecundity matrix
#' @param amax - number of age classes
#' @return spawning biomass per recruit
#' @export
getSSBPR2 <- function (Ms, mats, weights, fecun = 1, amax, R0 = 1, FMs = NULL, ns, season){
    ## account for seasonal natural mortalities
    if(inherits(Ms, "matrix")){
        Mtot <- apply(Ms, 1, sum)
    }else{
        Mtot <- Ms
    }
    if(is.null(FMs)){
        FMtot <- 0
        FMs <- matrix(0, nrow = amax, ncol = ns)
    }else if(inherits(FMs, "matrix")){
        FMtot <- apply(FMs, 1, sum)
    }else{
        FMtot <- FMs
    }
    ## first quarter/whole year
    NAA <- matrix(NA, nrow = amax, ncol = ns)
    NAA[1] <- R0
    MAA <- cumsum(Mtot)
    FAA <- cumsum(FMtot)
    ZAA <- MAA + FAA
    NAA[2:(amax  - 1),1] <- R0 * exp(-ZAA[1:(amax - 2)])
    NAA[amax,1] <- R0 * exp(-ZAA[amax - 1])/(1 - exp(-ZAA[amax]))
    ## other quarters
    if(ns > 1){
        for(s in 2:ns){
            ZAA <- Ms[,s] + FMs[,s]
            NAA[,s] <- NAA[,s-1] * exp(-ZAA)
        }
    }
    SBPR_per_season <- apply(NAA * mats * weights * fecun, 2, sum)
    ## season before spawning? season_before_spawning <- c(ns,(1:(ns-1)))[season]
    SBPR <- SBPR_per_season[season]
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
        ## beta <- log(5 * h) / (0.8 * R0)
        ## alpha <- exp(beta * R0)/SSBPR0
        ## rec <- alpha * SSB * exp(-beta * SSB)
        rec <- bp * SSB * exp(-beta * SSB)
    }else if(method == "average"){
        rec <- R0
    }else if(method == "hockey-stick"){
        rec <- ifelse(SSB > bp, R0, SSB * R0/bp)
    }else if(method == "bent-hyperbola"){  ## Watts-Bacon bent hyperbola
        rec <- beta * (SSB + sqrt(bp^2 + (gamma^2)/4) - sqrt((SSB-bp)^2 + (gamma^2)/4))
    }else print("Stock-recruitment method not known! Implemented methods: 'bevholt', 'ricker', 'average', and 'hockey-stick'.")

    return (rec)
}
