

#' @name genConvs
#' @description Get converged simulates from a resMSE object
#' @export
get.converged <- function(mse, convyears = "all", convhcrs = "all", out = 0, verbose = FALSE){

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


#' @name gen.noise
#' @export
gen.noise <- function(n, sd, rho=0, bias.cor = 0, mv=FALSE, dat=NULL){

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


#' @name est.depletion
#' @export
est.depletion <- function(dat, set=NULL, fmin = 0.0001,
                          fmax = 10, nrep = 100, verbose = TRUE,
                          method = "percentile",
                          tol = 0.0001, do.opt = TRUE){

    if(!any(names(dat) == "ref")) stop("Reference points are missing in dat. Use est.ref.levels.stochastic to estimate reference points.")

    if(is.null(set)){
        set <- check.set()
        nrep <- 1
    }
    refs <- dat$ref
    ny <- dat$ny
    ns <- dat$ns

    depl <- dat$depl
    depl.quant <- dat$depl.quant
    depl.prob <- dat$depl.prob

    if(depl.quant %in% c("Bmsy","Blim")){
        outopt <- 2
        blim <- dat$ref$Blim[ny]
    }else if(depl.quant %in% c("SSBmsy","SSBlim")){
        outopt <- 3
        blim <- dat$ref$SSBlim[ny]
    }else stop("depl.quant not implemented. Please use Bmsy, Blim,SSBmsy or SSBlim. Or implement others.")


    ## errors
    errs <- list()
    errs$eF <- lapply(as.list(1:nrep), function(x) gen.noise(ny, set$noiseF[1], set$noiseF[2], bias.cor = set$noiseF[3]))
    errs$eR <- lapply(as.list(1:nrep), function(x) gen.noise(ny, set$noiseR[1], set$noiseR[2], bias.cor = set$noiseR[3]))
    errs$eM <- lapply(as.list(1:nrep), function(x) gen.noise(ny, set$noiseM[1], set$noiseM[2], bias.cor = set$noiseM[3]))
    errs$eH <- lapply(as.list(1:nrep), function(x) gen.noise(ny, set$noiseH[1], set$noiseH[2], bias.cor = set$noiseH[3]))
    errs$eR0 <- lapply(as.list(1:nrep), function(x) gen.noise(ny, set$noiseR0[1], set$noiseR0[2], bias.cor = set$noiseR0[3]))
    errs$eMat <- lapply(as.list(1:nrep), function(x) gen.noise(ny, set$noiseMat[1], set$noiseMat[2], bias.cor = set$noiseMat[3]))
    errs$eSel <- lapply(as.list(1:nrep), function(x) gen.noise(ny, set$noiseSel[1], set$noiseSel[2], bias.cor = set$noiseSel[3]))
    errs$eW <- lapply(as.list(1:nrep), function(x) gen.noise(ny, set$noiseW[1], set$noiseW[2], bias.cor = set$noiseW[3]))

    frel <- dat$FM/max(dat$FM)

    fn <- function(logfabs, frel, depl, depl.prob, nrep, dat, set, errs, outopt, optFn=1){
        datx <- dat
        setx <- set
        fpat <- frel * exp(logfabs) / dat$ns
        datx$FM <- fpat
        dreal <- sapply(1:nrep, function(x){
            setx$eF <- errs$eF[[x]]
            setx$eR <- errs$eR[[x]]
            setx$eM <- errs$eM[[x]]
            setx$eH <- errs$eH[[x]]
            setx$eR0 <- errs$eR0[[x]]
            setx$eMat <- errs$eMat[[x]]
            setx$eSel <- errs$eSel[[x]]
            setx$eW <- errs$eW[[x]]
            initpop(datx, setx, out.opt = outopt)
        })
        if(method == "mean"){
            drealQ <- mean(dreal)
        }else if(method == "median"){
            drealQ <- quantile(dreal, probs = 0.5)
        }else if(method == "percentile"){
            drealQ <- quantile(dreal, probs = depl.prob)
        }
        if(optFn==1) return((drealQ - depl)^2)
        if(optFn==2) return(drealQ)
        if(optFn==3) return(dreal)
    }

    ## opt <- nlminb(log(5), fn, lower = log(fmin), upper = log(fmax), frel = frel, depl = depl,
    ##               depl.prob = depl.prob,
    ##               nrep = nrep, dat = dat, set=set, errs=errs, outopt = outopt, optFn = 1)
    if(do.opt){
        opt <- optimize(fn, c(log(fmin),log(fmax)), frel = frel, depl = depl, depl.prob = depl.prob,
                        nrep = nrep, dat = dat, set=set, errs=errs, outopt = outopt, optFn = 1, tol = tol)
        fabs <- exp(opt$minimum)
    }else{
        fabs <- max(apply(dat$FM,1,sum))
    }
    fvals <- frel * fabs / ns
    tmp <- fn(log(fabs), frel, depl, depl.prob, nrep, dat, set, errs, outopt=outopt, optFn = 2)
    dreal <- round(tmp,5)
    tmp <- fn(log(fabs), frel, depl, depl.prob, nrep, dat, set, errs, outopt=outopt, optFn = 3)
    risk <- round(mean(tmp * dat$ref[[depl.quant]][ny] < blim) * 100,1)
    fmsyfac <- round(fabs / dat$ref$Fmsy[ny],2)

    if(verbose){
        print(paste0("Target depletion level as ",depl.prob, "% quantile of ", depl, " ", depl.quant,
                     " -- Realised depletion level at ", dreal, " ", depl.quant,
                     " with absolute F equal to ",round(fabs,3)," (",fmsyfac," * Fmsy). This corresponds to a risk of ",
                     risk,"% (P(B[last] < B[lim])."))
    }

    dat$FM <- fvals

    return(dat)
}



#' @name est.productivity
#' @export
est.productivity <- function(dat, set= NULL,
                    ny = 100,
                    fmax = 10,
                    tsSplit = 8,
                    plot = TRUE){

    dat$ny <- ny
    ns <- dat$ns
    nt <- ny * ns

    ## noise
    if(is.null(set)) set <- check.set()
    set$noiseF <- c(0,0,0)
    set$noiseR <- c(0,0,0)
    set$noiseR0 <- c(0,0,0)
    set$noiseH <- c(0,0,0)
    set$noiseM <- c(0,0,0)
    set$noiseW <- c(0,0,0)
    set$noiseMat <- c(0,0,0)
    set$noiseSel <- c(0,0,0)
    set$noiseImp <- c(0,0,0)


    ##
    len1 <- len3 <- floor(ny/tsSplit)
    len2 <- ny - len1 - len3

    ## increasing effort
    dat$FM <- matrix(c(rep(0,len1),
                   seq(0, fmax, length.out = len2),
                rep(fmax,len3)) / ns, ncol=ns, nrow = ny)
    ## CHECK: how to estimate productivity with time variant M?
    dat$M <- matrix(dat$M[1,], ncol=ns, nrow=1)
    dat <- check.dat(dat)
    pop1 <- initpop(dat, set)
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
    dat$FM <- matrix(c(rep(fmax, len1),
                   seq(fmax, 0, length.out = len2),
                rep(0, len3)) / ns, ncol=ns, nrow=ny)
    dat$M <- matrix(dat$M[1,], ncol=ns, nrow=1)
    dat <- check.dat(dat)
    pop2 <- initpop(dat, set)
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



#' @name est.productivity
#' @export
est.productivity.stochastic <- function(dat, set= NULL,
                         fmax = 10,
                         nf = 1e3,
                         prob = c(0.1,0.9),
                         ncores = parallel::detectCores()-1,
                         plot = TRUE){

    amax <- dat$amax + 1
    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    asmax <- amax * ns
    ## noise
    if(is.null(set)) set <- check.set()
    nyref <- set$refYears
    nrep <- set$refN
    nyrefmsy <- set$refYearsMSY

    ##CHECK: set$noiseR <- c(dat$sigmaR, dat$rhoR, 1)
    dist <- NULL
    if(!(set$refMethod %in% c("mean","median"))){
        stop("'set$refMethod' not known! Has to be 'mean' or 'median'!")
    }

    ## errors (have to be re-used for estimation of Bmsy)
    errs <- vector("list", nrep)
    for(i in 1:nrep){
        errs[[i]] <- vector("list", 7)
        errs[[i]]$eF <- gen.noise(nyref, set$noiseF[1], set$noiseF[2], set$noiseF[3])
        errs[[i]]$eR <- gen.noise(nyref, set$noiseR[1], set$noiseR[2], set$noiseR[3])
        errs[[i]]$eM <- gen.noise(nyref, set$noiseM[1], set$noiseM[2], set$noiseM[3])
        errs[[i]]$eH <- gen.noise(nyref, set$noiseH[1], set$noiseH[2], set$noiseH[3])
        errs[[i]]$eW <- gen.noise(nyref, set$noiseW[1], set$noiseW[2], set$noiseW[3])
        errs[[i]]$eR0 <- gen.noise(nyref, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
        errs[[i]]$eMat <- gen.noise(nyref, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
        errs[[i]]$eSel <- gen.noise(nyref, set$noiseSel[1], set$noiseSel[2], set$noiseSel[3])
        errs[[i]]$eImp <- gen.noise(nyref, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
    }

    datx <- dat
    ## ##
    ## datx$yvec <- rep(1:nyref, each = ns)
    ## datx$svec <- rep(1:ns, each = nyref)
    ## datx$s1vec <- seq(1, nyref * ns, ns)
    ## datx$as2a <- rep(1:amax, each = ns)
    ## datx$as2s <- rep(1:ns, amax)
    ## datx$inds <- seq(1,asmax,ns)

    ## natural mortality
    ms <- NULL
    for(i in 1:ns){
        tmp0 <- unique(dat$M[,i])
        if(is.null(ms) || length(tmp0) == nrow(ms)){
            ms <- cbind(ms,tmp0)
        }else if(length(tmp0) == 1){
            ms <- cbind(ms,rep(tmp0, length.out = ns))
        }else stop("You are natural mortality (M) is time-variant but M does not vary consitently among seasons. Please review dat$M or contact the package maintainer.")
    }
    mtv <- nrow(ms)
    mind <- match(dat$M[,1], ms[,1])
    ## M selectivity
    if(length(dat$Msel) > 1){
        msels <- dat$Msel[!duplicated(dat$Msel)]
        mseltv <- length(msel)
    }else{
        msels <- dat$Msel[1]
        mseltv <- 1
    }
    if(mseltv > 1 && mseltv != mtv) stop("Both natural mortality over time (dat$M) and over age (dat$Msel) are time-variant, but do not have the same dimensions. This is not yet implemented, please let both vary equally or keep one of them constant.")
    alltv <- max(c(mtv, mseltv))
    ## selectivity
    if(length(dat$sel) > 1){
        sels <- dat$sel[!duplicated(dat$sel)]
        seltv <- length(sel)
    }else{
        sels <- dat$sel[1]
        seltv <- 1
    }
    if(seltv > 1 && alltv > 1 && seltv != alltv) stop("Both gear selectivity (dat$sel) and natural mortality (dat$M or dat$Msel) are time-variant, but do not have the same dimensions. This is not yet implemented, please let both vary equally or keep one of them constant.")
    alltv <- max(c(alltv,seltv))

    if(alltv > 1){
        if(mtv == alltv){
            mtv <- 1:mtv
        }else mtv <- rep(mtv, length.out = alltv)
        if(mseltv == alltv){
            mseltv <- 1:mseltv
        }else mseltv <- rep(mseltv, length.out = alltv)
        if(seltv == alltv){
            seltv <- 1:seltv
        }else seltv <- rep(seltv, length.out = alltv)
    }


    ##
    blims <- vector("list",alltv)
    means <- vector("list",alltv)
    meds <- vector("list",alltv)
    lo <- vector("list",alltv)
    up <- vector("list",alltv)
    for(i in 1:alltv){
        datx$M <- t(as.matrix(ms[mtv[i],]))
        datx$Msel <- msels[mseltv[i]]
        datx$sel <- sels[seltv[i]]
        fms <- seq(0, fmax, length.out = nf)
        tmp2 <- vector("list", nf)
        for(fx in 1:nf){
            tmp0 <- mclapply.all.os(as.list(1:nrep), function(x){
                setx <- c(set, errs[[x]])
                pop <- simpop(log(fms[fx]), datx, setx, out=0)
                tsb <- tail(pop$TSB,1)
                esb <- tail(pop$ESB,1)
                ssb <- tail(pop$SSB,1)
                cw <- tail(pop$CW,1)
                sp <- tail(pop$SP,1)
                return(c(TSB = tsb, SSB = ssb, ESB = esb, CW = cw, SP = sp))
            }, mc.cores = ncores)
            tmp1 <- do.call(rbind, tmp0)
            tmp2[[fx]] <- cbind(f = rep(fms[fx],nrep), tmp1)
        }

        ## est blim as fraction of B corresponding to 0.5 MSY (ICES WKBUT 2013, Cadrin 1999)
        bs <- do.call(rbind, lapply(tmp2, function(x) x[,2]))
        sps <- do.call(rbind, lapply(tmp2, function(x) x[,6]))
        blims[[i]] <- rep(NA, nrep)
        for(j in 1:nrep){
            msy <- max(sps[,j], na.rm=TRUE)
            blims[[i]][j] <- bs[which.min(abs(sps[,j] - msy/2)),i]
        }

        means[[i]] <- as.data.frame(do.call(rbind,lapply(tmp2,
                                                    function(x) apply(x,2, mean, na.rm=TRUE))))
        meds[[i]] <- as.data.frame(do.call(rbind,lapply(tmp2,
                                                   function(x) apply(x,2, median, na.rm=TRUE))))
        lo[[i]] <- as.data.frame(do.call(rbind,lapply(tmp2,
                                                 function(x) apply(x,2, quantile, prob=min(prob),
                                                                   na.rm=TRUE))))
        up[[i]] <- as.data.frame(do.call(rbind,lapply(tmp2,
                                                 function(x) apply(x,2, quantile, prob=max(prob),
                                                                   na.rm=TRUE))))
    }


    if(plot){
        cols <- rep(c("darkred","dodgerblue","darkgreen","darkorange","purple","gray","black","goldenrod"),100)
        plot(meds[[1]]$TSB, meds[[1]]$SP, ty = 'n',
             ylim = range(0,sapply(meds,function(x) x$SP),
                          sapply(lo,function(x) x$SP),
                          sapply(up,function(x) x$SP), na.rm =TRUE),
             xlim = range(0,sapply(meds,function(x) x$TSB),
                          sapply(lo,function(x) x$TSB),
                          sapply(up,function(x) x$TSB), na.rm =TRUE),
             xlab = "TSB", ylab = "SP")
        for(i in 1:alltv){
            if(i <= 3){
                polygon(c(lo[[i]]$TSB, rev(up[[i]]$TSB)), c(lo[[i]]$SP, rev(up[[i]]$SP)), border = NA,
                        col = rgb(t(col2rgb(cols[i])/255), alpha = 0.2))
            }
            lines(meds[[i]]$TSB, meds[[i]]$SP, col=cols[i], lwd=2)
        }
        legend("topright",
               legend = c("M = 0.2", "M = 0.3"),  ## CHECK: adjust
               col = cols[1:alltv],
               lwd = 1.5)
    }

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
baranov <- function(F, M, N){
    Z <- F + M
    return(F/Z * N * (1 - exp(-Z)))
}


#' @name predCatch
#'
#' @param seasons vector with season indices
#' @param ns number of seasons
#' @param h steepness
#'
#' @details get predicted catch for TAC period or difference between provided
#'     and predicted catch
predCatch <- function(logFM,
                      NAA, MAAy,
                      sel, weight,
                      seasons, ns, y, h2, asmax, mat, pzbm, spawning,
                      R0, SR, bp, recBeta, recGamma, eR,
                      indage0,
                      iaFM,
                      TAC = NULL,
                      out = 0){
    Ctmp <- 0
    NAAtmp <- NAA

    for(i in 1:length(seasons)){
        s <- seasons[i]
        FAA <- exp(logFM) * iaFM[s] * sel
        MAA <- MAAy[[s]] ## MAAy[,s]
        Ztmp <- FAA + MAA
        ## recruitment
        if(spawning[s] > 0 && i > 1){
            SSB0 <- get.ssb0(MAA, mat, weight,
                               fecun=1, asmax, ns, spawning,
                             R0, indage0, s)
            ## Survivors from previous season/year
            SSBtmp <- sum(NAAtmp * weight  * mat  * exp(-pzbm * Ztmp))
            rec <- spawning[s] * recfunc(h = h2, SPR0 = SSB0/R0, SSB = SSBtmp,
                             R0 = R0, method = SR, bp = bp,
                             beta = recBeta, gamma = recGamma) * eR
            rec[rec<0] <- 1e-10
            NAAtmp[indage0] <- rec
        }
        Ctmp <- Ctmp + sum(baranov(FAA, MAA, NAAtmp) * weight)
        ## print(paste0("TSB: ",round(sum(NAAtmp*weight))))
        ## print(paste0("Ctmp: ",round(Ctmp)))
        ## Aging
        NAAtmp <- NAAtmp * exp(-Ztmp)
        NAAtmp[asmax] <- NAAtmp[asmax] + NAAtmp[asmax-1]
        for(as in (asmax-1):2) NAAtmp[as] <- NAAtmp[as-1]
        NAAtmp[indage0] <- 0
    }

    if(out == 0){
        return(Ctmp)
    }else{
##         return((Ctmp - TAC)^2)  ## CHECK: needs if clause if Ctmp or TAC == 0?
##         return((log(Ctmp) - log(TAC))^2)  ## CHECK: needs if clause if Ctmp or TAC == 0?
         return(sqrt(mean(Ctmp - TAC)^2))
    }
}

#' @name get.f
#' @details get FM accounting for seasons
#' @export
get.f <- function(TAC,
                   NAA, MAA,
                   sel, weight,
                   seasons, ns, y, h, asmax, mat,
                   pzbm, spawning,
                   R0, SR, bp, recBeta,
                  recGamma, eR,
                  indage0,
                  iaFM,
                  lastFM = 0.01, upper = 100,
                  tac_cut_off = 0.01){

    ## browser()

    ## TAC
    ## iaFM
    ## predCatch(log(0.48), NAA, MAA, sel, weight,
    ##           seasons, ns, y, h, asmax, mat, pzbm, spawning,
    ##           R0, SR, bp, recBeta, recGamma, eR,
    ##           indage0,
    ##           iaFM,
    ##           TAC = TAC,
    ##           out = 0)

    if(TAC < tac_cut_off){
        return(0)
    }else{
        opt <- nlminb(start = log(lastFM), objective = predCatch,
                      NAA = NAA, MAAy = MAA,
                      sel = sel, weight = weight,
                      seasons = seasons, ns = ns, y = y,
                      h2 = h, asmax = asmax, mat = mat,
                      pzbm = pzbm, spawning = spawning,
                      R0 = R0, SR = SR, bp = bp, recBeta = recBeta,
                      recGamma = recGamma, eR = eR,
                      indage0 = indage0,
                      iaFM = iaFM,
                      TAC = TAC,
                      out = 1,
                      lower = -10, upper = log(upper),
                      control = list(rel.tol = 1e-3,
                                     iter.max = 1e4,
                                     eval.max = 1e4))

##    print(paste0("obj: ",round(opt$objective,2), "- fm: ",round(exp(opt$par),3)))
        return(exp(opt$par))
    }
}


#' @name getSel
#' @description Function to estimate selectivity ogive
#' @param pars - parameters, either Ls50 (length at 50% selectivity) and Ls95
#'     (length at 95% selectivity) for logistic or LFS, sdl, sdr
#' @param mids - midlengths
#' @param plba - probability of being in mids given age
#' @param type - logistic or dnormal
#'
getSel <- function(pars, mids, plba, type = "logistic", age = FALSE, ages = NULL){
    if(!age){
        if(type == "logistic"){
            Ls50 <- pars[["Ls50"]]
            Ls95 <- pars[["Ls95"]]
            stopifnot(all(is.numeric(c(Ls50,Ls95))))
            n <- max(c(length(Ls50),length(Ls95)))
            sel <- vector("list", n)
            for(i in 1:n){
                selL <- (1 /(1 + exp(-log(19)*(mids - Ls50[i])/(Ls95[i] - Ls50[i]))))
                dims <- dim(plba)
                selA <- matrix(NA, ncol = dims[3], nrow = dims[1])
                for(j in 1:dim(plba)[3]){
                    selA[,j] <- apply(t(plba[,,j]) * selL, 2, sum)
                }
                ##    selA <- apply(t(plba) * selL, 2, sum)
                ##    selA[1] <- 1e-9 # it should be zero for age 0
                sel[[i]] <- selA
            }
        }else if(type == "dnormal"){
            LFS <- pars[["LFS"]]
            sl <- pars[["sl"]]
            sr <- pars[["sr"]]
            stopifnot(all(is.numeric(c(LFS,sl,sr))))
            n <- max(c(length(LFS),length(sl),length(sr)))
            sel <- vector("list", n)
            for(i in 1:n){
                ind <- mids <= LFS
                selL <- rep(NA,length(mids))
                selL[ind] <- 2^(-((mids[ind] - LFS[i])/sl[i] *
                                  (mids[ind] - LFS[i])/sl[i]))
                selL[!ind] <- 2^(-((mids[!ind] - LFS[i])/sr[i] *
                                   (mids[!ind] - LFS[i])/sr[i]))
                dims <- dim(plba)
                selA <- matrix(NA, ncol = dims[3], nrow = dims[1])
                for(j in 1:dim(plba)[3]){
                    selA[,j] <- apply(t(plba[,,j]) * selL, 2, sum)
                }
                ##    selA <- apply(t(plba) * selL, 2, sum)
                ##    selA[1] <- 1e-9 # it should be zero for age 0
                sel[[i]] <- selA
            }

        }else{
            stop("Only 'logistic' and 'dnormal' implemented.")
        }
    }else{
        mids <- ages
        if(is.null(ages)) stop("No age classes provided! Use argument 'ages' for selection by age.")
        if(type == "logistic"){
            Ls50 <- pars[["Ls50"]]
            Ls95 <- pars[["Ls95"]]
            stopifnot(all(is.numeric(c(Ls50,Ls95))))
            n <- max(c(length(Ls50),length(Ls95)))
            sel <- vector("list", n)
            for(i in 1:n){
                sel[[i]] <- matrix((1 /(1 + exp(-log(19)*(mids - Ls50[i])/(Ls95[i] - Ls50[i])))),
                                   ncol = ncol(ages), nrow = nrow(ages))
            }
        }else if(type == "dnormal"){
            LFS <- pars[["LFS"]]
            sl <- pars[["sl"]]
            sr <- pars[["sr"]]
            stopifnot(all(is.numeric(c(LFS,sl,sr))))
            n <- max(c(length(LFS),length(sl),length(sr)))
            sel <- vector("list", n)
            for(i in 1:n){
                ind <- mids <= LFS
                selL <- rep(NA,length(mids))
                selL[ind] <- 2^(-((mids[ind] - LFS[i])/sl[i] *
                                  (mids[ind] - LFS[i])/sl[i]))
                selL[!ind] <- 2^(-((mids[!ind] - LFS[i])/sr[i] *
                                   (mids[!ind] - LFS[i])/sr[i]))
                sel[[i]] <- matrix(selL, ncol = ncol(ages), nrow = nrow(ages))
            }

        }else{
            stop("Only 'logistic' and 'dnormal' implemented.")
        }
    }
    return(sel)
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
#' @param Linf - Linf of vBGF
#' @param K - K of vBGF
#' @param mids - midlengths
getM <- function(Linf, K, mids, a = 0.55, b = 1.61, c = 1.44){
    n <- max(c(length(a),length(b),length(c)))
    maxM <- rep(NA, n)
    for(i in 1:n){
        selL <- exp(a[i] - b[i] * log(mids) + c[i] * log(Linf) + log(K))
        selL[mids < 10] <- exp(a[i] - b[i] * log(10) + c[i] * log(Linf) + log(K))
        selL <- round(selL, 3)
        maxM[i] <- max(selL)
    }
    return(maxM)
}


#' @name getMsel
#' @description Function to estimate selectivity of natural mortality
#' @param Linf - Linf of vBGF
#' @param K - K of vBGF
#' @param mids - midlengths
#' @param plba - probability of being in mids given age
getMsel <- function(Linf, K, mids, plba, a = 0.55, b = 1.61, c = 1.44, below10 = FALSE){
    if(is.null(below10)) below10 <- FALSE
    n <- max(c(length(a),length(b),length(c)))
    sel <- vector("list", n)
    for(i in 1:n){
        selL <- exp(a[i] - b[i] * log(mids) + c[i] * log(Linf) + log(K))
        if(!below10) selL[mids < 10] <- exp(a[i] - b[i] * log(10) + c[i] * log(Linf) + log(K))
        dims <- dim(plba)
        selA <- matrix(NA, ncol = dims[3], nrow = dims[1])
        for(j in 1:dim(plba)[3]){
            selA[,j] <- apply(t(plba[,,j]) * selL, 2, sum)
        }
        maxM <- max(selA)
        sel[[i]] <- selA/maxM
    }
    return(sel)
}




#' @name get.ssb0
#' @description Function to calculate SSB (F=0)
#' @param Z - total mortality
#' @param mat - maturity ogive
#' @param fecun - fecundity matrix
#' @param amax - number of age classes
#' @return spawning biomass per recruit
#' @export
get.ssb0 <- function (M, mat, weight, fecun = 1,
                       asmax, ns, spawning,
                      R0, indage0, season, FM=NULL){

    if(is.null(FM)){
        FM <- 0
    }
    ZAA <-  M + FM

    NAAS <- initdistR(M, FM=FM, ns, asmax, indage0, spawning, R0)
##    print(NAAS)

    ## SSB0 season dependent
    while(season > 1){
        NAAS <- NAAS * exp(-ZAA)
        NAAS[asmax] <- NAAS[asmax] + NAAS[asmax-1]
        for(as in (asmax-1):2) NAAS[as] <- NAAS[as-1]
        ## NAAS[1] <- 0  ## CHECK: ?
        NAAS[indage0] <- 0
        season <- season - 1
    }

    ## SSB0
    SBB0 <- sum(NAAS * mat * weight * fecun)

    return(SBB0)
}


#' @name recfunc
#' @description Function to calculate recruitment (Beverton - Holt)
#' @param h - steepness
#' @param R0 - recruitment in unfished population
#' @param SPR0 - spawning biomass produced by one recrut in its lifetime
#' @param SSB - spawning biomass
#' @param bp - breakpoint for hockey-stick SR
#' @param method - SR type
#'
#' @export
recfunc <- function(h, SPR0, SSB,  R0 = 1e6, method = "bevholt", bp = 0,
                    beta = 0, gamma = 0){

    if(method == "bevholt"){
        alpha <- SPR0 * (1-h)/(4*h)
        beta <- (5*h-1) / (4*h*R0)
        rec <- SSB / (alpha + beta * SSB)
    }else if(method == "ricker"){
        ## beta <- log(5 * h) / (0.8 * R0)
        ## alpha <- exp(beta * R0)/SPR0
        ## rec <- alpha * SSB * exp(-beta * SSB)
        rec <- bp * SSB * exp(-beta * SSB)
    }else if(method == "average"){
        rec <- rep(R0, length(SSB))
    }else if(method == "hockey-stick" || method == "segreg"){
        rec <- ifelse(SSB > bp, R0, SSB * R0/bp)
    }else if(method == "bent-hyperbola"){  ## Watts-Bacon bent hyperbola
        rec <- beta * (SSB + sqrt(bp^2 + (gamma^2)/4) - sqrt((SSB-bp)^2 + (gamma^2)/4))
    }else print("Stock-recruitment method not known! Implemented methods: 'bevholt', 'ricker', 'average', and 'hockey-stick' ('segreg').")

    return (rec)
}



#' @name initdistR
#' @export
initdistR <- function(M, FM=NULL, ns, asmax, indage0, spawning, R0=1){

    ## TODO: include here argument to get initdistR for specific season?

    if(is.null(FM)){
        FM <- 0
    }

    NAA2 <- NAA <- matrix(0, asmax, ns)
    NAA[indage0,] <- R0 * spawning
    ZAA <-  M + FM
    ## each season
    for(as in (indage0+1):asmax)
        NAA[as,] <- NAA[as-1,] * exp(-ZAA[as-1])
    ## only keep age groups present relative to end of year (last season)
    indx <- rep(NA,ns)
    for(s in 1:ns){
        indi <- seq(ns+2-s+indage0-1,asmax,ns)
        NAA2[indi,s] <- NAA[indi,s]
        if(asmax %in% indi) indx[s] <- 0
    }
    ##
    indi <- which(indx == 0)
    if(indi != 1)
        indx[1:(indi-1)] <- (indi-1):1
    if(indi != ns)
        indx[(indi+1):ns] <- (ns:(indi+1)) - 1

    ## keep last age group for every season
    indi <- which(NAA2[asmax,]==0)
    ##
    tmp <- NAA[asmax,] * exp(-indx*ZAA[asmax])
    NAA2[asmax,indi] <- tmp[indi]
    ## plus group correction
    NAA2[asmax,] <- NAA2[asmax,] / (1 - exp(-sum(ns * ZAA[asmax])))
    ## combine seasons
    NAAS <- rowSums(NAA2)
    ## remove recruits
    NAAS[indage0] <- 0

    return(NAAS)
}


#' @name remove.noise
#' @export
remove.noise <- function(set){

    ind <- grep("noise",names(set))
    n <- length(ind)
    for(i in 1:n){
        set[[ind[i]]] <- c(0,0,0)
    }

    return(set)
}



#' @name get.annual
#' @export
get.annual <- function (intime, vec, type = "mean"){
    anntime <- intime[which(intime%%1 == 0)]
    nanntime <- length(anntime)
    nstepvec <- rep(0, nanntime)
    floortime <- floor(intime)
    for (i in 1:nanntime) {
        nstepvec[i] <- sum(anntime[i] == floortime)
    }
    nsteps <- max(nstepvec)
    anntime <- anntime[which(nstepvec == max(nstepvec))]
    nanntime <- length(anntime)
    annvec <- rep(0, nanntime)
    for (i in 1:nanntime) {
        inds <- which(anntime[i] == floortime)
        if (is(type, "function")) {
            annvec[i] <- type(vec[inds])
        }
        else {
            if (type == "mean") {
                annvec[i] <- mean(vec[inds])
            }
            if (type == "sum") {
                annvec[i] <- sum(vec[inds])
            }
        }
    }
    return(list(anntime = anntime, annvec = annvec))
}




#' @name repar.sr
#' @description Reparameterise stock-recruitment relationship
#' @param dat
#'
#' @details Convert steepness and virbin biomass into other parameterisations.
#'     So far implemented: h and vb into breakpoint (bp) and R0 for segmented
#'     regression (hockey-stick)
#'
#'
#' @return spawning biomass per recruit
#' @export
repar.sr <- function (dat){

    if(is.null(dat$h)){
        stop(paste0("Reparameterisation of stock-recruitment parameters requires ",
                    "information about steepness (dat$h)."))
    }else{
        h <- dat$h
    }
    if(is.null(dat$vb)){
        stop(paste0("Reparameterisation of stock-recruitment parameters requires ",
                    "information about virgin biomass (dat$vb)."))
    }else{
        vb <- dat$vb
    }
    M <- dat$M[1,] * as.numeric(t(dat$Msel[[1]]))  ## CHECK: TVM: which M if time-variant?
    sel <-  as.numeric(t(dat$sel[[1]]))
    mat <- as.numeric(t(dat$mat))
    weight <- as.numeric(t(dat$weight))

    ## Estimate SPR0
    SPR0 <- get.ssb0(M = M, mat = mat, weight = weight, fecun = 1,
                     asmax = dat$asmax, ns = dat$ns, spawning = dat$spawning,
                     R0 = 1, indage0 = dat$indage0,
                     season = dat$ns, ## CHECK: which season?
                     FM=NULL)

    if(dat$SR == "hockey-stick" || dat$SR == "segreg"){
        a <- 5 * h/SPR0
        bp <- vb/(a * SPR0)
        R0 <- a * bp
    }else{
        stop(paste0("Only implemented for hockey-stick so far."))
    }

    dat$R0 <- R0
    dat$bp <- bp

    return(dat)
}


#' @name mclapply.windows
#' @description Alternative parallelisation for windows
#'
#' @importFrom parallel detectCores makeCluster clusterExport parLapply
#'
#' @details Reference: https://www.r-bloggers.com/2014/07/implementing-mclapply-on-windows-a-primer-on-embarrassingly-parallel-computation-on-multicore-systems-with-r/
#'
#' @export
mclapply.windows <- function(...,mc.cores = parallel::detectCores()-1) {
    ## Create a cluster
    size.of.list <- length(list(...)[[1]])
    cl <- parallel::makeCluster(spec = min(size.of.list, mc.cores) )

    ## Find out the names of the loaded packages
    loaded.package.names <- c(
        ## Base packages
        sessionInfo()$basePkgs,
        ## Additional packages
        names( sessionInfo()$otherPkgs ))

    tryCatch( {
       ## Copy over all of the objects within scope to all clusters
       this.env <- environment()
       while( identical( this.env, globalenv() ) == FALSE ) {
           parallel::clusterExport(cl,
                         ls(all.names=TRUE, env=this.env),
                         envir=this.env)
           this.env <- parent.env(environment())
       }
       parallel::clusterExport(cl,
                     ls(all.names=TRUE, env=globalenv()),
                     envir=globalenv())

       ## Load the libraries on all the clusters
       ## N.B. length(cl) returns the number of clusters
       parallel::parLapply(
                     cl,
                     1:length(cl),
                     function(xx){
           lapply(loaded.package.names, function(yy) {
               require(yy , character.only=TRUE)})
       })

       ## Run the lapply in parallel
       return( parLapply( cl, ...) )

    }, finally = {
       ## Stop the cluster
       stopCluster(cl)
    })
}


#' @name mclapply.all.os
#' @description mclapply comaptible with all OS
#'
#' @importFrom parallel mclapply
#'
#' @export
mclapply.all.os <- switch(
    Sys.info()[['sysname']],
   Windows = {iamse::mclapply.windows},
   Linux   = {parallel::mclapply},
   Darwin  = {parallel::mclapply}
)
