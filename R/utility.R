

#' genConvs
#' @description Get converged simulates from a resMSE object
#'
#' @param mse info
#' @param convyears info
#' @param convhcrs info
#' @param out info
#' @param verbose info
#'
#' @export
get.converged <- function(mse, convyears = "all", convhcrs = "all", out = 0, verbose = FALSE) {

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


#' sdconv
#'
#' @param mu info
#' @param sd info
#'
#' @export
sdconv <- function(mu, sd) (log(1 + ((sd^2)/(mu^2))))^0.5


#' muconv
#'
#' @param mu info
#' @param sd info
#'
#' @export
muconv <- function(mu, sd) log(mu) - 0.5 * log(1 + ((sd^2)/(mu^2)))


#' gen.noise
#'
#' @param n info
#' @param sd info
#' @param rho info
#' @param bias.cor info
#' @param mv info
#' @param dat info
#' @param by.asmax info
#' @param by.length info
#' @param hist info
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats rnorm rlnorm
#' @importFrom utils tail
#'
#' @export
gen.noise <- function(n,
                      sd,
                      rho = 0,
                      bias.cor = 0,
                      mv = FALSE,
                      dat = NULL,
                      by.asmax = FALSE,
                      by.length = FALSE,
                      hist = NULL) {

    if(mv){
        ## multivariate noise
        stopifnot(!is.null(dat))
        if(by.length){
            stopifnot(!is.null(dat$plba))
            amax <- dim(dat$plba)[2]
        }else{
            if(by.asmax){
                amax <- dat$asmax
            }else{
                amax <- dat$amax + 1
            }
        }
        Sigma <- matrix(NA, amax, amax)
        for(i in 1:amax) for(j in 1:amax) Sigma[i,j] = rho^abs(i - j) * sd^2
        res <- MASS::mvrnorm(n, rep(0,ncol(Sigma)), Sigma)
      if(bias.cor == 1){
        ## TODO: res is not a matrix!
        ## browser()
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

            if(!is.null(hist)){
                res[1] <- rho * tail(log(hist),1) + sqrt(1 - rho^2) * rnoise[1]
            }else{
                res[1] <- 0 ## rnoise[1]  ## HERE:
            }

            if(n > 1){
                for(i in 2:n) res[i] <- rho * res[i-1] + sqrt(1 - rho^2) * rnoise[i]
            }

            res <- exp(res)
            ## res <- res/mean(res)  ## TODO: this removes all noise, when n = 1!

        }else{
            res <- exp(rnoise)
        }
    }
    return(res)
}


#' get.errs
#'
#' @param dat info
#' @param set info
#' @param x info
#' @param hist info
#' @param rep info
#'
#'
#' @importFrom utils tail
#'
#' @export
get.errs <- function(dat, set, x, hist = NULL, rep = FALSE) {

    n <- length(x)
  nsurv <- length(dat$surveyTimes)

    if(!is.null(hist) && "errs" %in% names(hist)){
        if(rep){
            hist.errs <- hist$errs$rep
        }else{
            hist.errs <- hist$errs$time
        }
    }else{
        hist.errs <- list(eF = NULL,
                          eR = NULL,
                          eM = NULL,
                          eH = NULL,
                          eR0 = NULL,
                          eAlpha = NULL,
                          eBeta = NULL,
                          eMat = NULL,
                          eSel = NULL,
                          eW = NULL,
                          eImp = NULL,
                          eC = NULL)
        hist.errs$eI <- vector("list", nsurv)
        hist.errs$eCmvA <- NULL
        hist.errs$eCmvL <- NULL
        hist.errs$eImvA <- vector("list", nsurv)
        hist.errs$eImvL <- vector("list", nsurv)
        hist.errs$eE <- NULL
    }

    if(rep){
        errs.in <- set$errs$rep
        noise <- set$noise$rep
    }else{
        errs.in <- set$errs$time
        noise <- set$noise$time
    }
    eF <- errs.in$eF[x]
    eR <- errs.in$eR[x]
    eM <- errs.in$eM[x]
    eH <- errs.in$eH[x]
    eR0 <- errs.in$eR0[x]
    eAlpha <- errs.in$eAlpha[x]
    eBeta <- errs.in$eBeta[x]
    eMat <- errs.in$eMat[x]
    eSel <- errs.in$eSel[x]
    eW <- errs.in$eW[x]
    eImp <- errs.in$eImp[x]
    eC <- errs.in$eC[x]
    eI <- list()
    for(i in 1:nsurv){
        eI[[i]] <- errs.in$eI[[i]][x]
    }
    if(is.null(errs.in$eCmvA)){
        eCmvA <- NULL
    }else{
        if(inherits(errs.in$eCmvA, "matrix")){
            eCmvA <- try(errs.in$eCmvA[x,], silent = TRUE)
            if(inherits(eCmvA, "try-error")) eCmvA <- NA
        }else{
            eCmvA <- t(as.matrix(errs.in$eCmvA))
        }
    }
    if(is.null(errs.in$eCmvL)){
        eCmvL <- NULL
    }else{
        if(inherits(errs.in$eCmvL, "matrix")){
            eCmvL <- try(errs.in$eCmvL[x,], silent = TRUE)
            if(inherits(eCmvL, "try-error")) eCmvL <- NA
        }else{
            eCmvL <- t(as.matrix(errs.in$eCmvL))
        }
    }
    if(is.null(errs.in$eImvA)){
        eImvA <- NULL
    }else{
        eImvA <- list()
        for(i in 1:nsurv){
            if(inherits(errs.in$eImvA[[i]], "matrix")){
                eImvA[[i]] <- try(errs.in$eImvA[[i]][x,], silent = TRUE)
                if(inherits(eImvA[[i]], "try-error")) eImvA[[i]] <- NA
            }else{
                eImvA[[i]] <- t(as.matrix(errs.in$eImvA[[i]]))
            }
        }
    }
    if(is.null(errs.in$eImvL)){
        eImvL <- NULL
    }else{
        eImvL <- list()
        for(i in 1:nsurv){
            if(inherits(errs.in$eImvL[[i]], "matrix")){
                eImvL[[i]] <- try(errs.in$eImvL[[i]][x,], silent = TRUE)
                if(inherits(eImvL[[i]], "try-error")) eImvL[[i]] <- NA
            }else{
                eImvL[[i]] <- t(as.matrix(errs.in$eImvL[[i]]))
            }
        }
    }
    eE <- errs.in$eE[x]

    if(is.null(eF) || all(is.na(eF))) eF <- gen.noise(n, noise$F[1], noise$F[2],
                                                      bias.cor = noise$F[3],
                                                      hist = tail(hist.errs$eF,1))
    if(is.null(eR) || all(is.na(eR))) eR <- gen.noise(n, noise$R[1], noise$R[2],
                                                      bias.cor = noise$R[3],
                                                      hist = tail(hist.errs$eR,1))
    if(is.null(eM) || all(is.na(eM))) eM <- gen.noise(n, noise$M[1], noise$M[2],
                                                      bias.cor = noise$M[3],
                                                      hist = tail(hist.errs$eM,1))
    if(is.null(eH) || all(is.na(eH))) eH <- gen.noise(n, noise$H[1], noise$H[2],
                                                      bias.cor = noise$H[3],
                                                      hist = tail(hist.errs$eH,1))
    if(is.null(eR0) || all(is.na(eR0))) eR0 <- gen.noise(n, noise$R0[1], noise$R0[2],
                                                         bias.cor = noise$R0[3],
                                                         hist = tail(hist.errs$eR0,1))
    if(is.null(eAlpha) || all(is.na(eAlpha))) eAlpha <- gen.noise(n, noise$Alpha[1], noise$Alpha[2],
                                                                  bias.cor = noise$Alpha[3],
                                                                  hist = tail(hist.errs$eAlpha,1))
    if(is.null(eBeta) || all(is.na(eBeta))) eBeta <- gen.noise(n, noise$Beta[1], noise$Beta[2],
                                                               bias.cor = noise$Beta[3],
                                                               hist = tail(hist.errs$eBeta,1))
    if(is.null(eMat) || all(is.na(eMat))) eMat <- gen.noise(n, noise$Mat[1], noise$Mat[2],
                                                            bias.cor = noise$Mat[3],
                                                            hist = tail(hist.errs$eMat,1))
    if(is.null(eSel) || all(is.na(eSel))) eSel <- gen.noise(n, noise$Sel[1], noise$Sel[2],
                                                            bias.cor = noise$Sel[3],
                                                          , hist = tail(hist.errs$eSel,1))
    if(is.null(eW) || all(is.na(eW))) eW <- gen.noise(n, noise$W[1], noise$W[2],
                                                      bias.cor = noise$W[3],
                                                      hist = tail(hist.errs$eW,1))
    if(is.null(eImp) || all(is.na(eImp))) eImp <- gen.noise(n, noise$Imp[1], noise$Imp[2],
                                                            bias.cor = noise$Imp[3],
                                                            hist = tail(hist.errs$eImp,1))
    if(is.null(eC) || all(is.na(eC))) eC <- gen.noise(n, noise$C[1], noise$C[2],
                                                      bias.cor = noise$C[3],
                                                      hist = tail(hist.errs$eC,1))
    if(is.null(eI) || length(eI) == 0 || all(is.na(eI))){
        eI <- list()
        for(i in 1:nsurv){
            eI[[i]] <- gen.noise(n, noise$I[1], noise$I[2],
                                 bias.cor = noise$I[3],
                                 hist = tail(hist.errs$eI[[i]],1))
        }
    }
    if(is.null(eCmvA) || all(is.na(eCmvA)))
        eCmvA <- gen.noise(n, noise$CmvA[1], noise$CmvA[2],
                           bias.cor = noise$CmvA[3],
                           mv = TRUE, dat = dat,
                           by.asmax = FALSE,
                           hist = tail(hist.errs$eCmvA,1))
    if((is.null(eCmvL) || all(is.na(eCmvL))) &&
       !is.null(dat$plba))
        eCmvL <- gen.noise(n, noise$CmvL[1], noise$CmvL[2],
                           bias.cor = noise$CmvL[3],
                           mv = TRUE, dat = dat,
                           by.length = TRUE,
                           hist = tail(hist.errs$eCmvL,1))
    if(is.null(eImvA) || all(is.na(eImvA))){
        eImvA <- list()
        for(i in 1:nsurv){
            eImvA[[i]] <- gen.noise(n, noise$ImvA[1], noise$ImvA[2],
                                    bias.cor = noise$ImvA[3],
                                    mv = TRUE, dat = dat,
                                    by.asmax = FALSE,
                                    hist = tail(hist.errs$eImvA[[i]],1))
        }
    }

    if(is.null(eImvL) || all(is.na(eImvL)) &&
       !is.null(dat$plba)){
        eImvL <- list()
        for(i in 1:nsurv){
            eImvL[[i]] <- gen.noise(n, noise$ImvL[1], noise$ImvL[2],
                                    bias.cor = noise$ImvL[3],
                                    mv = TRUE, dat = dat,
                                    by.asmax = TRUE,
                                    hist = tail(hist.errs$eImvL[[i]],1))
        }
    }
    if(is.null(eE) || all(is.na(eE))) eE <- gen.noise(n, noise$E[1], noise$E[2],
                                                      bias.cor = noise$E[3],
                                                      hist = tail(hist.errs$eE,1))


    if(!is.null(hist) && "errs" %in% names(hist) && !rep){
        errs <- list(eF = c(hist.errs$eF, eF),
                     eR = c(hist.errs$eR, eR),
                     eM = c(hist.errs$eM, eM),
                     eH = c(hist.errs$eH, eH),
                     eR0 = c(hist.errs$eR0,eR0),
                     eAlpha = c(hist.errs$eAlpha, eAlpha),
                     eBeta = c(hist.errs$eBeta, eBeta),
                     eMat = c(hist.errs$eMat, eMat),
                     eSel = c(hist.errs$eSel, eSel),
                     eW = c(hist.errs$eW, eW),
                     eImp = c(hist.errs$eImp, eImp),
                     eC = c(hist.errs$eC, eC))
        errs$eI <- list()
        for(i in 1:nsurv){
            errs$eI[[i]] = c(hist.errs$eI[[i]], eI[[i]])
        }
        errs$eCmvA <- rbind(hist.errs$eCmvA, eCmvA)
        errs$eCmvL <- rbind(hist.errs$eCmvL, eCmvL)
        errs$eImvA <- list()
        for(i in 1:nsurv){
            errs$eImvA[[i]] = rbind(hist.errs$eImvA[[i]], eImvA[[i]])
        }
        errs$eImvL <- list()
        for(i in 1:nsurv){
            errs$eImvL[[i]] = rbind(hist.errs$eImvL[[i]], eImvL[[i]])
        }
        errs$eE <- c(hist.errs$eE, eE)
    }else{
        errs <- list(eF = eF,
                     eR = eR,
                     eM = eM,
                     eH = eH,
                     eAlpha = eAlpha,
                     eBeta = eBeta,
                     eR0 = eR0,
                     eMat = eMat,
                     eSel = eSel,
                     eW = eW,
                     eImp = eImp,
                     eC = eC,
                     eI = eI,
                     eCmvA = eCmvA,
                     eCmvL = eCmvL,
                     eImvA = eImvA,
                     eImvL = eImvL,
                     eE = eE)
    }

    return(errs)
}



#' Estimate depletion under alternative fishing mortalities
#'
#' `est.depletion()` evaluates the depletion level (typically biomass relative
#' to unfished biomass, \eqn{B/B_0}, or a related quantity) for a range of
#' fishing mortalities. The function simulates the IAMSE operating model over
#' a grid of fishing mortality values between `fmin` and `fmax`, using
#' `nrep` stochastic replicates, and summarises the resulting depletion
#' distribution according to the chosen `method`.
#'
#' Optionally, an optimisation step can be performed to identify the fishing
#' mortality that leads to a target depletion level, controlled by the
#' convergence tolerance `tol` and the flag `do.opt`.
#'
#' @param dat Data object as returned by [check.dat()], containing life-history
#'   and stock information for the operating model.
#' @param set Optional settings list as returned by [check.set()]. If `NULL`
#'   (default), internal defaults are used where possible. When supplied,
#'   `set` controls aspects such as the number of years, projection horizon,
#'   and noise structure.
#' @param fmin Numeric value giving the minimum fishing mortality (or
#'   F-multiplier) considered in the depletion calculation. Default is `0.0001`.
#' @param fmax Numeric value giving the maximum fishing mortality considered in
#'   the depletion calculation. Default is `10`.
#' @param nrep Integer specifying the number of stochastic replicates used for
#'   each fishing mortality level. Larger values give more stable depletion
#'   estimates at the cost of increased computation time. Default is `100`.
#' @param verbose Logical; if `TRUE`, print progress messages and basic
#'   diagnostics during the calculation. Default is `TRUE`.
#' @param method Character string indicating how the depletion distribution is
#'   summarised across replicates. The default `"percentile"` typically uses
#'   quantiles of the simulated depletion distribution. Other methods may be
#'   available depending on the internal implementation.
#' @param B0K Logical; if `TRUE`, treat the carrying capacity \eqn{K} as
#'   equivalent to unfished biomass \eqn{B_0} when computing depletion. If
#'   `FALSE` (default), \eqn{B_0} is used explicitly if available.
#' @param tol Numeric tolerance used in any optimisation or root-finding
#'   procedures (e.g. when solving for the fishing mortality associated with a
#'   target depletion level). Default is `0.0001`.
#' @param do.opt Logical; if `TRUE` (default), perform optimisation to identify
#'   fishing mortality values corresponding to target depletion levels. If
#'   `FALSE`, the function only evaluates depletion over the grid defined by
#'   `fmin` and `fmax` without attempting to optimise.
#'
#' @return An object (typically a list or data frame) summarising depletion
#'   as a function of fishing mortality, including summary statistics across
#'   replicates. The exact structure depends on the internal implementation.
#'   The object is returned invisibly if the function is mainly used for its
#'   side effects (e.g. storing results in `dat`).
#'
#' @export
est.depletion <- function(dat, set=NULL, fmin = 0.0001,
                          fmax = 10, nrep = 100, verbose = TRUE,
                          method = "percentile", B0K = FALSE,
                          tol = 0.0001, do.opt = TRUE) {

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

    if(B0K){
        outopt <- 6
        blim <- dat$ref$Blim[ny]
    }else{
        if(depl.quant %in% c("Bmsy","Blim")){
            outopt <- 2
            blim <- dat$ref$Blim[ny]
        }else if(depl.quant %in% c("SSBmsy","SSBlim")){
            outopt <- 3
            blim <- dat$ref$SSBlim[ny]
        }else stop("depl.quant not implemented. Please use Bmsy, Blim,SSBmsy or SSBlim. Or implement others.")
    }

    ## errors
    errs <- list()
    errs$time <- errs$rep <- vector("list", nrep)
    for(i in 1:nrep){
        errs$time[[i]] <- get.errs(dat, set, 1:ny)
        errs$rep[[i]] <- get.errs(dat, set, 1, rep = TRUE)
    }

    frel <- dat$FM/max(dat$FM)

    fn <- function(logfabs, frel, depl, depl.prob, nrep, dat,
                   set, errs, outopt, optFn=1){
        datx <- dat
        setx <- set
        fpat <- frel * exp(logfabs) / dat$ns
        datx$FM <- fpat
        dreal <- sapply(1:nrep, function(x){
            setx$errs <- list(time = errs$time[[x]],
                              rep = errs$rep[[x]])
            initpop(datx, setx, out.opt = outopt)
        })
        if(method == "mean"){
            drealQ <- mean(dreal)
        }else if(method == "median"){
            drealQ <- stats::quantile(dreal, probs = 0.5)
        }else if(method == "percentile"){
            drealQ <- stats::quantile(dreal, probs = depl.prob)
        }
        if(optFn==1) return((drealQ - depl)^2)
        if(optFn==2) return(drealQ)
        if(optFn==3) return(dreal)
    }


    ## opt <- nlminb(log(5), fn, lower = log(fmin), upper = log(fmax), frel = frel, depl = depl,
    ##               depl.prob = depl.prob,
    ##               nrep = nrep, dat = dat, set=set, errs=errs, outopt = outopt, optFn = 1)
    if(do.opt){
        opt <- stats::optimize(fn, c(log(fmin),log(fmax)), frel = frel, depl = depl,
                               depl.prob = depl.prob,
                               nrep = nrep, dat = dat, set=set, errs=errs,
                               outopt = outopt, optFn = 1, tol = tol)
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



#' Estimate productivity curve
#'
#' `est.productivity()` estimates a productivity relationship for the IAMSE
#' operating model, typically linking biomass (or depletion) to surplus
#' production or long-term yield. The operating model is simulated over a
#' range of fishing mortalities up to `fmax` for `ny` years, and summary
#' quantities are extracted to characterise the stock's productivity.
#'
#' The function can be used to construct productivity curves or clouds that
#' relate state variables (e.g. biomass, depletion) to production, which can
#' then be analysed or plotted for diagnostics and comparison across stocks.
#'
#' @param dat Data object as returned by [check.dat()], containing life-history
#'   and stock information for one or more stocks.
#' @param set Optional settings list as returned by [check.set()]. If `NULL`
#'   (default), internal defaults are used where possible. When supplied, `set`
#'   controls aspects such as the number of replicates, noise structure, and
#'   other operating-model options.
#' @param ny Integer giving the number of years to simulate when estimating the
#'   productivity relationship. Larger values allow the system to approach
#'   equilibrium or typical dynamic behaviour under each fishing level.
#'   Default is `100`.
#' @param fmax Numeric value giving the maximum fishing mortality (or
#'   F-multiplier) considered when constructing the productivity curve.
#'   Default is `10`.
#' @param nf Integer giving the number of fishing-mortality levels (or points
#'   along the F gradient) to consider between 0 and `fmax`. Default is
#'   `1e3` (1000 points), providing a relatively fine grid.
#' @param tsSplit Integer controlling how the simulated time series are split
#'   into segments when summarising productivity (for example, discarding
#'   an initial burn-in and using the remaining periods for analysis). The
#'   precise interpretation depends on the internal implementation. Default is
#'   `8`.
#' @param plot Logical; if `TRUE` (default), produce diagnostic plots of the
#'   productivity relationship (e.g. yield or surplus production vs. biomass
#'   or F). If `FALSE`, no plots are produced.
#'
#' @return An object (typically a list or data frame) containing productivity
#'   summaries as a function of fishing mortality and/or biomass for each
#'   stock. The exact structure depends on the internal implementation. The
#'   object is usually returned invisibly if the function is called primarily
#'   for its side effects (e.g. plotting).
#'
#' @export
est.productivity <- function(dat, set= NULL,
                             ny = 100,
                             fmax = 10,
                             nf = 1e3,
                             tsSplit = 8,
                             plot = TRUE) {

    dat$ny <- ny
    ns <- dat$ns
    nt <- ny * ns

    ## noise
    if(is.null(set)) set <- check.set()
    set$noise$time <- lapply(set$noise$time, function(x)
        c(0,0,0))
    set$noise$rep <- lapply(set$noise$rep, function(x)
        c(0,0,0))

    nyref <- set$refYears
    nrep <- set$refN
    nrep <- 1
    nyrefmsy <- set$refYearsMSY

    ## TODO: they are all 1!
    ## errors (have to be re-used for estimation of Bmsy)
    errs <- list()
    errs$time <- errs$rep <- vector("list", nrep)
    for(i in 1:nrep){
        errs$time[[i]] <- get.errs(dat, set, 1:nyref)
        errs$rep[[i]] <- get.errs(dat, set, 1, rep = TRUE)
    }

    datx <- dat

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
    blims <- rep(NA, alltv)
    prods <- vector("list", alltv)
    for(i in 1:alltv){
        datx$M <- t(as.matrix(ms[mtv[i],]))
        datx$Msel <- msels[mseltv[i]]
        datx$sel <- sels[seltv[i]]
        fms <- seq(0, fmax, length.out = nf)
        tmp2 <- vector("list", nf)
        for(fx in 1:nf){
            tmp0 <- lapply(1, function(x){
                setx <- set
                setx$errs <- list(time = errs$time[[x]],
                                  rep = errs$rep[[x]])
                logfm <- log(max(fms[fx], 1e-30))
                pop <- simpop(logfm, datx, setx, out=0)
                tsb <- tail(pop$TSB,1)
                esb <- tail(pop$ESB,1)
                ssb <- tail(pop$SSB,1)
                cw <- tail(pop$CW,1)
                sp <- tail(pop$SP,1)
                return(c(TSB = tsb, SSB = ssb, ESB = esb, CW = cw, SP = sp))
            })
            tmp1 <- do.call(rbind, tmp0)
            tmp2[[fx]] <- cbind(f = fms[fx], tmp1)
        }

        bs <- do.call(rbind, lapply(tmp2, function(x) x[,2]))
        sps <- do.call(rbind, lapply(tmp2, function(x) x[,6]))

        msy <- max(sps[,1], na.rm=TRUE)
        blims[i] <- bs[which.min(abs(sps[,1] - msy/2)),i]
        prods[[i]] <- as.data.frame(do.call(rbind, tmp2))
    }


    ## old
    if (FALSE) {
        ## increasing effort
        ##
        len1 <- len3 <- floor(ny/tsSplit)
        len2 <- ny - len1 - len3

        dat$FM <- matrix(c(rep(0,len1),
                           seq(0, fmax, length.out = len2),
                           rep(fmax,len3)) / ns, ncol=ns, nrow = ny)
        ## CHECK: how to estimate productivity with time variant M?
        datx$M <- matrix(datx$M[1,], ncol=ns, nrow=1)
        datx <- check.datx(datx)
        pop1 <- initpop(datx, set)
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
        datx$FM <- matrix(c(rep(fmax, len1),
                            seq(fmax, 0, length.out = len2),
                            rep(0, len3)) / ns, ncol=ns, nrow=ny)
        datx$M <- matrix(datx$M[1,], ncol=ns, nrow=1)
        datx <- check.datx(datx)
        pop2 <- initpop(datx, set)
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

    }

    if (plot) {

        dati <- prods[[1]]
        dati <- dati[dati$SP > 0, ]

        plot(dati$TSB, dati$SP, ty='n',
             ylab = "SP", xlab = "TSB")
        abline(h = 0, col= "grey70")
        lines(dati$TSB, dati$SP, ty='b',
              col = "dodgerblue2")

    }

    res <- list(ests = prods,
                blims = blims)
    return(res)
}



#' Estimate stochastic productivity relationship
#'
#' @param dat info
#' @param set info
#' @param fmax info
#' @param nf info
#' @param prob info
#' @param ncores info
#' @param plot info
#'
#' @return An object (typically a list or data frame) containing stochastic
#'   productivity summaries as a function of fishing mortality (and possibly
#'   biomass or depletion) for each stock, including the quantiles specified
#'   in \code{prob}. The exact structure depends on the internal implementation.
#'   The object is usually returned invisibly if the function is called
#'   primarily for its side effects (for example plotting).
#'
#' @importFrom parallel detectCores
#' @export
est.productivity.stochastic <- function(dat, set= NULL,
                                        fmax = 10,
                                        nf = 1e3,
                                        prob = c(0.1,0.9),
                                        ncores = 1,
                                        plot = TRUE) {

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
    errs <- list()
    errs$time <- errs$rep <- vector("list", nrep)
    for(i in 1:nrep){
        errs$time[[i]] <- get.errs(dat, set, 1:nyref)
        errs$rep[[i]] <- get.errs(dat, set, 1, rep = TRUE)
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
            tmp0 <- parallel::mclapply(as.list(1:nrep), function(x){
                setx <- set
                setx$errs <- list(time = errs$time[[x]],
                                  rep = errs$rep[[x]])
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
                                                        function(x) apply(x,2, stats::median, na.rm=TRUE))))
        lo[[i]] <- as.data.frame(do.call(rbind,lapply(tmp2,
                                                      function(x) apply(x,2, stats::quantile, prob=min(prob),
                                                                        na.rm=TRUE))))
        up[[i]] <- as.data.frame(do.call(rbind,lapply(tmp2,
                                                      function(x) apply(x,2, stats::quantile, prob=max(prob),
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
        ## legend("topright",
        ##        legend = c("M = 0.2", "M = 0.3"),  ## CHECK: adjust
        ##        col = cols[1:alltv],
        ##        lwd = 1.5)
    }

    res <- list(meds = meds,
                means = means,
                lo = lo,
                up = up,
                blims = blims)
    return(res)

}


#' fpat
#'
#' @param fmax info
#' @param fscen info
#'
#' @export
fpat <- function(fmax, fscen = 1) {
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



#' baranov
#'
#' @param F info
#' @param M info
#' @param N info
#'
#' @export
baranov <- function(F, M, N){
    Z <- F + M
    return(F/Z * N * (1 - exp(-Z)))
}


#' Predict catch over a TAC period
#'
#' `predCatch()` computes the predicted catch over a TAC period, given fishing
#' mortality, numbers at age, natural mortality, selectivity, and other
#' biological and stock–recruitment parameters. If a target TAC is supplied,
#' the function can also return the difference between the predicted catch and
#' the specified TAC.
#'
#' This is an internal helper used by higher-level IAMSE functions when
#' searching for fishing mortality levels that match a given TAC or when
#' projecting catches under a given management rule.
#'
#' @param logFM Array, matrix, or vector of log fishing mortality values
#'   (e.g. by age and season) for the TAC period.
#' @param NAA Numbers-at-age at the start of the TAC period (typically a
#'   vector or matrix indexed by age and possibly season).
#' @param MAA Natural mortality-at-age (and, if applicable, season) used in
#'   the projection.
#' @param sel Selectivity-at-age (and possibly season) applied to fishing
#'   mortality when computing catch.
#' @param weight Mean weight-at-age (and possibly season) used to convert
#'   numbers caught into catch biomass.
#' @param seasons Integer vector with the season indices included in the TAC
#'   period (e.g. quarters or months within a year).
#' @param ns Integer giving the total number of seasons per year.
#' @param y Integer index of the current year in the projection.
#' @param h2 Steepness (or related transformed parameter) of the
#'   stock–recruitment relationship, used when updating recruitment.
#' @param asmax Integer giving the maximum spawning age (or related age limit)
#'   used when computing spawning biomass and recruitment.
#' @param mat Maturity-at-age (and possibly season), used to calculate
#'   spawning biomass.
#' @param pzbm Additional parameter controlling the timing or proportion of
#'   biomass contributing to spawning (exact role depends on the internal
#'   implementation).
#' @param spawning Indicator (e.g. logical or integer vector) identifying the
#'   spawning season(s) within the year.
#' @param R0 Unfished recruitment level used in the stock–recruitment
#'   relationship.
#' @param SR Character string or code specifying the type of
#'   stock–recruitment function (e.g. Beverton–Holt, Ricker), as used in the
#'   operating model.
#' @param bp Additional parameter (or vector of parameters) for the
#'   stock–recruitment relationship (e.g. breakpoints or shape parameters).
#' @param recBeta Parameter (or vector) controlling the slope or shape of the
#'   recruitment function.
#' @param recGamma Parameter (or vector) controlling the curvature or shape of
#'   the recruitment function.
#' @param eR Stochastic recruitment deviation(s) applied to the expected
#'   recruitment (e.g. lognormal errors).
#' @param indage0 Integer index specifying the recruitment age (age at which
#'   new recruits enter the age structure).
#' @param TAC Optional numeric value giving the target total allowable catch
#'   for the TAC period. If `NULL` (default), only the predicted catch is
#'   returned. If non-`NULL`, the function can return the difference between
#'   predicted catch and `TAC`, depending on `out`.
#' @param out Integer flag controlling the type of output returned. The
#'   precise meaning depends on the internal implementation but typically
#'   distinguishes between returning the predicted catch, the difference
#'   between catch and `TAC`, or possibly additional diagnostics.
#'
#' @return A numeric value or vector containing the predicted catch over the
#'   TAC period, or the difference between predicted and target catch if
#'   `TAC` is supplied and `out` is set accordingly. The exact structure of
#'   the return value depends on the internal implementation. This function is
#'   intended for internal use within IAMSE.
#'
#' @keywords internal
#'
predCatch <- function(logFM,
                      NAA, MAA,
                      sel, weight,
                      seasons, ns, y, h2, asmax, mat, pzbm, spawning,
                      R0, SR, bp, recBeta, recGamma, eR,
                      indage0,
                      TAC = NULL,
                      out = 0) {
    Ctmp <- 0
    NAAtmp <- NAA

    for(s in seasons){
        FAA <- exp(logFM) * sel
        Ztmp <- FAA + MAA
        ## recruitment
        if(spawning[s] > 0 && s > 1){
            ## Survivors from previous season/year
            SSB0 <- get.ssb0(MAA, mat, weight,
                             fecun=1, asmax, ns, spawning,
                             R0, indage0, s)
            SSBtmp <- sum(NAA * weight  * mat  * exp(-pzbm * Ztmp))
            rec <- recfunc(h = h2, SSBPR0 = SSB0/R0, SSB = SSBtmp,
                           R0 = R0, method = SR, bp = bp,
                           beta = recBeta, gamma = recGamma)
            rec[rec<0] <- 1e-10
            NAAtmp[indage0] <- rec * spawning[s] * eR
        }
        Ctmp <- Ctmp + sum(baranov(FAA, MAA, NAAtmp) * weight)
        ## Aging
        NAAtmp <- NAAtmp * exp(-Ztmp)
        NAAtmp[asmax] <- NAAtmp[asmax] + NAAtmp[asmax-1]
        for(as in (asmax-1):2) NAAtmp[as] <- NAAtmp[as-1]
        NAAtmp[indage0] <- 0
    }

    if(out == 0){
        return(Ctmp)
    }else{
        return((TAC - Ctmp)^2)
    }
}



#' Solve for fishing mortality corresponding to a target TAC
#'
#' `get.f()` finds the fishing mortality (FM) that produces a given total
#' allowable catch (TAC) over a TAC period, accounting for seasonal structure,
#' natural mortality, selectivity, and recruitment. Internally, it repeatedly
#' calls the catch-prediction logic (e.g. [predCatch()]) and uses a numerical
#' root-finding or optimisation routine to match predicted catch to the target
#' `TAC`.
#'
#' This function is typically called internally by higher-level IAMSE functions
#' when implementing HCRs that are specified in terms of TAC (rather than
#' directly in terms of fishing mortality).
#'
#' @param TAC Numeric value giving the target total allowable catch for the TAC
#'   period (e.g. over the seasons specified in `seasons`).
#' @param NAA Numbers-at-age at the start of the TAC period (typically a vector
#'   or matrix indexed by age and possibly season).
#' @param MAA Natural mortality-at-age (and, if applicable, season) used in the
#'   projection.
#' @param sel Selectivity-at-age (and possibly season) applied to fishing
#'   mortality when computing catch.
#' @param weight Mean weight-at-age (and possibly season) used to convert
#'   numbers caught into catch biomass.
#' @param seasons Integer vector with the season indices included in the TAC
#'   period (e.g. quarters or months within a year).
#' @param ns Integer giving the total number of seasons per year.
#' @param y Integer index of the current year in the projection.
#' @param h Steepness (or related parameter) of the stock–recruitment
#'   relationship, used when updating recruitment along the projection.
#' @param asmax Integer giving the maximum spawning age (or related age limit)
#'   used when computing spawning biomass and recruitment.
#' @param mat Maturity-at-age (and possibly season), used to calculate
#'   spawning biomass.
#' @param pzbm Additional parameter controlling the timing or proportion of
#'   biomass contributing to spawning (exact role depends on the internal
#'   implementation).
#' @param spawning Indicator (e.g. logical or integer vector) identifying the
#'   spawning season(s) within the year.
#' @param R0 Unfished recruitment level used in the stock–recruitment
#'   relationship.
#' @param SR Character string or code specifying the type of
#'   stock–recruitment function (e.g. Beverton–Holt, Ricker), as used in the
#'   operating model.
#' @param bp Additional parameter (or vector of parameters) for the
#'   stock–recruitment relationship (e.g. breakpoints or shape parameters).
#' @param recBeta Parameter (or vector) controlling the slope or shape of the
#'   recruitment function.
#' @param recGamma Parameter (or vector) controlling the curvature or shape of
#'   the recruitment function.
#' @param eR Stochastic recruitment deviation(s) applied to the expected
#'   recruitment (e.g. lognormal errors).
#' @param indage0 Integer index specifying the recruitment age (age at which
#'   new recruits enter the age structure).
#' @param lastFM Numeric value giving the starting (or typical) fishing
#'   mortality used as an initial guess in the numerical search. Default is
#'   `0.01`.
#' @param upper Numeric value giving the upper bound (on the log scale) for
#'   the search over fishing mortality. Default is `log(100)`, corresponding
#'   to a relatively high maximum F.
#'
#' @details
#' The function typically:
#' \enumerate{
#'   \item Defines an objective function measuring the difference between the
#'         predicted catch (for a given FM) and the target `TAC`.
#'   \item Searches over FM (on a log scale) between a lower bound and
#'         `upper`, using `lastFM` as an initial value.
#'   \item Returns the FM that minimises this difference (or sets the root to
#'         zero), within numerical tolerance.
#' }
#' The exact optimisation method depends on the internal implementation (e.g.
#' use of [stats::optimize()] or [stats::uniroot()]).
#'
#' @return A numeric value giving the fishing mortality (or F-multiplier) that
#'   produces a predicted catch matching `TAC` as closely as possible under
#'   the given biological and seasonal assumptions.
#'
#' @export
get.f <- function(TAC,
                  NAA, MAA,
                  sel, weight,
                  seasons, ns, y, h, asmax, mat,
                  pzbm, spawning,
                  R0, SR, bp, recBeta,
                  recGamma, eR,
                  indage0,
                  lastFM = 0.01, upper = log(100)){

    opt <- nlminb(start = log(lastFM), objective = predCatch,
                  NAA = NAA, MAA = MAA,
                  sel = sel, weight = weight,
                  seasons = seasons, ns = ns, y = y,
                  h2 = h, asmax = asmax, mat = mat,
                  pzbm = pzbm, spawning = spawning,
                  R0 = R0, SR = SR, bp = bp, recBeta = recBeta,
                  recGamma = recGamma, eR = eR,
                  indage0 = indage0,
                  TAC = TAC,
                  out = 1,
                  lower = -10, upper = upper,
                  control = list(rel.tol = 1e-3))

    ##    print(paste0("obj: ",round(opt$objective,2), "- fm: ",round(exp(opt$par),3)))
    return(exp(opt$par))
}


#' getSel
#'
#' @description Function to estimate selectivity ogive
#'
#' @param L50 length at 50\% selectivity
#' @param L95 length at 95\% selectivity
#' @param mids midlengths
#' @param plba probability of being in mids given age
#'
getSel <- function(L50, L95, mids, plba) {
    n <- max(c(length(L50),length(L95)))
    sel <- vector("list", n)
    for(i in 1:n){
        selL <- (1 /(1 + exp(-log(19)*(mids - L50[i])/(L95[i] - L50[i]))))
        ## dims <- dim(plba)
        ## selA <- matrix(NA, ncol = dims[3], nrow = dims[1])
        ## for(j in 1:dim(plba)[3]){
        ##     selA[,j] <- apply(t(plba[,,j]) * selL, 2, sum)
        ## }
        selA <- apply(t(plba) * selL, 2, sum)
        ##        selA[1] <- 1e-9 # it should be zero for age 0
        sel[[i]] <- selA
    }
    return(sel)
}


#' getMat
#'
#' @description Function to estimate maturity at age
#'
#' @param Lm50 length at 50\% maturity
#' @param Lm95 length at 95\% maturity
#' @param mids midlengths
#' @param plba probability of being in mids given age
#'
getMat <- function(Lm50, Lm95, mids, plba) {
    ## maturity at length
    matL <- (1 /(1 + exp(-log(19)*(mids - Lm50)/(Lm95 - Lm50))))
    ## maturity at age
    ## dims <- dim(plba)
    ## matA <- matrix(NA, ncol = dims[3], nrow = dims[1])
    ## for(i in 1:dim(plba)[3]){
    ##     matA[,i] <- apply(t(plba[,,i]) * matL, 2, sum)
    ## }
    matA <- apply(t(plba) * matL, 2, sum)
    ##   matA <- c(1e-9,matA[-1])
    return(matA)
}


#' getM
#'
#' @description Function to estimate selectivity of natural mortality
#'
#' @param Linf Linf of vBGF
#' @param K K of vBGF
#' @param mids midlengths
#' @param a info
#' @param b info
#' @param c info
#'
getM <- function(Linf, K, mids, a = 0.55, b = 1.61, c = 1.44) {
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


#' getMsel
#'
#' @description Function to estimate selectivity of natural mortality
#'
#' @param Linf Linf of vBGF
#' @param K K of vBGF
#' @param mids midlengths
#' @param plba probability of being in mids given age
#' @param a info
#' @param b info
#' @param c info
#' @param scale logical; scale to have maximum 1 (default: true)
#' @param lunit info
#' @param useBelow10 info
#'
getMsel <- function(Linf, K, mids, plba, a = 0.55, b = 1.61, c = 1.44,
                    scale = TRUE, lunit = "cm", useBelow10 = FALSE) {
    if(lunit == "mm"){
        mids <- mids / 10
        Linf <- Linf / 10
    }else if(lunit != "cm"){
        stop("Can only work with lunit 'mm' or 'cm'!")
    }
    n <- max(c(length(a),length(b),length(c)))
    sel <- vector("list", n)
    for(i in 1:n){
        selL <- exp(a[i] - b[i] * log(mids) + c[i] * log(Linf) + log(K))
        if(!useBelow10) {
            selL[mids < 10] <- exp(a[i] - b[i] * log(10) + c[i] * log(Linf) + log(K))
        }

        ## dims <- dim(plba)
        ## selA <- matrix(NA, ncol = dims[3], nrow = dims[1])
        ## for(j in 1:dim(plba)[3]){
        ##     selA[,j] <- apply(t(plba[,,j]) * selL, 2, sum)
        ## }
        selA <- apply(t(plba) * selL, 2, sum)
        if(scale){
            maxM <- max(selA)
        }else{
            maxM <- 1
        }
        sel[[i]] <- selA/maxM
    }
    return(sel)
}




#' Calculate spawning stock biomass at F = 0 (SSB0)
#'
#' `get.ssb0()` computes spawning stock biomass under zero fishing mortality
#' (often referred to as SSB0), based on natural mortality, maturity, weight,
#' and fecundity-at-age, accounting for seasonal structure and the timing of
#' spawning. Optionally, a non-zero fishing mortality (`FM`) can be supplied to
#' obtain spawning biomass under that F instead.
#'
#' @param M Natural mortality-at-age (and, if applicable, season), typically
#'   given as a vector or matrix of instantaneous natural mortality rates.
#' @param mat Maturity ogive, giving the proportion mature at each age (and
#'   season, if relevant).
#' @param weight Mean weight-at-age (and possibly season), used to convert
#'   numbers of mature fish into spawning biomass.
#' @param fecun Fecundity-at-age (and possibly season). Defaults to `1`, in
#'   which case fecundity is assumed to be proportional to weight and maturity
#'   only.
#' @param asmax Integer giving the maximum age class (or the number of age
#'   classes) used in the calculation.
#' @param ns Integer giving the number of seasons per year.
#' @param spawning Indicator (e.g. logical or integer vector of length `ns`)
#'   identifying the spawning season(s) within the year.
#' @param R0 Unfished recruitment level used to scale the numbers-at-age under
#'   F = 0. Depending on the implementation, this may be used to convert
#'   per-recruit quantities into absolute biomass.
#' @param indage0 Integer index specifying the recruitment age (age at which
#'   new recruits enter the age structure).
#' @param season Integer or vector specifying the season(s) for which SSB is
#'   evaluated; typically aligned with `spawning` and the seasonal structure.
#' @param FM Optional fishing mortality-at-age (and season). If `NULL`
#'   (default), fishing mortality is set to zero and SSB0 is computed. If
#'   non-`NULL`, the function returns spawning biomass under the combined
#'   mortality `M + FM`.
#'
#' @return A numeric value giving spawning stock biomass under F = 0 (SSB0),
#'   or under the supplied `FM` if not `NULL`. The exact scaling (per recruit
#'   vs. absolute biomass) depends on how `R0` is used internally.
#'
#' @export
get.ssb0 <- function (M, mat, weight, fecun = 1,
                      asmax, ns, spawning,
                      R0, indage0, season, FM=NULL) {

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
        season <- season - 1
    }

    ## SSB0
    SBB0 <- sum(NAAS * mat * weight * fecun)

    return(SBB0)
}


#' recfunc
#'
#' @description Function to calculate recruitment (Beverton - Holt)
#'
#' @param h steepness
#' @param SSBPR0 spawning biomass produced by one recrut in its lifetime
#' @param SSB spawning biomass
#' @param R0 recruitment in unfished population
#' @param method SR type
#' @param bp breakpoint for hockey-stick SR
#' @param beta info
#' @param gamma info
#' @param alpha info
#'
#' @export
recfunc <- function(h, SSBPR0, SSB,  R0 = 1e6, method = "bevholt", bp = 0,
                    beta = 0, gamma = 0, alpha = 0) {

    if(method == "bevholt"){
        alpha <- SSBPR0 * (1-h)/(4*h)
        beta <- (5*h-1) / (4*h*R0)
        rec <- SSB / (alpha + beta * SSB)
    }else if(method == "ricker"){
        ## beta <- log(5 * h) / (0.8 * R0)
        ## alpha <- exp(beta * R0)/SSBPR0
        ## rec <- alpha * SSB * exp(-beta * SSB)
        rec <- bp * SSB * exp(-beta * SSB)
    }else if(method == "average"){
        rec <- rep(R0, length(SSB))
    }else if(method == "hockey-stick"){
        rec <- ifelse(SSB > bp, R0, SSB * R0/bp)
    }else if(method == "bent-hyperbola"){  ## Watts-Bacon bent hyperbola
        rec <- beta * (SSB + sqrt(bp^2 + (gamma^2)/4) - sqrt((SSB-bp)^2 + (gamma^2)/4))
    }else if(method == "bevholt2"){
        rec <- SSB / (alpha + beta * SSB)
    }else print("Stock-recruitment method not known! Implemented methods: 'bevholt', 'ricker', 'average', and 'hockey-stick'.")

    return (rec)
}



#' initdistR
#'
#' @param M info
#' @param FM info
#' @param ns info
#' @param asmax info
#' @param indage0 info
#' @param spawning info
#' @param R0 info
#'
#' @export
initdistR <- function(M, FM=NULL, ns, asmax, indage0, spawning, R0=1){

    if(is.null(FM)){
        FM <- 0
    }

    NAA2 <- NAA <- matrix(0, asmax, ns)
    NAA[indage0,] <- R0 * spawning
    ZAA <- M + FM

    ## each season
    for(as in (indage0+1):asmax)
        NAA[as,] <- NAA[as-1,] * exp(-ZAA[as-1])
    ## only keep age groups present relative to end of year (last season)
    for(s in 1:ns){
        indi <- seq(ns+2-s+indage0-1,asmax,ns)  ## NEW:
        ## indi <- seq(s+indage0-1,asmax,ns)
        NAA2[indi,s] <- NAA[indi,s]
    }
    ## keep last age group for every season
    NAA2[asmax,] <- NAA2[asmax,] + NAA[asmax,] * exp(-ZAA[asmax]) - NAA2[asmax,] + NAA[asmax,] - NAA[asmax,] * exp(-ZAA[asmax])
    ## better than which:
    ## indi <- which(NAA2[asmax,]==0)
    ## NAA2[asmax,indi] <- NAA[asmax,indi] * exp(-ZAA[asmax])

    ## plus group correction
    NAA2[asmax,] <- NAA2[asmax,] / (1 - exp(-sum(ZAA[(asmax-ns+1):asmax])))
    ## NAA2[(asmax-ns+1):asmax,] <- NAA2[(asmax-ns+1):asmax,] / (1 - exp(-sum(ZAA[(asmax-ns+1):asmax])))
    ## combine seasons
    NAAS <- rowSums(NAA2)
    ## remove recruits
    NAAS[indage0] <- 0

    return(NAAS)
}


check.par <- function(x, n){
    if(length(x) != n) c(x, rep(x, n-length(x))) else x
}


##' get.leslie.r
##'
##' @description Function to estimate intrinsic growth rate r and generation
##'     time
##'
##' @param dat iamse data set
##'
##' @details Function adopted from Henning Winker's spmpriors package
##'
##' @export
get.leslie.r <- function(dat){

    MAA <- dat$M[nrow(dat$M),] * dat$Msel[[1]]
    FM <- 0
    ns <- dat$ns
    asmax <- dat$asmax
    indage0 <- dat$indage0
    spawning <- dat$spawning
    R0 <- dat$R0
    weight <- dat$weight
    mat <- dat$mat
    fecun <- dat$fecun
    season <- 1

    ## mean spawner weight at age
    swt <- weight * mat

    ZAA <- MAA + FM
    NAAS <- initdistR(MAA, FM = FM, ns, asmax, indage0, spawning, R0)
    while(season > 1){
        NAAS <- NAAS * exp(-ZAA)
        NAAS[asmax] <- NAAS[asmax] + NAAS[asmax - 1]
        for(as in (asmax - 1):2) NAAS[as] <- NAAS[as - 1]
        season <- season - 1
    }
    NAAS[1] <- R0
    NAAS <- NAAS / R0

    ##  compute unfished Spawning biomass per recruit (SBR0)
    ssb0 <- get.ssb0(MAA,
                     mat, weight, fecun,
                     asmax, ns, spawning,
                     R0, indage0, season)
    ssbpr0 <- ssb0 / R0
    ## ssbpr0 <- sum(NAAS * weight * mat) ## same

    R <- recfunc(dat$h, ssbpr0, 1, dat$R0, dat$SR, dat$bp,
                 dat$recBeta, dat$recGamma, dat$recAlpha)

    ## Make Leslie matrix
    lm <- mat.or.vec(dat$amax+1, dat$amax+1)
    lm[1,] <- R * as.numeric(by(swt, dat$as2a, mean))

    ## fill rest of Matrix with Survival
    maa.annual <- as.numeric(by(MAA, dat$as2a, sum))
    for(i in 2:dat$amax){
        lm[i,(i-1)] <- exp(-maa.annual[i-1])
    }

    ## Net reproductive rate
    nrr <- sum(NAAS * swt)

    ## Intrinsic rate of population increase and generation time
    res <- list(r = log(as.numeric(eigen(lm)$values[1])),
                gt = sum(as.vector(t(dat$ages)) * NAAS * swt) / nrr)

    return(res)
}
