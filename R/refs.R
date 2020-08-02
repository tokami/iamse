## estimate reference levels for OM



#' @name simpopR
#' @export
simpopR <- function(FM, dat, set){

    ny <- set$refYears
    ns <- dat$nseasons
    nt <- ny * ns
    ## errors
    eF <- exp(rnorm(ny, 0, set$sigmaF) - set$sigmaF^2/2)
    eR <- genDevs(ny, set$sigmaR, set$rhoR)
    eM <- exp(rnorm(ny, 0, set$sigmaM) - set$sigmaM^2/2)
    eH <- exp(rnorm(ny, 0, set$sigmaH) - set$sigmaH^2/2)
    eMat <- exp(rnorm(ny, 0, set$sigmaMat) - set$sigmaMat^2/2)
    eR0 <- exp(rnorm(ny, 0, set$sigmaR0) - set$sigmaR0^2/2)
    eImp <- exp(rnorm(ny, 0, set$sigmaImp) - set$sigmaImp^2/2)
    ## some variables
    amax <- dat$amax + 1
    weight <- dat$weight
    weights <- dat$weights
    weightF <- dat$weightF
    weightFs <- dat$weightFs
    R0 <- dat$R0
    pzbm <- dat$pzbm
    sels <- dat$sels
    ## containers
    NAA <- rep(0, amax)
    TSB <- TSB1plus <- SSB <- CW <- matrix(0, nrow=ny, ncol=ns)
    ## initialize starting values
    NAA[1] <- exp(dat$initN[1]) * R0
    for(a in 2:amax)
        NAA[a] <- NAA[a-1] * exp(-(dat$M[a-1] + FM * dat$sel[a-1])) * exp(dat$initN[a])
    ## loop
    for(y in 1:ny){
        ## adding noise
        FAA <- FM / ns * sels * eF[y]
        MAA <- dat$Ms * eM[y]
        ZAA <- MAA + FAA
        maty <- dat$mats * eMat[y]
        hy <- dat$h * eH[y]
        R0y <- dat$R0 * eR0[y]
        ## recruitment
        SSBtemp <- sum(NAA * weights[,1] * maty[,1] * exp(-pzbm * ZAA[,1])) ## pre-recruitment mort
        SSBPR0 <- getSSBPR0(dat$M * eM[y], dat$mat * eMat[y], 1, amax)
        rec <- recfunc(h = hy, SSBPR0 = SSBPR0, SSB = SSBtemp,
                       R0 = R0y, method = dat$SR, bp = dat$bp)
        rec[rec<0] <- 1e-10
        NAA[1] <- rec * eR[y]
        ## seasons
        for(s in 1:ns){
            CAA <- baranov(FAA[,s], MAA[,s], NAA)
            CW[y,s] <- sum(CAA * weightFs[,s])
            ## TSB
            TSB[y,s] <- sum(NAA * weights[,s])
            TSB1plus[y,s] <- sum(NAA[-1] * weights[-1,s])
            ## SSB
            SSB[y,s] <- sum(NAA * weights[,s] * maty[,s] * exp(-pzbm * ZAA[,s]))
            NAA <- Ntemp <- NAA * exp(-ZAA[,s])
            if(s == ns){
                NAA[amax] <- Ntemp[amax] + Ntemp[amax-1]
                for(a in 2:(amax-1)) NAA[a] <- Ntemp[a-1]
            }
        }
    }

    ## surplus production (for reflevs)
    TSBy <- TSB[,1]
    TSB1plusy <- TSB1plus[,1]
    SSBy <- SSB[,1]
    CWy <- apply(CW, 1, sum)
    SP <- rep(NA, ny)
    for(y in 1:(ny-1)){
        SP[y] <- TSBy[y+1] - TSBy[y] + CWy[y]
    }

    ## return
    return(list(CW=CWy,
                TSB=TSBy,
                TSB1plus=TSB1plusy,
                SP=SP,
                SSB=SSBy))
}


#' @name estRef
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @export
#'
estRef <- function(dat, set=NULL, fvec = seq(0,5,0.1),
                   ncores=parallel::detectCores()-1,
                   ref = c("Fmsy","Bmsy","MSY","B0"),
                   plot = FALSE){

    ny <- dat$ny
    ns <- dat$nseasons
    nt <- ny * ns

    if(is.null(set)) set <- checkSet()
    ## Remove variability
    set$sigmaF <- 0
    set$sigmaR <- 0
    set$rhoR <- 0
    set$sigmaM <- 0
    set$sigmaH <- 0
    set$sigmaMat <- 0
    set$sigmaR0 <- 0
    set$sigmaImp <- 0

    if(any(ref %in% c("Fmsy","Bmsy","MSY"))){
        res0 <- parallel::mclapply(as.list(fvec),
                                   function(x){
                                       with(simpop(log(x), dat, set, out = 0),
                                            c(x,
                                              head(tail(TSB,2),1),
                                              head(tail(SP,2),1),
                                              head(tail(CW,2),1)))
                                   },
                                   mc.cores = ncores)

        ## refs
        res <- do.call(cbind, res0)
        msy <- res[4, which.max(res[3,])]
        fmsy <- res[1, which.max(res[3,])]
        bmsy <- res[2, which.max(res[3,])]
        binf <- res[2, which.max(res[2,])]

        if(plot){
            plot(fvec, res[3,], ty='b')
        }
    }


    ## B0
    if(any(ref == "B0")){
        unfished <- simpop(log(1e-20), dat, set, out = 0)
        b0 <- tail(unfished$TSB,1) ## CHECK: why TSB1plus?
    }


    refs <- list()
    if(any(ref == "Fmsy")) refs$Fmsy = fmsy
    if(any(ref == "MSY")) refs$MSY = msy
    if(any(ref == "Bmsy")) refs$Bmsy = bmsy
    if(any(ref == "B0")) refs$B0 = b0

    dat$ref <- refs

    ## return
    return(dat)
}


#' @name estRefStoch
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @export
#'
estRefStoch <- function(dat, set=NULL, fvec = seq(0,5,0.1),
                   ncores=parallel::detectCores()-1,
                   ref = c("Fmsy","Bmsy","MSY","B0"),
                   plot = FALSE){

    ny <- dat$ny
    ns <- dat$nseasons
    nt <- ny * ns
    nyref <- set$refYears
    nrep <- set$refN

    if(is.null(set)) set <- checkSet()
    dist <- NULL

    ## errors (have to be re-used for estimation of Bmsy)
    errs <- vector("list", nrep)
    for(i in 1:nrep){
        errs[[i]]$eR <- genDevs(nyref, set$sigmaR, set$rhoR)
        errs[[i]]$eF <- exp(rnorm(nyref, 0, set$sigmaF) - set$sigmaF^2/2)
        errs[[i]]$eM <- exp(rnorm(nyref, 0, set$sigmaM) - set$sigmaM^2/2)
        errs[[i]]$eR0 <- exp(rnorm(nyref, 0, set$sigmaR0) - set$sigmaR0^2/2)
        errs[[i]]$eH <- exp(rnorm(nyref, 0, set$sigmaH) - set$sigmaH^2/2)
        errs[[i]]$eMat <- exp(rnorm(nyref, 0, set$sigmaMat) - set$sigmaMat^2/2)
        errs[[i]]$eImp <- exp(rnorm(nyref, 0, set$sigmaImp) - set$sigmaImp^2/2)
    }


    if(any(ref %in% c("Fmsy","Bmsy","MSY"))){

        ## Fmsy & MSY
        res <- parallel::mclapply(as.list(1:nrep), function(x){
            setx <- c(set, errs[[x]])
            opt <- optimise(function(x) unlist(simpop(x, dat, setx, out=1)),
                            log(c(0.001,10)), maximum = TRUE)
            c(exp(opt$maximum), opt$objective)
        }, mc.cores = ncores)
        dist <- do.call(rbind, res)

        ## Bmsy
        bmsy <- sapply(1:nrep, function(x){
            setx <- c(set, errs[[x]])
            tail(simpop(log(dist[x,1]), dat, setx, out = 0)$TSB, 1)
            })
        dist <- cbind(dist, bmsy)

        colnames(dist) <- c("Fmsy","MSY","Bmsy")
    }

    ## B0
    if(any(ref == "B0")){
        b0 <- sapply(1:nrep, function(x){
            setx <- c(set, errs[[x]])
            tail(simpop(log(1e-20), dat, setx, out = 0)$TSB, 1)
            })
        dist <- cbind(dist, b0)
        colnames(dist) <- c(colnames(dist)[-ncol(dist)], "B0")
    }

    ## remove runs where long-term SP is smaller or equal to 0
    nami <- colnames(dist)
    if("MSY" %in% nami){
        dist2 <- dist[dist[,2] > 0,]
    }else
        dist2 <- dist
    meds <- apply(dist2,2, mean)


    refs <- list()
    refs$dist <- dist
    if(any(ref == "Fmsy")) refs$Fmsy = as.numeric(meds[which(nami == "Fmsy")])
    if(any(ref == "MSY")) refs$MSY = as.numeric(meds[which(nami == "MSY")])
    if(any(ref == "Bmsy")) refs$Bmsy = as.numeric(meds[which(nami == "Bmsy")])
    if(any(ref == "B0")) refs$B0 = as.numeric(meds[which(nami == "B0")])

    dat$ref <- refs

    ## return
    return(dat)
}




#' @name estRefBent
#'
#' @importFrom parallel detectCores
#' @importFrom parallel mclapply
#'
#' @export
#'
estRefBent <- function(dat, set=NULL, fvec = seq(0,5,0.1),
                   ncores=parallel::detectCores()-1,
                   ref = c("Fmsy","Bmsy","MSY","B0"),
                   plot = FALSE){
    ny <- dat$ny
    ns <- dat$nseasons
    nt <- ny * ns
    nyref <- set$refYears
    nrep <- set$refN

    if(is.null(set)) set <- checkSet()
    ## Remove variability
    set$sigmaF <- 0
    set$sigmaR <- 0
    set$rhoR <- 0
    set$sigmaM <- 0
    set$sigmaH <- 0
    set$sigmaMat <- 0
    set$sigmaR0 <- 0
    set$sigmaImp <- 0

    browser()

    fvec <- seq(0,1,0.001)


    ssbprs <- sapply(fvec, sbr, dat = dat, out=0)
    yprs <- sapply(fvec, sbr, dat = dat, out=1)

    bp <- dat$bp
    beta <- dat$recBeta
    gamma <- dat$recGamma
    K <- sqrt(bp^2 + gamma^2/4)
    ssbe <- (2*K/(ssbprs*beta) - 2*bp - 2*K) / ((1/(ssbprs^2*beta^2)) - 2/(ssbprs*beta))

    ssbe <- (2*K - 2 *bp - 2*K) / ((beta * ssbprs)^(-2) - 2 * (beta*ssbprs)^-1)


    ye <- yprs * ssbe / ssbprs

    plot(ye)

    plot(ye, ty='l')

    2 * beta
    1/ssbprs

    ye

    ssbprs


    ## B0
    if(any(ref == "B0")){
        unfished <- simpop(log(1e-20), dat, set, out = 0)
        b0 <- tail(unfished$TSB,1) ## CHECK: why TSB1plus?
    }


    refs <- list()
    if(any(ref == "Fmsy")) refs$Fmsy = fmsy
    if(any(ref == "MSY")) refs$MSY = msy
    if(any(ref == "Bmsy")) refs$Bmsy = bmsy
    if(any(ref == "B0")) refs$B0 = b0

    dat$ref <- refs

    ## return
    return(dat)
}


#' @name sbr
#' @export
sbr <- function(FM, dat, out=0){
    ## some variables
    amax <- dat$amax + 1
    ## initialize starting values
    NAA <- rep(0, amax)
    NAA[1] <- 1
    for(a in 2:amax)
        NAA[a] <- NAA[a-1] * exp(-(dat$M[a-1] + FM * dat$sel[a-1]))

    SSBPR <- sum(NAA * dat$weight * dat$mat * dat$fecun)
    YPR <- sum(baranov(FM*dat$sel, dat$M, NAA) * dat$weightFs)

    if(out == 0) res <- SSBPR else  res <- YPR

    return(res)
}
