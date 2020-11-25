## observation model


## make observations with given CV for catch and index (only one index time series so far)
#' @name obsmod
#' @export
obsmod <- function(dat, hist, set, years = NULL){
    ## parameters
    CW <- hist$CW
    q <- dat$q
    ## indices
    if(is.null(years)){
        idx <- 1:length(CW)
    }else{
        idx <- years
    }
    ny <- length(idx)
    ns <- dat$ns
    nt <- ny * ns
    nsC <- set$catchSeasons
    seasonStart <- seq(0,1-1/ns,1/ns)
    if("obsC" %in% names(hist$obs)) ny <- 1
    Ms <- dat$M / ns
    ## errors
    eC <- set$eC
    eI <- set$eI
    if(is.null(eC)) eC <- exp(rnorm(ny, 0, set$CVC) - set$CVC^2/2)
    if(is.null(eI)) eI <- exp(rnorm(ny, 0, set$CVI) - set$CVI^2/2)

    nsurv <- length(set$surveyTimes)
    if(length(q) < nsurv) q <- rep(q, nsurv)

    if("obsC" %in% names(hist$obs)){
        ## adding observations
        ## -----------------------
        y <- nrow(hist$TSB)
        ## catch
        if(ns > 1){
            if(nsC == 1){
                if(inherits(CW[y,],"matrix")) ctemp <- apply(CW[y,], 1, sum) else ctemp <- sum(CW[y,])
                obsC <- c(hist$obs$obsC, ctemp * eC[y-dat$ny])
                timeC <- c(hist$obs$timeC, tail(hist$obs$timeC,1)+1)
            }else if(nsC == ns){
                obsC <- c(hist$obs$obsC, as.numeric(t(CW[y,] * eC[y-dat$ny])))
                timeC <- c(hist$obs$timeC, tail(hist$obs$timeC, nsC) + 1)
            }else{
                stop("Seasons of catch observations and operating model do not match. This aggregation is not yet implemented!")
            }
        }else{
            if(nsC > 1) writeLines("Set dat$ns to > 1 for seasonal catches. Generating annual catches!")
            obsC <- c(hist$obs$obsC, CW[y,] * eC[y-dat$ny])
            timeC <- c(hist$obs$timeC, tail(hist$obs$timeC,1)+1)
        }

        ## survey indices
        obsI <- vector("list", nsurv)
        timeI <- vector("list", nsurv)
        for(i in 1:nsurv){
            ## closest season
            tmp <- seasonStart[seasonStart < set$surveyTimes[i]]
            idxS <- which.min((tmp - set$surveyTimes[i])^2)
            ## correct NAA for survey timing dependent on seasons in opmod
            Z <- hist$FAA[,y,idxS] + Ms
            surveyTime <- set$surveyTimes[i] - seasonStart[idxS]
            naa <- exp(log(hist$NAA[,y,idxS]) - Z * surveyTime)
            if(inherits(naa, "matrix")){
                esb <- apply(naa, 2, function(x) sum(x * dat$weightF * dat$sel))
            }else{
                esb <- sum(naa * dat$weightF * dat$sel)
            }
            obsI[[i]] <- c(hist$obs$obsI[[i]], q[i] * esb * eI[y-dat$ny])
            timeI[[i]] <- c(hist$obs$timeI[[i]], floor(tail(hist$obs$timeI[[i]],1)) + 1 + set$surveyTimes[i])
        }

     }else{
        ## initalizing observations
        ## -----------------------
        ## catch
        if(ns > 1){
            if(nsC == 1){
                obsC <- apply(CW[idx,], 1, sum) * eC
                timeC <- 1:ny
            }else if(nsC == ns){
                obsC <- as.numeric(t(CW[idx,] * eC))
                timeC <- rep(1:ny, each = ns) + rep(seasonStart, ny)
            }else{
                stop("Catch observation seasons and operating model seasons do not match. Not yet implemented!")
            }
        }else{
            if(nsC > 1) writeLines("Set dat$ns to > 1 for seasonal catches. Generating annual catches!")
            obsC <- CW[idx,] * eC
            timeC <- 1:ny
        }
        ## survey indices
        obsI <- vector("list", nsurv)
        timeI <- vector("list", nsurv)
        for(i in 1:nsurv){
            ## closest season
            tmp <- seasonStart[seasonStart < set$surveyTimes[i]]
            idxS <- which.min((tmp - set$surveyTimes[i])^2)
            ## correct NAA for survey timing dependent on seasons in opmod
            Z <- hist$FAA[,idx,idxS] + Ms
            surveyTime <- set$surveyTimes[i] - seasonStart[idxS]
            naa <- exp(log(hist$NAA[,idx,idxS]) - Z * surveyTime)
            esb <- apply(naa, 2, function(x) sum(x * dat$weightF * dat$sel))
            obsI[[i]] <- q[i] * esb * eI
            timeI[[i]] <- 1:ny + set$surveyTimes[i]
        }
    }

    ## output indices
    out <- hist
    obs <- list()
    obs$obsC <- obsC
    obs$timeC <- timeC
    obs$obsI <- obsI
    obs$timeI <- timeI
    ## obs$obsIsd <- sdI
    out$obs <- obs
    if("refs" %in% names(hist)) out$refs <- hist$refs
    return(out)
}
