## observation model


## make observations with given CV for catch and index (only one index time series so far)
#' @name obsmod
#' @export
obsmod <- function(specdat, hist, set, years = NULL){
    ## parameters
    ESB <- hist$ESB
    CW <- hist$CW
    q <- specdat$q
    ## indices
    if(is.null(years)){
        idx <- 1:length(CW)
    }else{
        idx <- years
    }
    ny <- length(idx)
    if("obsC" %in% names(hist$obs)) ny <- 1
    ## errors
    eC <- set$eC
    eI <- set$eI
    if(is.null(eC)) eC <- rnorm(ny, 0, set$CVC) - set$CVC^2/2
    if(is.null(eI)) eI <- rnorm(ny, 0, set$CVI) - set$CVI^2/2

    nsurv <- length(set$surveyTimes)
    if(length(q) < nsurv) q <- rep(q, nsurv)

    ## adding observation
    if("obsC" %in% names(hist$obs)){
        y <- length(hist$ESB)
        tmp1 <- length(hist$obs$obsC)
        tmp <- y - tmp1
        if(tmp > 0){
            ## add catch obs
            obsC <- c(hist$obs$obsC, exp(log(CW[y]) + eC[y-specdat$ny]))
            timeC <- c(hist$obs$timeC, tail(hist$obs$timeC,1)+1)
            ## add survey obs
            Z <- hist$FAA[y,] + specdat$M
            obsI <- vector("list", nsurv)
            timeI <- vector("list", nsurv)
            for(i in 1:nsurv){
                naa <- exp(log(hist$NAA[y,]) - Z * set$surveyTimes[i])
                if(inherits(naa, "matrix")){
                    esb <- apply(naa, 1, function(x) sum(x * specdat$weightF * specdat$sel))
                }else{
                    esb <- sum(naa * specdat$weightF * specdat$sel)
                }
                obsI[[i]] <- c(hist$obs$obsI[[i]], exp(log(q[i]) + log(esb) + eI[y-specdat$ny]))
                timeI[[i]] <- c(hist$obs$timeI[[i]], floor(tail(hist$obs$timeI[[i]],1)) + 1 + set$surveyTimes[i])
                ## sdI[[1]] <- c(hist$obs$obsIsd[[i]], set$CVI)
            }
        }else return(hist)
    ## initalizing observations
    }else{
        ## catch
        obsC <- exp(log(CW[idx]) + eC)
        timeC <- 1:ny
        ## survey indices
        obsI <- vector("list", nsurv)
        timeI <- vector("list", nsurv)
        Z <- hist$FAA[idx,] + specdat$M
        for(i in 1:nsurv){
            naa <- exp(log(hist$NAA[idx,]) - Z * set$surveyTimes[i])
            esb <- apply(naa, 1, function(x) sum(x * specdat$weightF * specdat$sel))
            obsI[[i]] <- exp(log(q[i]) + log(esb) + eI)
            timeI[[i]] <- 1:ny + set$surveyTimes[i]
            ## sdI[[i]] <- rep(set$CVI, length(obsI1))
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
