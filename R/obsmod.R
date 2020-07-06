## observation model


## make observations with given CV for catch and index (only one index time series so far)
#' @name obsmod
#' @export
obsmod <- function(specdat, hist, set, years = NULL){
    ## parameters
    TSB <- hist$TSB
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

    ## adding observation
    if("obsC" %in% names(hist$obs)){
        y <- length(hist$TSB)
        tmp1 <- length(hist$obs$obsC)
        tmp <- y - tmp1
        if(tmp > 0){
            obsC <- c(hist$obs$obsC, exp(log(CW[y]) + eC[y-specdat$ny]))
            obsI <- list()
            obsI[[1]] <- c(hist$obs$obsI[[1]], exp(log(q[1]) + log(TSB[y]) + eI[y-specdat$ny]))
            timeC <- c(hist$obs$timeC, tail(hist$obs$timeC,1)+1)
            timeI <- list()
            timeI[[1]] <- c(hist$obs$timeI[[1]], tail(hist$obs$timeI[[1]],1)+1)
            sdI <- list()
            sdI[[1]] <- c(hist$obs$obsIsd[[1]], set$CVI)
        }else return(hist)
    ## initalizing observations
    }else{
        ## catch
        obsC <- exp(log(CW[idx]) + eC)
        timeC <- 1:ny
        ## index time 1
        obsI1 <- exp(log(q[1]) + log(TSB[idx]) + eI)
        obsI <- list(obsI1)
        timeI1 <- 1:ny
        timeI <- list(timeI1)
        sdI1 <- rep(set$CVI, length(obsI1))
        sdI <- list(sdI1)
    }

    ## output indices
    out <- hist
    obs <- list()
    obs$obsC <- obsC
    obs$timeC <- timeC
    obs$obsI <- obsI
    obs$timeI <- timeI
    obs$obsIsd <- sdI
    out$obs <- obs
    if("refs" %in% names(hist)) out$refs <- hist$refs
    return(out)
}
