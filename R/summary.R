

#' @name sumMSE
#'
#' @export
#'
sumMSE <- function(mse){
    ##
    quants <- c("TSB","SSB","ESB","CW","TACs","FM")
    nquants <- length(quants)
    nmse <- length(mse)
    resList <- vector("list",nmse)
    for(i in 1:nmse){
        msei <- mse[[i]]
        res <- vector("list",length(nquants))
        for(j in 1:nquants){
            quant <- quants[j]
            if(quant %in% c("TSB","SSB","ESB")){
                tmp <- do.call(rbind,lapply(msei, function(x) apply(x[[quant]],1,mean)))
            }else if(quant %in% c("CW","FM")){
                tmp <- do.call(rbind,lapply(msei, function(x) apply(x[[quant]],1,sum)))
            }else if(quant == "TACs"){
                tmp <- do.call(rbind,lapply(msei, function(x) x[[quant]]))
            }
            res[[j]] <- apply(tmp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=TRUE)
        }
        names(res) <- quants
        resList[[i]] <- res
    }
    return(resList)
}
