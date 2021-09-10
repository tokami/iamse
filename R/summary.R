

#' @name summary.mse
#'
#' @export
#'
summary.mse <- function(mse, yearly = TRUE){

    ##
    quants <- c("TSB","SSB","ESB","CW","TSBfinal","TACs","FM","TSBfinalSea")
    nquants <- length(quants)
    nmse <- length(mse)
    resList <- vector("list",nmse)
    for(i in 1:nmse){
        msei <- mse[[i]]
        res <- vector("list",length(nquants))
        for(j in 1:nquants){
            quant <- quants[j]
            if(quant %in% c("TSB","SSB","ESB","TSBfinalSea")){
                if(yearly && inherits(msei[[1]][["FM"]], "matrix")){
                    tmp <- do.call(rbind,lapply(msei, function(x) apply(x[[quant]],1,mean)))
                }else{
                    tmp <- do.call(rbind,lapply(msei, function(x) as.vector(t(x[[quant]]))))
                }
            }else if(quant %in% c("CW","FM")){
                if(yearly && inherits(msei[[1]][["FM"]], "matrix")){
                    tmp <- do.call(rbind,lapply(msei, function(x) apply(x[[quant]],1,sum)))
                }else{
                    tmp <- do.call(rbind,lapply(msei, function(x) as.vector(t(x[[quant]]))))
                }
            }else if(quant %in% c("TSBfinal","TACs")){
                ## CHECK: possible to extract sub-yearly?
                tmp <- do.call(rbind,lapply(msei, function(x) x[[quant]]))
            }
            res[[j]] <- apply(tmp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=TRUE)
        }
        names(res) <- quants
        resList[[i]] <- res
    }

    return(resList)
}
