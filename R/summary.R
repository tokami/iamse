#' Summarise an IAMSE object
#'
#' This is the `summary()` method for objects of class `"iamse"`. It provides a
#' compact overview of the main components of an IAMSE object, such as the
#' number of stocks, harvest control rules (HCRs), simulation years and
#' replicates, and possibly key performance metrics or reference points,
#' depending on what is stored in `object`.
#'
#' Typically, users will not call `summary.iamse()` directly, but will instead
#' use [summary()] on an `"iamse"` object:
#' \preformatted{
#'   summary(x)
#' }
#'
#' @param object An object of class `"iamse"`, usually returned by IAMSE
#'   helper functions (e.g. a wrapper around data, settings, and/or results).
#' @param ... Further arguments passed to or from other methods. Currently not
#'   used, but included for method compatibility.
#'
#' @return An object of class `"summary.iamse"` (or similar), which is
#'   typically printed to the console. The function is usually called for its
#'   side effects (printing a readable summary), and the return value is
#'   returned invisibly.
#'
#' @export
#' @method summary iamse
summary.iamse <- function(object, ...){
    ##
    quants <- c("TSB","SSB","ESB","CW","TSBfinal","SSBfinal","TACs","FM")
    nquants <- length(quants)
    nmse <- length(object)
    resList <- vector("list",nmse)

    for(i in 1:nmse){
        msei <- object[[i]]
        res <- vector("list",length(nquants))
        for(j in 1:nquants){
            quant <- quants[j]
            if(length(msei) == 0){
                res[[j]] <- NA
            }else{
                if(quant %in% c("TSB","SSB","ESB")){
                    tmp <- do.call(rbind,lapply(msei,
                                                function(x) apply(x[[quant]],1,mean)))
                }else if(quant %in% c("CW","FM")){
                    tmp <- do.call(rbind,lapply(msei,
                                                function(x) apply(x[[quant]],1,sum)))
                }else if(quant %in% c("TSBfinal","SSBfinal","TACs")){
                    tmp <- do.call(rbind,lapply(msei, function(x) x[[quant]]))
                }
                res[[j]] <- apply(tmp, 2, quantile, probs=c(0.025,0.5,0.975), na.rm=TRUE)
            }
        }
        names(res) <- quants
        resList[[i]] <- res
    }

    return(resList)
}
