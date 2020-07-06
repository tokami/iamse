#' @name estMets
#' @export
estMets <- function(mse, specdat, set){
    ##
    mets <- c("BBmsyFL","avCatch","avCatchLast5y","BBmsyLowest",
              "PBBlim","PBBlimFirst5y","PBBlimLast5y","CatchCV", "avRelCatch",
              "avRelCatchLast5y","converged","avCatchFirst5y",
              "avRelCatchFirst5y", "BBmsy5y","TimeRecov",
              "AAVC","AAVB")
    nmets <- length(mets)
    nmse <- length(mse)
    nrep <- set$nrep
    nysim <- set$nysim
    finalYear <- specdat$ny + set$nysim
    fifthYear <- specdat$ny + 5
    last5Years <- (finalYear - 4) : finalYear
    first5Years <- (specdat$ny + 1) : (specdat$ny + 5)
    simYears <- (specdat$ny + 1) : finalYear
    resList <- vector("list",nmse)
    ##
    for(i in 1:nmse){
        msei <- mse[[i]]
        res <- matrix(NA, nrow=nmets, ncol=3)
        ## "BBmsyFL"
        tmp <- unlist(lapply(msei, function(x) x$TSB[finalYear]/x$refs$Bmsy))
        res[1,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "avCatch"   ### also implement that for rel catch (to ref hcr)
        tmp <- unlist(lapply(msei, function(x) x$CW[simYears]))
        res[2,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "avCatchLast5y"
        tmp <- unlist(lapply(msei, function(x) x$CW[last5Years]))
        res[3,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "BBmsyLowest"
        tmp <- unlist(lapply(lapply(msei, function(x) x$TSB[simYears]/x$refs$Bmsy),
                             function(y) min(y, na.rm=TRUE)))
        res[4,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlim"
        tmp <- unlist(lapply(msei, function(x) mean(x$TSB[simYears]/(0.3*x$refs$Bmsy) < 1)))
        res[5,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlimFirst5y"
        tmp <- unlist(lapply(msei, function(x) mean(x$TSB[first5Years]/(0.5*x$refs$Bmsy) < 1))) ## HERE:
        res[6,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlimLast5y"
        tmp <- unlist(lapply(msei, function(x) mean(x$TSB[last5Years]/(0.3*x$refs$Bmsy) < 1)))
        res[7,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "CatchCV"
        tmp <- unlist(lapply(msei, function(x) sd(x$CW[simYears])/mean(x$CW[simYears])))
        res[8,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "avRelCatch"
        if("refFmsy" %in% set$hcr){
            tmpRef <- lapply(resMSE[[which(set$hcr == "refFmsy")]], function(x) x$CW[simYears])
            tmp0 <- lapply(msei, function(x) x$CW[simYears])
            tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
            res[9,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        }else{
            res[9,] <- rep(NA,3)
        }
        ## "avCRelatchLast5y"
        if("refFmsy" %in% set$hcr){
            tmpRef <- lapply(resMSE[[which(set$hcr == "refFmsy")]], function(x) x$CW[last5Years])
            tmp0 <- lapply(msei, function(x) x$CW[last5Years])
            tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
            res[10,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        }else{
            res[10,] <- rep(NA,3)
        }
        ## "converged"
        tmp <- unlist(lapply(msei, function(x)
            length(which(x$tacs$conv == FALSE))))
        res[11,] <- c(NA, (1 - (sum(tmp) / (nysim * nrep))) * 100, NA)
        ## "avCatchFirst5y"
        tmp <- unlist(lapply(msei, function(x) x$CW[first5Years]))
        res[12,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "avRelCatchFirst5y"
        if("refFmsy" %in% set$hcr){
            tmpRef <- lapply(resMSE[[which(set$hcr == "refFmsy")]], function(x) x$CW[first5Years])
            tmp0 <- lapply(msei, function(x) x$CW[first5Years])
            tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
            res[13,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        }else{
            res[13,] <- rep(NA,3)
        }
        ## "BBmsy5y"
        tmp <- unlist(lapply(msei, function(x) x$TSB[fifthYear]/x$refs$Bmsy))
        res[14,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "TimeRecov"
        tmp <- unlist(lapply(lapply(msei, function(x) x$TSB[simYears]/x$refs$Bmsy),
                             function(y) min(which(y > 1),na.rm=TRUE)))
        tmp[is.infinite(tmp)] <- NA  ## some sims might never recover...
        res[15,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE)
        ## "AAVC"
        tmp <- unlist(lapply(lapply(msei,
                                    function(x) (((x$CW[simYears] - x$CW[simYears+1])/
                                                 x$CW[simYears+1])^2)^0.5),
                             mean ,na.rm=TRUE))
        tmp[is.infinite(tmp)] <- NA
        res[16,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE)
        ## "AAVB"
        tmp <- unlist(lapply(lapply(msei,
                                    function(x) (((x$TSB[simYears] - x$TSB[simYears+1])/
                                                 x$TSB[simYears+1])^2)^0.5),
                             mean ,na.rm=TRUE))
        tmp[is.infinite(tmp)] <- NA
        res[17,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE)
        ## names
        rownames(res) <- mets
        colnames(res) <- c("2.5%","50%","97.5%")
        resList[[i]] <- res
    }
    names(resList) <- set$mss
    return(resList)
}
