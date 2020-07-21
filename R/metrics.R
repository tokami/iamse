#' @name estMets
#' @export
estMets <- function(mse, dat){

    refs <- dat$ref
    ##
    mets <- c("BBmsyFL","avCatch","avCatchLast5y","BBmsyLowest",
              "PBBlim","PBBlimFirst5y","PBBlimLast5y","CatchCV", "avRelCatch",
              "avRelCatchLast5y","converged","avCatchFirst5y",
              "avRelCatchFirst5y", "BBmsy5y","TimeRecov",
              "AAVC","AAVB")
    nmets <- length(mets)
    hcrs <- names(mse)
    nhcr <- length(mse)
    nrep <- length(mse[[1]])
    nysim <- nrow(mse[[1]][[1]]$tacs)
    dims <- dim(mse[[1]][[1]]$TSB)
    ny <- dims[1] - nysim
    ns <- dims[2]

    finalYear <- ny + nysim
    fifthYear <- ny + 5
    last5Years <- (finalYear - 4) : finalYear
    first5Years <- (ny + 1) : (ny + 5)
    simYears <- (ny + 1) : finalYear

    resList <- vector("list",nhcr)
    ##
    for(i in 1:nhcr){
        msei <- mse[[i]]
        res <- matrix(NA, nrow=nmets, ncol=3)
        ## "BBmsyFL"
        tmp <- unlist(lapply(msei, function(x) x$TSBfinal[finalYear]/refs$Bmsy))
        res[1,] <- c(NA,mean(tmp),NA) ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "avCatch"   ### also implement that for rel catch (to ref hcr)
        tmp <- unlist(lapply(msei, function(x) apply(x$CW,1,sum)[simYears]))
        res[2,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "avCatchLast5y"
        tmp <- unlist(lapply(msei, function(x) apply(x$CW,1,sum)[last5Years]))
        res[3,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "BBmsyLowest"
        tmp <- unlist(lapply(lapply(msei, function(x) x$TSBfinal[simYears]/refs$Bmsy),
                             function(y) min(y, na.rm=TRUE)))
        res[4,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlim"
        tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[simYears]/(0.15*refs$B0) < 1)))
        tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
        res[5,] <- c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2])
        ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlimFirst5y"
        tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[first5Years]/(0.15*refs$B0) < 1)))
        tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
        res[6,] <- c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2])
##        res[6,] <- c(sd(tmp),mean(tmp),NA) ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlimLast5y"
        tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last5Years]/(0.15*refs$B0) < 1)))
        tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
        res[7,] <- c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2])
##        res[7,] <- c(sd(tmp),mean(tmp),NA) ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "CatchCV"
        tmp <- unlist(lapply(msei, function(x) sd(apply(x$CW,1,sum)[simYears])/mean(x$CW[simYears])))
        res[8,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "avRelCatch"
        if("refFmsy" %in% hcrs){
            tmpRef <- lapply(resMSE[[which(hcrs == "refFmsy")]], function(x) apply(x$CW,1,mean)[simYears])
            tmp0 <- lapply(msei, function(x) apply(x$CW,1,sum)[simYears])
            tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
            res[9,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        }else{
            res[9,] <- rep(NA,3)
        }
        ## "avCRelatchLast5y"
        if("refFmsy" %in% hcrs){
            tmpRef <- lapply(resMSE[[which(hcrs == "refFmsy")]], function(x) apply(x$CW,1,mean)[last5Years])
            tmp0 <- lapply(msei, function(x) apply(x$CW,1,sum)[last5Years])
            tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
            res[10,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        }else{
            res[10,] <- rep(NA,3)
        }
        ## "converged"
        tmp <- unlist(lapply(msei, function(x)
            length(which(x$tacs$conv == TRUE))
##            length(which(x$tacs$id != "cc"))
            ))
        res[11,] <- c(NA, mean(tmp/nysim) * 100, NA)
        ## "avCatchFirst5y"
        tmp <- unlist(lapply(msei, function(x) apply(x$CW,1,sum)[first5Years]))
        res[12,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "avRelCatchFirst5y"
        if("refFmsy" %in% hcrs){
            tmpRef <- lapply(resMSE[[which(hcrs == "refFmsy")]], function(x) apply(x$CW,1,sum)[first5Years])
            tmp0 <- lapply(msei, function(x) apply(x$CW,1,mean)[first5Years])
            tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
            res[13,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        }else{
            res[13,] <- rep(NA,3)
        }
        ## "BBmsy5y"
        tmp <- unlist(lapply(msei, function(x) x$TSBfinal[fifthYear]/refs$Bmsy))
        res[14,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "TimeRecov"
        tmp <- unlist(lapply(lapply(msei, function(x) x$TSBfinal[simYears]/refs$Bmsy),
                             function(y) min(which(y > 1),na.rm=TRUE)))
        tmp[is.infinite(tmp)] <- NA  ## some sims might never recover...
        res[15,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE)
        ## "AAVC"
        tmp <- unlist(lapply(lapply(msei,
                                    function(x) (((apply(x$CW,1,sum)[simYears] -
                                                   apply(x$CW,1,sum)[simYears+1])/
                                                 apply(x$CW,1,sum)[simYears+1])^2)^0.5),
                             mean ,na.rm=TRUE))
        tmp[is.infinite(tmp)] <- NA
        res[16,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE)
        ## "AAVB"
        tmp <- unlist(lapply(lapply(msei,
                                    function(x) (((x$TSBfinal[simYears] -
                                                   x$TSBfinal[simYears+1])/
                                                 x$TSBfinal[simYears+1])^2)^0.5),
                             mean ,na.rm=TRUE))
        tmp[is.infinite(tmp)] <- NA
        res[17,] <- quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE)
        ## names
        rownames(res) <- mets
        colnames(res) <- c("2.5%","50%","97.5%")
        resList[[i]] <- round(res,2)
    }
    names(resList) <- hcrs
    return(resList)
}
