


#' @name estConsMets
#'
#' @description Estimate consistency of metrics across different sample sizes
#'
#' @export
estConsMets <- function(mse, dat, mets = "all",
                        sampleSizes = c(0.1,0.25,0.5,0.75,1)){

    nsamp <- length(sampleSizes)
    nhcr <- length(mse)
    ## nreps <- sapply(mse, length)
    nreps <- sapply(mse, length)

    res <- vector("list", nsamp)
    for(samp in 1:nsamp){
        sampi <- floor(sampleSizes[samp] * nreps)
        ind <- sapply(1:nhcr, function(x) sample(1:nreps[x], sampi[x]))
        msei <- lapply(1:nhcr, function(x) mse[[x]][ind[[x]]])
        res[[samp]] <- estMets(msei, dat = dat, mets = mets)
    }

    return(res)
}




#' @name estMets
#' @export
estMets <- function(mse, dat, mets = "all"){

    nhcr <- length(mse)
    nrep <- length(mse[[1]])
    nquant <- length(mse[[1]][[1]])
    nysim <- nrow(mse[[1]][[1]]$tacs)
    dims <- dim(mse[[1]][[1]]$CW)
    ny <- dims[1] - nysim
    ns <- dims[2]

    finalYear <- ny + nysim
    fifthYear <- ny + 5
    last5Years <- (finalYear - 4) : finalYear
    first5Years <- (ny + 1) : (ny + 5)
    simYears <- (ny + 1) : finalYear
    last10Years <- (finalYear - 9) : finalYear
    first10Years <- (ny + 1) : (ny + 10)
    last15Years <- (finalYear - 14) : finalYear
    first15Years <- (ny + 1) : (ny + 15)

    hcrs <- names(mse)
    reffmsyInd <- which(hcrs == "refFmsy")
    refyield <- lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[simYears]) ## HERE:

    metsAll <- c("CMSY","CMSYmean","PBBlim","AAVC",
                 "CMSYST","PBBlimST",
                 "CMSYLT","PBBlimLT",
                 "CMSYMaxAge","PBBlimMaxAge",
                 "BBmsy","BBmsyLT",
                 ## OLDER:
                 "BBmsyFL","avCatch",
                 "avCatchFirst5y","avCatchLast5y","BBmsyLowest",
                 "PBBlim2",
                 "PBBlimFirst5y","PBBlimLast5y","CatchCV", "avRelCatch",
                 "avRelCatchLast5y","converged",
                 "avRelCatchFirst5y", "BBmsy5y","TimeRecov",
                 "AAVCFirst5y","AAVB")
    if(mets[1] == "all") mets <- metsAll
    if(any(which(!mets %in% metsAll))) writeLines(paste0("Metric ",
                                                         paste0(mets[which(!mets %in% metsAll)], collapse=", "),
                                                                " not defined and thus omitted."))
    mets <- metsAll[which(metsAll %in% mets)]
    nmets <- length(mets)

    ## TODO: make part of mse
    refs <- dat$ref
    if(is.null(refs)) stop("Reference levels missing.")
    ##
    res2 <- vector("list", nhcr)
    for(hcr in 1:nhcr){
        nrep <- length(mse[[hcr]])
        msei <- mse[[hcr]]
        res <- NULL
        ## CMSY
        if(any(mets == "CMSY")){
            ## tmp <- unlist(lapply(msei, function(x) mean(apply(x$CW,1,sum)[last5Years] / refs$MSY)))## HERE:
            ## res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
            tmp <- unlist(lapply(msei, function(x) median(apply(x$CW,1,sum)[simYears] / refs$MSY)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        if(any(mets == "CMSYmean")){
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x) mean(apply(msei[[x]]$CW,1,sum)[simYears] / refyield[[indi[x]]]))
##                res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
                res <- rbind(res, c(NA, mean(tmp,na.rm=TRUE), NA))
            }else writeLines("CMSYmean could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlim"
        if(any(mets == "PBBlim")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[simYears] / refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ## "AAVC"
        if(any(mets == "AAVC")){
            tmp <- unlist(lapply(lapply(msei,
                                        function(x) (((apply(x$CW,1,sum)[simYears] -
                                                       apply(x$CW,1,sum)[simYears+1])/
                                                      apply(x$CW,1,sum)[simYears+1])^2)^0.5),
                                 mean ,na.rm=TRUE))
            tmp[is.infinite(tmp)] <- NA
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE))
        }
        ## CMSYST
        if(any(mets == "CMSYST")){
            tmp <- unlist(lapply(msei, function(x) median(apply(x$CW,1,sum)[first5Years] / refs$MSY)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "PBBlimST"
        if(any(mets == "PBBlimST")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[first5Years] / refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ## CMSYLT
        if(any(mets == "CMSYLT")){
            tmp <- unlist(lapply(msei, function(x) median(apply(x$CW,1,sum)[last15Years] / refs$MSY)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "PBBlimLT"
        if(any(mets == "PBBlimLT")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last15Years] / refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ## CMSYMaxAge
        if(any(mets == "CMSYMaxAge")){
            tmp <- unlist(lapply(msei, function(x) median(apply(x$CW,1,sum)[min(simYears):(min(simYears)+dat$amax)] /
                                                   refs$MSY)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "PBBlimMaxAge"
        if(any(mets == "PBBlimMaxAge")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[min(simYears):(min(simYears)+dat$amax)]/
                                                        refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ## "BBmsy"
        if(any(mets == "BBmsy")){
            tmp <- unlist(lapply(msei, function(x) median(x$TSBfinal[simYears] / refs$Bmsy)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "BBmsyLT"
        if(any(mets == "BBmsyLT")){
            tmp <- unlist(lapply(msei, function(x) median(x$TSBfinal[last5Years] / refs$Bmsy)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## OLDER:
        ## "BBmsyFL"
        if(any(mets == "BBmsyFL")){
            tmp <- unlist(lapply(msei, function(x) x$TSBfinal[finalYear]/refs$Bmsy))
            res <- rbind(res, c(NA,mean(tmp),NA)) ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        }
        ## "avCatch"   ### also implement that for rel catch (to ref hcr)
        if(any(mets == "avCatch")){
            tmp <- unlist(lapply(msei, function(x) apply(x$CW,1,sum)[simYears]))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "avCatchFirst5y"
        if(any(mets == "avCatchFirst5y")){
            tmp <- unlist(lapply(msei, function(x) apply(x$CW,1,sum)[first5Years]))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "avCatchLast5y"
        if(any(mets == "avCatchLast5y")){
            tmp <- unlist(lapply(msei, function(x) apply(x$CW,1,sum)[last5Years]))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "BBmsyLowest"
        if(any(mets == "BBmsyLowest")){
            tmp <- unlist(lapply(lapply(msei, function(x) x$TSBfinal[simYears]/refs$Bmsy),
                                 function(y) min(y, na.rm=TRUE)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## ## "PBBlim"
        ## if(any(mets == "PBBlim")){
        ##     tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[simYears]/refs$Blim < 1)))
        ##     tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
        ##     res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        ## }
        ## "PBBlim2"
        if(any(mets == "PBBlim2")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSB[simYears,1]/refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlimFirst5y"
        if(any(mets == "PBBlimFirst5y")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[first5Years]/refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ##        res[6,] <- c(sd(tmp),mean(tmp),NA) ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlimLast5y"
        if(any(mets == "PBBlimLast5y")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last5Years]/refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ##        res[7,] <- c(sd(tmp),mean(tmp),NA) ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "CatchCV"
        if(any(mets == "CatchCV")){
            tmp <- unlist(lapply(msei, function(x) sd(apply(x$CW,1,sum)[simYears])/mean(x$CW[simYears])))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "avRelCatch"
        if(any(mets == "avRelCatch")){
            if("refFmsy" %in% hcrs){
                tmpRef <- lapply(resMSE[[which(hcrs == "refFmsy")]], function(x) apply(x$CW,1,mean)[simYears])
                tmp0 <- lapply(msei, function(x) apply(x$CW,1,sum)[simYears])
                tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
                res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
            }else{
                res <- rbind(res, rep(NA,3))
            }
        }
        ## "avCRelatchLast5y"
        if(any(mets == "avCRelatchLast5y")){
            if("refFmsy" %in% hcrs){
                tmpRef <- lapply(resMSE[[which(hcrs == "refFmsy")]], function(x) apply(x$CW,1,mean)[last5Years])
                tmp0 <- lapply(msei, function(x) apply(x$CW,1,sum)[last5Years])
                tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
                res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
            }else{
                res <- rbind(res, rep(NA,3))
            }
        }
        ## "converged"
        if(any(mets == "converged")){
            tmp <- unlist(lapply(msei, function(x)
                length(which(x$tacs$conv == TRUE))
                ))
            res <- rbind(res, c(NA, mean(tmp/nysim) * 100, NA))
        }
        ## "avRelCatchFirst5y"
        if(any(mets == "avRelCatchFirst5y")){
            if("refFmsy" %in% hcrs){
                tmpRef <- lapply(resMSE[[which(hcrs == "refFmsy")]], function(x) apply(x$CW,1,sum)[first5Years])
                tmp0 <- lapply(msei, function(x) apply(x$CW,1,mean)[first5Years])
                tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
                res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
            }else{
                res <- rbind(res, rep(NA,3))
            }
        }
        ## "BBmsy5y"
        if(any(mets == "BBmsy5y")){
            tmp <- unlist(lapply(msei, function(x) x$TSBfinal[fifthYear]/refs$Bmsy))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "TimeRecov"
        if(any(mets == "TimeRecov")){
            tmp <- unlist(lapply(lapply(msei, function(x) x$TSBfinal[simYears]/refs$Bmsy),
                                 function(y) min(which(y > 1),na.rm=TRUE)))
            tmp[is.infinite(tmp)] <- NA  ## some sims might never recover...
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE))
        }
        if(any(mets == "AAVCFirst5y")){
            tmp <- unlist(lapply(lapply(msei,
                                        function(x) (((apply(x$CW,1,sum)[first5Years] -
                                                       apply(x$CW,1,sum)[first5Years+1])/
                                                      apply(x$CW,1,sum)[first5Years+1])^2)^0.5),
                                 mean ,na.rm=TRUE))
            tmp[is.infinite(tmp)] <- NA
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE))
        }
        ## "AAVB"
        if(any(mets == "AAVB")){
            tmp <- unlist(lapply(lapply(msei,
                                        function(x) (((x$TSBfinal[simYears] -
                                                       x$TSBfinal[simYears+1])/
                                                      x$TSBfinal[simYears+1])^2)^0.5),
                                 mean ,na.rm=TRUE))
            tmp[is.infinite(tmp)] <- NA
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975),na.rm=TRUE))
        }

        ## names
        rownames(res) <- mets
        colnames(res) <- c("2.5%","50%","97.5%")
        res2[[hcr]] <- round(res, 3)
    }
    names(res2) <- hcrs


    class(res2) <- "met"


    return(res2)
}
