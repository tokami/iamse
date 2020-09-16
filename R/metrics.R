


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
    tacID <- names(mse)
    tacID2 <- sapply(strsplit(as.character(tacID), "-"), "[[", 1)

    res <- vector("list", nsamp)
    for(samp in 1:nsamp){
        sampi <- floor(sampleSizes[samp] * nreps)
        sampi[tacID2 %in% c("noF", "refFmsy")] <- nreps[tacID2 %in% c("noF", "refFmsy")] ## otherwise C/MSY not possible
        ind <- sapply(1:nhcr, function(x) sample(1:nreps[x], sampi[x]))
        msei <- lapply(1:nhcr, function(x) mse[[x]][ind[[x]]])
        names(msei) <- names(mse)
        names(mse[[1]])
        names(msei[[1]])
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
    first5Years <- (ny + 1) : (ny + 10)
    simYears <- (ny + 1) : finalYear
    last10Years <- (finalYear - 9) : finalYear
    first10Years <- (ny + 1) : (ny + 10)
    last15Years <- (finalYear - 14) : finalYear
    first15Years <- (ny + 1) : (ny + 15)

    hcrs <- names(mse)
    reffmsyInd <- which(hcrs == "refFmsy")
    refyield <- list(
        "simYears" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[simYears]),
        "first5Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[first5Years]),
        "last15Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[last15Years]),
        "last5Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[last5Years]),
        "last10Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[last10Years]))

    refF <- list(
        "finalYear" = lapply(mse[[reffmsyInd]], function(x) apply(x$FM,1,sum)[finalYear]))
    refB <- list(
        "finalYear" = lapply(mse[[reffmsyInd]], function(x) x$TSBfinal[finalYear]))

    metsAll <- c("CMSY","CMSYmean","avCatch",
                 "CMSYlast5","CMSYlast10",
                 "PBBlim","PBBlimlast5","PBBlimlast10",
                 "PBBlim1","PBBlim3",
                 "AAVC",
                 "CMSYST","PBBlimST",
                 "CMSYLT","PBBlimLT","AAVC2",
                 "CMSY2","CMSYST2","CMSYLT2",
                 "CMSYMaxAge","PBBlimMaxAge",
                 "BBmsy","BBmsyLT",
                 "BBmsyFL","FFmsyFL",
                 ## OLDER:
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
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[simYears] / refyield[[1]][[indi[x]]]))
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                tmp2 <- try(wilcox.test(as.numeric(tmp),
                                   alternative="two.sided",
                                   correct=TRUE,
                                   conf.int=TRUE,
                                   conf.level=0.95), silent=TRUE)
                if(hcrs[hcr] == "noF"){
                    res <- rbind(res, c(0,
                                        0,
                                        0,
                                        sei,
                                        ni))
                }else if(hcrs[hcr] == "refFmsy"){
                    res <- rbind(res, c(1,
                                        1,
                                        1,
                                        sei,
                                        ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSY could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlim3"
        if(any(mets == "PBBlim1")){
            ry <- vector("list",nysim)
            for(y in 1:nysim){
                yi <- ((ny+1):finalYear)[y]
                tmp <- sapply(msei, function(x) x$TSBfinal[yi] / refs$Blim < 1)
                tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
                ry[[y]] <- c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2])
            }
            tmp <- do.call(rbind,ry)
            res <- rbind(res, c(NA,tmp[nysim,],NA, NA, NA))
        }
        ## "PBBlim" # used for probHCR, but not averaged over years first!
        if(any(mets == "PBBlim")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[simYears] / refs$Blim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "PBBlim3"
        if(any(mets == "PBBlim3")){
            ry <- vector("list",nysim)
            for(y in 1:nysim){
                yi <- ((ny+1):finalYear)[y]
                tmp <- sapply(msei, function(x) x$TSBfinal[yi] / refs$Blim < 1)
                tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
                ry[[y]] <- c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2])
            }
            tmp <- do.call(rbind,ry)
            res <- rbind(res, c(NA,max(tmp[,2]),NA, NA, NA))
        }
        ## "AAVC"
        if(any(mets == "AAVC")){
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[simYears] -
                                                     apply(x$CW,1,sum)[simYears-1]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[simYears], na.rm=TRUE)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp2 <- try(wilcox.test(as.numeric(tmp),
                                    alternative="two.sided",
                                    correct=TRUE,
                                    conf.int=TRUE,
                                    conf.level=0.95), silent=TRUE)
            if(hcrs[hcr] == "noF"){
                res <- rbind(res, c(0,
                                    0,
                                    0,
                                    sei,
                                    ni))
            }else{
                res <- rbind(res, c(tmp2$conf.int[1],
                                    tmp2$estimate,
                                    tmp2$conf.int[2],
                                    sei,
                                    ni))
            }
        }
        ## CMSYST
        if(any(mets == "CMSYST")){
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[first5Years] / refyield[[2]][[indi[x]]]))
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                tmp2 <- try(wilcox.test(as.numeric(tmp),
                                   alternative="two.sided",
                                   correct=TRUE,
                                   conf.int=TRUE,
                                   conf.level=0.95), silent=TRUE)
                if(hcrs[hcr] == "noF"){
                    res <- rbind(res, c(0,
                                        0,
                                        0,
                                        sei,
                                        ni))
                }else if(hcrs[hcr] == "refFmsy"){
                    res <- rbind(res, c(1,
                                        1,
                                        1,
                                        sei,
                                        ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYST could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlimST"
        if(any(mets == "PBBlimST")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[first5Years] / refs$Blim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## CMSYLT
        if(any(mets == "CMSYLT")){
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[last15Years] / refyield[[3]][[indi[x]]]))
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                tmp2 <- try(wilcox.test(as.numeric(tmp),
                                   alternative="two.sided",
                                   correct=TRUE,
                                   conf.int=TRUE,
                                   conf.level=0.95), silent=TRUE)
                if(hcrs[hcr] == "noF"){
                    res <- rbind(res, c(0,
                                        0,
                                        0,
                                        sei,
                                        ni))
                }else if(hcrs[hcr] == "refFmsy"){
                    res <- rbind(res, c(1,
                                        1,
                                        1,
                                        sei,
                                        ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlimLT"
        if(any(mets == "PBBlimLT")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last15Years] / refs$Blim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## CMSYlast5
        if(any(mets == "CMSYlast5")){
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[last5Years] / refyield[["last5Years"]][[indi[x]]]))
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                tmp2 <- try(wilcox.test(as.numeric(tmp),
                                   alternative="two.sided",
                                   correct=TRUE,
                                   conf.int=TRUE,
                                   conf.level=0.95), silent=TRUE)
                if(hcrs[hcr] == "noF"){
                    res <- rbind(res, c(0,
                                        0,
                                        0,
                                        sei,
                                        ni))
                }else if(hcrs[hcr] == "refFmsy"){
                    res <- rbind(res, c(1,
                                        1,
                                        1,
                                        sei,
                                        ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## CMSYlast10
        if(any(mets == "CMSYlast10")){
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[last10Years] / refyield[["last10Years"]][[indi[x]]]))
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                tmp2 <- try(wilcox.test(as.numeric(tmp),
                                   alternative="two.sided",
                                   correct=TRUE,
                                   conf.int=TRUE,
                                   conf.level=0.95), silent=TRUE)
                if(hcrs[hcr] == "noF"){
                    res <- rbind(res, c(0,
                                        0,
                                        0,
                                        sei,
                                        ni))
                }else if(hcrs[hcr] == "refFmsy"){
                    res <- rbind(res, c(1,
                                        1,
                                        1,
                                        sei,
                                        ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
            print(                                        tmp2$estimate)
        }
        ## "PBBlimlast5"
        if(any(mets == "PBBlimlast5")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last5Years] / refs$Blim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "PBBlimlast10"
        if(any(mets == "PBBlimlast10")){
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last10Years] / refs$Blim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "AAVC2"
        if(any(mets == "AAVC2")){
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[simYears] -
                                                     apply(x$CW,1,sum)[simYears-1]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[simYears], na.rm=TRUE)))
            mu <- mean(tmp)
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            if(hcrs[hcr] == "noF"){
                res <- rbind(res, c(0,
                                    0,
                                    0,
                                    sei,
                                    ni))
            }else{
                res <- rbind(res, c(mu - qnorm(1-0.05/2)*sei,
                                    mu,
                                    mu + qnorm(1-0.05/2)*sei,
                                    sei,
                                    ni))
            }
        }
        ## CMSY2 (mean)
        if(any(mets == "CMSY2")){
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x) mean(apply(msei[[x]]$CW,1,sum)[simYears] /
                                                               refyield[[1]][[indi[x]]]))
                mu <- mean(tmp,na.rm=TRUE)
                vari <- var(tmp,na.rm=TRUE)
                sei <- sqrt(vari/length(tmp))
                res <- rbind(res, c(mu - qt(0.975,df=length(tmp)-1)*sei,
                                    mu,
                                    mu + qt(0.975,df=length(tmp)-1)*sei,
                                    sei,
                                    length(tmp)))
            }else writeLines("CMSY2 could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## CMSYST2 (mean)
        if(any(mets == "CMSYST2")){
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x) mean(apply(msei[[x]]$CW,1,sum)[first5Years] /
                                                               refyield[[2]][[indi[x]]]))
                mu <- mean(tmp,na.rm=TRUE)
                vari <- var(tmp,na.rm=TRUE)
                sei <- sqrt(vari/length(tmp))
                res <- rbind(res, c(mu - qt(0.975,df=length(tmp)-1)*sei,
                                    mu,
                                    mu + qt(0.975,df=length(tmp)-1)*sei,
                                    sei,
                                    length(tmp)))
            }else writeLines("CMSYST2 could not be estimated, because no rule 'refFmsy' not found.")
        }
        if(any(mets == "CMSYLT2")){
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x) mean(apply(msei[[x]]$CW,1,sum)[last15Years] /
                                                               refyield[[3]][[indi[x]]]))
                mu <- mean(tmp,na.rm=TRUE)
                vari <- var(tmp,na.rm=TRUE)
                sei <- sqrt(vari/length(tmp))
                res <- rbind(res, c(mu - qt(0.975,df=length(tmp)-1)*sei,
                                    mu,
                                    mu + qt(0.975,df=length(tmp)-1)*sei,
                                    sei,
                                    length(tmp)))
            }else writeLines("CMSYLT2 could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## OLDER (not used in probHCR):
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
        if(any(mets == "CMSYold")){
            ## tmp <- unlist(lapply(msei, function(x) mean(apply(x$CW,1,sum)[last5Years] / refs$MSY)))## HERE:
            ## res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
            tmp <- unlist(lapply(msei, function(x) median(apply(x$CW,1,sum)[simYears] / refs$MSY)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## OLDER:
        ## "BBmsyFL"
        if(any(mets == "BBmsyFL")){
            indi <- as.numeric(names(msei))
            tmp <- sapply(1:length(msei), function(x) msei[[x]]$TSBfinal[finalYear]/
                                                      refB[["finalYear"]][[indi[x]]]) ##refs$Bmsy)
            res <- rbind(res, c(quantile(tmp, probs = c(0.25, 0.5, 0.75), na.rm=TRUE),NA,NA))
        }
        if(any(mets == "FFmsyFL")){
            indi <- as.numeric(names(msei))
            tmp <- sapply(1:length(msei), function(x) apply(msei[[x]]$FM,1,sum)[finalYear]/
                                                      refF[["finalYear"]][[indi[x]]]) ##refs$Fmsy)
            res <- rbind(res, c(quantile(tmp, probs = c(0.25, 0.5, 0.75), na.rm=TRUE),NA,NA))
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
        colnames(res) <- c("loCI","mu","upCI","se","n")
        res2[[hcr]] <- round(res, 3)
    }
    names(res2) <- hcrs


    class(res2) <- "met"


    return(res2)
}
