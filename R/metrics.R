


#' est.cons.mets
#'
#' @description Estimate consistency of metrics across different sample sizes
#'
#' @export
est.cons.mets <- function(mse, dat, mets = "all",
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
        res[[samp]] <- est.metrics(msei, dat = dat, mets = mets)
    }

    return(res)
}




#' est.metrics
#' @export
est.metrics <- function(mse, dat, set = NULL,
                        mets = "all",
                        nysim = NULL){

    nysim0 <- nysim


    nhcr <- length(mse)
    nrep <- length(mse[[1]])
    nquant <- length(mse[[1]][[1]])
    ## IMPROVE:
    if(!is.null(set)){
        nysim <- set$nysim
    }else if(any(names(mse[[1]][[1]]$tacs) == "assessInt") &&
       any(!names(mse) %in% c("noF","refFmsy")) &&
       !is.na(mse[[which(!names(mse) %in%
                                c("noF","refFmsy"))[1]]][[1]]$tacs$assessInt[1])){
        assessInt <- mse[[which(!names(mse) %in%
                                c("noF","refFmsy"))[1]]][[1]]$tacs$assessInt[1]
        nysim <- nrow(mse[[1]][[1]]$tacs) * assessInt
    }else{
        nysim <- nrow(mse[[1]][[1]]$tacs)
    }

    if(!is.null(nysim0)) nysim <- nysim0
    dims <- dim(mse[[1]][[1]]$CW)
    ny <- dat$ny
    ns <- dat$ns
    ## ny <- dims[1] - nysim
    ## ns <- dims[2]
    amax <- dat$amax

    finalYear <- ny + nysim
    fifthYear <- ny + 5
    tenthYear <- ny + 10
    last5Years <- (finalYear - 4) : finalYear
    first5Years <- (ny + 1) : (ny + 5)
    simYears <- (ny + 1) : finalYear
    last10Years <- (finalYear - 9) : finalYear
    first10Years <- (ny + 1) : (ny + 10)
    last15Years <- (finalYear - 14) : finalYear
    last35Years <- (finalYear - 34) : finalYear
    first15Years <- (ny + 1) : (ny + 15)
    last5AmaxYears <- ifelse(amax>5,amax-4,1):amax + ny
    amaxYears <- 1:amax + ny
    amaxYears2 <- amax + ny
    amaxYears3 <- 1:amax + ny ## floor(amax/2):amax + ny
    amaxYears4 <- floor(amax/2) + ny
    ## for stand. metrics
    newYears <- 1:nysim + ny
    newYears_st <- 1:5 + ny
    newYears_lt <- (nysim-10):nysim + ny
    ## newYears_lt <- (nysim-24):nysim + ny ## HERE:
    ## newYears_lt <- (nysim-4):nysim + ny     ## HERE:
    newYears_y1 <- ny + 1
    newYears_y2 <- ny + 2
    ## for kobe plot
    newYears2 <- round(nysim - nysim/2) + ny
    newYears2_st <- 1 + ny
    newYears2_lt <- nysim + ny
    newYears2_mt <- round(nysim - nysim/2) + ny  + (-2:2)
    amaxYears5 <- 3:amax + ny

    hcrs <- names(mse)
    reffmsyInd <- which(hcrs == "refFmsy")
    if(length(reffmsyInd) > 0){
        refyield <- list(
            "simYears" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[simYears]),
            "first5Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[first5Years]),
            "first10Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[first10Years]),
            "last15Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[last15Years]),
            "last35Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[last35Years]),
            "last5Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[last5Years]),
            "last10Years" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[last10Years]),
            "last5AmaxYears" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[last5AmaxYears]),
            "amaxYears" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[amaxYears]),
            "newYears" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[newYears]),
            "newYears_st" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[newYears_st]),
            "newYears_lt" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[newYears_lt]),
            "amaxYears5" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[amaxYears5]),
            "newYears_y1" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[newYears_y1]),   "newYears_y2" = lapply(mse[[reffmsyInd]], function(x) apply(x$CW,1,sum)[newYears_y2])
        )
    refF <- list(
        "finalYear" = lapply(mse[[reffmsyInd]], function(x) apply(x$FM,1,sum)[finalYear]),
        "tenthYear" = lapply(mse[[reffmsyInd]], function(x) apply(x$FM,1,sum)[tenthYear]))
    refB <- list(
        "finalYear" = lapply(mse[[reffmsyInd]], function(x) x$TSBfinal[finalYear]),
        "newYears" = lapply(mse[[reffmsyInd]], function(x) x$ESBfinal[newYears]),
        "tenthYear" = lapply(mse[[reffmsyInd]], function(x) x$TSBfinal[tenthYear]))

    }


    metsAll <- c("CMSY",
                 "PBBlim",
                 "PBBlimSSB",
                 "AAVC",
                 "AAVB",
                 "CMSYmean","avCatch",
                 "CMSYfirst10",
                 "CMSYlast5","CMSYlast10",
                 "PBBlimfirst10","PBBlimfirst10SSB",
                 "PBBlimlast5","PBBlimlast10","PBBlimlast10SSB",
                 "PBBlim1","PBBlim3",
                 "CMSYST","PBBlimST","PBBlim1ST","PBBlim3ST",
                 "CMSYLT","PBBlimLT","PBBlim1LT","PBBlim3LT",
                 "CMSYamax","PBBlimamax","AAVCamax",
                 "CMSYamax2","PBBlimamax2",
                 "BBmsyamax",
                 "FFmsyamax",
                 "AAVC2",
                 "CMSY2","CMSYST2","CMSYLT2",
                 "CMSYMaxAge","PBBlimMaxAge",
                 "BBmsy","BBmsyLT",
                 "BBmsyFL","FFmsyFL",
                 "AAVCfirst10",
                 "AAVBfirst10",
                 ## OLDER:
                 "avCatchFirst5y","avCatchLast5y","BBmsyLowest",
                 "PBBlim2",
                 "PBBlimFirst5y","PBBlimFirst10y","PBBlimLast5y","CatchCV", "avRelCatch",
                 "avRelCatchLast5y","converged",
                 "avRelCatchFirst5y", "BBmsy5y","TimeRecov",
                 "AAVCFirst5y",
                 "CMSYnew","PBBlimnew","AAVCnew","AAVBnew",
                 "CMSYnew_st","PBBlimnew_st","AAVCnew_st","AAVBnew_st",
                 "CMSYnew_lt","PBBlimnew_lt","AAVCnew_lt","AAVBnew_lt",
                 "Cnew","Cnew_st","Cnew_lt",
                 "BBmsynew","FFmsynew",
                 "BBmsynew_st","FFmsynew_st",
                 "BBmsynew_lt","FFmsynew_lt",
                 "CMSYamax5","PBBlimamax5","AAVCamax5",
                 "CMSYnew_y1","PBBlimnew_y1",
                 "CMSYnew_y2","PBBlimnew_y2",
                 "timeBmsy",
                 "t2risk5",
                 "t2bmsy"
                 )
    if(mets[1] == "all") mets <- metsAll
    if(any(which(!mets %in% metsAll))) writeLines(paste0("Metric ",
                                                         paste0(mets[which(!mets %in% metsAll)], collapse=", "),
                                                         " not defined and thus omitted."))

    mets <- metsAll[which(metsAll %in% mets)]
    nmets <- length(mets)

    ## TODO: make part of mse
    refs <- dat$ref
    if(is.null(refs)) stop("Reference levels missing.")
    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    nyall <- ny + nysim
    ntall <- nyall * ns
    nrefs <- nrow(refs)
    if(nrefs < ntall){
        refs <- rbind(refs, refs[rep((nrefs-ns+1):nrefs, ntall-nt),])
    }
    ##
    res2 <- vector("list", nhcr)
    for(hcr in 1:nhcr){
        metsUsed <- NULL
        nrep <- length(mse[[hcr]])
        print(names(mse)[hcr])
        msei <- mse[[hcr]]
        res <- NULL
        if(length(msei) == 0) break()
        ## CMSY
        if(any(mets == "CMSY")){
            metsUsed <- c(metsUsed, "CMSY")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[simYears] / refyield[["simYears"]][[indi[x]]]))
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
        ## "PBBlim" # used for probHCR, but not averaged over years first!
        if(any(mets == "PBBlim")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlim")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[simYears] / refs$Blim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp),
                             n = length(tmp),
                             conf.level = 0.95,
                             correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        if(any(mets == "PBBlimSSB")){
            if(is.null(refs$SSBlim)) stop("There is no Blim in dat$ref! Please add a SSBlim!")
            metsUsed <- c(metsUsed, "PBBlimSSB")
            tmp <- unlist(lapply(msei, function(x) mean(x$SSBfinal[simYears] / refs$SSBlim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "AAVC"
        if(any(mets == "AAVC")){
            metsUsed <- c(metsUsed, "AAVC")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[simYears[-1]] -
                                                     apply(x$CW,1,sum)[simYears[-length(simYears)]]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[simYears[-1]], na.rm=TRUE)))
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
        ## "PBBlim1"
        if(any(mets == "PBBlim1")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlim1")
            ry <- vector("list",nysim)
            for(y in 1:nysim){
                yi <- ((ny+1):finalYear)[y]
                tmp <- sapply(msei, function(x) x$TSBfinal[yi] / refs$Blim[yi] < 1)
                tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
                ry[[y]] <- c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2])
            }
            tmp <- do.call(rbind,ry)
            res <- rbind(res, c(NA,tmp[nysim,],NA, NA, NA))
        }
        ## "PBBlim3"
        if(any(mets == "PBBlim3")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlim3")
            ry <- vector("list",nysim)
            for(y in 1:nysim){
                yi <- ((ny+1):finalYear)[y]
                tmp <- sapply(msei, function(x) x$TSBfinal[yi] / refs$Blim[yi] < 1)
                tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
                ry[[y]] <- c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2])
            }
            tmp <- do.call(rbind,ry)
            res <- rbind(res, c(max(tmp[,1]),max(tmp[,2]),max(tmp[,3]), NA, NA))
        }
        ## CMSYST
        if(any(mets == "CMSYST")){
            metsUsed <- c(metsUsed, "CMSYST")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[first5Years] / refyield[["first5Years"]][[indi[x]]]))
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
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimST")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[first5Years] / refs$Blim[first5Years] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        if(any(mets == "PBBlim1ST")){
            metsUsed <- c(metsUsed, "PBBlim1ST")
            ry <- vector("list",nysim)
            for(y in 1:nysim){
                yi <- ((ny+1):finalYear)[y]
                tmp <- sapply(msei, function(x) x$TSBfinal[yi] / refs$Blim < 1)
                tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
                ry[[y]] <- c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2])
            }
            tmp <- do.call(rbind,ry)
            res <- rbind(res, c(mean(tmp[1:5,1]),mean(tmp[1:5,2]),mean(tmp[1:5,3]), NA, NA))
        }
        if(any(mets == "PBBlim3ST")){
            metsUsed <- c(metsUsed, "PBBlim3ST")
            ry <- vector("list",nysim)
            for(y in 1:nysim){
                yi <- ((ny+1):finalYear)[y]
                tmp <- sapply(msei, function(x) x$TSBfinal[yi] / refs$Blim < 1)
                tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
                ry[[y]] <- c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2])
            }
            tmp <- do.call(rbind,ry)
            res <- rbind(res, c(max(tmp[1:5,1]),max(tmp[1:5,2]),max(tmp[1:5,3]), NA, NA))
        }
        ## CMSYLT
        if(any(mets == "CMSYLT")){
            metsUsed <- c(metsUsed, "CMSYLT")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[last15Years] / refyield[["last15Years"]][[indi[x]]]))
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
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimLT")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last15Years] / refs$Blim[last15Years] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        if(any(mets == "PBBlim1LT")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlim1LT")
            ry <- vector("list",nysim)
            for(y in 1:nysim){
                yi <- ((ny+1):finalYear)[y]
                tmp <- sapply(msei, function(x) x$TSBfinal[yi] / refs$Blim < 1)
                tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
                ry[[y]] <- c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2])
            }
            tmp <- do.call(rbind,ry)
            res <- rbind(res, c(mean(tmp[6:40,1]),mean(tmp[6:40,2]),mean(tmp[6:40,3]), NA, NA))
        }
        if(any(mets == "PBBlim3LT")){
            metsUsed <- c(metsUsed, "PBBlim3LT")
            ry <- vector("list",nysim)
            for(y in 1:nysim){
                yi <- ((ny+1):finalYear)[y]
                tmp <- sapply(msei, function(x) x$TSBfinal[yi] / refs$Blim < 1)
                tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
                ry[[y]] <- c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2])
            }
            tmp <- do.call(rbind,ry)
            res <- rbind(res, c(max(tmp[6:40,1]),max(tmp[6:40,2]),max(tmp[6:40,3]), NA, NA))
        }
        ## CMSYamax
        if(any(mets == "CMSYamax")){
            metsUsed <- c(metsUsed, "CMSYamax")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[amaxYears] / refyield[["amaxYears"]][[indi[x]]]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlimamax"
        if(any(mets == "PBBlimamax")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimamax")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[amaxYears] / refs$Blim[amaxYears] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "AAVCamax"
        if(any(mets == "AAVCamax")){
            metsUsed <- c(metsUsed, "AAVCamax")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[amaxYears[-1]] -
                                                     apply(x$CW,1,sum)[amaxYears[-length(amaxYears)]]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[amaxYears[-1]], na.rm=TRUE)))
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
        if(any(mets == "CMSYamax2")){
            metsUsed <- c(metsUsed, "CMSYamax2")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[amaxYears3] / refyield[["amaxYears"]][[indi[x]]]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlimamax2"
        if(any(mets == "PBBlimamax2")){
            metsUsed <- c(metsUsed, "PBBlimamax2")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[amaxYears4] / refs$Blim[amaxYears4] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "BBmsy"
        if(any(mets == "BBmsyamax")){
            metsUsed <- c(metsUsed, "BBmsyamax")
            tmp <- unlist(lapply(msei, function(x) median(x$TSBfinal[amaxYears2] / refs$Bmsy[amaxYears2])))
            res <- rbind(res, quantile(tmp, probs = c(0.2, 0.5, 0.8), na.rm=TRUE))
        }
        ## "FFmsy"
        if(any(mets == "FFmsyamax")){
            metsUsed <- c(metsUsed, "FFmsyamax")
            tmp <- unlist(lapply(msei, function(x) median(apply(x$FM,1,sum)[amaxYears2] / refs$Fmsy[amaxYears2])))
            res <- rbind(res, quantile(tmp, probs = c(0.2, 0.5, 0.8), na.rm=TRUE))
        }
        ## CMSYlast5
        if(any(mets == "CMSYlast5")){
            metsUsed <- c(metsUsed, "CMSYlast5")
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
            metsUsed <- c(metsUsed, "CMSYlast10")
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
        }
        ## CMSYfirst10
        if(any(mets == "CMSYfirst10")){
            metsUsed <- c(metsUsed, "CMSYfirst10")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x)
                    median(apply(msei[[x]]$CW,1,sum)[first10Years] /
                           refyield[["first10Years"]][[indi[x]]]))
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
        ## "PBBlimlast5"
        if(any(mets == "PBBlimlast5")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimlast5")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last5Years] / refs$Blim[last5Years] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "PBBlimlast10"
        if(any(mets == "PBBlimlast10")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimlast10")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[last10Years] / refs$Blim[last10Years] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "PBBlimlast10"
        if(any(mets == "PBBlimlast10SSB")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimlast10SSB")
            tmp <- unlist(lapply(msei, function(x) mean(x$SSBfinal[last10Years] / refs$SSBlim[last10Years] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "PBBlimfirst10"
        if(any(mets == "PBBlimfirst10")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimfirst10")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[first10Years] / refs$Blim[first10Years] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "PBBlimfirst10"
        if(any(mets == "PBBlimfirst10SSB")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimfirst10SSB")
            tmp <- unlist(lapply(msei, function(x) mean(x$SSBfinal[first10Years] / refs$SSBlim[first10Years] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "AAVC2"
        if(any(mets == "AAVC2")){
            metsUsed <- c(metsUsed, "AAVC2")
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
            metsUsed <- c(metsUsed, "CMSY2")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x) mean(apply(msei[[x]]$CW,1,sum)[simYears] /
                                                               refyield[["simYears"]][[indi[x]]]))
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
            metsUsed <- c(metsUsed, "CMSYST2")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x) mean(apply(msei[[x]]$CW,1,sum)[first5Years] /
                                                               refyield[["first5Years"]][[indi[x]]]))
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
            metsUsed <- c(metsUsed, "CMSYLT2")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- sapply(1:length(msei), function(x) mean(apply(msei[[x]]$CW,1,sum)[last15Years] /
                                                               refyield[["last15Years"]][[indi[x]]]))
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
        ## "AAVCfirst10"
        if(any(mets == "AAVCfirst10")){
            metsUsed <- c(metsUsed, "AAVCfirst10")
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[first10Years[-1]] -
                                                     apply(x$CW,1,sum)[first10Years[-length(first10Years)]]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[first10Years[-1]], na.rm=TRUE)))
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
        ## "AAVBfirst10"
        if(any(mets == "AAVBfirst10")){
            metsUsed <- c(metsUsed, "AAVBfirst10")
            tmp <- sapply(msei, function(x) (sum(abs(x$TSBfinal[first10Years[-1]] -
                                                     x$TSBfinal[first10Years[-length(first10Years)]]),na.rm=TRUE)/
                                             sum(x$TSBfinal[first10Years[-1]], na.rm=TRUE)))
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
        ## "AAVB"
        if(any(mets == "AAVB")){
            metsUsed <- c(metsUsed, "AAVB")
            tmp <- sapply(msei, function(x) (sum(abs(x$TSBfinal[simYears[-1]] -
                                                     x$TSBfinal[simYears[-length(simYears)]]),na.rm=TRUE)/
                                             sum(x$TSBfinal[simYears[-1]], na.rm=TRUE)))
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
        ## OLDER (not used in probHCR):
        ## CMSYMaxAge
        if(any(mets == "CMSYMaxAge")){
            metsUsed <- c(metsUsed, "CMSYMaxAge")
            tmp <- unlist(lapply(msei, function(x) median(apply(x$CW,1,sum)[min(simYears):(min(simYears)+dat$amax)] /
                                                   refs$MSY)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "PBBlimMaxAge"
        if(any(mets == "PBBlimMaxAge")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimMaxAge")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[min(simYears):(min(simYears)+dat$amax)]/
                                                        refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ## "BBmsy"
        if(any(mets == "BBmsy")){
            metsUsed <- c(metsUsed, "BBmsy")
            tmp <- unlist(lapply(msei, function(x) median(x$TSBfinal[simYears] / refs$Bmsy)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## "BBmsyLT"
        if(any(mets == "BBmsyLT")){
            metsUsed <- c(metsUsed, "BBmsyLT")
            tmp <- unlist(lapply(msei, function(x) median(x$TSBfinal[last5Years] / refs$Bmsy)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        if(any(mets == "CMSYold")){
            metsUsed <- c(metsUsed, "CMSYold")
            ## tmp <- unlist(lapply(msei, function(x) mean(apply(x$CW,1,sum)[last5Years] / refs$MSY)))## HERE:
            ## res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
            tmp <- unlist(lapply(msei, function(x) median(apply(x$CW,1,sum)[simYears] / refs$MSY)))
            res <- rbind(res, quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE))
        }
        ## OLDER:
        ## "BBmsyFL"
        if(any(mets == "BBmsyFL")){
            metsUsed <- c(metsUsed, "BBmsyFL")
            indi <- as.numeric(names(msei))
            tmp <- sapply(1:length(msei), function(x) msei[[x]]$TSBfinal[finalYear]/
                                                      refB[["finalYear"]][[indi[x]]]) ##refs$Bmsy)
            res <- rbind(res, c(quantile(tmp, probs = c(0.25, 0.5, 0.75), na.rm=TRUE),NA,NA))
        }
        if(any(mets == "FFmsyFL")){
            metsUsed <- c(metsUsed, "FFmsyFL")
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
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSB[simYears,1]/refs$Blim < 1)))
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1],tmp$estimate,tmp$conf.int[2]))
        }
        ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlimFirst5y"
        if(any(mets == "PBBlimFirst5y")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimFirst5y")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[first5Years]/refs$Blim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95,
                             correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate,
                                tmp$conf.int[2], sei, ni))
        }
        if(any(mets == "PBBlimFirst10y")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimFirst10y")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[first10Years]/refs$Blim < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95,
                             correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate,
                                tmp$conf.int[2], sei, ni))
        }
        ##        res[6,] <- c(sd(tmp),mean(tmp),NA) ##quantile(tmp, probs = c(0.025, 0.5, 0.975), na.rm=TRUE)
        ## "PBBlimLast5y"
        if(any(mets == "PBBlimLast5y")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
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
            metsUsed <- c(metsUsed, "avRelCatch")
            if("refFmsy" %in% hcrs){
                tmpRef <- lapply(resMSE[[which(hcrs == "refFmsy")]], function(x) apply(x$CW,1,mean)[simYears])
                tmp0 <- lapply(msei, function(x) apply(x$CW,1,sum)[simYears])
                tmp <- unlist(lapply(1:nrep, function(x) tmp0[[x]] / tmpRef[[x]]))
                mu <- mean(tmp,na.rm=TRUE)
                vari <- var(tmp,na.rm=TRUE)
                sei <- sqrt(vari/length(tmp))
                res <- rbind(res, c(quantile(tmp, probs = c(0.025, 0.5, 0.975),
                                             na.rm=TRUE),
                                    sei,
                                    length(tmp)))
            }else{
                res <- rbind(res, rep(NA,5))
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
        ## NEW:
        ## CMSYamax
        if(any(mets == "CMSYnew")){
            metsUsed <- c(metsUsed, "CMSYnew")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[newYears] / refyield[["newYears"]][[indi[x]]]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        if(any(mets == "Cnew")){
            metsUsed <- c(metsUsed, "Cnew")
            indi <- as.numeric(names(msei))
            tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                apply(msei[[x]]$CW,1,sum)[newYears]))
            meani <- mean(tmp)
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
            }else{
                res <- rbind(res, c(tmp2$conf.int[1],
                                    tmp2$estimate,
                                    tmp2$conf.int[2],
                                    sei,
                                    ni))
            }
        }
        ## "PBBlimamax"
        if(any(mets == "PBBlimnew")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimnew")
            ## Average risk over years per rep
            tmp <- unlist(lapply(msei, function(x)
                mean((x$TSBfinal[newYears] / refs$Blim[newYears]) < 1)))

            ## Average risk over reps per year
            ## tmp <- sapply(newYears, function(x) mean(sapply(msei, function(y) y$TSBfinal[x] / refs$Blim[x] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95,
                             correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate,
                                tmp$conf.int[2], sei, ni))

            ## try B relative to B of Fmsy ref rule
            ##     indi <- as.numeric(names(msei))
            ##     tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
            ##         msei[[x]]$ESBfinal[newYears] / refB[["newYears"]][[indi[x]]]))
            ##     meani <- mean(tmp)
            ##     vari <- var(tmp)
            ##     ni <- length(tmp)
            ##     sei <- sqrt(vari/ni)
            ##     erri <- 1.96 * (sqrt(vari)/sqrt(ni))
            ##     tmp2 <- try(wilcox.test(as.numeric(tmp),
            ##                        alternative="two.sided",
            ##                        correct=TRUE,
            ##                        conf.int=TRUE,
            ##                        conf.level=0.95), silent=TRUE)
            ## res <- rbind(res, c(tmp2$conf.int[1],
            ##                     tmp2$estimate,
            ##                     tmp2$conf.int[2],
            ##                     sei,
            ##                     ni))
        }
        ## "AAVCamax"
        if(any(mets == "AAVCnew")){
            metsUsed <- c(metsUsed, "AAVCnew")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[newYears[-1]] -
                                                     apply(x$CW,1,sum)[newYears[-length(newYears)]]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[newYears[-1]], na.rm=TRUE)))
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

        if(any(mets == "AAVBnew")){
            metsUsed <- c(metsUsed, "AAVBnew")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(x$TSBfinal[newYears[-1]] -
                                                     x$TSBfinal[newYears[-length(newYears)]]),na.rm=TRUE)/
                                             sum(x$TSBfinal[newYears[-1]], na.rm=TRUE)))
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

        ## Short-term

        ## CMSYamax
        if(any(mets == "CMSYnew_st")){
            metsUsed <- c(metsUsed, "CMSYnew_st")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[newYears_st] / refyield[["newYears_st"]][[indi[x]]]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        if(any(mets == "Cnew_st")){
            metsUsed <- c(metsUsed, "Cnew_st")
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[newYears_st]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
        }
        ## "PBBlimamax"
        if(any(mets == "PBBlimnew_st")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimnew_st")
            ## Average risk over years per rep
            tmp <- unlist(lapply(msei, function(x)
                mean((x$TSBfinal[newYears_st] / refs$Blim[newYears_st]) < 1)))

            ## Average risk over reps per year
            ## tmp <- sapply(newYears, function(x) mean(sapply(msei, function(y) y$TSBfinal[x] / refs$Blim[x] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95,
                             correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate,
                                tmp$conf.int[2], sei, ni))

            ## try B relative to B of Fmsy ref rule
            ##     indi <- as.numeric(names(msei))
            ##     tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
            ##         msei[[x]]$ESBfinal[newYears] / refB[["newYears"]][[indi[x]]]))
            ##     meani <- mean(tmp)
            ##     vari <- var(tmp)
            ##     ni <- length(tmp)
            ##     sei <- sqrt(vari/ni)
            ##     erri <- 1.96 * (sqrt(vari)/sqrt(ni))
            ##     tmp2 <- try(wilcox.test(as.numeric(tmp),
            ##                        alternative="two.sided",
            ##                        correct=TRUE,
            ##                        conf.int=TRUE,
            ##                        conf.level=0.95), silent=TRUE)
            ## res <- rbind(res, c(tmp2$conf.int[1],
            ##                     tmp2$estimate,
            ##                     tmp2$conf.int[2],
            ##                     sei,
            ##                     ni))
        }
        ## "AAVCamax"
        if(any(mets == "AAVCnew_st")){
            metsUsed <- c(metsUsed, "AAVCnew_st")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[newYears_st[-1]] -
                                                     apply(x$CW,1,sum)[newYears_st[-length(newYears_st)]]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[newYears_st[-1]], na.rm=TRUE)))
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
if(any(mets == "AAVBnew_st")){
            metsUsed <- c(metsUsed, "AAVBnew_st")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(x$TSBfinal[newYears_st[-1]] -
                                                     x$TSBfinal[newYears_st[-length(newYears_st)]]),na.rm=TRUE)/
                                             sum(x$TSBfinal[newYears_st[-1]], na.rm=TRUE)))
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


        ## long-term
        ## CMSYamax
        if(any(mets == "CMSYnew_lt")){
            metsUsed <- c(metsUsed, "CMSYnew_lt")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[newYears_lt] / refyield[["newYears_lt"]][[indi[x]]]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        if(any(mets == "Cnew_lt")){
            metsUsed <- c(metsUsed, "Cnew_lt")
            indi <- as.numeric(names(msei))
            tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                apply(msei[[x]]$CW,1,sum)[newYears_lt]))
            meani <- mean(tmp)
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
            }else{
                res <- rbind(res, c(tmp2$conf.int[1],
                                    tmp2$estimate,
                                    tmp2$conf.int[2],
                                    sei,
                                    ni))
            }
        }
        ## "PBBlimamax"
        if(any(mets == "PBBlimnew_lt")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimnew_lt")
            ## Average risk over years per rep
            tmp <- unlist(lapply(msei, function(x)
                mean((x$TSBfinal[newYears_lt] / refs$Blim[newYears_lt]) < 1)))

            ## Average risk over reps per year
            ## tmp <- sapply(newYears, function(x) mean(sapply(msei, function(y) y$TSBfinal[x] / refs$Blim[x] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95,
                             correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate,
                                tmp$conf.int[2], sei, ni))

            ## try B relative to B of Fmsy ref rule
            ##     indi <- as.numeric(names(msei))
            ##     tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
            ##         msei[[x]]$ESBfinal[newYears] / refB[["newYears"]][[indi[x]]]))
            ##     meani <- mean(tmp)
            ##     vari <- var(tmp)
            ##     ni <- length(tmp)
            ##     sei <- sqrt(vari/ni)
            ##     erri <- 1.96 * (sqrt(vari)/sqrt(ni))
            ##     tmp2 <- try(wilcox.test(as.numeric(tmp),
            ##                        alternative="two.sided",
            ##                        correct=TRUE,
            ##                        conf.int=TRUE,
            ##                        conf.level=0.95), silent=TRUE)
            ## res <- rbind(res, c(tmp2$conf.int[1],
            ##                     tmp2$estimate,
            ##                     tmp2$conf.int[2],
            ##                     sei,
            ##                     ni))
        }
        ## "AAVCamax"
        if(any(mets == "AAVCnew_lt")){
            metsUsed <- c(metsUsed, "AAVCnew_lt")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[newYears_lt[-1]] -
                                                     apply(x$CW,1,sum)[newYears_lt[-length(newYears_lt)]]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[newYears_lt[-1]], na.rm=TRUE)))
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

if(any(mets == "AAVBnew_lt")){
            metsUsed <- c(metsUsed, "AAVBnew_lt")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(x$TSBfinal[newYears_lt[-1]] -
                                                     x$TSBfinal[newYears_lt[-length(newYears_lt)]]),na.rm=TRUE)/
                                             sum(x$TSBfinal[newYears_lt[-1]], na.rm=TRUE)))
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


        ## first year

        ## CMSY
        if(any(mets == "CMSYnew_y1")){
            metsUsed <- c(metsUsed, "CMSYnew_y1")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[newYears_y1]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlimamax"
        if(any(mets == "PBBlimnew_y1")){
            metsUsed <- c(metsUsed, "PBBlimnew_y1")
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")

            ## Average risk over years per rep
            ## tmp <- unlist(lapply(msei, function(x)
            ##     mean((x$TSBfinal[newYears_y1] / refs$Blim[newYears_y1]) < 1)))

            ## Average risk over reps per year
            tmp <- sapply(newYears_y1, function(x) mean(sapply(msei, function(y) y$TSBfinal[x] / refs$Blim[x] < 1)))
            ## vari <- var(tmp)
            ## ni <- length(tmp)
            ## sei <- sqrt(vari/ni)
            ## tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95,
            ##                  correct = FALSE)
            ## res <- rbind(res, c(tmp$conf.int[1], tmp$estimate,
            ##                     tmp$conf.int[2], sei, ni))

            res <- rbind(res, c(NA, tmp,
                                NA, NA, NA))

            ## try B relative to B of Fmsy ref rule
            ##     indi <- as.numeric(names(msei))
            ##     tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
            ##         msei[[x]]$ESBfinal[newYears] / refB[["newYears"]][[indi[x]]]))
            ##     meani <- mean(tmp)
            ##     vari <- var(tmp)
            ##     ni <- length(tmp)
            ##     sei <- sqrt(vari/ni)
            ##     erri <- 1.96 * (sqrt(vari)/sqrt(ni))
            ##     tmp2 <- try(wilcox.test(as.numeric(tmp),
            ##                        alternative="two.sided",
            ##                        correct=TRUE,
            ##                        conf.int=TRUE,
            ##                        conf.level=0.95), silent=TRUE)
            ## res <- rbind(res, c(tmp2$conf.int[1],
            ##                     tmp2$estimate,
            ##                     tmp2$conf.int[2],
            ##                     sei,
            ##                     ni))
        }

        ## CMSY
        if(any(mets == "CMSYnew_y2")){
            metsUsed <- c(metsUsed, "CMSYnew_y2")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[newYears_y2] / refyield[["newYears_y2"]][[indi[x]]]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlimamax"
        if(any(mets == "PBBlimnew_y2")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimnew_y2")
            ## Average risk over years per rep
            tmp <- unlist(lapply(msei, function(x)
                mean((x$TSBfinal[newYears_y2] / refs$Blim[newYears_y2]) < 1)))

            ## Average risk over reps per year
            ## tmp <- sapply(newYears, function(x) mean(sapply(msei, function(y) y$TSBfinal[x] / refs$Blim[x] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95,
                             correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate,
                                tmp$conf.int[2], sei, ni))

            ## try B relative to B of Fmsy ref rule
            ##     indi <- as.numeric(names(msei))
            ##     tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
            ##         msei[[x]]$ESBfinal[newYears] / refB[["newYears"]][[indi[x]]]))
            ##     meani <- mean(tmp)
            ##     vari <- var(tmp)
            ##     ni <- length(tmp)
            ##     sei <- sqrt(vari/ni)
            ##     erri <- 1.96 * (sqrt(vari)/sqrt(ni))
            ##     tmp2 <- try(wilcox.test(as.numeric(tmp),
            ##                        alternative="two.sided",
            ##                        correct=TRUE,
            ##                        conf.int=TRUE,
            ##                        conf.level=0.95), silent=TRUE)
            ## res <- rbind(res, c(tmp2$conf.int[1],
            ##                     tmp2$estimate,
            ##                     tmp2$conf.int[2],
            ##                     sei,
            ##                     ni))
        }


        ## for kobe

        ## "BBmsy"
        if(any(mets == "BBmsynew")){
            metsUsed <- c(metsUsed, "BBmsynew")
            tmp <- unlist(lapply(msei, function(x) median(x$TSBfinal[newYears2_mt] / refs$Bmsy[newYears2_mt])))
            res <- rbind(res, quantile(tmp, probs = c(0.2, 0.5, 0.8), na.rm=TRUE))
        }
        ## "FFmsy"
        if(any(mets == "FFmsynew")){
            metsUsed <- c(metsUsed, "FFmsynew")
            tmp <- unlist(lapply(msei, function(x) median(apply(x$FM,1,sum)[newYears2_mt] / refs$Fmsy[newYears2_mt])))
            res <- rbind(res, quantile(tmp, probs = c(0.2, 0.5, 0.8), na.rm=TRUE))
        }

        ## short-term
        ## "BBmsy"
        if(any(mets == "BBmsynew_st")){
            metsUsed <- c(metsUsed, "BBmsynew_st")
            tmp <- unlist(lapply(msei, function(x) median(x$TSBfinal[newYears_st] / refs$Bmsy[newYears_st])))
            res <- rbind(res, quantile(tmp, probs = c(0.2, 0.5, 0.8), na.rm=TRUE))
        }
        ## "FFmsy"
        if(any(mets == "FFmsynew_st")){
            metsUsed <- c(metsUsed, "FFmsynew_st")
            tmp <- unlist(lapply(msei, function(x) median(apply(x$FM,1,sum)[newYears_st] / refs$Fmsy[newYears_st])))
            res <- rbind(res, quantile(tmp, probs = c(0.2, 0.5, 0.8), na.rm=TRUE))
        }

        ## long-term
        ## "BBmsy"
        if(any(mets == "BBmsynew_lt")){
            metsUsed <- c(metsUsed, "BBmsynew_lt")
            tmp <- unlist(lapply(msei, function(x) median(x$TSBfinal[newYears_lt] / refs$Bmsy[newYears_lt])))
            res <- rbind(res, quantile(tmp, probs = c(0.2, 0.5, 0.8), na.rm=TRUE))
        }
        ## "FFmsy"
        if(any(mets == "FFmsynew_lt")){
            metsUsed <- c(metsUsed, "FFmsynew_lt")
            tmp <- unlist(lapply(msei, function(x) median(apply(x$FM,1,sum)[newYears_lt] / refs$Fmsy[newYears_lt])))
            res <- rbind(res, quantile(tmp, probs = c(0.2, 0.5, 0.8), na.rm=TRUE))
        }

        ## CMSYamax5
        if(any(mets == "CMSYamax5")){
            metsUsed <- c(metsUsed, "CMSYamax5")
            if(length(reffmsyInd) > 0){
                indi <- as.numeric(names(msei))
                tmp <- unlist(lapply(as.list(1:length(msei)), function(x)
                    apply(msei[[x]]$CW,1,sum)[amaxYears5] / refyield[["amaxYears5"]][[indi[x]]]))
                meani <- mean(tmp)
                vari <- var(tmp)
                ni <- length(tmp)
                sei <- sqrt(vari/ni)
                erri <- 1.96 * (sqrt(vari)/sqrt(ni))
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
                ## }else if(hcrs[hcr] == "refFmsy"){
                ##     res <- rbind(res, c(1,
                ##                         1,
                ##                         1,
                ##                         sei,
                ##                         ni))
                }else{
                    res <- rbind(res, c(tmp2$conf.int[1],
                                        tmp2$estimate,
                                        tmp2$conf.int[2],
                                        sei,
                                        ni))
                }
            }else writeLines("CMSYLT could not be estimated, because no rule 'refFmsy' not found.")
        }
        ## "PBBlimamax5"
        if(any(mets == "PBBlimamax5")){
            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "PBBlimamax5")
            tmp <- unlist(lapply(msei, function(x) mean(x$TSBfinal[amaxYears5] / refs$Blim[amaxYears5] < 1)))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp <- prop.test(sum(tmp), n = length(tmp), conf.level = 0.95, correct = FALSE)
            res <- rbind(res, c(tmp$conf.int[1], tmp$estimate, tmp$conf.int[2], sei, ni))
        }
        ## "AAVCamax5"
        if(any(mets == "AAVCamax5")){
            metsUsed <- c(metsUsed, "AAVCamax5")
            ## tmp <- unlist(lapply(lapply(msei,
            ##                             function(x) (((apply(x$CW,1,sum)[simYears] -
            ##                                            apply(x$CW,1,sum)[simYears+1])/
            ##                                           apply(x$CW,1,sum)[simYears+1])^2)^0.5),
            ##                      median ,na.rm=TRUE))
            tmp <- sapply(msei, function(x) (sum(abs(apply(x$CW,1,sum)[amaxYears5[-1]] -
                                                     apply(x$CW,1,sum)[amaxYears5[-length(amaxYears5)]]),na.rm=TRUE)/
                                             sum(apply(x$CW,1,sum)[amaxYears5[-1]], na.rm=TRUE)))
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

        if(any(mets == "timeBmsy")){
            if(is.null(refs$Bmsy)) stop("There is no Bmsy in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "timeBmsy")
            tmp <- unlist(lapply(msei, function(x){
                indi <- which(x$TSBfinal[simYears] > refs$Bmsy[1])
                if(length(indi) > 0) min(indi) else tail(nysim,1)
            }))
            vari <- var(tmp)
            ni <- length(tmp)
            sei <- sqrt(vari/ni)
            tmp2 <- try(wilcox.test(as.numeric(tmp),
                                    alternative="two.sided",
                                    correct=TRUE,
                                    conf.int=TRUE,
                                    conf.level=0.95), silent=TRUE)
            res <- rbind(res, c(tmp2$conf.int[1], tmp2$estimate, tmp2$conf.int[2],
                                sei, ni))
        }

        if(any(mets == "t2risk5")){

            if(is.null(refs$Blim)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "t2risk5")
            ## risk by year
            ## x <- msei[[1]]
            tmp <- lapply(msei, function(x) x$TSBfinal[newYears] / refs$Blim[newYears])
            tmp2 <- do.call(rbind, tmp)
            tmp3 <- apply(tmp2, 2, function(x) mean(x < 1))
            est <- min(which(tmp3 < 0.05))

            ## uncertainty for rebuilding time via bootstrapping
            rebtimes <- replicate(500, {
                indi <- sample(1:nrow(tmp2), nrow(tmp2), replace = TRUE)
                tmp3u <- apply(tmp2[indi,], 2, function(x) mean(x < 1))
                min(which(tmp3u < 0.05))
            })

            unc <- quantile(rebtimes, c(0.025, 0.5, 0.975))

            unc[is.infinite(unc)] <- NA

            res <- rbind(res, c(unc[1], unc[2],
                                unc[3], NA, NA))

        }

        if(any(mets == "t2bmsy")){

            if(is.null(refs$Bmsy)) stop("There is no Blim in dat$ref! Please add a Blim!")
            metsUsed <- c(metsUsed, "t2bmsy")
            ## risk by year
            ## x <- msei[[1]]
            tmp <- lapply(msei, function(x) x$TSBfinal[newYears] / refs$Bmsy[newYears])
            tmp2 <- do.call(rbind, tmp)
            tmp3 <- apply(tmp2, 2, function(x) mean(x < 1))
            est <- min(which(tmp3 < 0.5))

            ## uncertainty for rebuilding time via bootstrapping
            rebtimes <- replicate(500, {
                indi <- sample(1:nrow(tmp2), nrow(tmp2), replace = TRUE)
                tmp3u <- apply(tmp2[indi,], 2, function(x) mean(x < 1))
                min(which(tmp3u < 0.5))
            })

            unc <- quantile(rebtimes, c(0.025, 0.5, 0.975))

            unc[is.infinite(unc)] <- NA

            res <- rbind(res, c(unc[1], unc[2],
                                unc[3], NA, NA))

        }

        ## names
        rownames(res) <- metsUsed
        colnames(res) <- c("loCI","mu","upCI","se","n")
        res2[[hcr]] <- round(res, 3)
    }
    names(res2) <- hcrs


    class(res2) <- "met"


    return(res2)
}
