##' @name runMSE
##'
##' @importFrom parallel mclapply
##' @importFrom parallel detectCores
##'
##' @export
runMSE <- function(dat, set, ncores=parallel::detectCores()-1, verbose=TRUE){
    if(ncores > 1) verbose <- FALSE

    ## define constant catch (resort HCR if something not converging)
    defHCRconscat()

    ## Variables
    hcrs <- set$hcr
    hcrs2 <- sapply(hcrs, function(x) unlist(strsplit(as.character(x), "-"))[1])
    nhcrs <- length(hcrs)
    nysim <- set$nysim
    nrep <- set$nrep
    ny <- dat$ny
    nyall <- ny + nysim
    ns <- dat$ns
    nt <- nyall * ns

    ## Checks
    ## Reference points provided
    if(!any(names(dat) == "ref")) stop("Reference levels have to be part of dat. Use estRef to estimate them.")
    refs <- dat$ref
    ## Natural mortality
    ## per year
    if(length(dat$M) == ny){
        dat$M <- c(dat$M, rep(tail(dat$M,1), nysim))
    }else if(length(dat$M) < nyall){
        if(verbose) warning("Length of vector with annual natural mortality (M) is shorter than sum of historical and projection period. Using last provided natural mortality value for missing years.")
        dat$M <- c(dat$M, rep(tail(dat$M,1), nyall - length(dat$M)))
    }
    ## per season
    if(length(dat$Ms) == ny * ns){
        dat$Ms <- c(dat$Ms, rep(tail(dat$Ms,1), nysim * ns))
    }else if(length(dat$M) < (ny * ns + nysim)){
        if(verbose) warning("Length of vector with seasonal natural mortality (Ms) is shorter than sum of historical and projection period. Using last provided natural mortality value for missing years.")
        dat$Ms <- c(dat$Ms, rep(tail(dat$Ms,ns), nt - length(dat$Ms)))
    }
    ## per length
    if(length(dat$Msels) == ny * ns){
        dat$Msels <- c(dat$Msels, rep(tail(dat$Msels,1), nysim * ns))
    }else if(length(dat$Msels) != 1 && length(dat$Msels) < nt){
        if(verbose) warning("Length of list with relative natural mortality at age (Msels) is shorter than sum of historical and projection period. Using last provided values for missing years.")
        dat$Msels <- c(dat$Msels, rep(tail(dat$Msels,1), nt - length(dat$Msels)))
    }
    dat$Msel <- lapply(dat$Msels, rowMeans)
    ## Indices
    dat$yvec <- rep(1:nyall, each = ns)
    dat$svec <- rep(1:ns, each = nyall)
    dat$s1vec <- seq(1, nt, ns)

    ## parallel loop
    res <- parallel::mclapply(as.list(1:nrep), function(x){

        if(verbose) writeLines(paste0("Running replicate: ", x))

        ## set seed
        if(is.numeric(set$seed)) set.seed(set$seed + x)

        ## pop list with errors
        pop <- initPop(dat, set)
        ## add reference levels
        pop$refs <- refs
        popList <- vector("list", nhcrs)
        for(i in 1:nhcrs){
            popList[[i]] <- pop
        }

        repList <- vector("list", nhcrs)
        popListx <- popList
        setx <- set
        ## errors
        setx$eF <- genNoise(nysim, set$noiseF[1], set$noiseF[2], set$noiseF[3])
        setx$eR <- genNoise(nysim, set$noiseR[1], set$noiseR[2], set$noiseR[3])
        setx$eM <- genNoise(nysim, set$noiseM[1], set$noiseM[2], set$noiseM[3])
        setx$eH <- genNoise(nysim, set$noiseH[1], set$noiseH[2], set$noiseH[3])
        setx$eR0 <- genNoise(nysim, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
        setx$eMat <- genNoise(nysim, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
        setx$eSel <- genNoise(nysim, set$noiseSel[1], set$noiseSel[2], set$noiseSel[3])
        setx$eImp <- genNoise(nysim, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
        setx$eC <- genNoise(nysim, set$noiseC[1], set$noiseC[2], set$noiseC[3])
        setx$eI <- list()
        for(i in 1:length(set$surveyTimes)){
           setx$eI[[i]] <- genNoise(nysim, set$noiseI[1], set$noiseI[2], set$noiseI[3])
        }

        ## loop
        for(i in 1:nhcrs){
            hcri <- hcrs[i]
            poptmp <- popListx[[i]]
            poptmp$tacs <- NULL
            for(y in 1:nysim){
                poptmp <- advancePop(dat = dat,
                                     hist = poptmp,
                                     set = setx,
                                     hcr = hcri,
                                     year = y)
            }
            popListx[[i]] <- poptmp
            gc()
        }
        repList[[x]] <- popListx
    }, mc.cores = ncores)


    ## Debugging printing
    if(any(sapply(res, length) != nhcrs)){
        ind <- which(sapply(res, length) != nhcrs)[1]
        writeLines(paste0("Info about failed replicate ",ind,": length = ",
                          length(res[[ind]]), " value = ", res[[ind]], " names = ",
                          names(res[[ind]]),
                          " length(res) = ", length(res)))
        writeLines("Warning messages: ")
        warnings()
        stop(paste0("Replicate ",ind," does not have the correct length."))
    }


    ## sort res of reps for each ms together
    resList <- vector("list", nhcrs)

    res2 <- lapply(1:nhcrs, function(x){
        tmp <- vector("list",nrep)
        for(i in 1:nrep){
            tmp[[i]] <- res[[i]][[x]]
        }
        names(tmp) <- 1:nrep
        resList[[x]] <- tmp
    })
    names(res2) <- hcrs

    class(res2) <- "mse"

    return(res2)
}
