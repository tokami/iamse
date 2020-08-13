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

    if(!any(names(dat) == "ref")) stop("Reference levels have to be part of dat. Use estRef to estimate them.")
    refs <- dat$ref

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
            if(hcrs2[i] %in% c("refFmsy","consF")){
                ## Reference rule Fmsy
                for(y in 1:nysim){
                    poptmp$tacs <- gettacs(tacs = poptmp$tacs, id = hcri, TAC=NA)
                    poptmp <- advancePop(dat = dat, hist = poptmp, set = setx,
                                         tacs = poptmp$tacs)
                }
            }else{
                ## Any other HCR
                for(y in 1:nysim){
                    poptmp$tacs <- estTAC(inp = poptmp$inp, hcr = hcri, tacs = poptmp$tacs)
                    poptmp <- advancePop(dat = dat, hist = poptmp, set = setx,
                                         tacs = poptmp$tacs)
                }
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
