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
        setx$eF <- rlnorm(nysim, muconv(1,set$sigmaF), sdconv(1,set$sigmaF))
        setx$eR <- genDevs(nysim, set$sigmaR, set$rhoR)
        setx$eM <- rlnorm(nysim, muconv(1,set$sigmaM), sdconv(1,set$sigmaM))
        setx$eH <- rlnorm(nysim, muconv(1,set$sigmaH), sdconv(1,set$sigmaH))
        setx$eR0 <- rlnorm(nysim, muconv(1,set$sigmaR0), sdconv(1,set$sigmaR0))
        setx$eMat <- rlnorm(nysim, muconv(1,set$sigmaMat), sdconv(1,set$sigmaMat))
        setx$eImp <- rlnorm(nysim, muconv(1,set$sigmaImp), sdconv(1,set$sigmaImp))
        setx$eC <- rlnorm(nysim, muconv(1,set$CVC), sdconv(1,set$CVC))
        setx$eI <- list()
        for(i in 1:nsurv){
           setx$eI[[i]] <- rlnorm(nysim, muconv(1,set$CVI), sdconv(1,set$CVI))
        }

        ## loop
        for(i in 1:nhcrs){
            hcri <- hcrs[i]
            poptmp <- popListx[[i]]
            poptmp$tacs <- NULL
            if(hcri == "refFmsy"){
                ## Reference rule Fmsy
                for(y in 1:nysim){
                    poptmp$tacs <- gettacs(tacs = poptmp$tacs, id = "refFmsy", TAC=NA)
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
