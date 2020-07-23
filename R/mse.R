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
        setx$eC <- exp(rnorm(nysim, 0, set$CVC) - set$CVC^2/2)
        setx$eI <- exp(rnorm(nysim, 0, set$CVI) - set$CVI^2/2)
        setx$eF <- exp(rnorm(nysim, 0, set$sigmaF) - set$sigmaF^2/2)
        setx$eR <- genDevs(nysim, set$sigmaR, set$rhoR)
        setx$eM <- exp(rnorm(nysim, 0, set$sigmaM) - set$sigmaM^2/2)
        setx$eH <- exp(rnorm(nysim, 0, set$sigmaH) - set$sigmaH^2/2)
        setx$eR0 <- exp(rnorm(nysim, 0, set$sigmaR0) - set$sigmaR0^2/2)
        setx$eMat <- exp(rnorm(nysim, 0, set$sigmaMat) - set$sigmaMat^2/2)
        setx$eImp <- exp(rnorm(nysim, 0, set$sigmaImp) - set$sigmaImp^2/2)
        setx$eC <- exp(rnorm(nysim, 0, set$CVC) - set$CVC^2/2)
        setx$eI <- exp(rnorm(nysim, 0, set$CVI) - set$CVI^2/2)


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
        print(res[[ind]])
        stop("Some error in res loop.")
    }


    ## sort res of reps for each ms together
    resList <- vector("list", nhcrs)

    res2 <- lapply(1:nhcrs, function(x){
        tmp <- vector("list",nrep)
        for(i in 1:nrep){
            tmp[[i]] <- res[[i]][[x]]
        }
        resList[[x]] <- tmp
    })
    names(res2) <- hcrs

    return(res2)
}
