##' @name runMSE
##'
##' @importFrom parallel mclapply
##'
##' @export
runMSE <- function(dat, set, ref, ncores=detectCores()-1, verbose=TRUE){

    if(ncores > 1) verbose <- FALSE

    ## Variables
    hcrs <- set$hcr
    nhcrs <- length(hcrs)
    nysim <- set$nysim
    nrep <- set$nrep

    ## parallel loop
    res <- parallel::mclapply(as.list(1:nrep), function(x){

        if(verbose) writeLines(paste0("Running replicate: ", x))

        ## pop list with errors
        pop <- initPop(dat, set)
        ## add reference levels
        pop$refs <- ref$refs
        ## add observation noise
        pop <- obsmod(specdat = dat, hist = pop, set = set,
                      years = (dat$ny-(set$nyhist-1)):dat$ny)
        popList <- vector("list", nhcrs)
        for(i in 1:nhcrs){
            popList[[i]] <- pop
        }

        repList <- vector("list", nhcrs)
        popListx <- popList
        setx <- set
        ## errors
        setx$eC <- rnorm(nysim, 0, set$CVC) - set$CVC^2/2
        setx$eI <- rnorm(nysim, 0, set$CVI) - set$CVI^2/2
        setx$eF <- rnorm(nysim, 0, set$sigmaF) - set$sigmaF^2/2
        setx$eR <- genDevs(nysim, set$sigmaR, set$rhoR)
        setx$eM <- rnorm(nysim, 0, set$sigmaM) - set$sigmaM^2/2
        setx$eH <- rnorm(nysim, 0, set$sigmaH) - set$sigmaH^2/2
        setx$eMat <- rnorm(nysim, 0, set$sigmaMat) - set$sigmaMat^2/2
        ## loop
        for(i in 1:nhcrs){
            hcri <- hcrs[i]
            poptmp <- popListx[[i]]
            poptmp$tacs <- NULL
            if(hcri == "refFmsy"){
                ## Reference rule Fmsy
                for(y in 1:nysim){
                    poptmp$tacs <- estTAC(inp = NA, hcr = hcri, stab = setx$stab,
                                          tacs = poptmp$tacs, hist = poptmp)
                    poptmp <- advancePop(specdat = dat, hist = poptmp, set = setx,
                                         tacs = poptmp$tacs)
                    poptmp <- obsmod(specdat = dat, hist = poptmp, set = setx)
                }

            }else{
                ## Any other HCR
                for(y in 1:nysim){
                    inp <- poptmp$obs
                    poptmp$tacs <- estTAC(inp = inp, hcr = hcri, stab = setx$stab,
                                          tacs = poptmp$tacs, hist = poptmp)
                    poptmp <- advancePop(specdat = dat, hist = poptmp, set = setx,
                                         tacs = poptmp$tacs)
                    poptmp <- obsmod(specdat = dat, hist = poptmp, set = setx)
                }
            }
            popListx[[i]] <- poptmp
            gc()
        }
        repList[[x]] <- popListx
    }, mc.cores = ncores)


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
