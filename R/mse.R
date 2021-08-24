##' @name run.mse
##'
##' @importFrom parallel mclapply
##' @importFrom parallel detectCores
##'
##' @export
run.mse <- function(dat, set, ncores=parallel::detectCores()-1, verbose=TRUE, dbg = 0){

    if(ncores > 1) verbose <- FALSE
    ## define constant catch (resort HCR if something not converging)
    def.hcr.conscat()

    ## Variables
    hcrs <- set$hcr
    hcrs2 <- sapply(hcrs, function(x) unlist(strsplit(as.character(x), "-"))[1])
    nhcrs <- length(hcrs)
    nysim <- set$nysim
    nrep <- set$nrep
    ny <- dat$ny
    ns <- dat$ns
    nt <- ny * ns
    nyall <- ny + nysim
    ntall <- nyall * ns

    ## Checks
    ## --------------------

    ## Reference points provided
    if(!any(names(dat) == "ref")) stop("Reference levels have to be part of dat. Use est.ref.levels/est.ref.levels.stochastic to estimate them.")
    refs <- dat$ref

    nrefs <- nrow(refs)
    if(nrefs < ntall){
        refs <- rbind(refs, refs[rep((nrefs-ns+1):nrefs, ntall-nt),])
    }

    ## Natural mortality
    ## --------------------
    if(!inherits(dat$M, "matrix")) stop(paste0("dat$M is not a matrix. It should be a matrix with ",ns," columns and at least ",ny," rows."))
    if(ncol(dat$M) != ns) stop(paste0("dat$M has ",ncol(dat$M), "columns. It should have ",ns," columns."))
    ## per season
    if(nrow(dat$M) == ny){
        dat$M <- rbind(dat$M, matrix(rep(tail(dat$M,ns), nyall - ny), nrow = nyall - ny, ncol = ns))
        if(verbose) writeLines("No M provided for projection period. Using M in the last historical year for the projection period.")
    }else if(nrow(dat$M) < nyall){
        if(verbose) warning("Length of vector with natural mortality (dat$M) is shorter than the historical period. Using last provided natural mortality value for missing years and projection period.")
        dat$M <- rbind(dat$M, matrix(rep(tail(dat$M,ns), ny - nrow(dat$M)), nrow = ny-nrow(dat$M), ncol=ns))
    }
    ## per age
    if(length(dat$Msel) == nt){
        dat$Msel <- c(dat$Msel, rep(tail(dat$Msel,1), ntall - nt))
        if(verbose) writeLines("No Msel provided for projection period. Using Msel in the last historical year for the projection period.")
    }else if(length(dat$Msel) != 1 && length(dat$Msel) < nt){
        if(verbose) warning("Length of list with relative natural mortality at age (dat$Msel) is shorter than the historical period. Using last provided values for missing time step and projection period.")
        dat$Msel <- c(dat$Msel, rep(tail(dat$Msel,1), nt - length(dat$Msel)))
    }

    ## Selectivity
    ## --------------------
    ## per age
    if(length(dat$sel) == nt){
        dat$sel <- c(dat$sel, rep(tail(dat$sel,1), ntall - nt))
    }else if(length(dat$sel) != 1 && length(dat$sel) < nt){
        if(verbose) warning("Length of list with selectivity at age (dat$sel) is shorter than the historical period. Using last provided values for missing years and projection period.")
        dat$sel <- c(dat$sel, rep(tail(dat$sel,1), nt - length(dat$sel)))
    }


    ## Overwritting indices (accounting for projection period)
    dat$yvec <- rep(1:nyall, each = ns)
    dat$svec <- rep(1:ns, each = nyall)
    dat$s1vec <- seq(1, ntall, ns)

    ## parallel loop
    res <- parallel::mclapply(as.list(1:nrep), function(x){

        if(verbose) writeLines(paste0("Running replicate: ", x))

        ## set seed
        if(is.numeric(set$seed)) set.seed(set$seed + x)

        ## pop list with errors
        pop <- initpop(dat, set, dbg = dbg)
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
        setx$eF <- gen.noise(nysim, set$noiseF[1], set$noiseF[2], set$noiseF[3])
        setx$eR <- gen.noise(nysim, set$noiseR[1], set$noiseR[2], set$noiseR[3])
        setx$eM <- gen.noise(nysim, set$noiseM[1], set$noiseM[2], set$noiseM[3])
        setx$eH <- gen.noise(nysim, set$noiseH[1], set$noiseH[2], set$noiseH[3])
        setx$eW <- gen.noise(nysim, set$noiseW[1], set$noiseW[2], set$noiseW[3])
        setx$eR0 <- gen.noise(nysim, set$noiseR0[1], set$noiseR0[2], set$noiseR0[3])
        setx$eMat <- gen.noise(nysim, set$noiseMat[1], set$noiseMat[2], set$noiseMat[3])
        setx$eSel <- gen.noise(nysim, set$noiseSel[1], set$noiseSel[2], set$noiseSel[3])
        setx$eImp <- gen.noise(nysim, set$noiseImp[1], set$noiseImp[2], set$noiseImp[3])
        setx$eC <- gen.noise(nysim, set$noiseC[1], set$noiseC[2], set$noiseC[3])
        setx$eI <- list()
        for(i in 1:length(dat$surveyTimes)){
           setx$eI[[i]] <- gen.noise(nysim, set$noiseI[1], set$noiseI[2], set$noiseI[3])
        }
        setx$eE <- gen.noise(nysim, set$noiseE[1], set$noiseE[2], set$noiseE[3])

        ## loop
        for(i in 1:nhcrs){
            hcri <- hcrs[i]
            poptmp <- popListx[[i]]
            poptmp$tacs <- NULL
            for(y in 1:nysim){
                poptmp <- advancepop(dat = dat,
                                     hist = poptmp,
                                     set = setx,
                                     hcr = hcri,
                                     year = y,
                                     verbose = verbose,
                                     dbg = dbg)
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

    class(res2) <- "iamse"

    return(res2)
}
