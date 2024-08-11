##' @name run.mse
##'
##' @importFrom parallel mclapply
##' @importFrom parallel detectCores
##'
##' @export
run.mse <- function(dat, set,
                    fest = NULL,
                    ncores=parallel::detectCores()-1,
                    verbose=TRUE){

    if(ncores > 1) verbose <- FALSE

    ## define constant catch (resort HCR if something not converging)
    def.hcr.conscat(set = set)

    ## Variables
    hcrs <- set$hcr
    hcrs2 <- sapply(hcrs,
                    function(x)
                        unlist(strsplit(as.character(x), "-"))[1])
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

    if(!is.null(fest)){
        if(!inherits(fest, "matrix")){
            fest <- t(as.matrix(fest))
        }
        if(nrow(fest) < nrep) stop("Provided fishing mortality has to have as many rows as requested replicates (set$nrep)! Please check.")
        if(any(is.null(set$errs)) ||
           length(set$errs$time) < nrep) stop("No errors provided in set or number of erros do not correspond to number of replicates (set$nrep)! Please check.")
        fhist.flag <- TRUE
        errs <- set$errs
        set$errs <- NULL
    }else{
        fhist.flag <- FALSE
    }

    ## Natural mortality
    ## --------------------
    if(!inherits(dat$M, "matrix")) stop(paste0("dat$M is not a matrix. It should be a matrix with ",ns," columns and at least ",ny," rows."))
    if(ncol(dat$M) != ns) stop(paste0("dat$M has ",ncol(dat$M), "columns. It should have ",ns," columns."))
    ## per season
    if(nrow(dat$M) == ny){
        dat$M <- rbind(dat$M, matrix(rep(tail(dat$M,1), nyall - ny),
                                     nrow = nyall - ny, ncol = ns))
        if(verbose) writeLines("No M provided for projection period. Using M in the last historical year for the projection period.")
    }else if(nrow(dat$M) < nyall){
        if(verbose) warning("Length of vector with natural mortality (dat$M) is shorter than the historical period. Using last provided natural mortality value for missing years and projection period.")
        dat$M <- rbind(dat$M, matrix(rep(tail(dat$M,ns), ny - nrow(dat$M)),
                                     nrow = ny-nrow(dat$M), ncol=ns))
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
    if(ncores > 1){
        res <- parallel::mclapply(as.list(1:nrep), function(x){

            if(verbose) writeLines(paste0("Running replicate: ", x))

            ## set seed
            if(is.numeric(set$seed)) set.seed(set$seed + x)

            setx <- set
            datx <- dat
            ## errors
            if(fhist.flag){
                datx$FM <- matrix(fest[x,]/datx$ns,
                                  datx$ny, datx$ns)
                setx$errs <- list()
                setx$errs$time <- errs$time[[x]]
                setx$errs$rep <- errs$rep[[x]]
            }

            ## pop list with errors
            pop <- initpop(datx, setx)
            ## add reference levels
            pop$refs <- refs
            popList <- vector("list", nhcrs)
            for(i in 1:nhcrs){
                popList[[i]] <- pop
            }

            repList <- vector("list", nhcrs)
            popListx <- popList

            ## setx$errs <- list()
            ## setx$errs$time <- get.errs(datx, setx, (ny+1):(ny+nysim), pop)
            ## setx$errs$rep <- pop$errs$rep

            ## loop
            for(i in 1:nhcrs){
                hcri <- hcrs[i]
                poptmp <- popListx[[i]]
                poptmp$tacs <- NULL
                for(y in 1:nysim){
                    poptmp <- advancepop(dat = datx,
                                         hist = poptmp,
                                         set = setx,
                                         hcr = hcri,
                                         year = y,
                                         verbose = verbose)
                }
                popListx[[i]] <- poptmp
                gc()
            }
            ## repList[[x]] <- popListx
            return(popListx)
        }, mc.cores = ncores)

    }else{

        res <- vector("list", nrep)
        for(x in 1:nrep){

            if(verbose) writeLines(paste0("Running replicate: ", x))

            ## set seed
            if(is.numeric(set$seed)) set.seed(set$seed + x)

            setx <- set
            datx <- dat
            ## errors
            if(fhist.flag){
                datx$FM <- matrix(fest[x,]/datx$ns,
                                  datx$ny, datx$ns)
                setx$errs <- list()
                setx$errs$time <- errs$time[[x]]
                setx$errs$rep <- errs$rep[[x]]
            }

            ## pop list with errors
            pop <- initpop(datx, setx)
            ## add reference levels
            pop$refs <- refs
            popList <- vector("list", nhcrs)
            for(i in 1:nhcrs){
                popList[[i]] <- pop
            }

            repList <- vector("list", nhcrs)
            popListx <- popList

            ## setx$errs <- list()
            ## setx$errs$time <- get.errs(datx, setx, (ny+1):(ny+nysim), popList[[i]])
            ## setx$errs$rep <- pop$errs$rep

            ## loop
            for(i in 1:nhcrs){
                hcri <- hcrs[i]
                poptmp <- popListx[[i]]
                poptmp$tacs <- NULL
                for(y in 1:nysim){
                    poptmp <- advancepop(dat = datx,
                                         hist = poptmp,
                                         set = setx,
                                         hcr = hcri,
                                         year = y,
                                         verbose = verbose)
                }
                popListx[[i]] <- poptmp
                gc()
            }
            ## repList[[x]] <- popListx
            res[[x]] <- popListx
        }
    }



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
