#' check.set
#' @export
check.set <- function(set = NULL){

    ## empty list
    if(is.null(set)) set <- list()

    ## noise vectors (SD, rho, bias correction)
    if(is.null(set$noise)){
        set$noise <- list()
    }

    ## for noise over time
    if(is.null(set$noise$time)){
        set$noise$time <- list()
    }
    if(is.null(set$noise$time$R)){
        set$noise$time$R <- c(0,0,0)
    }else if(length(set$noise$time$R) != 3){
        writeLines("'set$noise$time$R' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$F)){
        set$noise$time$F <- c(0,0,0)
    }else if(length(set$noise$time$F) != 3){
        writeLines("'set$noise$time$F' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$M)){
        set$noise$time$M <- c(0,0,0)
    }else if(length(set$noise$time$M) != 3){
        writeLines("'set$noise$time$M' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$H)){
        set$noise$time$H <- c(0,0,0)
    }else if(length(set$noise$time$H) != 3){
        writeLines("'set$noise$time$H' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$R0)){
        set$noise$time$R0 <- c(0,0,0)
    }else if(length(set$noise$time$R0) != 3){
        writeLines("'set$noise$time$R0' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$Alpha)){
        set$noise$time$Alpha <- c(0,0,0)
    }else if(length(set$noise$time$Alpha) != 3){
        writeLines("'set$noise$time$Alpha' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$Beta)){
        set$noise$time$Beta <- c(0,0,0)
    }else if(length(set$noise$time$Beta) != 3){
        writeLines("'set$noise$time$Beta' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$Mat)){
        set$noise$time$Mat <- c(0,0,0)
    }else if(length(set$noise$time$Mat) != 3){
        writeLines("'set$noise$time$Mat' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$Sel)){
        set$noise$time$Sel <- c(0,0,0)
    }else if(length(set$noise$time$Sel) != 3){
        writeLines("'set$noise$time$Sel' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$W)){
        set$noise$time$W <- c(0,0,0)
    }else if(length(set$noise$time$W) != 3){
        writeLines("'set$noise$time$W' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$Imp)){
        set$noise$time$Imp <- c(0,0,0)
    }else if(length(set$noise$time$Imp) != 3){
        writeLines("'set$noise$time$Imp' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$C)){
        set$noise$time$C <- c(0,0,0)
    }else if(length(set$noise$time$C) != 3){
        writeLines("'set$noise$time$C' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set[['noiseI']])){
        set$noise$time$I <- c(0,0,0)
    }else if(length(set$noise$time$I) != 3){
        writeLines("'set$noise$time$I' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$CmvA)){
        set$noise$time$CmvA <- c(0,0,0)
    }else if(length(set$noise$time$CmvA) != 3){
        writeLines("'set$noise$time$CmvA' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$CmvL)){
        set$noise$time$CmvL <- c(0,0,0)
    }else if(length(set$noise$time$CmvL) != 3){
        writeLines("'set$noise$time$CmvL' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$ImvA)){
        set$noise$time$ImvA <- c(0,0,0)
    }else if(length(set$noise$time$ImvA) != 3){
        writeLines("'set$noise$time$ImvA' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$ImvL)){
        set$noise$time$ImvL <- c(0,0,0)
    }else if(length(set$noise$time$ImvL) != 3){
        writeLines("'set$noise$time$ImvL' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$time$E)){
        set$noise$time$E <- c(0,0,0)
    }else if(length(set$noise$time$E) != 3){
        writeLines("'set$noise$time$E' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }

    ## noise between replicates
    if(is.null(set$noise$rep)){
        set$noise$rep <- list()
    }
    if(is.null(set$noise$rep$R)){
        set$noise$rep$R <- c(0,0,0)
    }else if(length(set$noise$rep$R) != 3){
        writeLines("'set$noise$rep$R' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$F)){
        set$noise$rep$F <- c(0,0,0)
    }else if(length(set$noise$rep$F) != 3){
        writeLines("'set$noise$rep$F' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$M)){
        set$noise$rep$M <- c(0,0,0)
    }else if(length(set$noise$rep$M) != 3){
        writeLines("'set$noise$rep$M' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$H)){
        set$noise$rep$H <- c(0,0,0)
    }else if(length(set$noise$rep$H) != 3){
        writeLines("'set$noise$rep$H' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$R0)){
        set$noise$rep$R0 <- c(0,0,0)
    }else if(length(set$noise$rep$R0) != 3){
        writeLines("'set$noise$rep$R0' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$Alpha)){
        set$noise$rep$Alpha <- c(0,0,0)
    }else if(length(set$noise$rep$Alpha) != 3){
        writeLines("'set$noise$rep$Alpha' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$Beta)){
        set$noise$rep$Beta <- c(0,0,0)
    }else if(length(set$noise$rep$Beta) != 3){
        writeLines("'set$noise$rep$Beta' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$Mat)){
        set$noise$rep$Mat <- c(0,0,0)
    }else if(length(set$noise$rep$Mat) != 3){
        writeLines("'set$noise$rep$Mat' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$Sel)){
        set$noise$rep$Sel <- c(0,0,0)
    }else if(length(set$noise$rep$Sel) != 3){
        writeLines("'set$noise$rep$Sel' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$W)){
        set$noise$rep$W <- c(0,0,0)
    }else if(length(set$noise$rep$W) != 3){
        writeLines("'set$noise$rep$W' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$Imp)){
        set$noise$rep$Imp <- c(0,0,0)
    }else if(length(set$noise$rep$Imp) != 3){
        writeLines("'set$noise$rep$Imp' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$C)){
        set$noise$rep$C <- c(0,0,0)
    }else if(length(set$noise$rep$C) != 3){
        writeLines("'set$noise$rep$C' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set[['noiseI']])){
        set$noise$rep$I <- c(0,0,0)
    }else if(length(set$noise$rep$I) != 3){
        writeLines("'set$noise$rep$I' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$CmvA)){
        set$noise$rep$CmvA <- c(0,0,0)
    }else if(length(set$noise$rep$CmvA) != 3){
        writeLines("'set$noise$rep$CmvA' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$CmvL)){
        set$noise$rep$CmvL <- c(0,0,0)
    }else if(length(set$noise$rep$CmvL) != 3){
        writeLines("'set$noise$rep$CmvL' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$ImvA)){
        set$noise$rep$ImvA <- c(0,0,0)
    }else if(length(set$noise$rep$ImvA) != 3){
        writeLines("'set$noise$rep$ImvA' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$ImvL)){
        set$noise$rep$ImvL <- c(0,0,0)
    }else if(length(set$noise$rep$ImvL) != 3){
        writeLines("'set$noise$rep$ImvL' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noise$rep$E)){
        set$noise$rep$E <- c(0,0,0)
    }else if(length(set$noise$rep$E) != 3){
        writeLines("'set$noise$rep$E' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }

    ## maximum F for baranov solution for F given TAC
    if(is.null(set$maxF)) set$maxF <- 5

    ## for estimation of ref levels
    if(is.null(set$refN)) set$refN <- 1e3
    if(is.null(set$refYears)){
        set$refYears <- 300
    }else if(set$refYears < 100){
        writeLines("set$refYears has to be at least 100 years, as the median surplus production over the last 50 years are used for reference estimation. Setting to 100.")
        set$refYears <- 100
    }
    if(is.null(set$refYearsMSY)){
        set$refYearsMSY <- 100
    }else if(set$refYearsMSY >= set$refYears){
        writeLines("set$refYearsMSY cannot be longer than 'set$refYears'. Setting to half of 'set$refYears'.")
        set$refYearsMSY <- floor(set$refYears / 2)
    }
    if(is.null(set$refMethod)){
        set$refMethod <- "mean"
    }else if(set$refYears < 100){
        writeLines("'set$refMethod' has to be mean or median. Setting to 'mean'.")
        set$refMethod <- "mean"
    }

    ## number of years available to assessment method
    if(is.null(set$nysim)) set$nysim <- 35
    if(is.null(set$nrep)) set$nrep <- 50

    ## Assessment
    if(is.null(set$assessmentTiming)) set$assessmentTiming <- 1
    if(is.null(set$assessmentInterval)) set$assessmentInterval <- 1
    if(is.null(set$assessmentIntYear)) set$assessmentIntYear <- 0

    ## HCR
    if(is.null(set$hcr)) set$hcr <- NULL ## c(def.hcr.ref(),def.hcr.ref(consF = "fmsy"))
    ## define constant catch rule
    ## def.hcr.conscat()
    if(is.null(set$stab)) set$stab <- FALSE

    ## burn in period
    if(is.null(set$burnin)) set$burnin <- 5e2

    ## seed
    if(is.null(set$seed)) set$seed <- NA

    ## SP type used for estimation of reference levels (and porduction curve)
    ## SP based on TSB (spType == 0) or on ESB (spType == 1)
    if(is.null(set$spType)) set$spType <- 0

    ## Record states at end of year or beginning?
    if(is.null(set$recordLast)) set$recordLast <- 0

    if(is.null(set$surveyBeforeRec)) set$surveyBeforeRec <- FALSE


    ## return
    return(set)
}
