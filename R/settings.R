#' @name check.set
#' @export
check.set <- function(set = NULL){

    ## empty list
    if(is.null(set)) set <- list()

    ## noise vectors (SD, rho, bias correction)
    if(is.null(set$noiseR)){
        set$noiseR <- c(0,0,0)
    }else if(length(set$noiseR) != 3){
        writeLines("'set$noiseR' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseF)){
        set$noiseF <- c(0,0,0)
    }else if(length(set$noiseF) != 3){
        writeLines("'set$noiseF' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseM)){
        set$noiseM <- c(0,0,0)
    }else if(length(set$noiseM) != 3){
        writeLines("'set$noiseM' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseH)){
        set$noiseH <- c(0,0,0)
    }else if(length(set$noiseH) != 3){
        writeLines("'set$noiseH' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseR0)){
        set$noiseR0 <- c(0,0,0)
    }else if(length(set$noiseR0) != 3){
        writeLines("'set$noiseR0' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseMat)){
        set$noiseMat <- c(0,0,0)
    }else if(length(set$noiseMat) != 3){
        writeLines("'set$noiseMat' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseSel)){
        set$noiseSel <- c(0,0,0)
    }else if(length(set$noiseSel) != 3){
        writeLines("'set$noiseSel' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseW)){
        set$noiseW <- c(0,0,0)
    }else if(length(set$noiseW) != 3){
        writeLines("'set$noiseW' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseImp)){
        set$noiseImp <- c(0,0,0)
    }else if(length(set$noiseImp) != 3){
        writeLines("'set$noiseImp' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseC)){
        set$noiseC <- c(0,0,0)
    }else if(length(set$noiseC) != 3){
        writeLines("'set$noiseC' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set[['noiseI']])){
        set$noiseI <- c(0,0,0)
    }else if(length(set$noiseI) != 3){
        writeLines("'set$noiseI' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseCmv)){
        set$noiseCmv <- c(0,0,0)
    }else if(length(set$noiseCmv) != 3){
        writeLines("'set$noiseCmv' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseImv)){
        set$noiseImv <- c(0,0,0)
    }else if(length(set$noiseImv) != 3){
        writeLines("'set$noiseImv' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor. Setting to c(0,0,0)!")
    }
    if(is.null(set$noiseE)){
        set$noiseE <- c(0,0,0)
    }else if(length(set$noiseE) != 3){
        writeLines("'set$noiseE' needs to be a vector with 3 values corresponding to: sd, rho, bias.cor (see gen.noise). Setting to c(0,0,0)!")
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
    if(is.null(set$hcr)) set$hcr <- c(def.hcr.ref(),def.hcr.ref(consF = "fmsy"))
    ## define constant catch rule
    def.hcr.conscat()
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

    ## return
    return(set)
}
