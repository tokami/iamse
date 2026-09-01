## -------------------------------------------------------------------------
## Conditioning an iamse operating model on a fitted assessment model
## -------------------------------------------------------------------------


#' @name initdist.ad
#'
#' @title AD-safe initial age distribution
#'
#' @description Vectorised, automatic-differentiation-safe equivalent of
#'   [initdistR()]. `initdistR()` assigns into a plain numeric matrix, which
#'   drops the `advector` class and therefore cannot be taped by `RTMB`. This
#'   version returns the identical result from a closed-form expression.
#'
#' @details `initdistR()` builds `NAA[a, s] = R0 * spawning[s] *
#'   exp(-sum(ZAA[indage0:(a-1)]))`, keeps only the age-season combinations
#'   present at the end of the year, applies a plus-group correction to the
#'   last age and sums over seasons. Only the exponentiated cumulative
#'   mortality depends on `FM`, so the season mask can be precomputed with
#'   [initdist.weights()] and the remainder is a simple product.
#'
#' @param M natural mortality at age-season (length `asmax`)
#' @param FM fishing mortality at age-season (length `asmax`)
#' @param w season weights from [initdist.weights()]
#' @param asmax number of age-season groups
#' @param indage0 index of the recruitment age-season group
#' @param ns number of seasons
#' @param R0 unfished recruitment
#'
#' @return Numeric (or `advector`) vector of length `asmax`.
#'
#' @keywords internal
initdist.ad <- function(M, FM, w, asmax, indage0, ns, R0 = 1){
    ZAA <- M + FM
    ## mortality below the recruitment age-season group is masked out so that
    ## the cumulative sum starts at indage0
    pad <- as.numeric(seq_len(asmax) >= indage0)
    Zp <- ZAA * pad
    ## survival from indage0 up to the start of each age-season group
    cum <- exp(-(cumsum(Zp) - Zp))
    ## plus group correction, applied to the last age-season group only
    plus <- 1 / (1 - exp(-sum(ZAA[(asmax - ns + 1):asmax])))
    lastm <- as.numeric(seq_len(asmax) == asmax)
    cum * (w * R0) * (1 + lastm * (plus - 1))
}


#' @name initdist.weights
#'
#' @title Season weights for the initial age distribution
#'
#' @description Parameter-independent part of [initdist.ad()]: for every
#'   age-season group, the sum of the spawning proportions of the seasons in
#'   which that group is present at the end of the year.
#'
#' @param ns number of seasons
#' @param asmax number of age-season groups
#' @param indage0 index of the recruitment age-season group
#' @param spawning spawning proportion per season
#'
#' @return Numeric vector of length `asmax`.
#'
#' @keywords internal
initdist.weights <- function(ns, asmax, indage0, spawning){
    w <- rep(0, asmax)
    for(s in 1:ns){
        indi <- seq(ns + 2 - s + indage0 - 1, asmax, ns)
        w[indi] <- w[indi] + spawning[s]
    }
    ## last age-season group is kept for every season (see initdistR)
    w[asmax] <- sum(spawning)
    ## recruits are removed from the initial distribution (see initdistR)
    w[indage0] <- 0
    w
}


#' @name prep.hist
#'
#' @title Precompute the parameter-independent parts of the historic OM
#'
#' @description Collects everything [sim.hist()] needs that does not depend on
#'   fishing mortality, so that it can be computed once outside the automatic
#'   differentiation tape. This includes all noise-scaled life-history
#'   quantities per year and, crucially, the unfished spawning biomass per
#'   recruit (`SSBPR0`) per year and season: [get.ssb0()] is always evaluated
#'   at `FM = 0`, so it is independent of the quantity being estimated.
#'
#' @details The layout mirrors [initpop()] (`R/opmod.R`) exactly, including
#'   its use of `mean(errs$time$eX * errs$rep$eX)` during the burn-in, its use
#'   of `errs$time$eAlpha[1]` / `errs$time$eBeta[1]` (year 1, not year `y`) for
#'   the stock-recruitment shape parameters, and its use of `dat$beta` rather
#'   than `dat$recBeta` for the `beta` argument of [recfunc()]. The latter two
#'   are quirks of `initpop()`; they are reproduced here so that [sim.hist()]
#'   and [initpop()] agree exactly. Neither is used by the Beverton-Holt
#'   stock-recruitment relationship.
#'
#' @param dat stock data list, as returned by [check.dat()]
#' @param set settings list, as returned by [check.set()]
#' @param errs list with `time` and `rep` components, as returned by
#'   [get.errs()]
#'
#' @return A list of precomputed quantities consumed by [sim.hist()].
#'
#' @keywords internal
prep.hist <- function(dat, set, errs){

    ny <- dat$ny
    ns <- dat$ns
    asmax <- dat$asmax
    amax <- dat$amax + 1
    indage0 <- dat$indage0
    as2s <- dat$as2s

    weight <- dat$weight
    if(!inherits(weight, "matrix")) weight <- matrix(weight, asmax, ny)
    weightF <- dat$weightF
    if(!inherits(weightF, "matrix")) weightF <- matrix(weightF, asmax, ny)

    R0 <- check.par(dat$R0, ny)
    h <- check.par(dat$h, ny)
    alpha <- check.par(dat$recAlpha, ny)
    ## NOTE: initpop() overwrites the recruitment beta with dat$beta (the
    ## survey hyperstability parameter) before the main loop; mirrored here.
    beta <- dat$beta

    mselFlag <- inherits(dat$Msel, "list") && length(dat$Msel) > 1
    selFlag <- inherits(dat$sel, "list") && length(dat$sel) > 1

    et <- errs$time
    er <- errs$rep

    ## ---- per-year quantities ------------------------------------------
    MAA <- sely <- maty <- weighty <- weightFy <- vector("list", ny)
    hy <- R0y <- alphay <- betay <- eRy <- rep(NA_real_, ny)
    for(y in 1:ny){
        selInd <- ifelse(selFlag, y, 1)
        mselInd <- ifelse(mselFlag, y, 1)
        sely[[y]] <- as.numeric(t(dat$sel[[selInd]])) * et$eSel[y] * er$eSel
        MAA[[y]] <- dat$M[y, as2s] * as.numeric(t(dat$Msel[[mselInd]])) *
            et$eM[y] * er$eM
        maty[[y]] <- as.numeric(t(dat$mat)) * et$eMat[y] * er$eMat
        weighty[[y]] <- as.numeric(t(weight[, y])) * et$eW[y] * er$eW
        weightFy[[y]] <- as.numeric(t(weightF[, y])) * et$eW[y] * er$eW
        hy[y] <- h[y] * et$eH[y] * er$eH
        R0y[y] <- R0[y] * et$eR0[y] * er$eR0
        alphay[y] <- ifelse(length(alpha) == 0, NA_real_,
                            alpha[y] * et$eAlpha[1] * er$eAlpha)
        betay[y] <- ifelse(length(beta) == 0, NA_real_,
                           beta[min(y, length(beta))] * et$eBeta[1] * er$eBeta)
        eRy[y] <- et$eR[y] * er$eR
    }

    ## ---- unfished spawning biomass per recruit, per year and season ----
    ## get.ssb0() is called with FM = 0 in initpop(), so this does not depend
    ## on the fishing mortality being estimated.
    ssbpr0 <- matrix(NA_real_, ny, ns)
    for(y in 1:ny){
        for(s in 1:ns){
            if(dat$spawning[s] > 0){
                r0i <- R0y[y] * eRy[y]
                ssbpr0[y, s] <- get.ssb0(MAA[[y]], maty[[y]], weighty[[y]],
                                         fecun = 1, asmax, ns, dat$spawning,
                                         indage0 = indage0, R0 = r0i,
                                         season = s, FM = 0) / r0i
            }
        }
    }

    ## ---- burn-in quantities -------------------------------------------
    eFbi <- mean(et$eF * er$eF)
    eMbi <- mean(et$eM * er$eM)
    eRbi <- mean(et$eR * er$eR)
    eR0bi <- mean(et$eR0 * er$eR0)
    bi <- list(
        eF = eFbi,
        sel = as.numeric(t(dat$sel[[1]])) * et$eSel[1] * er$eSel,
        mat = as.numeric(t(dat$mat)) * et$eMat[1] * er$eMat,
        weight = as.numeric(t(weight[, 1])) * et$eW[1] * er$eW,
        h = h[1] * et$eH[1] * er$eH,
        R0 = R0[1] * eR0bi,
        eR = eRbi,
        alpha = if(length(alpha) == 0) NA_real_ else alpha[1] * et$eAlpha[1] * er$eAlpha,
        beta = if(length(beta) == 0) NA_real_ else beta[1] * et$eBeta[1] * er$eBeta)
    bi$M <- dat$M[1, as2s] * as.numeric(t(dat$Msel[[1]])) * eMbi
    bi$ssbpr0 <- rep(NA_real_, ns)
    for(s in 1:ns){
        if(dat$spawning[s] > 0){
            r0i <- bi$R0 * bi$eR
            bi$ssbpr0[s] <- get.ssb0(bi$M, bi$mat, bi$weight, fecun = 1,
                                     asmax, ns, dat$spawning,
                                     indage0 = indage0, R0 = r0i,
                                     season = s, FM = 0) / r0i
        }
    }

    ## initpop() clamps negative recruitment. For the stock-recruitment
    ## relationships below, recruitment is guaranteed non-negative as long as
    ## steepness stays below one, so the (tape-expanding) clamp can be skipped.
    clamp <- !(dat$SR %in% c("bevholt", "average", "hockey-stick", "ricker")) ||
        any(c(hy, bi$h) >= 1, na.rm = TRUE)

    list(ny = ny, ns = ns, asmax = asmax, indage0 = indage0,
         spawning = dat$spawning, pzbm = dat$pzbm, SR = dat$SR,
         bp = dat$bp, recGamma = dat$recGamma,
         burnin = if(is.null(set$burnin)) 0 else set$burnin,
         MAA = MAA, sely = sely, maty = maty,
         weighty = weighty, weightFy = weightFy,
         hy = hy, R0y = R0y, alphay = alphay, betay = betay, eRy = eRy,
         ssbpr0 = ssbpr0, bi = bi, clamp = clamp,
         eF = et$eF * er$eF,
         idw = initdist.weights(ns, asmax, indage0, dat$spawning))
}


#' @name sim.hist
#'
#' @title Simulate the historic period of the operating model
#'
#' @description Automatic-differentiation-safe forward run of the
#'   age/season-structured operating model over the historic period, given an
#'   annual fishing mortality vector. Reproduces [initpop()] (`R/opmod.R`)
#'   exactly for `set$recordLast == 1`, but returns only the aggregated
#'   quantities needed for conditioning and can be taped by `RTMB`.
#'
#' @param pre precomputed quantities from [prep.hist()]
#' @param FM numeric (or `advector`) vector of length `pre$ny` with **annual**
#'   fishing mortality. Seasonal fishing mortality is `FM[y] / pre$ns`, the
#'   same convention [run.mse()] uses for its `fest` argument.
#'
#' @return A list with `esb` (exploitable stock biomass at the end of each
#'   year), `cw` (annual catch in weight), `ssb` (spawning stock biomass at
#'   the end of each year) and `rec` (annual recruitment).
#'
#' @keywords internal
sim.hist <- function(pre, FM){

    ny <- pre$ny
    ns <- pre$ns
    asmax <- pre$asmax
    indage0 <- pre$indage0
    spawning <- pre$spawning
    pzbm <- pre$pzbm
    bi <- pre$bi

    ## results are accumulated in lists rather than numeric vectors: assigning
    ## an advector into a numeric vector strips the AD class once the package
    ## is byte compiled (see initdist.ad)
    esb <- ssb <- cw <- rec <- vector("list", ny)

    ## ---- initial age distribution and burn-in --------------------------
    Fbi <- FM[1] / ns * bi$sel * bi$eF
    Zbi <- bi$M + Fbi
    NAAS <- initdist.ad(bi$M, Fbi, pre$idw, asmax, indage0, ns,
                        bi$R0 * bi$eR)
    if(is.numeric(pre$burnin) && pre$burnin > 0){
        for(y in 1:pre$burnin){
            for(s in 1:ns){
                if(spawning[s] > 0){
                    ssbi <- sum(NAAS * bi$weight * bi$mat * exp(-pzbm * Zbi))
                    recbi <- spawning[s] *
                        recfunc(h = bi$h, SSBPR0 = bi$ssbpr0[s], SSB = ssbi,
                                R0 = bi$R0, method = pre$SR, bp = pre$bp,
                                beta = bi$beta, gamma = pre$recGamma,
                                alpha = bi$alpha) * bi$eR
                    ## initpop() clamps negative recruitment (can occur when
                    ## noise pushes steepness above 1)
                    if(pre$clamp) recbi <- max(recbi, 1e-10)
                    NAAS[indage0] <- recbi
                }
                NAAS <- NAAS * exp(-Zbi)
                NAAS <- age.naas(NAAS, asmax, indage0)
            }
        }
    }

    ## ---- historic period -----------------------------------------------
    for(y in 1:ny){
        sely <- pre$sely[[y]]
        maty <- pre$maty[[y]]
        weighty <- pre$weighty[[y]]
        weightFy <- pre$weightFy[[y]]
        MAA <- pre$MAA[[y]]
        FAA <- FM[y] / ns * pre$eF[y] * sely
        ZAA <- MAA + FAA
        R0y <- pre$R0y[y]
        cwy <- 0
        recy <- 0
        esby <- 0
        ssby <- 0
        for(s in 1:ns){
            if(spawning[s] > 0){
                ssbi <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))
                reci <- spawning[s] *
                    recfunc(h = pre$hy[y], SSBPR0 = pre$ssbpr0[y, s],
                            SSB = ssbi, R0 = R0y, method = pre$SR,
                            bp = pre$bp, beta = pre$betay[y],
                            gamma = pre$recGamma,
                            alpha = pre$alphay[y]) * pre$eRy[y]
                if(pre$clamp) reci <- max(reci, 1e-10)
                recy <- recy + reci
                NAAS[indage0] <- reci
            }
            cwy <- cwy + sum(baranov(FAA, MAA, NAAS) * weightFy)
            NAAS <- NAAS * exp(-ZAA)
            if(s == ns){
                esby <- sum(NAAS * weighty * sely)
                ssby <- sum(NAAS * weighty * maty * exp(-pzbm * ZAA))
            }
            NAAS <- age.naas(NAAS, asmax, indage0)
        }
        cw[[y]] <- cwy
        rec[[y]] <- recy
        esb[[y]] <- esby
        ssb[[y]] <- ssby
    }

    list(esb = do.call("c", esb), cw = do.call("c", cw),
         ssb = do.call("c", ssb), rec = do.call("c", rec))
}


#' @name age.naas
#'
#' @title Advance numbers-at-age by one season
#'
#' @description Continuous ageing step of the operating model: every
#'   age-season group moves up one group, the last two are pooled into the
#'   plus group and the recruitment group is emptied. Identical to the
#'   in-place loop in [initpop()] but written as a single concatenation so
#'   that it can be taped by `RTMB` and stays cheap in tape size.
#'
#' @param NAAS numbers per age-season group
#' @param asmax number of age-season groups
#' @param indage0 index of the recruitment age-season group
#'
#' @return Vector of length `asmax`.
#'
#' @keywords internal
age.naas <- function(NAAS, asmax, indage0){
    out <- c(NAAS[1], NAAS[1:(asmax-2)], NAAS[asmax] + NAAS[asmax-1])
    out[indage0] <- 0
    out
}


#' @name annual.mean
#'
#' @title Aggregate a sub-annual time series to annual values
#'
#' @description Local replacement for the unexported `spict:::annual()`. For
#'   every calendar year covered by `time`, returns the mean (or another
#'   summary) of the values falling within that year.
#'
#' @param time numeric vector of times
#' @param value numeric vector of values, same length as `time`
#' @param fun summary function, default [mean()]
#'
#' @return A list with `time` (integer years) and `value`.
#'
#' @keywords internal
annual.mean <- function(time, value, fun = mean){
    keep <- !is.na(value) & !is.na(time)
    yr <- floor(time[keep])
    agg <- tapply(value[keep], yr, fun)
    list(time = as.numeric(names(agg)), value = as.numeric(agg))
}


#' @name check.spict.dat
#'
#' @title Check that a stock list and a SPiCT fit can be aligned
#'
#' @description Verifies that an `iamse` stock list (`dat`) and settings list
#'   (`set`) are compatible with a fitted `spict` object, and returns the
#'   assessment years the operating model's historic period maps onto.
#'
#' @details The operating model's historic period is aligned with the **last**
#'   `dat$ny` years of the assessment. This is the only alignment the
#'   conditioning supports, and it is checked rather than assumed: the
#'   prototype this replaces cut the assessment series with
#'   `x[(length(x) - dat$ny + 1):length(x)]`, which silently *drops* elements
#'   through negative indexing when `dat$ny` is larger than the assessment.
#'
#' @param dat stock data list, as returned by [check.dat()]
#' @param set settings list, as returned by [check.set()]
#' @param fit a fitted `spict` object
#' @param verbose print messages about the alignment
#'
#' @return Invisibly, a list with `years` (the assessment years used),
#'   `ind` (their positions in the assessment's annual series) and `nya` (the
#'   number of assessment years available).
#'
#' @export
check.spict.dat <- function(dat, set = NULL, fit, verbose = TRUE){

    if(!inherits(fit, "spictcls")) stop("'fit' is not a fitted spict object.")
    if(is.null(fit$opt)) stop("'fit' does not contain an optimisation result. Did fit.spict() fail?")
    if(is.null(set)) set <- check.set()

    ny <- dat$ny
    ns <- dat$ns

    ## ---- catch observations must be annual -----------------------------
    dtc <- unique(fit$inp$dtc)
    if(length(dtc) != 1 || !isTRUE(all.equal(dtc, 1)))
        stop(paste0("The spict fit uses catch observations over intervals of ",
                    paste(signif(dtc, 3), collapse = ", "),
                    " year(s). Conditioning compares annual catches and ",
                    "requires annual catch observations (inp$dtc == 1)."))

    yrsA <- unique(floor(fit$inp$timeC))
    if(length(yrsA) != length(fit$inp$timeC))
        stop("The spict fit has more than one catch observation per year.")
    nya <- length(yrsA)

    ## ---- historic period must fit into the assessment ------------------
    if(ny > nya)
        stop(paste0("dat$ny (", ny, ") is larger than the number of years in ",
                    "the spict fit (", nya, "). Shorten the historic period ",
                    "of the operating model or extend the assessment."))
    ind <- (nya - ny + 1):nya
    if(verbose && ny < nya)
        message("Aligning the ", ny, " historic years of the operating model ",
                "with the last ", ny, " of the ", nya,
                " assessment years (", yrsA[ind[1]], "-", yrsA[nya], ").")

    ## ---- annual dimensions of the stock list ---------------------------
    chk <- function(n, what)
        if(!(n %in% c(1, ny)))
            stop(paste0(what, " has length/dimension ", n,
                        ", which is neither 1 nor dat$ny (", ny, ")."))
    if(!inherits(dat$M, "matrix") || ncol(dat$M) != ns)
        stop(paste0("dat$M has to be a matrix with ", ns, " columns."))
    if(nrow(dat$M) < ny)
        stop(paste0("dat$M has ", nrow(dat$M), " rows but dat$ny is ", ny, "."))
    if(!inherits(dat$FM, "matrix") || any(dim(dat$FM) != c(ny, ns)))
        stop(paste0("dat$FM has to be a matrix with dimensions ", ny, " x ", ns, "."))
    chk(length(dat$sel), "dat$sel")
    chk(length(dat$Msel), "dat$Msel")
    if(inherits(dat$weight, "matrix")) chk(ncol(dat$weight), "dat$weight")
    if(inherits(dat$weightF, "matrix")) chk(ncol(dat$weightF), "dat$weightF")

    ## ---- reference points ----------------------------------------------
    if(!"ref" %in% names(dat))
        stop("Reference levels have to be part of dat. Use est.ref.levels() or est.ref.levels.stochastic().")
    for(r in c("ESBmsy", "MSY", "Fmsy"))
        if(!r %in% names(dat$ref) || !is.finite(dat$ref[[r]][1]))
            stop(paste0("dat$ref$", r, " is missing or not finite."))

    ## ---- the OM's ESB has to be an end-of-year quantity ------------------
    if(!identical(as.numeric(set$recordLast), 1))
        stop("set$recordLast has to be 1: the operating model's exploitable biomass must be recorded at the end of the year to be comparable with the biomass of a spict fit.")

    invisible(list(years = yrsA[ind], ind = ind, nya = nya))
}


#' @name get.spict.targets
#'
#' @title Extract annual conditioning targets from a SPiCT fit
#'
#' @description Pulls the quantities the operating model is conditioned on out
#'   of a fitted `spict` object, on an annual time step aligned with the
#'   operating model's historic period.
#'
#' @details Relative biomass is taken at the **last sub-annual time step of
#'   each year**, which is the point a `spict` fit is directly comparable with
#'   the operating model's `ESBfinal` (recorded after the decay of the last
#'   season, see [initpop()]). Relative fishing mortality is a within-year
#'   mean and is used only to start the optimiser and for diagnostics.
#'
#' @param fit a fitted `spict` object
#' @param dat stock data list
#' @param set settings list
#' @param verbose print messages about the alignment
#'
#' @return A list with `bbmsy`, `ffmsy`, `cmsy` (all of length `dat$ny`),
#'   the assessment `years`, and the `spict` quantities `sdf`, `MSY`, `Fmsy`
#'   and `Bmsy`.
#'
#' @export
get.spict.targets <- function(fit, dat, set = NULL, verbose = TRUE){

    ali <- check.spict.dat(dat, set, fit, verbose = verbose)

    gp <- function(nm, exp. = TRUE){
        v <- try(spict::get.par(nm, fit, exp = exp.), silent = TRUE)
        if(inherits(v, "try-error")) return(NA_real_)
        as.numeric(v[, 2])
    }

    ind <- fit$inp$indest
    tt <- fit$inp$time[ind]

    ## relative biomass at the end of each year
    bb <- spict::get.par("logBBmsy", fit, exp = TRUE)[ind, 2]
    ieoy <- sapply(ali$years, function(y){
        w <- which(tt < y + 1 - 1e-8)
        if(length(w) == 0) NA_integer_ else max(w)
    })
    if(any(is.na(ieoy)))
        stop("Could not locate end-of-year biomass estimates for all assessment years.")
    bbmsy <- as.numeric(bb[ieoy])

    ## relative fishing mortality as a within-year mean
    ff <- spict::get.par("logFFmsy", fit, exp = TRUE)[ind, 2]
    fa <- annual.mean(tt, as.numeric(ff))
    ffmsy <- fa$value[match(ali$years, fa$time)]

    ## predicted catch relative to MSY
    msy <- gp("logMSY")
    if(!is.finite(msy)) msy <- gp("logMSYd")
    cp <- spict::get.par("logCpred", fit, exp = TRUE)[, 2]
    cp <- as.numeric(cp[fit$inp$timeCpred %in% fit$inp$timeC])
    cmsy <- cp[ali$ind] / msy

    list(bbmsy = bbmsy, ffmsy = ffmsy, cmsy = cmsy,
         years = ali$years, nya = ali$nya,
         sdf = gp("logsdf"), MSY = msy,
         Fmsy = { x <- gp("logFmsy"); if(is.finite(x)) x else gp("logFmsyd") },
         Bmsy = { x <- gp("logBmsy"); if(is.finite(x)) x else gp("logBmsyd") })
}


#' @name zero.noise
#'
#' @title Set all noise standard deviations to zero
#'
#' @param set settings list
#'
#' @return `set` with every noise standard deviation set to zero.
#'
#' @keywords internal
zero.noise <- function(set){
    for(lev in c("time", "rep")){
        for(nm in names(set$noise[[lev]])){
            v <- set$noise[[lev]][[nm]]
            if(is.numeric(v) && length(v) >= 1){
                v[1] <- 0
                set$noise[[lev]][[nm]] <- v
            }
        }
    }
    set
}


#' @name condition.on.spict
#'
#' @title Condition an operating model on a SPiCT assessment
#'
#' @description Finds the historic fishing mortality that makes the
#'   exploitable stock biomass of the `iamse` operating model track the
#'   biomass of a fitted `spict` assessment as closely as possible.
#'   Selectivity and all life-history parameters are taken as given; fishing
#'   mortality is the only quantity estimated.
#'
#' @details The operating model's exploitable biomass at the end of each year
#'   (`ESBfinal`, see [initpop()]) relative to `dat$ref$ESBmsy` is fitted to
#'   the assessment's relative biomass. Catch is simulated and returned, but
#'   is **not** part of the objective: for a given operating model, fishing
#'   mortality determines biomass and catch jointly, and the age-structured
#'   production is not the assessment's surplus-production curve, so the two
#'   cannot in general both be matched. The returned catch trajectory is the
#'   honest diagnostic of how far apart the two models are.
#'
#'   Two parameterisations of fishing mortality are available. `method = "rw"`
#'   estimates one fishing mortality per year, penalised by a random walk
#'   whose standard deviation is taken from the assessment (`logsdf`).
#'   `method = "scale"` keeps the shape of the assessment's fishing mortality
#'   and estimates a single level parameter; it is the robust fallback when
#'   the free fit is unstable.
#'
#'   With `nrep > 1` a separate fishing mortality trajectory is estimated for
#'   each replicate, conditional on that replicate's process noise (mainly
#'   recruitment deviations). Feed the result straight into [run.mse()]:
#'
#'   \preformatted{
#'   res <- condition.on.spict(dat, set, fit, nrep = set$nrep)
#'   set$noise$time$F <- c(0, 0, 1)   ## F is already conditioned
#'   set$errs <- res$errs
#'   mse <- run.mse(dat, set, fest = res$FM)
#'   }
#'
#' @param dat stock data list, as returned by [check.dat()], including
#'   reference levels from [est.ref.levels.stochastic()]
#' @param set settings list, as returned by [check.set()]
#' @param fit a fitted `spict` object
#' @param nrep number of replicates. `nrep = 1` (default) uses noise-free
#'   deviations and gives a single deterministic fit. `nrep > 1` draws one set
#'   of deviations per replicate and returns them in `errs`.
#' @param method `"rw"` (default) or `"scale"`, see Details
#' @param target `"bbmsy"` (default) matches relative biomass as it is;
#'   `"shape"` matches only the shape of the biomass trajectory by dividing
#'   both series by their geometric mean. Use `"shape"` when the operating
#'   model and the assessment disagree strongly about `Bmsy / B0`, which
#'   `condition.on.spict()` reports.
#' @param sdB observation standard deviation of log relative biomass
#' @param sdF standard deviation of the random walk on log fishing mortality.
#'   Defaults to the assessment's `logsdf`.
#' @param sdF1 standard deviation of the log-scale prior anchoring the first
#'   year's fishing mortality at the assessment's `F/Fmsy * dat$ref$Fmsy`
#' @param est.sdB estimate `sdB`? If `TRUE`, fishing mortality becomes a
#'   random effect and `sdB` is estimated by maximum marginal likelihood.
#' @param rescale.R0 scale `dat$R0` so that the operating model's `MSY`
#'   matches the assessment's. Exact for a Beverton-Holt operating model, and
#'   it changes neither fishing mortality nor any relative quantity, but
#'   reference levels have to be re-estimated afterwards.
#' @param maxtries number of restarts from jittered starting values for
#'   replicates that do not converge
#' @param gr.tol convergence tolerance on the maximum absolute gradient
#' @param verbose print progress
#'
#' @return An object of class `iamse.cond`: a list with
#'   \describe{
#'     \item{FM}{`nrep x dat$ny` matrix of annual fishing mortality, ready for
#'       the `fest` argument of [run.mse()]}
#'     \item{errs}{list of per-replicate deviations, for `set$errs`}
#'     \item{esb, cw, ssb, rec}{`nrep x dat$ny` matrices of the conditioned
#'       operating model}
#'     \item{obs}{the assessment targets, see [get.spict.targets()]}
#'     \item{conv}{per-replicate convergence information}
#'     \item{dat}{`dat`, with `R0` rescaled if `rescale.R0 = TRUE`}
#'     \item{bmsy.frac}{`Bmsy` relative to the unfished level in the
#'       operating model (`ESBmsy/ESB0`) and in the assessment (`Bmsy/K`).
#'       A large difference means the two models disagree about where `Bmsy`
#'       sits, and matching `B/Bmsy` then transfers the assessment's perceived
#'       status rather than its biomass}
#'   }
#'
#' @seealso [get.spict.targets()], [check.spict.dat()], [run.mse()]
#'
#' @export
condition.on.spict <- function(dat, set = NULL, fit, nrep = 1,
                               method = c("rw", "scale"),
                               target = c("bbmsy", "shape"),
                               sdB = 0.1, sdF = NULL, sdF1 = 0.5,
                               est.sdB = FALSE, rescale.R0 = FALSE,
                               maxtries = 3, gr.tol = 1e-3,
                               verbose = TRUE){

    if(!requireNamespace("RTMB", quietly = TRUE))
        stop("Package 'RTMB' is required for condition.on.spict(). Please install it.")
    if(!requireNamespace("spict", quietly = TRUE))
        stop("Package 'spict' is required for condition.on.spict(). Please install it.")

    method <- match.arg(method)
    target <- match.arg(target)
    if(is.null(set)) set <- check.set()

    obs <- get.spict.targets(fit, dat, set, verbose = verbose)
    ny <- dat$ny
    esbmsy <- dat$ref$ESBmsy[1]
    fmsy <- dat$ref$Fmsy[1]

    if(is.null(sdF)) sdF <- obs$sdf
    if(!is.finite(sdF) || sdF <= 0)
        stop("Could not take a random walk standard deviation from the spict fit (logsdf). Please supply 'sdF'.")

    ## warn about F noise being layered on top of a conditioned trajectory
    if(isTRUE(set$noise$time$F[1] > 0))
        warning("set$noise$time$F[1] is ", signif(set$noise$time$F[1], 3),
                ". The historic fishing mortality returned here is already ",
                "conditioned on the assessment; set set$noise$time$F <- c(0, 0, 1) ",
                "before run.mse() to avoid adding further noise on top of it.")

    ## ---- deviations per replicate ---------------------------------------
    seti <- if(nrep == 1) zero.noise(set) else set
    errs <- list(time = vector("list", nrep), rep = vector("list", nrep))
    for(i in 1:nrep){
        if(is.numeric(set$seed)) set.seed(set$seed + i)
        seti$errs <- list()
        errs$time[[i]] <- get.errs(dat, seti, 1:ny)
        errs$rep[[i]] <- get.errs(dat, seti, 1, rep = TRUE)
    }

    ## ---- targets ---------------------------------------------------------
    obsB <- obs$bbmsy
    logObs <- log(obsB)
    if(target == "shape") logObs <- logObs - mean(logObs)

    ## starting values: the assessment's F/Fmsy on the operating model's scale
    logF0 <- log(pmax(obs$ffmsy * fmsy, 1e-6))
    logFpat <- logF0 - mean(logF0)

    FM <- esb <- cw <- ssb <- rec <- matrix(NA_real_, nrep, ny)
    conv <- data.frame(rep = 1:nrep, converged = FALSE,
                       objective = NA_real_, maxgrad = NA_real_,
                       tries = NA_integer_, sdB = NA_real_)

    for(i in 1:nrep){

        pre <- prep.hist(dat, set, list(time = errs$time[[i]], rep = errs$rep[[i]]))

        jnll <- function(p){
            ret <- 0
            if(method == "rw"){
                logF <- p$logF
                ret <- ret - sum(RTMB::dnorm(logF[-1], logF[-ny], sdF, log = TRUE))
                ret <- ret - RTMB::dnorm(logF[1], logF0[1], sdF1, log = TRUE)
            }else{
                logF <- logFpat + p$logScaleF
            }
            fm <- exp(logF)
            sim <- sim.hist(pre, fm)
            logPred <- log(sim$esb / esbmsy)
            if(target == "shape") logPred <- logPred - mean(logPred)
            ret <- ret - sum(RTMB::dnorm(logObs, logPred, exp(p$logSdB), log = TRUE))
            RTMB::ADREPORT(sim$esb)
            RTMB::ADREPORT(sim$cw)
            RTMB::ADREPORT(sim$ssb)
            RTMB::ADREPORT(sim$rec)
            RTMB::ADREPORT(fm)
            ret
        }

        for(k in 1:maxtries){
            jit <- if(k == 1) 0 else stats::rnorm(ny, 0, 0.2 * (k - 1))
            pars <- if(method == "rw")
                        list(logF = logF0 + jit, logSdB = log(sdB))
                    else
                        list(logScaleF = mean(logF0) + jit[1], logSdB = log(sdB))
            map <- list()
            if(!est.sdB) map$logSdB <- factor(NA)
            rnd <- if(est.sdB && method == "rw") "logF" else NULL

            obj <- try(RTMB::MakeADFun(jnll, pars, map = map, random = rnd,
                                       silent = TRUE), silent = TRUE)
            if(inherits(obj, "try-error")) next
            opt <- try(stats::nlminb(obj$par, obj$fn, obj$gr,
                                     control = list(iter.max = 1e3, eval.max = 1e3)),
                       silent = TRUE)
            if(inherits(opt, "try-error")) next
            mg <- max(abs(obj$gr(opt$par)))
            ok <- opt$convergence == 0 && is.finite(mg) && mg < gr.tol
            if(ok && !is.null(rnd)){
                sdr <- try(RTMB::sdreport(obj), silent = TRUE)
                ok <- !inherits(sdr, "try-error") && isTRUE(sdr$pdHess)
            }
            if(ok) break
        }

        conv$tries[i] <- k
        conv$converged[i] <- ok
        if(!ok) next
        conv$objective[i] <- opt$objective
        conv$maxgrad[i] <- mg
        conv$sdB[i] <- exp(if(est.sdB) opt$par[["logSdB"]] else log(sdB))

        sim <- sim.hist(pre, if(method == "rw")
                                 exp(obj$env$last.par.best[names(obj$env$last.par.best) == "logF"])
                             else
                                 exp(logFpat + obj$env$last.par.best[["logScaleF"]]))
        FM[i, ] <- if(method == "rw")
                       exp(obj$env$last.par.best[names(obj$env$last.par.best) == "logF"])
                   else
                       exp(logFpat + obj$env$last.par.best[["logScaleF"]])
        esb[i, ] <- sim$esb
        cw[i, ] <- sim$cw
        ssb[i, ] <- sim$ssb
        rec[i, ] <- sim$rec

        if(verbose) cat("\r  conditioned replicate ", i, " / ", nrep, "   ", sep = "")
    }
    if(verbose && nrep > 0) cat("\n")

    nok <- sum(conv$converged)
    if(nok < nrep)
        stop(paste0("Only ", nok, " of ", nrep, " replicates converged. ",
                    "Try method = \"scale\", a larger 'sdB', more 'maxtries', ",
                    "or check that the operating model can reach the biomass ",
                    "levels of the assessment."))

    ## ---- how far apart are the two models about where Bmsy sits? --------
    ## The two models normalise by their own Bmsy, so a large difference in
    ## Bmsy/B0 means matching B/Bmsy transfers the assessment's perception of
    ## stock status into the operating model rather than its biomass as such.
    omBmsy <- if(all(c("ESBmsy", "ESB0") %in% names(dat$ref)) &&
                 is.finite(dat$ref$ESB0[1]) && dat$ref$ESB0[1] > 0)
                  dat$ref$ESBmsy[1] / dat$ref$ESB0[1] else NA_real_
    spBmsy <- {
        k <- try(as.numeric(spict::get.par("logK", fit, exp = TRUE)[, 2]), silent = TRUE)
        if(inherits(k, "try-error") || !is.finite(k) || !is.finite(obs$Bmsy))
            NA_real_ else obs$Bmsy / k
    }
    lvl <- c(om = omBmsy, spict = spBmsy)
    if(verbose && all(is.finite(lvl))){
        message("Bmsy relative to the unfished level: ", signif(omBmsy, 3),
                " in the operating model (ESBmsy/ESB0) vs ", signif(spBmsy, 3),
                " in the assessment (Bmsy/K).",
                if(abs(log(omBmsy / spBmsy)) > log(1.5))
                    " These differ substantially, so B/Bmsy is not a common currency: the conditioned operating model reproduces the assessment's perceived status, not its absolute biomass. Compare the catch diagnostic, and consider target = \"shape\"." else "")
    }

    ## ---- optional rescaling of R0 ----------------------------------------
    if(rescale.R0){
        if(!is.finite(obs$MSY)) stop("The spict fit does not provide a finite MSY, cannot rescale R0.")
        fac <- obs$MSY / dat$ref$MSY[1]
        dat$R0 <- dat$R0 * fac
        message("dat$R0 scaled by ", signif(fac, 4),
                " so that the operating model's MSY matches the assessment's. ",
                "Reference levels have to be re-estimated (est.ref.levels.stochastic()); ",
                "fishing mortality and all relative quantities are unchanged.")
    }

    res <- list(FM = FM, errs = errs, esb = esb, cw = cw, ssb = ssb, rec = rec,
                obs = obs, conv = conv, dat = dat,
                method = method, target = target,
                sdF = sdF, bmsy.frac = lvl,
                esbmsy = esbmsy, fmsy = fmsy, msy = dat$ref$MSY[1])
    class(res) <- c("iamse.cond", "list")
    res
}


#' @name plotiamse.conditioning
#'
#' @title Diagnostic plots for a conditioned operating model
#'
#' @description Six panels comparing the conditioned operating model with the
#'   assessment it was conditioned on: fishing mortality (absolute and
#'   relative to `Fmsy`), exploitable biomass (absolute and relative to
#'   `ESBmsy`), catch relative to `MSY`, and recruitment.
#'
#'   The catch panel is the one to look at: catch is **not** part of the
#'   objective, so it shows how far the age-structured production of the
#'   operating model is from the assessment's surplus-production curve.
#'
#' @param res an object of class `iamse.cond` from [condition.on.spict()]
#' @param years optional vector of years for the x axis
#' @param col colour of the operating model replicates
#' @param col.assess colour of the assessment
#'
#' @return Invisibly `NULL`, called for its side effect.
#'
#' @export
plotiamse.conditioning <- function(res, years = NULL,
                                   col = "dodgerblue2",
                                   col.assess = "darkorange2"){

    stopifnot(inherits(res, "iamse.cond"))
    ny <- ncol(res$FM)
    if(is.null(years)) years <- res$obs$years
    nrep <- nrow(res$FM)

    panel <- function(om, ass, ylab, hline = NA){
        rng <- range(c(om, ass), na.rm = TRUE)
        plot(years, om[1, ], type = "n", ylim = c(0.9, 1.1) * rng,
             xlab = "Year", ylab = ylab)
        if(is.finite(hline)) graphics::abline(h = hline, col = "grey70", lty = 2)
        for(i in 1:nrep)
            graphics::lines(years, om[i, ],
                            col = grDevices::adjustcolor(col, if(nrep > 1) 0.4 else 1))
        if(nrep > 1)
            graphics::lines(years, apply(om, 2, stats::median), col = col, lwd = 2)
        if(!all(is.na(ass)))
            graphics::lines(years, ass, col = col.assess, lwd = 2, lty = 2)
        graphics::box(lwd = 1.5)
    }

    op <- graphics::par(mfrow = c(3, 2), mar = c(4, 4, 1, 1),
                        oma = c(0, 0, 2.5, 0))
    on.exit(graphics::par(op))

    panel(res$FM, res$obs$ffmsy * res$fmsy, "F")
    panel(res$FM / res$fmsy, res$obs$ffmsy, "F / Fmsy", hline = 1)
    panel(res$esb, rep(NA_real_, ny), "ESB")
    panel(res$esb / res$esbmsy, res$obs$bbmsy, "ESB / ESBmsy   (fitted)", hline = 1)
    panel(res$cw / res$msy, res$obs$cmsy, "Catch / MSY   (not fitted)", hline = 1)
    panel(res$rec, rep(NA_real_, ny), "Recruitment")

    graphics::mtext("operating model (solid) vs assessment (dashed)",
                    outer = TRUE, side = 3, line = 0.8, cex = 0.9)
    invisible(NULL)
}


#' @name print.iamse.cond
#'
#' @title Print a conditioned operating model
#'
#' @param x an object of class `iamse.cond`
#' @param ... ignored
#'
#' @return Invisibly `x`.
#'
#' @export
print.iamse.cond <- function(x, ...){
    ny <- ncol(x$FM)
    cat("Operating model conditioned on a spict assessment\n")
    cat("  replicates      : ", nrow(x$FM), " (all converged)\n", sep = "")
    cat("  historic years  : ", ny, " (", x$obs$years[1], "-",
        x$obs$years[ny], " of ", x$obs$nya, " assessment years)\n", sep = "")
    cat("  method / target : ", x$method, " / ", x$target, "\n", sep = "")
    cat("  F random walk sd: ", signif(x$sdF, 3), "\n", sep = "")
    fit.rmse <- sqrt(mean((log(x$obs$bbmsy) -
                           log(apply(x$esb, 2, stats::median) / x$esbmsy))^2))
    cat("  fitted   log rmse ESB/ESBmsy : ", signif(fit.rmse, 3), "\n", sep = "")
    cw.rmse <- sqrt(mean((log(x$obs$cmsy) -
                          log(apply(x$cw, 2, stats::median) / x$msy))^2))
    cat("  unfitted log rmse C/MSY      : ", signif(cw.rmse, 3), "\n", sep = "")
    cat("  Bmsy / unfished : ", signif(x$bmsy.frac[["om"]], 3),
        " (operating model) vs ", signif(x$bmsy.frac[["spict"]], 3),
        " (assessment)\n", sep = "")
    cat("  F/Fmsy last year: ", signif(stats::median(x$FM[, ny]) / x$fmsy, 3),
        " (operating model) vs ", signif(x$obs$ffmsy[ny], 3),
        " (assessment)\n", sep = "")
    invisible(x)
}
