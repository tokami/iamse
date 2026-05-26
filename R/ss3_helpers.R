

#' Resolve an SS3 executable path
#'
#' Finds an SS3 executable in a user-supplied directory, or downloads one
#' using `r4ss::get_ss3_exe()` when not available.
#'
#' @param exe_path Character path to an SS3 executable or a directory where
#'   an executable should be found/downloaded.
#'
#' @return Character path to an executable SS3 binary.
#' @keywords internal
#' @noRd
get_ss3_exe_path <- function(exe_path = NULL) {

  exe_path <- path.expand(exe_path)

  if (dir.exists(exe_path) && grep("ss3", dir(exe_path))) {

    exe_path <- file.path(exe_path, "ss3") ## TODO: only for linux?

  } else {
    dir.create(exe_path, FALSE, TRUE)
    exe_path <- r4ss::get_ss3_exe(dir = exe_path)
  }

  if (.Platform$OS.type != "windows") {
    Sys.chmod(exe_path, mode = "0755")
  }

  if (!is.null(exe_path)) return(exe_path)

  ## TODO: not implemented yet:
  opt <- getOption("iamse.ss3_exe")
  if (!is.null(opt)) return(opt)

  stop(
    "No SS3 executable supplied. Set `exe_path=` or option ",
    "`options(iamse.ss3_exe = '/path/to/ss3')`."
  )
}


#' Copy a bundled SS3 template
#'
#' Copies one of the package SS3 templates from `inst/extdata/ss3` into a
#' working directory.
#'
#' @param path Destination directory.
#' @param template Template folder name inside bundled SS3 templates.
#'
#' @return Invisibly returns `path`.
#' @keywords internal
#' @noRd
copy_ss3_template <- function(path, template = "simple_small") {

  src <- system.file("extdata", "ss3", template, package = "iamse")

  if (src == "") {
    stop("Could not find bundled SS3 template '", template, "'.")
  }

  dir.create(path, recursive = TRUE, showWarnings = FALSE)

  files <- list.files(src, full.names = TRUE, all.files = TRUE, no.. = TRUE)

  ## Do not copy old cached R objects
  files <- files[!grepl("\\.rds$", files, ignore.case = TRUE)]

  ok <- file.copy(
    from = files,
    to = path,
    recursive = TRUE,
    overwrite = TRUE
  )

  if (!all(ok)) {
    stop("Failed to copy one or more SS3 template files to: ", path)
  }

  invisible(path)
}


#' Read SS3 template input files
#'
#' Reads SS3 input files from a template directory and applies small
#' compatibility fixes for missing tag-related slots.
#'
#' @param path Directory containing SS3 input files.
#'
#' @return A list returned by `r4ss::SS_read()` with patched `dat`/`ctl` slots.
#' @keywords internal
#' @noRd
load_ss3_template_inputs <- function(path) {

  inputs <- r4ss::SS_read(
    dir = path,
    verbose = FALSE
  )

  ## Workaround for r4ss no-tag templates
  if (is.null(inputs$dat$do_tags) || length(inputs$dat$do_tags) == 0) {
    inputs$dat$do_tags <- 0
  }
  if (is.null(inputs$dat$N_tag_groups) || length(inputs$dat$N_tag_groups) == 0) {
    inputs$dat$N_tag_groups <- 0
  }
  if (is.null(inputs$dat$N_recap_events) || length(inputs$dat$N_recap_events) == 0) {
    inputs$dat$N_recap_events <- 0
  }
  if (is.null(inputs$ctl$TG_custom) || length(inputs$ctl$TG_custom) == 0) {
    inputs$ctl$TG_custom <- 0
  }

  inputs
}



#' Update SS3 inputs from iamse observations
#'
#' Fills and modifies SS3 `dat`, `ctl`, and `fore` objects using iamse-style
#' observations and simulation data.
#'
#' @param inputs List from `r4ss::SS_read()` containing `dat`, `ctl`, and `fore`.
#' @param obs Observation list (catch, index, and optional length-composition).
#' @param sim_dat iamse data list used to derive biology and fishery settings.
#' @param catch_se Catch standard error used in SS3 catch table.
#' @param cpue_se_log Optional log-scale SE for CPUE; defaults to
#'   `sim_dat$sd_i` when available.
#' @param cpue_month Month assigned to CPUE observations.
#' @param fix_biology Logical; whether to set and fix key biology parameters.
#' @param fix_selectivity Logical; whether to set and fix selectivity settings.
#' @param fix_q Logical; whether to set and fix survey catchability.
#' @param use_lencomp Logical; include length compositions if `TRUE`.
#' @param ny_lencomp Number of years to sample for length comps, or "all".
#'
#' @return Updated SS3 input list with modified `dat`, `ctl`, and `fore`.
#' @keywords internal
#' @noRd
update_ss3_from_obs <- function(inputs,
                                obs,
                                sim_dat,
                                catch_se = 0.01,
                                cpue_se_log = NULL,
                                cpue_month = 6,
                                fix_biology = TRUE,
                                fix_selectivity = TRUE,
                                fix_q = TRUE,
                                use_lencomp = FALSE,
                                ny_lencomp = 5) {

  ss_dat <- inputs$dat
  ctl <- inputs$ctl
  fore <- inputs$fore


  ## ---------------- helper ----------------
  set_par <- function(tab, pattern, init, lo = NULL, hi = NULL,
                      prior = init, phase = -99,
                      pr_sd = 99, pr_type = 0) {
    ii <- grep(pattern, rownames(tab))
    if (length(ii) == 0 || !is.finite(init)) return(tab)

    if (!is.null(lo)) tab[ii, "LO"] <- lo
    if (!is.null(hi)) tab[ii, "HI"] <- hi

    tab[ii, "INIT"] <- init
    tab[ii, "PRIOR"] <- prior
    tab[ii, "PR_SD"] <- pr_sd
    tab[ii, "PR_type"] <- pr_type
    tab[ii, "PHASE"] <- phase

    tab
  }

  ## ---------------- checks ----------------

  if (is.null(obs$timeC) || is.null(obs$obsC)) {
    stop("obs must contain 'timeC' and 'obsC'.")
  }

  if (is.null(obs$timeI) || is.null(obs$obsI) ||
        length(obs$timeI) < 1 || length(obs$obsI) < 1) {
    stop("obs must contain at least one survey index in 'timeI[[1]]' and 'obsI[[1]]'.")
  }

  years <- sort(unique(as.integer(obs$timeC)))
  if (length(years) == 0 || anyNA(years)) {
    stop("Could not determine valid years from obs$timeC.")
  }

  ## 1955 - 2024
  ss_dat$styr <- min(years)
  ss_dat$endyr <- max(years)

  ## ---------------- catch ----------------

  keep_c <- !is.na(obs$timeC) & !is.na(obs$obsC)

  ss_dat$catch <- rbind(
    data.frame(
      year = -999,
      seas = 1,
      fleet = 1,
      catch = 1e-20,
      catch_se = catch_se
    ),
    data.frame(
      year = as.integer(obs$timeC[keep_c]),
      seas = 1,
      fleet = 1,
      catch = as.numeric(obs$obsC[keep_c]),
      catch_se = catch_se
    )
  )

  ## ---------------- survey index ----------------

  ind_year <- as.integer(obs$timeI[[1]])
  ind_obs <- as.numeric(obs$obsI[[1]])

  if (all(is.na(ind_obs))) stop("all obsI are NA! check!")

  keep_i <- !is.na(ind_year) & !is.na(ind_obs) &
    ind_year >= ss_dat$styr & ind_year <= ss_dat$endyr

  if (is.null(cpue_se_log)) {
    cpue_se_log <- if (!is.null(sim_dat$sd_i)) sim_dat$sd_i else 0.2
  }

  ## gives error if specified here
  ## ss_dat$CPUEinfo <- data.frame(
  ##   fleet = 2,
  ##   units = 1,
  ##   errtype = 0,
  ##   SD_report = 1,
  ##   stringsAsFactors = FALSE
  ## )
  ## rownames(ss_dat$CPUEinfo) <- "SURVEY1"

  ss_dat$CPUE <- data.frame(
    year = ind_year[keep_i],
    month = cpue_month,
    index = 1,
    obs = ind_obs[keep_i],
    se_log = cpue_se_log
  )

  ## ---------------- optional length comps ----------------

  if (!isTRUE(use_lencomp)) {
    ss_dat$use_lencomp <- 0
    ss_dat$len_info <- NULL
    ss_dat$lencomp <- NULL
  } else {

    ss_dat$use_lencomp <- 1

    # observed length comp matrix
    obsL <- obs$obsCL

    # convert column names like "01", "02", ... to numeric
    obs_bins <- as.numeric(colnames(obsL))

    ## number of years specified
    if (ny_lencomp == "all") {
      row_ind <- 1:nrow(obsL)
    } else {
      row_ind <- sort(sample(1:nrow(obsL),
                             as.integer(ny_lencomp)))  ## TODO: not ideal if they are X different every run
    }

    obsL_prop <- obsL[row_ind,] / rowSums(obsL[row_ind,])

    # keep only bins used in SS3
    keep <- obs_bins %in% ss_dat$lbin_vector
    obsL_keep <- obsL_prop[, keep, drop = FALSE]

    # rename columns to SS3 format
    colnames(obsL_keep) <- paste0("f", obs_bins[keep])

    years_L <- as.numeric(rownames(obsL_keep))

    ## formula used by Lucas
    Nsamp <- rowSums(obsL[row_ind,]) / sum(obsL[row_ind,]) * 200

    lencomp_new <- data.frame(
      year  = years_L,
      month = 1,
      fleet = 1,
      sex   = 1,
      part  = 0,
      Nsamp = Nsamp,
      obsL_keep,
      check.names = FALSE
    )

    ## zeros for males
    male_cols <- paste0("m", ss_dat$lbin_vector)

    for (col in male_cols) {
      lencomp_new[[col]] <- 0
    }

    lencomp_new <- lencomp_new[, names(ss_dat$lencomp)]
    ss_dat$lencomp <- lencomp_new

  }

  ## ---------------- no age data / no tags ----------------

  ## Works:
  ## ss_dat$N_agebins <- 0
  ## ss_dat$agebin_vector <- NULL
  ## ss_dat$N_ageerror_definitions <- 0
  ## ss_dat$ageerror <- NULL
  ## ss_dat$age_info <- NULL
  ## ss_dat$agecomp <- NULL
  ## ss_dat$use_MeanSize_at_Age_obs <- 0
  ## ss_dat$MeanSize_at_Age_obs <- NULL


  ## trying:
  ss_dat$N_agebins <- 32 ## TODO make sim_dat dependent
  ss_dat$agebin_vector <- 0:31
  ss_dat$N_ageerror_definitions <- 1
  tmp <- as.data.frame(rbind(rep(-1, ss_dat$N_agebins + 1),
                             rep(0.001, ss_dat$N_agebins + 1)))
  colnames(tmp) <- paste0("age",c(ss_dat$agebin_vector, max(ss_dat$agebin_vector) + 1))
  ss_dat$ageerror <- tmp
  tmp <- as.data.frame(rbind(c(-1,0.001,0,0,0,0,0.001)))
  colnames(tmp) <- c("mintailcomp","addtocomp",
                     "combine_M_F","CompressBins",
                     "CompError","ParmSelect","minsamplesize")
  rownames(tmp) <- "Fishery"
  ss_dat$age_info <- tmp
  ss_dat$agecomp <- NULL
  ss_dat$use_MeanSize_at_Age_obs <- 0
  ss_dat$MeanSize_at_Age_obs <- NULL

  ss_dat$Lbin_method <- 3


  ## tags
  ss_dat$do_tags <- 0
  ss_dat$N_tag_groups <- 0
  ss_dat$N_recap_events <- 0

  ctl$TG_custom <- 0

  ## ---------------- length bins ----------------

  ## Species-appropriate population length bins
  lmax_pop <- ceiling(sim_dat$Linf * 1.4)

  ss_dat$lbin_method <- 2
  ss_dat$binwidth <- sim_dat$binwidth
  ss_dat$lbin_vector_pop <- seq(
    from = 0,
    to = lmax_pop,
    by = sim_dat$binwidth
  )
  ss_dat$N_lbinspop <- length(ss_dat$lbin_vector_pop)
  ss_dat$minimum_size <- 0
  ss_dat$maximum_size <- lmax_pop


  ## ---------------- biology ----------------

  ## general settings
  ctl$time_vary_auto_generation[] <- 1
  ctl$Exp_Decay <- 0
  ctl$Growth_Age_for_L1 <- 0.09
  ctl$Growth_Age_for_L2 <- 999


  ss_dat$Nages <- sim_dat$amax

  if (isTRUE(fix_biology) && !is.null(ctl$MG_parms)) {

    ## do not estimate male pars
    ii_male <- grep("_Mal_", rownames(ctl$MG_parms))

    if (length(ii_male) > 0) {
      ctl$MG_parms[ii_male, "PHASE"] <- -99
    }

    M_value <- mean(sum(sim_dat$M[1,]) * sim_dat$Msel[[1]])

    ## per-age (annual) mean length from OM
    la_base <- if (is.list(sim_dat$LA)) sim_dat$LA[[1]] else if (is.matrix(sim_dat$LA)) sim_dat$LA[1, ] else as.numeric(sim_dat$LA)
    ns_val  <- if (!is.null(sim_dat$ns) && sim_dat$ns >= 1L) sim_dat$ns else 1L
    la_annual <- la_base[seq(1L, length(la_base), by = ns_val)]

    ## L_at_Amin from youngest age mean length (was previously hardcoded 0)
    L_at_Amin <- if (length(la_annual) > 0 && is.finite(la_annual[1]) && la_annual[1] > 0) {
      la_annual[1]
    } else {
      min(la_base[is.finite(la_base) & la_base > 0], na.rm = TRUE)
    }
    L_at_Amax <- sim_dat$Linf

    ## maturity-at-length from OM's dat$mat (can exceed 1 for old ages)
    if (!is.null(sim_dat$mat) && length(sim_dat$mat) > 0) {
      mat_base   <- as.numeric(sim_dat$mat)
      mat_annual <- mat_base[seq(1L, length(mat_base), by = ns_val)]
      max_mat    <- max(mat_annual[is.finite(mat_annual)], na.rm = TRUE)
      if (!is.finite(max_mat) || max_mat <= 0) max_mat <- 1

      ## normalize to [0,1] and find Mat50 / L95 via interpolation on length axis
      mat_norm <- mat_annual / max_mat
      ok_mat   <- is.finite(mat_norm) & is.finite(la_annual) & mat_norm > 0.05 & mat_norm < 0.99
      if (sum(ok_mat) >= 2L) {
        oo     <- order(la_annual[ok_mat])
        la_ok  <- la_annual[ok_mat][oo]
        mn_ok  <- mat_norm[ok_mat][oo]
        Mat50_om <- tryCatch(stats::approx(mn_ok, la_ok, xout = 0.5,  rule = 2)$y, error = function(e) NA_real_)
        Lm95_om  <- tryCatch(stats::approx(mn_ok, la_ok, xout = 0.95, rule = 2)$y, error = function(e) NA_real_)
      } else {
        Mat50_om <- NA_real_
        Lm95_om  <- NA_real_
      }
      if (is.finite(Mat50_om) && is.finite(Lm95_om) && Lm95_om > Mat50_om) {
        mat_slope <- -log(19) / (Lm95_om - Mat50_om)
        Lm50_used <- Mat50_om
      } else {
        mat_slope <- -log(19) / (sim_dat$Lm95 - sim_dat$Lm50)
        Lm50_used <- sim_dat$Lm50
      }
      Eggs_alpha_new <- sim_dat$lwa * max_mat
    } else {
      max_mat        <- 1
      mat_slope      <- -log(19) / (sim_dat$Lm95 - sim_dat$Lm50)
      Lm50_used      <- sim_dat$Lm50
      Eggs_alpha_new <- sim_dat$lwa
    }

    ## CV for age-length relationship — honour dat$CVlen if set
    cv_len <- if (!is.null(sim_dat$CVlen) && is.finite(sim_dat$CVlen)) sim_dat$CVlen else 0.1

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "NatM_p_1_Fem", M_value,
      lo = 0.001, hi = max(2, M_value * 3),
      phase = -1, pr_type = 3, pr_sd = 10
    )

    min_pop_len <- min(ss_dat$lbin_vector_pop)
    max_pop_len <- max(ss_dat$lbin_vector_pop)

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "L_at_Amin_Fem",
      L_at_Amin,
      lo = min_pop_len,
      hi = max_pop_len,
      phase = -3, pr_sd = 10
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "L_at_Amax_Fem",
      sim_dat$Linf,
      lo = min_pop_len,
      hi = max_pop_len,
      phase = -2, pr_sd = 10
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "VonBert_K_Fem", sim_dat$K,
      lo = 0.001, hi = max(2, sim_dat$K * 5),
      phase = -3, pr_sd = 0.05
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "CV_young_Fem", cv_len,
      lo = 0.001, hi = 5,
      phase = -4, pr_sd = 0.5
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "CV_old_Fem", cv_len,
      lo = 0.001, hi = 5,
      phase = -4, pr_sd = 0.5
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Wtlen_1_Fem", sim_dat$lwa,
      lo = 0,
      hi = 3,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Wtlen_2_Fem", sim_dat$lwb,
      lo = 2,
      hi = 4,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Mat50", Lm50_used,
      lo = 0.0001,
      hi = 1000,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Mat_slope", mat_slope,
      lo = -2,
      hi = 4,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Eggs_alpha", Eggs_alpha_new,
      lo = 0,
      hi = max(3, Eggs_alpha_new * 5),
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Eggs_beta", sim_dat$lwb,
      lo = 0,
      hi = 10,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "NatM_p_1_Mal", M_value,
      lo = 0.001, hi = max(2, M_value * 3),
      phase = -1, pr_type = 3, pr_sd = 10
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "L_at_Amin_Mal",
      L_at_Amin,
      lo = min_pop_len,
      hi = max_pop_len,
      phase = -1, pr_sd = 10
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "L_at_Amax_Mal",
      sim_dat$Linf,
      lo = min_pop_len,
      hi = max_pop_len,
      phase = -1, pr_sd = 10
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "VonBert_K_Mal", sim_dat$K,
      lo = 0.001, hi = max(2, sim_dat$K * 5),
      phase = -3, pr_sd = 0.05
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "CV_young_Mal", cv_len,
      lo = 0.001, hi = 5,
      phase = -4, pr_sd = 0.5
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "CV_old_Mal", cv_len,
      lo = 0.001, hi = 5,
      phase = -4, pr_sd = 0.5
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Wtlen_1_Mal", sim_dat$lwa,
      lo = 0,
      hi = 3,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Wtlen_2_Mal", sim_dat$lwb,
      lo = 2,
      hi = 4,
      phase = -99
    )

  }

  ctl$First_Mature_Age <- 0
  ctl$fecundity_option <- 2

  ## ---------------- stock-recruitment ----------------

  ctl$SR_function <- 3
  ctl$MainRdevYrFirst <- ss_dat$styr
  ctl$MainRdevYrLast <- ss_dat$endyr
  ctl$F_ballpark_year <- ss_dat$styr


  if (!is.null(ctl$SR_parms)) {
    rn <- rownames(ctl$SR_parms)

    ii <- grep("R0|LnR0|log_R0|SR_LN\\(R0\\)", rn)
    if (length(ii) > 0 && !is.null(sim_dat$R0)) {
      ## iamse dat$R0 is female R0 (SSB uses /2, recruitment multiplies by 2).
      ## SS3 SR_LN(R0) is total R0 for both sexes → multiply by 2.
      lnR0_total <- log(sim_dat$R0 * 2)
      ctl$SR_parms[ii[1], "LO"] <- 0.0001
      ctl$SR_parms[ii[1], "HI"] <- 20
      ctl$SR_parms[ii[1], "INIT"] <- lnR0_total
      ctl$SR_parms[ii[1], "PRIOR"] <- lnR0_total
      ctl$SR_parms[ii[1], "PR_SD"] <- 99
      ctl$SR_parms[ii[1], "PR_type"] <- 0
      ctl$SR_parms[ii[1], "PHASE"] <- 1
    }

    ii <- grep("steep|Steep", rn)
    if (length(ii) > 0 && !is.null(sim_dat$h)) {
      ctl$SR_parms[ii[1], "INIT"] <- sim_dat$h
      ctl$SR_parms[ii[1], "PRIOR"] <- sim_dat$h
      ctl$SR_parms[ii[1], "PR_SD"] <- 0.24
      ctl$SR_parms[ii[1], "PR_type"] <- 3
      ctl$SR_parms[ii[1], "PHASE"] <- -1
    }

    ii <- grep("sigmaR", rn)
    if (length(ii) > 0 && !is.null(sim_dat$sd_r)) {
      ctl$SR_parms[ii[1], "INIT"] <- sim_dat$sd_r
      ctl$SR_parms[ii[1], "PRIOR"] <- sim_dat$sd_r
      ctl$SR_parms[ii[1], "PR_SD"] <- 99
      ctl$SR_parms[ii[1], "PR_type"] <- 0
      ctl$SR_parms[ii[1], "PHASE"] <- -6
    }

    ii <- grep("regime", rn)
    if (length(ii) > 0) {
      ctl$SR_parms[ii[1], "LO"] <- -5
      ctl$SR_parms[ii[1], "HI"] <- 5
      ctl$SR_parms[ii[1], "INIT"] <- 0
      ctl$SR_parms[ii[1], "PRIOR"] <- 0
      ctl$SR_parms[ii[1], "PR_SD"] <- 99
      ctl$SR_parms[ii[1], "PR_type"] <- 0
      ctl$SR_parms[ii[1], "PHASE"] <- -99
    }

    ii <- grep("autocorr", rn)
    if (length(ii) > 0) {
      ctl$SR_parms[ii[1], "LO"] <- 0
      ctl$SR_parms[ii[1], "HI"] <- 2
      ctl$SR_parms[ii[1], "INIT"] <- 0
      ctl$SR_parms[ii[1], "PRIOR"] <- 1
      ctl$SR_parms[ii[1], "PR_SD"] <- 99
      ctl$SR_parms[ii[1], "PR_type"] <- 0
      ctl$SR_parms[ii[1], "PHASE"] <- -99
    }
  }

  ctl$do_recdev <- 3
  ctl$recdev_phase <- 1
  ctl$Fcast_recr_phase <- -4

  ctl$last_early_yr_nobias_adj <- min(years)
  ctl$first_yr_fullbias_adj <- min(years)
  ctl$last_yr_fullbias_adj <- max(years)
  ctl$first_recent_yr_nobias_adj <- max(years)

  ctl$F_ballpark <- 0.03
  ctl$F_ballpark_year <- -1
  ctl$maxF <- 4

  ## Init F
  ctl$init_F["InitF_seas_1_flt_1Fishery", ] <- c(
    LO = 0,
    HI = 10,
    INIT = 1e-20,
    PRIOR = 1,
    PR_SD = 999,
    PR_type = 0,
    PHASE = -1
  )

  ctl$maxlambdaphase <- 15

  ## ctl$lambdas ## okay

  ctl$more_stddev_reporting <- 0



  ## ---------------- selectivity ----------------

  if (isTRUE(fix_selectivity) && !is.null(ctl$size_selex_parms)) {

    ## logistic pattern
    ctl$size_selex_types[1] <- 1

    L50 <- sim_dat$Ls50
    L95 <- sim_dat$Ls95

    sel_width <- log(19) / (L95 - L50)
    ## sel_width <- (L95 - L50)

    extra <- grep("^SizeSel_P_[3-6]_Fishery\\(1\\)$", rownames(ctl$size_selex_parms))
    if (length(extra) > 0) {
      ctl$size_selex_parms <- ctl$size_selex_parms[-extra,]
    }

    ctl$size_selex_parms <- set_par(
      ctl$size_selex_parms, "SizeSel_P_1_Fishery\\(1\\)",
      sim_dat$Ls50,
      lo = 0.001,
      hi = max(100, sim_dat$Ls50 * 3),
      phase = -99
    )

    ctl$size_selex_parms <- set_par(
      ctl$size_selex_parms, "SizeSel_P_2_Fishery\\(1\\)",
      sel_width,
      lo = 0,
      hi = 100,
      phase = -99
    )

  }


  ## Remove stale age selectivity parameters, if present
  if (!is.null(ctl$age_selex_parms)) {
    ctl$age_selex_parms <- NULL
  }

  ## ---------------- survey Q ----------------

  ## ctl$Q_options ## okay

  if (!is.null(ctl$Q_parms)) {

    iq <- grep("LnQ", rownames(ctl$Q_parms))

    ctl$Q_parms[iq, "LO"] <- -15
    ctl$Q_parms[iq, "HI"] <- 15
    ctl$Q_parms[iq, "INIT"] <- log(sim_dat$q)
    ctl$Q_parms[iq, "PRIOR"] <- log(sim_dat$q)
    ctl$Q_parms[iq, "PR_SD"] <- 1
    ctl$Q_parms[iq, "PR_type"] <- 0
    ctl$Q_parms[iq, "PHASE"] <- -1

  }


  ## ---------------- forecast ----------------

  fore$FirstYear_for_caps_and_allocations <- ss_dat$endyr + 1
  fore$Ydecl <- ss_dat$endyr + 1
  fore$Yinit <- ss_dat$endyr + 1
  ## Disable SS3-internal control rule — the R HCR code handles any hockey
  ## stick adjustment. Without this, the template's default ControlRuleMethod=1
  ## fires when depletion < 0.4, reducing F to near zero for depleted stocks.
  fore$ControlRuleMethod <- 0
  fore$Flimitfraction <- 1.0


  ## ---------------- return ----------------

  ## TODO causes error in fitting (maybe fix later)
  ## ss_dat$Nsexes <- 1

  ## lengthbin 4-38

  inputs$dat <- ss_dat
  inputs$ctl <- ctl
  inputs$fore <- fore

  inputs
}


#' Get length-composition columns from an SS3 template
#'
#' @param template Template name bundled with the package.
#'
#' @return Character vector of length-composition column names.
#' @keywords internal
#' @noRd
get_ss3_lencomp_columns <- function(template = "simple_small") {
  tmp_dir <- tempfile("ss3_template_")
  dir.create(tmp_dir)

  copy_ss3_template(path = tmp_dir, template = template)
  inputs <- r4ss::SS_read(dir = tmp_dir, verbose = FALSE)

  meta_cols <- intersect(
    c("year", "month", "fleet", "sex", "part", "Nsamp"),
    names(inputs$dat$lencomp)
  )
  comp_cols <- setdiff(names(inputs$dat$lencomp), meta_cols)

  comp_cols
}


#' Prepare an SS3 run directory from iamse data
#'
#' Copies a template, updates SS3 input files from observations, and writes the
#' resulting files to disk.
#'
#' @inheritParams update_ss3_from_obs
#' @param path Destination directory for SS3 input files.
#' @param dat iamse data list used as simulation input (`sim_dat`).
#'
#' @return Invisibly returns `path`.
#' @keywords internal
#' @noRd
prepare_ss3_run <- function(path,
                            obs,
                            dat,
                            template = "simple_small",
                            catch_se = 0.01,
                            cpue_se_log = NULL,
                            cpue_month = 6,
                            fix_biology = TRUE,
                            fix_selectivity = TRUE,
                            fix_q = TRUE,
                            use_lencomp = FALSE,
                            ny_lencomp = 5) {

  copy_ss3_template(path = path, template = template)

  ## Read the copied SS3 input files
  inputs <- load_ss3_template_inputs(path = path)

  ## Modify files
  inputs <- update_ss3_from_obs(inputs, obs, dat,
                                catch_se = catch_se,
                                cpue_se_log = cpue_se_log,
                                cpue_month = cpue_month,
                                fix_biology = fix_biology,
                                fix_selectivity = fix_selectivity,
                                fix_q = fix_q,
                                use_lencomp = use_lencomp,
                                ny_lencomp = ny_lencomp
                                )

  ## Keep these explicit after modification too
  if (is.null(inputs$dat$do_tags) || length(inputs$dat$do_tags) == 0) {
    inputs$dat$do_tags <- 0
  }
  if (is.null(inputs$dat$N_tag_groups) || length(inputs$dat$N_tag_groups) == 0) {
    inputs$dat$N_tag_groups <- 0
  }
  if (is.null(inputs$dat$N_recap_events) || length(inputs$dat$N_recap_events) == 0) {
    inputs$dat$N_recap_events <- 0
  }
  if (is.null(inputs$ctl$TG_custom) || length(inputs$ctl$TG_custom) == 0) {
    inputs$ctl$TG_custom <- 0
  }

  r4ss::SS_write(inputs, dir = path, overwrite = TRUE)

  invisible(path)
}


#' Summarize basic SS3 fit diagnostics
#'
#' @param replist SS3 output list from `r4ss::SS_output()`.
#'
#' @return Invisibly returns a list of selected convergence/diagnostic fields.
#' @keywords internal
#' @noRd
check_ss3_fit <- function(replist) {

  out <- list(
    SS_version = replist$SS_version,
    Final_phase = replist$Final_phase,
    N_iterations = replist$N_iterations,
    maximum_gradient_component = replist$maximum_gradient_component,
    Nwarnings = replist$Nwarnings,
    warnings = replist$warnings,
    log_det_hessian = replist$log_det_hessian,
    parameters_with_highest_gradients =
      replist$parameters_with_highest_gradients
  )

  print(out)

  invisible(out)
}


#' Read an SS3 assessment directory
#'
#' Convenience wrapper around `r4ss::SS_output()` with basic path checks.
#'
#' @param ss3_dir Directory containing SS3 output files.
#' @param forecast Passed to `r4ss::SS_output()`.
#' @param verbose Passed to `r4ss::SS_output()`.
#' @param printstats Passed to `r4ss::SS_output()`.
#'
#' @return SS3 output list as returned by `r4ss::SS_output()`.
#' @keywords internal
#' @noRd
read_ss3_assessment <- function(ss3_dir,
                                forecast = TRUE,
                                verbose = FALSE,
                                printstats = FALSE) {

  ss3_dir <- path.expand(ss3_dir)

  if (!dir.exists(ss3_dir)) {
    stop("SS3 directory does not exist: ", ss3_dir)
  }

  if (!file.exists(file.path(ss3_dir, "Report.sso"))) {
    stop("Could not find Report.sso in: ", ss3_dir)
  }

  r4ss::SS_output(
    dir = ss3_dir,
    forecast = forecast,
    verbose = verbose,
    printstats = printstats
  )
}



#' Condition iamse on an SS3 assessment
#'
#' Reads a completed SS3 assessment and conditions iamse `dat`/`set` objects to
#' match SS3 trajectories and key quantities (growth, mortality, selectivity,
#' recruitment variation, and reference metrics).
#'
#' @param ss3_dir Directory containing SS3 outputs (must include `Report.sso`).
#' @param ns Optional target number of seasons for returned iamse processes.
#' @param dat Optional iamse data list. Defaults to `iamse::check.dat()`.
#' @param set Optional iamse settings list. Defaults to `iamse::check.set()` with
#'   zeroed process/replicate noise and reduced reference-point settings.
#' @param show_plots Logical; draw trajectory comparison plots.
#' @param verbose Logical; print informative messages.
#'
#' @return Invisibly returns a list with conditioned objects and diagnostics:
#'   `dat_cond`, `set_cond`, `pop`, `replist`, `extracted`, `report`,
#'   `rec_map`, and `trajectories`.
#'
#' @examples
#' \dontrun{
#' out <- condition_on_ss3("path/to/ss3_run")
#' str(out$dat_cond)
#' }
#' @export
condition_on_ss3 <- function(ss3_dir,
                             ns = NULL,
                             dat = NULL,
                             set = NULL,
                             use_init_naa = TRUE,
                             show_plots = TRUE,
                             verbose = TRUE) {

  if (is.null(dat)) dat <- iamse::check.dat()
  if (is.null(set)) {
    set <- iamse::check.set()

    set$noise$time$R <- c(0, 0, 0)
    set$noise$time$F <- c(0, 0, 0)
    set$noise$time$C <- c(0, 0, 0)
    set$noise$time$I <- c(0, 0, 0)
    set$noise$time$E <- c(0, 0, 0)
    set$noise$rep$M <- c(0, 0, 0)
    set$noise$rep$H <- c(0, 0, 0)
    set$noise$rep$R0 <- c(0, 0, 0)
    set$strict_ss3 <- FALSE
    set$catch_eq <- "baranov"
    set$recordLast <- 1

    ## some speed-up
    set$refN <- 100
    set$refYears <- 50
    set$refYearsMSY <- 10
    set$burnin <- 50
  }

  out <- .condi_internal(
    dat = dat,
    set = set,
    ss3_dir = ss3_dir,
    condition_fields = c("growth", "mat_weight", "selectivity",
                         "M", "sr", "q", "timing", "ref_ss3", "ref_iamse"),
    interp_method = "functional_first",
    mode = "trajectory_conditioning",
    strict = FALSE,
    forecast = TRUE,
    verbose = verbose
  )

  dat_cond <- out$dat
  years_hist <- sort(unique(subset(out$replist$timeseries, Era == "TIME" & Yr <= out$replist$endyr)$Yr))[seq_len(dat_cond$ny)]

  if (!is.null(ns)) {
    if (!is.numeric(ns) || length(ns) != 1 || ns < 1) stop("`ns` must be a positive integer or NULL")
    dat_cond <- .expand_processes_to_ns(dat_cond, ns_target = as.integer(ns))
    dat_cond$nyI <- dat_cond$ny
    dat_cond$surveyTimes <- rep(0.5, 1)
    dat_cond$nyC <- dat_cond$ny
    ## dat_cond$plba <- NULL  ## TODO: why?
  }

  if (is.list(dat_cond$Msel) && length(dat_cond$Msel) >= 1) {
    for (i in seq_along(dat_cond$Msel)) {
      mi <- as.numeric(dat_cond$Msel[[i]])
      if (length(mi) == dat_cond$asmax && all(is.na(mi))) {
        dat_cond$Msel[[i]] <- rep(1, dat_cond$asmax)
      }
    }
  }

  age_cols <- as.character(0:dat_cond$amax)
  naa_tab <- out$replist$natage_annual_2_with_fishery
  init_naa <- NULL
  if (!is.null(naa_tab) && nrow(naa_tab) > 0) {
    naa_row <- subset(naa_tab, Yr == min(naa_tab$Yr, na.rm = TRUE) &
                                 Sex == min(naa_tab$Sex, na.rm = TRUE))
    if (nrow(naa_row) > 0) {
      init_naa <- as.numeric(naa_row[1, age_cols])
    }
  }

  wt <- out$replist$wtatage
  if (!is.null(wt) && nrow(wt) > 0) {
    wf <- subset(wt, fleet == 1)
    if (nrow(wf) > 0) {
      sx <- min(wf$sex, na.rm = TRUE)
      wf <- subset(wf, sex == sx)
      w_mat <- matrix(NA_real_, nrow = dat_cond$asmax, ncol = dat_cond$ny)
      for (iy in seq_len(dat_cond$ny)) {
        y <- years_hist[iy]
        row_y <- subset(wf, year == y)
        if (nrow(row_y) < 1) row_y <- wf[which.min(abs(wf$year - y)), , drop = FALSE]
        wa <- as.numeric(row_y[1, age_cols])
        w_mat[, iy] <- rep(wa, each = dat_cond$ns)
      }
      dat_cond$weight <- w_mat
      dat_cond$weightF <- w_mat

      ## Re-derive maturity from fleet=-2 / fleet=1 wtatage ratio so that
      ## dat$weight × dat$mat = wtatage_fleet=-2 (SS3 spawning weight).
      ## This corrects endgrowth-derived or prior maturity that may not equal
      ## the fleet=-2/fleet=1 ratio (e.g. when mature fish weigh more than
      ## the fleet-1 weight at the same age, common in lobster models).
      wm <- subset(wt, fleet == -2)
      if (nrow(wm) > 0) {
        sx_m <- min(wm$sex, na.rm = TRUE)
        wm <- subset(wm, sex == sx_m)
        years_wm <- sort(unique(wm$year))
        na_ages <- length(age_cols)
        ratios <- matrix(NA_real_, nrow = length(years_wm), ncol = na_ages)
        for (iy in seq_along(years_wm)) {
          y <- years_wm[iy]
          row_m <- subset(wm, year == y)
          if (nrow(row_m) < 1) row_m <- wm[which.min(abs(wm$year - y)), , drop = FALSE]
          row_f2 <- subset(wf, year == y)
          if (nrow(row_f2) < 1) row_f2 <- wf[which.min(abs(wf$year - y)), , drop = FALSE]
          vm <- as.numeric(row_m[1, age_cols])
          vf <- as.numeric(row_f2[1, age_cols])
          ok <- is.finite(vm) & is.finite(vf) & vf > 0
          r <- rep(0, na_ages)
          r[ok] <- vm[ok] / vf[ok]
          r[r < 0] <- 0
          ## ratio > 1 is valid: mature fish (e.g. egg-bearing lobster females)
          ## weigh more than the fleet-1 weight at the same age, so
          ## dat$weight × dat$mat correctly reaches wtatage_fleet=-2
          ratios[iy, ] <- r
        }
        mat_age <- colMeans(ratios, na.rm = TRUE)
        dat_cond$mat <- rep(mat_age, each = dat_cond$ns)
      }
    }
  }

  rec_map <- .set_iamse_recruitment_errors_from_ss3(
    set = set,
    replist = out$replist,
    years = years_hist,
    method = "main_recdev",
    use_bias_correction = FALSE,
    verbose = FALSE
  )

  set_cond <- rec_map$set
  set_cond$strict_ss3 <- FALSE
  set_cond$catch_eq <- "baranov"
  if (use_init_naa && !is.null(init_naa)) set_cond$init_naa <- init_naa * 2

  ## for females only
  dat_cond$R0 <- dat_cond$R0 / 2

  dat_cond <- est.ref.levels.stochastic(dat_cond, set_cond)
  ## Blim for now:
  dat_cond$refdist$Blim <- 0.2 * dat_cond$refdist$B0
  dat_cond$ref$Blim <- 0.2 * dat_cond$ref$B0

  if (!is.null(dat_cond$sd_r)) set_cond$noise$time$R[1] <- dat_cond$sd_r
  if (!is.null(dat_cond$sd_f)) set_cond$noise$time$F[1] <- dat_cond$sd_f

  pop <- initpop(dat_cond, set_cond)
  tr <- .calc_traj_metrics(pop, out$replist, out$extracted,
                           ss3_ssb_source = "SpawnBio")


  if (isTRUE(show_plots)) .plot_traj(tr)

  dat_cond$conditioning_ss3_raw$trajectories <- tr

  invisible(list(
    dat = dat_cond,
    set = set_cond,
    pop = pop,
    replist = out$replist,
    extracted = out$extracted,
    report = out$report,
    rec_map = rec_map,
    trajectories = tr
  ))
}


#' Internal SS3-to-iamse conditioning engine
#'
#' Applies field-level conditioning of iamse `dat` from SS3 outputs in either
#' parameter- or trajectory-conditioning modes.
#'
#' @param dat iamse data list.
#' @param set iamse settings list.
#' @param ss3_dir Directory containing SS3 outputs.
#' @param condition_fields Character vector specifying which components to
#'   condition.
#' @param interp_method Interpolation method for remapping age/length processes.
#' @param mode Conditioning mode.
#' @param strict Logical; stop on missing conditioning targets when `TRUE`.
#' @param forecast Logical; passed to SS3 reader.
#' @param verbose Logical; print informative messages.
#'
#' @return A list with elements `dat`, `report`, `extracted`, and `replist`.
#' @keywords internal
#' @noRd
.condi_internal <- function(dat,
                            set,
                            ss3_dir,
                            condition_fields = c("growth", "mat_weight",
                                                 "selectivity",
                                                 "M", "sr", "q",
                                                 "timing", "ref_ss3",
                                                 "ref_iamse"),
                            interp_method = c("functional_first",
                                              "mono_spline",
                                              "linear"),
                            mode = c("parameter_conditioning",
                                     "trajectory_conditioning"),
                            strict = FALSE,
                            forecast = TRUE,
                            verbose = TRUE) {

  replist <- read_ss3_assessment(
    ss3_dir = ss3_dir,
    forecast = forecast,
    verbose = FALSE,
    printstats = FALSE
  )

  interp_method <- match.arg(interp_method)
  mode <- match.arg(mode)
  tar <- .extract_ss3_targets(replist, verbose = verbose)

  report <- data.frame(
    field = character(),
    old = numeric(),
    new = numeric(),
    method = character(),
    stringsAsFactors = FALSE
  )

  add_report <- function(field, old, new, method) {
    report <<- rbind(
      report,
      data.frame(field = field, old = as.numeric(old), new = as.numeric(new), method = method,
                 stringsAsFactors = FALSE)
    )
  }

  first_num <- function(x) {
    if (is.null(x) || length(x) < 1) return(NA_real_)
    if (is.list(x)) return(first_num(x[[1]]))
    as.numeric(x)[1]
  }

  recycle_like <- function(template, values, ny = NULL, asmax = NULL) {
    if (is.list(template)) {
      out <- template
      for (i in seq_along(out)) {
        n <- length(out[[i]])
        out[[i]] <- values[seq_len(min(length(values), n))]
        if (length(out[[i]]) < n) out[[i]] <- c(out[[i]], rep(tail(out[[i]], 1), n - length(out[[i]])))
      }
      return(out)
    }

    if (is.matrix(template)) {
      return(matrix(values[seq_len(asmax)], nrow = ny, ncol = asmax, byrow = TRUE))
    }

    values[seq_len(length(template))]
  }

  interp_curve <- function(x_old, y_old, x_new, method = "mono_spline") {
    keep <- is.finite(x_old) & is.finite(y_old)
    x_old <- as.numeric(x_old[keep])
    y_old <- as.numeric(y_old[keep])
    if (length(x_old) < 2) return(rep(first_num(y_old), length(x_new)))

    oo <- order(x_old)
    x_old <- x_old[oo]
    y_old <- y_old[oo]

    if (method == "linear") {
      return(stats::approx(x = x_old, y = y_old, xout = x_new, rule = 2)$y)
    }

    fun <- try(stats::splinefun(x_old, y_old, method = "monoH.FC"), silent = TRUE)
    if (inherits(fun, "try-error")) {
      return(stats::approx(x = x_old, y = y_old, xout = x_new, rule = 2)$y)
    }
    as.numeric(fun(x_new))
  }

  if (mode == "trajectory_conditioning") {
    ts <- subset(replist$timeseries, Era == "TIME")
    ts <- subset(ts, Yr <= replist$endyr)
    y_start <- min(ts$Yr, na.rm = TRUE)
    y_end <- max(ts$Yr, na.rm = TRUE)
    ny_new <- as.integer(y_end - y_start + 1)
    ns_new <- if (!is.null(replist$N_seasons) && is.finite(replist$N_seasons)) {
      as.integer(replist$N_seasons)
    } else {
      max(1L, length(unique(ts$Seas)))
    }

    old_ny <- dat$ny
    old_ns <- dat$ns
    old_amax <- dat$amax
    old_ages <- as.numeric(dat$ages)
    old_as2a <- dat$as2a
    amax_ss3 <- .infer_ss3_amax(replist)
    dat <- .rebuild_age_time_structure(dat, ny_new = ny_new, ns_new = ns_new, amax_new = amax_ss3)
    add_report("ny", old_ny, dat$ny, "trajectory_from_ss3")
    add_report("ns", old_ns, dat$ns, "trajectory_from_ss3")
    add_report("amax", old_amax, dat$amax, ifelse(is.na(amax_ss3), "trajectory_from_ss3_keep", "trajectory_from_ss3"))

    if (!is.null(dat$nyC) && is.finite(dat$nyC)) dat$nyC <- min(dat$nyC, dat$ny)
    if (!is.null(dat$nyE) && is.finite(dat$nyE)) dat$nyE <- min(dat$nyE, dat$ny)
    if (!is.null(dat$nyI) && length(dat$nyI) > 0) dat$nyI <- pmin(dat$nyI, dat$ny)

    old_M <- dat$M
    old_FM <- dat$FM
    dat$M <- matrix(first_num(old_M), nrow = dat$ny, ncol = dat$ns)
    dat$FM <- matrix(first_num(old_FM), nrow = dat$ny, ncol = dat$ns)

    if (!is.null(dat$LA)) {
      la_old <- if (is.list(dat$LA)) dat$LA[[1]] else as.numeric(dat$LA)
      la_new <- interp_curve(old_ages, la_old, dat$ages,
                             method = ifelse(interp_method == "linear", "linear", "mono_spline"))
      dat$LA <- la_new
    }
    if (!is.null(dat$weight)) {
      w_old <- if (is.list(dat$weight)) dat$weight[[1]] else as.numeric(dat$weight)
      w_new <- interp_curve(old_ages, w_old, dat$ages,
                            method = ifelse(interp_method == "linear", "linear", "mono_spline"))
      dat$weight <- w_new
    }
    if (!is.null(dat$weightF)) dat$weightF <- dat$weight
    if (!is.null(dat$mat)) {
      m_old <- if (is.list(dat$mat)) dat$mat[[1]] else as.numeric(dat$mat)
      m_new <- interp_curve(old_ages, m_old, dat$ages,
                            method = ifelse(interp_method == "linear", "linear", "mono_spline"))
      m_new[m_new < 0] <- 0
      m_new[m_new > 1] <- 1
      dat$mat <- m_new
    }
    if (!is.null(dat$Msel)) {
      ms_old <- if (is.list(dat$Msel)) dat$Msel[[1]] else as.numeric(dat$Msel)
      ms_new <- interp_curve(old_ages, ms_old, dat$ages,
                             method = ifelse(interp_method == "linear", "linear", "mono_spline"))
      ms_new[ms_new < 0] <- 0
      dat$Msel <- list(ms_new)
    }
    if (!is.null(dat$sel)) {
      s_old <- if (is.list(dat$sel)) dat$sel[[1]] else as.numeric(dat$sel)
      s_new <- interp_curve(old_ages, s_old, dat$ages,
                            method = ifelse(interp_method == "linear", "linear", "mono_spline"))
      s_new[s_new < 0] <- 0
      s_new[s_new > 1] <- 1
      dat$sel <- list(s_new)
    }

    if (!is.null(dat$selI) && length(dat$selI) > 0) {
      dat$selI <- lapply(dat$selI, function(v) {
        v <- as.numeric(v)
        age_vals <- as.numeric(stats::aggregate(v, by = list(old_as2a), FUN = mean)[, 2])
        rep(age_vals, each = dat$ns)
      })
    }

    if ("mat_weight" %in% condition_fields) {
      age_classes <- 0:dat$amax
      age_cols <- as.character(age_classes)
      years_unique <- sort(unique(ts$Yr))

      ## 1) Prefer endgrowth-derived annual vectors (as in MSEtool)
      used_mat_weight <- FALSE
      if (!is.null(replist$endgrowth) && nrow(replist$endgrowth) > 0) {
        eg <- replist$endgrowth
        sex_col <- if ("Sex" %in% names(eg)) "Sex" else NULL
        age_col <- if ("Age_Beg" %in% names(eg)) "Age_Beg" else if ("Age" %in% names(eg)) "Age" else NULL
        wt_col <- if ("Wt_Beg" %in% names(eg)) "Wt_Beg" else NULL
        mat_age_col <- if ("Age_Mat" %in% names(eg)) "Age_Mat" else NULL
        mat_len_col <- if ("Len_Mat" %in% names(eg)) "Len_Mat" else NULL

        if (!is.null(sex_col) && !is.null(age_col) && !is.null(wt_col)) {
          sx <- min(eg[[sex_col]], na.rm = TRUE)
          egf <- eg[eg[[sex_col]] == sx, , drop = FALSE]
          egf <- egf[order(egf[[age_col]]), , drop = FALSE]
          wa <- as.numeric(egf[[wt_col]])
          aa <- as.numeric(egf[[age_col]])
          if (length(wa) > 1 && length(aa) == length(wa)) {
            ## optional age-index alignment check against SS3 wtatage
            aa_shift <- aa
            if (!is.null(replist$wtatage) && nrow(replist$wtatage) > 0) {
              wt_ref <- subset(replist$wtatage, fleet == 1)
              if (nrow(wt_ref) > 0 && "sex" %in% names(wt_ref) && "year" %in% names(wt_ref)) {
                sx_w <- min(wt_ref$sex, na.rm = TRUE)
                wt_ref <- subset(wt_ref, sex == sx_w)
                yref <- min(wt_ref$year, na.rm = TRUE)
                wt_ref <- subset(wt_ref, year == yref)
                if (nrow(wt_ref) > 0) {
                  age_cols_w <- intersect(as.character(0:dat$amax), names(wt_ref))
                  if (length(age_cols_w) > 1) {
                    w_ref <- as.numeric(wt_ref[1, age_cols_w])
                    a_ref <- as.numeric(age_cols_w)
                    rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
                    w0 <- stats::approx(x = aa, y = wa, xout = a_ref, rule = 2)$y
                    w1 <- stats::approx(x = aa + 1, y = wa, xout = a_ref, rule = 2)$y
                    if (is.finite(rmse(w1, w_ref)) && rmse(w1, w_ref) + 1e-12 < rmse(w0, w_ref)) {
                      aa_shift <- aa + 1
                    }
                  }
                }
              }
            }

            wa_age <- stats::approx(x = aa_shift, y = wa, xout = age_classes, rule = 2)$y
            w_mat <- matrix(rep(rep(wa_age, each = dat$ns), dat$ny), nrow = dat$asmax, ncol = dat$ny)
            old_w <- dat$weight
            dat$weight <- w_mat
            dat$weightF <- w_mat
            add_report("weight", first_num(old_w), first_num(dat$weight), "ss3_endgrowth")

            if (!is.null(mat_age_col) && !is.null(mat_len_col)) {
              ma <- as.numeric(egf[[mat_age_col]])
              ml <- as.numeric(egf[[mat_len_col]])
              if (length(ma) == length(wa) && length(ml) == length(wa)) {
                mata <- ma * ml
                mata[!is.finite(mata)] <- 0
                mata[mata < 0] <- 0
                mata[mata > 1] <- 1
                mat_age <- stats::approx(x = aa_shift, y = mata, xout = age_classes, rule = 2)$y
                mat_full <- rep(mat_age, each = dat$ns)
                old_m <- dat$mat
                dat$mat <- mat_full
                add_report("mat", first_num(old_m), first_num(dat$mat), "ss3_endgrowth")
              }
            }
            used_mat_weight <- TRUE
          }
        }
      }

      ## 2) Fallback to wtatage/ratio approach
      if (!used_mat_weight && !is.null(replist$wtatage)) {
        wt <- replist$wtatage
        age_cols <- intersect(as.character(0:dat$amax), names(wt))

        if (length(age_cols) > 0) {
          ## weight-at-age from fleet 1 (female)
          wf <- subset(wt, fleet == 1)
          if (nrow(wf) > 0) {
            sx <- min(wf$sex, na.rm = TRUE)
            wf <- subset(wf, sex == sx)
            w_mat <- matrix(NA_real_, nrow = dat$asmax, ncol = dat$ny)

            for (iy in seq_len(dat$ny)) {
              y <- years_unique[iy]
              row_y <- subset(wf, year == y)
              if (nrow(row_y) < 1) row_y <- wf[which.min(abs(wf$year - y)), , drop = FALSE]
              wa <- as.numeric(row_y[1, age_cols])
              wa_exp <- rep(wa, each = dat$ns)
              w_mat[, iy] <- wa_exp
            }

            old_w <- dat$weight
            dat$weight <- w_mat
            dat$weightF <- w_mat
            add_report("weight", first_num(old_w), first_num(dat$weight), "ss3_wtatage")
          }

          ## maturity-at-age from fleet -2 / fleet 1 ratio where available
          wm <- subset(wt, fleet == -2)
          wf <- subset(wt, fleet == 1)
          if (nrow(wm) > 0 && nrow(wf) > 0) {
            sxm <- min(wm$sex, na.rm = TRUE)
            sxf <- min(wf$sex, na.rm = TRUE)
            yref <- max(years_unique)
            row_m <- subset(wm, year == yref & sex == sxm)
            row_f <- subset(wf, year == yref & sex == sxf)
            if (nrow(row_m) < 1) row_m <- wm[which.min(abs(wm$year - yref)), , drop = FALSE]
            if (nrow(row_f) < 1) row_f <- wf[which.min(abs(wf$year - yref)), , drop = FALSE]

            va <- as.numeric(row_m[1, age_cols])
            vb <- as.numeric(row_f[1, age_cols])
            ma <- rep(0, length(va))
            ok <- is.finite(va) & is.finite(vb) & vb > 0
            ma[ok] <- va[ok] / vb[ok]
            ma[ma < 0] <- 0
            ma[ma > 1] <- 1
            ma_exp <- rep(ma, each = dat$ns)

            old_m <- dat$mat
            dat$mat <- ma_exp
            add_report("mat", first_num(old_m), first_num(dat$mat), "ss3_wtatage_ratio")
          }
        }
      }
    }

    years_ss3 <- ts$Yr
    f_year <- rep(NA_real_, length(unique(years_ss3)))
    years_unique <- sort(unique(years_ss3))

    if (nrow(tar$annF_series) > 0) {
      f_year <- tar$annF_series$value[match(years_unique, tar$annF_series$year)]
    } else if ("F:_1" %in% names(ts)) {
      f_year <- as.numeric(ts$`F:_1`[match(years_unique, ts$Yr)])
    }
    f_year[!is.finite(f_year)] <- stats::median(f_year[is.finite(f_year)], na.rm = TRUE)
    if (all(!is.finite(f_year))) f_year <- rep(0, length(years_unique))

    fm_mat <- matrix(rep(f_year / dat$ns, each = dat$ns), nrow = dat$ny, ncol = dat$ns, byrow = TRUE)
    dat$FM <- fm_mat
    add_report("FM", first_num(old_FM), first_num(dat$FM), "ss3_annF_to_seasonal")
  }

  apply_scalar <- function(field, value, method = "direct") {
    if (!is.finite(value)) return(invisible(FALSE))
    old <- dat[[field]]
    dat[[field]] <<- value
    add_report(field, ifelse(length(old) > 0, old[1], NA_real_), value, method)
    TRUE
  }

  if ("sr" %in% condition_fields) {
    ok_r0 <- apply_scalar("R0", tar$R0)
    ok_h <- apply_scalar("h", tar$h)
    if (is.finite(tar$sigmaR)) apply_scalar("sd_r", tar$sigmaR)
    if (isTRUE(strict) && (!ok_r0 || !ok_h)) {
      stop("Could not condition SR parameters (R0/h) from SS3.")
    }
  }

  if ("q" %in% condition_fields) {
    ok_q <- apply_scalar("q", tar$q)
    if (isTRUE(strict) && !ok_q) stop("Could not condition q from SS3.")
  }

  if ("timing" %in% condition_fields) {
    if (is.finite(tar$pzbm)) {
      old <- dat$pzbm
      dat$pzbm <- tar$pzbm
      add_report("pzbm", first_num(old), first_num(dat$pzbm), "ss3_spawn_timing")
    } else if (isTRUE(verbose)) {
      message("Could not derive pzbm from SS3 timing fields; keeping existing dat$pzbm.")
    }
  }

  if ("M" %in% condition_fields) {
    if (is.finite(tar$M)) {
      old <- dat$M
      if (mode == "trajectory_conditioning") {
        dat$M[,] <- tar$M / dat$ns
        add_report("M", old[1, 1], dat$M[1, 1], "annual_M_to_seasonal")
      } else {
        dat$M[,] <- tar$M
        add_report("M", old[1, 1], tar$M, "scalar_to_matrix")
      }
    } else if (isTRUE(strict)) {
      stop("Could not condition natural mortality from SS3.")
    }
  }

  if ("growth" %in% condition_fields) {

    apply_scalar("Linf", tar$Linf)
    apply_scalar("K", tar$K)
    apply_scalar("t0", tar$t0)
    apply_scalar("Lm50", tar$Lm50)
    apply_scalar("lwa", tar$lwa)
    apply_scalar("lwb", tar$lwb)

    if (is.finite(tar$Lm50) && is.finite(tar$mat_slope) && tar$mat_slope != 0) {
      lm95 <- tar$Lm50 - log(19) / tar$mat_slope
      apply_scalar("Lm95", lm95, "mat50_slope")
    }

    if (is.null(dat$t0) || is.na(dat$t0)) dat$t0 <- 0

    if (is.finite(dat$Linf) && is.finite(dat$K) && length(dat$ages) > 0) {
      old_la <- dat$LA
      if (is.null(old_la)) {
        old_la <- dat$Linf * (1 - exp(-dat$K * (as.vector(t(dat$ages)) - dat$t0)))
      }
      la_vec <- dat$Linf * (1 - exp(-dat$K * (dat$ages - dat$t0)))
      dat$LA <- recycle_like(old_la, la_vec, ny = dat$ny, asmax = dat$asmax)
      add_report("LA", first_num(old_la), first_num(dat$LA), "vonbert_functional")
    }


    if(!any(names(dat) == "binwidth")){
      if(is.null(dat$LA)){
        dat$binwidth <- NULL
      }else dat$binwidth <- 1
    }
    binwidth <- dat$binwidth

    if (!is.null(dat$plba)) {
      plba_old <- dat$plba
      plba_age <- sapply(seq_len(ncol(plba_old)), function(j) {
        as.numeric(stats::aggregate(plba_old[, j], by = list(old_as2a), FUN = mean)[, 2])
      })
      dat$plba <- plba_age[rep(seq_len(nrow(plba_age)), each = dat$ns), , drop = FALSE]
    } else {
      if (!is.null(dat$LA)) {
        mids <- seq((dat$binwidth/2), dat$Linf * 1.2, by = dat$binwidth)
        dat$mids <- mids
        highs <- mids + binwidth/2
        lows <- mids - binwidth/2
        lbprobs <- function(mnl, sdl){
          return(stats::pnorm(highs, mnl, sdl) - stats::pnorm(lows, mnl, sdl))
        }
        vlprobs <- Vectorize(lbprobs, vectorize.args=c("mnl","sdl"))
        if(!any(names(dat) == "CVlen")){
          dat$CVlen <- 0.1
        }
        LA2 <- array(dat$LA, dim = c(dat$asmax,1))
        plba <- apply(apply(LA2, c(1), function(x) vlprobs(x, x * dat$CVlen)), c(1), t)
        ## for(i in 1:dim(plba)[3]){
        ##     tmp <- rowSums(plba[,,i])
        ##     plba[,,i] <- plba[,,i] / tmp
        ##     plba[tmp == 0,,i] <- 0
        ## }
        tmp <- rowSums(plba)
        dat$plba <- plba / tmp
        ##            plba[tmp == 0] <- 0
        ## plba <- t(vlprobs(LA, LA * dat$CVlen))
        ## plba <- plba / rowSums(plba)
      }
    }

    if (is.null(dat$pabl)) {
      dat$pabl <- t(dat$plba)
    }

    if (!("mat_weight" %in% condition_fields) && !is.null(dat$LA) && is.finite(dat$lwa) && is.finite(dat$lwb)) {
      old_w <- dat$weight
      la_base <- if (is.list(dat$LA)) dat$LA[[1]] else if (is.matrix(dat$LA)) dat$LA[1, ] else dat$LA
      w_vec <- dat$lwa * (la_base ^ dat$lwb)
      dat$weight <- recycle_like(old_w, w_vec, ny = dat$ny, asmax = dat$asmax)
      add_report("weight", first_num(old_w), first_num(dat$weight), "lw_functional")
    }
  }

  if ("selectivity" %in% condition_fields) {
    sel_from_ss3_age <- NULL
    sel_from_ss3_len <- NULL
    sel_factor_used <- NA_character_

    ## 1) Prefer age selectivity (Asel2)
    if (!is.null(replist$ageselex) && nrow(replist$ageselex) > 0) {
      asel <- subset(replist$ageselex, Fleet == 1)
      if (nrow(asel) > 0) {
        if ("Sex" %in% names(asel)) asel <- subset(asel, Sex == min(Sex, na.rm = TRUE))
        if ("Seas" %in% names(asel)) asel <- subset(asel, Seas == min(Seas, na.rm = TRUE))

        age_cols <- grep("^[0-9]+$", names(asel), value = TRUE)
        if (length(age_cols) > 1) {
          aa <- as.numeric(age_cols)
          age_target <- 0:dat$amax

          pick_factor <- function(tab) {
            if (!("Factor" %in% names(tab))) return(list(tab = tab, factor = "Asel"))
            fac_priority <- c("Asel2", "Asel", "sel*wt", "dead*wt", "sel_nums", "dead_nums")
            fac_all <- unique(as.character(tab$Factor))
            fac_try <- unique(c(fac_priority[fac_priority %in% fac_all], fac_all))
            best_fac <- NA_character_
            best_rng <- -Inf
            for (fac in fac_try) {
              subf <- subset(tab, Factor == fac)
              if (nrow(subf) < 1) next
              ss <- as.numeric(subf[1, age_cols])
              ss <- ss[is.finite(ss)]
              if (length(ss) < 2) next
              rr <- diff(range(ss))
              if (is.finite(rr) && rr > best_rng) {
                best_rng <- rr
                best_fac <- fac
              }
            }
            if (!is.finite(best_rng) || best_rng < 0.05) {
              if ("Asel2" %in% fac_all) best_fac <- "Asel2" else if ("Asel" %in% fac_all) best_fac <- "Asel" else best_fac <- fac_all[1]
            }
            list(tab = subset(tab, Factor == best_fac), factor = best_fac)
          }

          if (mode == "trajectory_conditioning" && "Yr" %in% names(asel)) {
            ts_now <- subset(replist$timeseries, Era == "TIME" & Yr <= replist$endyr)
            years_unique <- sort(unique(ts_now$Yr))
            if (length(years_unique) > 0) {
              sel_by_year <- vector("list", length(years_unique))
              for (iy in seq_along(years_unique)) {
                y <- years_unique[iy]
                row_y <- subset(asel, Yr == y)
                if (nrow(row_y) < 1) row_y <- asel[which.min(abs(asel$Yr - y)), , drop = FALSE]
                pf <- pick_factor(row_y)
                row_y <- pf$tab
                if (is.na(sel_factor_used)) sel_factor_used <- pf$factor
                row_y <- row_y[1, , drop = FALSE]
                ss <- as.numeric(row_y[1, age_cols])
                s_age <- stats::approx(x = aa, y = ss, xout = age_target, rule = 2)$y
                s_age[s_age < 0] <- 0
                s_age[s_age > 1] <- 1
                sel_by_year[[iy]] <- rep(s_age, each = dat$ns)
              }
              old_sel <- dat$sel
              dat$sel <- sel_by_year
              add_report("sel", first_num(old_sel), first_num(dat$sel), "ss3_ageselex_timevarying")
              sel_from_ss3_age <- sel_by_year[[length(sel_by_year)]]
            }
          }

          if (is.null(sel_from_ss3_age)) {
            y_last <- max(asel$Yr, na.rm = TRUE)
            row_last <- subset(asel, Yr == y_last)
            if (nrow(row_last) > 0) {
              pf <- pick_factor(row_last)
              row_last <- pf$tab
              sel_factor_used <- pf$factor
              row_last <- row_last[1, , drop = FALSE]
              ss <- as.numeric(row_last[1, age_cols])
              s_age <- stats::approx(x = aa, y = ss, xout = age_target, rule = 2)$y
              sel_from_ss3_age <- rep(s_age, each = dat$ns)
              sel_from_ss3_age[sel_from_ss3_age < 0] <- 0
              sel_from_ss3_age[sel_from_ss3_age > 1] <- 1
            }
          }
        }
      }
    }

    ## 2) Fallback to length selectivity mapped with LA
    if (is.null(sel_from_ss3_age) && !is.null(replist$sizeselex) && nrow(replist$sizeselex) > 0 && !is.null(dat$LA)) {
      ss_sel <- subset(replist$sizeselex, Factor == "Lsel" & Fleet == 1)
      if (nrow(ss_sel) > 0) {
        y_last <- max(ss_sel$Yr, na.rm = TRUE)
        ss_sel <- subset(ss_sel, Yr == y_last)
        if ("Sex" %in% names(ss_sel)) ss_sel <- subset(ss_sel, Sex == min(Sex, na.rm = TRUE))
        ss_sel <- ss_sel[1, , drop = FALSE]
        num_cols <- grep("^-?[0-9]+\\.?[0-9]*$", names(ss_sel), value = TRUE)
        if (length(num_cols) > 1) {
          l_mid <- as.numeric(num_cols)
          l_sel <- as.numeric(ss_sel[1, num_cols])
          oo <- order(l_mid)
          l_mid <- l_mid[oo]
          l_sel <- l_sel[oo]
          la_base <- if (is.list(dat$LA)) dat$LA[[1]] else if (is.matrix(dat$LA)) dat$LA[1, ] else dat$LA
          sel_from_ss3_len <- interp_curve(l_mid, l_sel, as.numeric(la_base), method = ifelse(interp_method == "linear", "linear", "mono_spline"))
          sel_from_ss3_len[sel_from_ss3_len < 0] <- 0
          sel_from_ss3_len[sel_from_ss3_len > 1] <- 1
        }
      }
    }

    if (!is.null(sel_from_ss3_age) || !is.null(sel_from_ss3_len)) {
      sel_vec <- if (!is.null(sel_from_ss3_age)) sel_from_ss3_age else sel_from_ss3_len
      method_lab <- if (!is.null(sel_from_ss3_age)) {
        if (!is.na(sel_factor_used)) paste0("ss3_ageselex_", sel_factor_used) else "ss3_ageselex"
      } else {
        "ss3_sizeselex_interp"
      }
      if (!(mode == "trajectory_conditioning" && is.list(dat$sel) && length(dat$sel) == dat$ny && !is.null(sel_from_ss3_age))) {
        old_sel <- dat$sel
        dat$sel <- recycle_like(old_sel, sel_vec, ny = dat$ny, asmax = dat$asmax)
        add_report("sel", first_num(old_sel), first_num(dat$sel), method_lab)
      }

      la_base <- if (is.list(dat$LA)) dat$LA[[1]] else if (is.matrix(dat$LA)) dat$LA[1, ] else dat$LA
      if (!is.null(la_base) && length(la_base) == length(sel_vec)) {
        s_ok <- is.finite(sel_vec)
        s_rng <- if (any(s_ok)) diff(range(sel_vec[s_ok])) else NA_real_
        has_l50 <- any(s_ok & sel_vec <= 0.5) && any(s_ok & sel_vec >= 0.5)
        has_l95 <- any(s_ok & sel_vec <= 0.95) && any(s_ok & sel_vec >= 0.95)

        ## avoid nonsensical Ls50/Ls95 when age selectivity is effectively flat
        if (is.finite(s_rng) && s_rng > 0.1 && has_l50 && has_l95) {
          l50_i <- which.min(abs(sel_vec - 0.5))
          l95_i <- which.min(abs(sel_vec - 0.95))
          if (is.finite(la_base[l50_i])) apply_scalar("Ls50", la_base[l50_i], method_lab)
          if (is.finite(la_base[l95_i])) apply_scalar("Ls95", la_base[l95_i], method_lab)
        } else if (is.finite(tar$Ls50) && is.finite(tar$Ls95)) {
          ## preserve SS3 length-based selectivity proxies when available
          apply_scalar("Ls50", tar$Ls50, "ss3_target_fallback")
          apply_scalar("Ls95", tar$Ls95, "ss3_target_fallback")
        }
      }

    } else if (is.finite(tar$Ls50) && is.finite(tar$Ls95) && !is.null(dat$LA)) {
      old_sel <- dat$sel
      k <- log(19) / (tar$Ls95 - tar$Ls50)
      la_base <- if (is.list(dat$LA)) dat$LA[[1]] else if (is.matrix(dat$LA)) dat$LA[1, ] else dat$LA
      sel_age <- 1 / (1 + exp(-k * (la_base - tar$Ls50)))
      sel_age[!is.finite(sel_age)] <- 0
      dat$sel <- recycle_like(old_sel, sel_age, ny = dat$ny, asmax = dat$asmax)
      add_report("sel", first_num(old_sel), first_num(dat$sel), "logistic_functional")
      apply_scalar("Ls50", tar$Ls50)
      apply_scalar("Ls95", tar$Ls95)
    } else {
      if (isTRUE(strict)) stop("Could not condition selectivity from SS3.")
      if (isTRUE(verbose)) message("Skipping selectivity conditioning: missing SS3 selectivity inputs.")
    }

    ## assumption for now: survey selectivity == fishery selectivity
    ## TODO: improve
    dat$selI <- dat$sel

  }

  if ("ref_ss3" %in% condition_fields) {
    dat$ref_ss3 <- tar$ref_ss3
    if (!is.null(dat$ref) && is.data.frame(dat$ref)) {
      dat$ref_ss3$ny <- nrow(dat$ref)
    }
  }


  ## simplify matrices
  if (inherits(dat$weight, "matrix")) {
    if (all(apply(dat$weight, 1, sd) == 0)) {
      dat$weight <- dat$weight[,1]
    }
  }
  if (inherits(dat$weightF, "matrix")) {
    if (all(apply(dat$weightF, 1, sd) == 0)) {
      dat$weightF <- dat$weightF[,1]
    }
  }
  if (inherits(dat$mat, "matrix")) {
    if (all(apply(dat$mat, 1, sd) == 0)) {
      dat$mat <- dat$mat[,1]
    }
  }
  if (inherits(dat$sel, "list")) {
    if (all(apply(do.call(cbind,dat$sel), 1, sd) == 0)) {
      dat$sel <- dat$sel[1]
    }
  }
  if (inherits(dat$selI, "list")) {
    if (all(apply(do.call(cbind,dat$selI), 1, sd) == 0)) {
      dat$selI <- dat$selI[1]
    }
  }

  dat$conditioning_ss3_raw <- tar$raw

  out <- list(
    dat = dat,
    report = report,
    extracted = tar,
    replist = replist)
  out
}


#' Extract conditioning targets from an SS3 report object
#'
#' Pulls scalar parameters, derived quantities, time series, and auxiliary
#' tables used to condition iamse objects.
#'
#' @param replist SS3 output list from `r4ss::SS_output()`.
#' @param verbose Logical; emit messages when values are missing.
#'
#' @return Named list of extracted targets and raw source tables.
#' @keywords internal
#' @noRd
.extract_ss3_targets <- function(replist, verbose = TRUE) {

  if (is.null(replist$parameters) || is.null(replist$timeseries)) {
    stop("replist must contain at least 'parameters' and 'timeseries'.")
  }

  p <- replist$parameters
  ts <- replist$timeseries
  dq <- replist$derived_quants

  get_used <- function(pattern, log_scale = FALSE) {
    ii <- grep(pattern, p$Label, ignore.case = TRUE)
    if (length(ii) < 1) return(NA_real_)
    val <- as.numeric(p$Value[ii[1]])
    if (!is.finite(val)) val <- as.numeric(p$Value[ii[1]])
    if (!is.finite(val)) return(NA_real_)
    if (isTRUE(log_scale)) return(exp(val))
    val
  }

  get_dq <- function(label) {
    if (is.null(dq) || !("Label" %in% names(dq))) return(NA_real_)
    ii <- which(dq$Label == label)
    if (length(ii) < 1) return(NA_real_)
    as.numeric(dq$Value[ii[1]])
  }

  get_last_dq_series <- function(prefix) {
    if (is.null(dq) || !("Label" %in% names(dq))) return(NA_real_)
    ii <- grep(paste0("^", prefix, "_[0-9]+$"), dq$Label)
    if (length(ii) < 1) return(NA_real_)
    yrs <- as.integer(sub(paste0("^", prefix, "_"), "", dq$Label[ii]))
    i_last <- ii[which.max(yrs)]
    as.numeric(dq$Value[i_last])
  }

  get_dq_series <- function(prefix) {
    if (is.null(dq) || !("Label" %in% names(dq))) {
      return(data.frame(year = integer(), value = numeric()))
    }
    ii <- grep(paste0("^", prefix, "_[0-9]+$"), dq$Label)
    if (length(ii) < 1) {
      return(data.frame(year = integer(), value = numeric()))
    }
    yrs <- as.integer(sub(paste0("^", prefix, "_"), "", dq$Label[ii]))
    vals <- as.numeric(dq$Value[ii])
    oo <- order(yrs)
    data.frame(year = yrs[oo], value = vals[oo])
  }

  out <- list(
    R0 = get_used("SR_LN\\(R0\\)|LnR0|log_R0|R0", log_scale = TRUE),
    h = get_used("SR_BH_steep|steep"),
    sigmaR = get_used("SR_sigmaR|sigmaR"),
    q = get_used("LnQ", log_scale = TRUE),
    M = get_used("NatM.*Fem|NatM_p_1_Fem|NatM_uniform_Fem"),
    Linf = get_used("L_at_Amax_Fem"),
    K = get_used("VonBert_K_Fem|K_Fem"),
    t0 = get_used("VonBert_t0_Fem|t0_Fem"),  ## TODO: where to get t0
    Lm50 = get_used("Mat50|Mat50%"),
    mat_slope = get_used("Mat_slope"),
    lwa = get_used("Wtlen_1_Fem"),
    lwb = get_used("Wtlen_2_Fem"),
    Ls50 = get_used("SizeSel_P_1_Fishery|Size_inflection_Fishery"),
    Ls95_width = get_used("SizeSel_P_2_Fishery|Size_95%width_Fishery"),
    ref_ss3 = list(
      SSB_MSY = get_dq("SSB_MSY"),
      annF_MSY = get_dq("annF_MSY"),
      Dead_Catch_MSY = get_dq("Dead_Catch_MSY"),
      SSB_Virgin = get_dq("SSB_Virgin"),
      Bratio_last = get_last_dq_series("Bratio")
    ),
    annF_series = get_dq_series("annF"),
    Bratio_series = get_dq_series("Bratio")
  )

  ## derive iamse-style pzbm from SS3 spawning timing
  spawn_seas <- suppressWarnings(as.integer(replist$Spawn_seas[1]))
  spawn_tis <- suppressWarnings(as.numeric(replist$Spawn_timing_in_season[1]))
  spawn_month <- suppressWarnings(as.numeric(replist$Spawn_month[1]))
  seas_dur <- suppressWarnings(as.numeric(replist$seasdurations))

  pzbm <- NA_real_
  if (length(seas_dur) > 0 && all(is.finite(seas_dur)) &&
        is.finite(spawn_seas) && spawn_seas >= 1 && spawn_seas <= length(seas_dur) &&
          is.finite(spawn_tis)) {
    before <- if (spawn_seas > 1) sum(seas_dur[seq_len(spawn_seas - 1)]) else 0
    within <- seas_dur[spawn_seas] * spawn_tis
    pzbm <- before + within
    if (is.finite(pzbm)) {
      pzbm <- max(0, min(1, pzbm))
    }
  } else if (is.finite(spawn_month)) {
    ## fallback: approximate spawning timing from month if within-season timing unavailable
    pzbm <- (spawn_month - 1) / 12
    pzbm <- max(0, min(1, pzbm))
  }
  out$pzbm <- pzbm

  if (is.finite(out$Ls50) && is.finite(out$Ls95_width) && out$Ls95_width > 0) {
    out$Ls95 <- out$Ls50 + out$Ls95_width
  } else {
    out$Ls95 <- NA_real_
  }

  out$raw <- list(
    parameters = p,
    timeseries = ts,
    derived_quants = dq,
    endgrowth = replist$endgrowth,
    wtatage = replist$wtatage,
    ageselex = replist$ageselex,
    sizeselex = replist$sizeselex,
    recruit = replist$recruit
  )

  if (isTRUE(verbose)) {
    miss <- names(out)[vapply(out, function(x) is.numeric(x) && !is.finite(x), logical(1))]
    if (length(miss) > 0) {
      message("Missing/non-finite extracted values: ", paste(miss, collapse = ", "))
    }
  }

  out
}



#' Map SS3 recruitment deviations to iamse recruitment errors
#'
#' Derives `set$errs$time$eR` from SS3 recruit outputs or `Main_RecrDev_YYYY`
#' parameters and maps them to requested years.
#'
#' @param set iamse settings list.
#' @param replist SS3 output list.
#' @param dat Optional iamse data list used when `years` is not supplied.
#' @param years Optional integer vector of years to map to.
#' @param method Mapping source: predicted/expected recruitment ratio or
#'   `Main_RecrDev` parameters.
#' @param use_bias_correction Logical; apply `exp(-0.5*sigmaR^2)` correction.
#' @param verbose Logical; print summary mapping message.
#'
#' @return A list with updated `set` and a `map` data frame.
#' @keywords internal
#' @noRd
.set_iamse_recruitment_errors_from_ss3 <- function(set,
                                                   replist,
                                                   dat = NULL,
                                                   years = NULL,
                                                   method = c("pred_exp", "main_recdev"),
                                                   use_bias_correction = FALSE,
                                                   verbose = TRUE) {

  method <- match.arg(method)

  if (is.null(replist$parameters) || !("Label" %in% names(replist$parameters))) {
    stop("replist$parameters with SS3 labels is required.")
  }

  if (is.null(years)) {
    if (is.null(dat) || is.null(dat$ny)) {
      stop("Provide either 'years' or a 'dat' object with dat$ny.")
    }
    years <- seq_len(dat$ny)
  }
  years <- as.integer(years)

  rec_years <- NULL
  rec_dev <- NULL
  eR <- NULL

  if (method == "pred_exp" && !is.null(replist$recruit)) {
    rtab <- replist$recruit
    yr_col <- if ("Yr" %in% names(rtab)) "Yr" else if ("year" %in% names(rtab)) "year" else NULL
    pr_col <- if ("pred_recr" %in% names(rtab)) "pred_recr" else NULL
    ex_col <- if ("exp_recr" %in% names(rtab)) "exp_recr" else NULL
    if (!is.null(yr_col) && !is.null(pr_col) && !is.null(ex_col)) {
      rec_years <- as.integer(rtab[[yr_col]])
      pr <- as.numeric(rtab[[pr_col]])
      ex <- as.numeric(rtab[[ex_col]])
      eR <- pr / ex
      rec_dev <- log(eR)
    }
  }

  if (is.null(eR) || all(!is.finite(eR))) {
    p <- replist$parameters
    ii <- grep("^Main_RecrDev_[0-9]+$", p$Label)
    if (length(ii) < 1) {
      stop("No SS3 recruitment deviations found (labels like Main_RecrDev_YYYY).")
    }
    rec_years <- as.integer(sub("^Main_RecrDev_", "", p$Label[ii]))
    rec_dev <- as.numeric(p$Value[ii])
    bad <- !is.finite(rec_dev)
    if (any(bad)) rec_dev[bad] <- as.numeric(p$Value[ii][bad])
    eR <- exp(rec_dev)
  }

  if (isTRUE(use_bias_correction)) {
    sig <- NA_real_
    if (!is.null(replist$parameters) && "Label" %in% names(replist$parameters)) {
      ip <- grep("SR_sigmaR|sigmaR", replist$parameters$Label)
      if (length(ip) > 0) sig <- as.numeric(replist$parameters$Value[ip[1]])
    }
    if (!is.finite(sig)) sig <- as.numeric(stats::sd(rec_dev[is.finite(rec_dev)]))
    if (is.finite(sig) && sig > 0) {
      eR <- eR * exp(-0.5 * sig^2)
    }
  }

  map_idx <- match(years, rec_years)
  eR_mapped <- eR[map_idx]

  if (all(!is.finite(eR_mapped))) {
    stop("Could not map any SS3 recdev years to requested years.")
  }

  if (any(!is.finite(eR_mapped))) {
    med <- stats::median(eR[is.finite(eR)], na.rm = TRUE)
    if (!is.finite(med)) med <- 1
    eR_mapped[!is.finite(eR_mapped)] <- med
  }

  if (is.null(set$errs) || !is.list(set$errs)) set$errs <- list()
  if (is.null(set$errs$time) || !is.list(set$errs$time)) set$errs$time <- list()
  if (is.null(set$errs$rep) || !is.list(set$errs$rep)) set$errs$rep <- list()

  set$errs$time$eR <- as.numeric(eR_mapped)
  set$errs$rep$eR <- 1

  if (isTRUE(verbose)) {
    message(
      "Mapped SS3 recdevs to iamse eR for ", length(years),
      " years (", min(years), "-", max(years), ")."
    )
  }

  list(
    set = set,
    map = data.frame(
      year = years,
      ss3_recdev = rec_dev[map_idx],
      eR = set$errs$time$eR
    )
  )
}



#' Infer SS3 maximum age from report tables
#'
#' Infers iamse-style `amax` from numeric age columns in SS3 tables, with
#' fallback to `Nages`/`nages` fields.
#'
#' @param replist SS3 output list.
#'
#' @return Integer `amax` when inferred; `NA_integer_` otherwise.
#' @keywords internal
#' @noRd
.infer_ss3_amax <- function(replist) {
  extract_age_cols <- function(tab) {
    if (is.null(tab) || !is.data.frame(tab) || ncol(tab) < 1) return(integer())
    nms <- names(tab)
    keep <- grepl("^[0-9]+$", nms)
    if (!any(keep)) return(integer())
    ages <- suppressWarnings(as.integer(nms[keep]))
    ages[is.finite(ages)]
  }

  age_candidates <- c(
    extract_age_cols(replist$natage_annual_2_with_fishery),
    extract_age_cols(replist$wtatage),
    extract_age_cols(replist$ageselex)
  )

  if (length(age_candidates) > 0) {
    age_candidates <- age_candidates[is.finite(age_candidates)]
    if (length(age_candidates) > 0) {
      amax <- max(age_candidates, na.rm = TRUE)
      if (is.finite(amax) && amax >= 1) return(as.integer(amax))
    }
  }

  nages_candidates <- c(replist$Nages, replist$nages)
  nages_candidates <- suppressWarnings(as.integer(nages_candidates))
  nages_candidates <- nages_candidates[is.finite(nages_candidates) & nages_candidates >= 2]
  if (length(nages_candidates) > 0) {
    return(as.integer(max(nages_candidates, na.rm = TRUE) - 1L))
  }

  warning("Could not infer SS3 maximum age from age columns or nages. Keeping existing dat$amax.")
  NA_integer_
}



#' Rebuild iamse age-time indexing vectors
#'
#' Recomputes dimensions and index vectors (`ages`, `as2a`, `as2s`, `s1avec`,
#' `yvec`, `svec`, `s1vec`, etc.) after changing years/seasons and optionally
#' maximum age.
#'
#' @param dat iamse data list.
#' @param ny_new New number of years.
#' @param ns_new New number of seasons.
#' @param amax_new Optional new maximum age to set before rebuilding indices.
#'
#' @return Updated `dat` list.
#' @keywords internal
#' @noRd
.rebuild_age_time_structure <- function(dat, ny_new, ns_new, amax_new = NULL) {
  if (!is.null(amax_new) && is.finite(amax_new) && amax_new >= 1) {
    dat$amax <- as.integer(amax_new)
  }

  amax_age <- dat$amax + 1
  asmax_new <- amax_age * ns_new
  age_idx <- 0:(amax_age - 1)
  season_frac <- ((1:ns_new) - 0.5) / ns_new
  ages_mat <- outer(age_idx, season_frac, "+")

  dat$ny <- ny_new
  dat$ns <- ns_new
  dat$asmax <- asmax_new
  dat$ages <- as.numeric(t(ages_mat))
  dat$as2a <- rep(1:amax_age, each = ns_new)
  dat$as2s <- rep(1:ns_new, times = amax_age)
  dat$s1avec <- seq(1, asmax_new, by = ns_new)
  dat$indage0 <- 1 ## which(dat$as2a == 1)  ## TODO error with ns = 4
  dat$yvec <- rep(1:ny_new, each = ns_new)
  dat$svec <- rep(1:ns_new, times = ny_new)
  dat$s1vec <- seq(1, ny_new * ns_new, by = ns_new)

  dat$spawning <- .resample_spawning(dat$spawning, ns_new)

  if (!is.null(dat$catchSeasons)) dat$catchSeasons <- 1 ## annual catches
  if (!is.null(dat$effortSeasons)) dat$effortSeasons <- 1 ## annual effort

  dat
}


#' Resample spawning timing to a new number of seasons
#'
#' Transfers spawning fractions from an existing seasonal partition to a new
#' seasonal grid while preserving total spawning probability.
#'
#' @param old_sp Numeric vector of old seasonal spawning fractions.
#' @param ns_new Target number of seasons.
#'
#' @return Numeric vector of length `ns_new` summing to 1 (or zeros when
#'   spawning is absent).
#' @keywords internal
#' @noRd
.resample_spawning <- function(old_sp, ns_new) {

  out <- rep(0, ns_new)

  if (is.null(old_sp) || !any(old_sp > 0, na.rm = TRUE)) {
    return(out)
  }

  old_sp[is.na(old_sp)] <- 0

  ## Special case:
  ## annual spawning = 1 converted to subannual time steps
  ## means all spawning occurs in the first time step
  if (length(old_sp) == 1 && ns_new > 1) {
    out[1] <- 1
    return(out)
  }

  old_sp <- old_sp / sum(old_sp)

  if (ns_new == 1) {
    return(1)
  }

  ns_old <- length(old_sp)

  old_breaks <- seq(0, 1, length.out = ns_old + 1)
  new_breaks <- seq(0, 1, length.out = ns_new + 1)

  for (i in seq_len(ns_old)) {
    for (j in seq_len(ns_new)) {

      overlap <- max(
        0,
        min(old_breaks[i + 1], new_breaks[j + 1]) -
          max(old_breaks[i], new_breaks[j])
      )

      if (overlap > 0) {
        out[j] <- out[j] + old_sp[i] * overlap / (1 / ns_old)
      }
    }
  }

  out / sum(out)
}

#' Expand age-based iamse processes to a target seasonal resolution
#'
#' Rebuilds indexing and rescales/repeats age-specific vectors and matrices so
#' they are consistent with a new `ns`.
#'
#' @param dat iamse data list.
#' @param ns_target Target number of seasons.
#'
#' @return Updated `dat` list at the requested seasonal resolution.
#' @keywords internal
#' @noRd
.expand_processes_to_ns <- function(dat, ns_target) {
  dat_old <- dat
  s1_old <- dat_old$s1avec
  dat <- .rebuild_age_time_structure(dat_old, ny_new = dat_old$ny, ns_new = ns_target)

  if (is.matrix(dat$M)) {
    m_ann <- dat$M[, 1]
  } else {
    m_ann <- as.numeric(dat$M)
    if (length(m_ann) == 1) m_ann <- rep(m_ann, dat$ny)
  }
  dat$M <- matrix(rep(m_ann / ns_target, each = ns_target), nrow = dat$ny, ncol = ns_target, byrow = TRUE)

  if (is.matrix(dat$FM)) {
    f_ann <- rowSums(dat$FM)
  } else {
    f_ann <- as.numeric(dat$FM)
    if (length(f_ann) == 1) f_ann <- rep(f_ann, dat$ny)
  }
  dat$FM <- matrix(rep(f_ann / ns_target, each = ns_target), nrow = dat$ny, ncol = ns_target, byrow = TRUE)

  expand_age_vec <- function(v) {
    v <- as.numeric(v)
    rep(v, each = ns_target)
  }

  if (!is.null(dat$LA)) {
    la_age <- dat_old$LA[s1_old]
    dat$LA <- expand_age_vec(la_age)
  }

  if (!is.null(dat$mat)) {
    mat_age <- as.numeric(dat_old$mat[s1_old])
    dat$mat <- expand_age_vec(mat_age)
  }

  if (!is.null(dat$weight)) {
    if (is.matrix(dat$weight)) {
      w_age <- dat_old$weight[s1_old, , drop = FALSE]
      w_new <- matrix(NA_real_, nrow = dat$asmax, ncol = dat$ny)
      for (iy in seq_len(dat$ny)) {
        w_new[, iy] <- expand_age_vec(w_age[, iy])
      }
      dat$weight <- w_new
      dat$weightF <- w_new
    } else {
      w_age <- dat_old$weight[s1_old]
      dat$weight <- matrix(rep(expand_age_vec(w_age), dat$ny), nrow = dat$asmax, ncol = dat$ny)
      dat$weightF <- dat$weight
    }
  }

  if (!is.null(dat$Msel) && is.list(dat$Msel) && length(dat$Msel) >= 1) {
    ms_age <- as.numeric(dat_old$Msel[[1]][s1_old])
    dat$Msel <- list(expand_age_vec(ms_age))
  }

  if (!is.null(dat$sel) && is.list(dat$sel) && length(dat$sel) >= 1) {
    dat$sel <- lapply(dat$sel, function(v) {
      v_age <- as.numeric(v[s1_old])
      expand_age_vec(v_age)
    })
  }

  dat
}


#' Calculate iamse vs SS3 trajectory comparison metrics
#'
#' Builds aligned catch, fishing mortality, and SSB trajectories and computes
#' simple RMSE and bias diagnostics.
#'
#' @param pop iamse population output object.
#' @param replist SS3 output list.
#' @param extracted List of extracted SS3 targets.
#' @param ss3_ssb_source Which SS3 SSB field to use.
#'
#' @return A list with aligned series and summary metric table.
#' @keywords internal
#' @noRd
.calc_traj_metrics <- function(pop, replist, extracted,
                               ss3_ssb_source = c("auto", "SpawnBio", "mature_bio")) {
  ss3_ssb_source <- match.arg(ss3_ssb_source)
  ts <- subset(replist$timeseries, Era == "TIME" & Yr <= replist$endyr)
  ss_years <- ts$Yr
  ss_ssb <- NULL
  ss_ssb_field <- NULL
  if (ss3_ssb_source == "SpawnBio") {
    ss_ssb <- ts$SpawnBio
    ss_ssb_field <- "SpawnBio"
  } else if (ss3_ssb_source == "mature_bio") {
    ss_ssb <- ts$mature_bio
    ss_ssb_field <- "mature_bio"
  } else {
    if ("SpawnBio" %in% names(ts)) {
      ss_ssb <- ts$SpawnBio
      ss_ssb_field <- "SpawnBio"
    } else if ("mature_bio" %in% names(ts)) {
      ss_ssb <- ts$mature_bio
      ss_ssb_field <- "mature_bio"
    } else {
      stop("No SS3 SSB field available in timeseries")
    }
  }
  ss_catch <- ts$`dead(B):_1`
  if (nrow(extracted$annF_series) > 0) {
    ss_f <- extracted$annF_series$value[match(ss_years, extracted$annF_series$year)]
  } else {
    ss_f <- ts$`F:_1`
  }

  ny_cmp <- min(length(ss_years), nrow(pop$SSB), nrow(pop$FM), nrow(pop$CW))
  yrs <- ss_years[seq_len(ny_cmp)]
  ia_ssb <- pop$SSB[seq_len(ny_cmp),1]
  ia_f <- rowSums(pop$FM[seq_len(ny_cmp), , drop = FALSE])
  ia_c <- rowSums(pop$CW[seq_len(ny_cmp), , drop = FALSE])
  ss_ssb <- ss_ssb[seq_len(ny_cmp)]
  ss_f <- ss_f[seq_len(ny_cmp)]
  ss_catch <- ss_catch[seq_len(ny_cmp)]

  rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
  bias <- function(x, y) mean(x - y, na.rm = TRUE)

  list(
    years = yrs,
    iamse = list(catch = ia_c, F = ia_f, SSB = ia_ssb),
    ss3 = list(catch = ss_catch, F = ss_f, SSB = ss_ssb),
    meta = list(ss3_ssb_field = ss_ssb_field),
    metrics = data.frame(
      series = c("catch", "F", "SSB"),
      rmse = c(rmse(ia_c, ss_catch), rmse(ia_f, ss_f), rmse(ia_ssb, ss_ssb)),
      bias = c(bias(ia_c, ss_catch), bias(ia_f, ss_f), bias(ia_ssb, ss_ssb))
    )
  )
}


#' Plot iamse and SS3 trajectory comparisons
#'
#' @param tr Trajectory object produced by `.calc_traj_metrics()`.
#'
#' @return Invisibly returns `NULL`; called for side-effect plotting.
#' @keywords internal
#' @noRd
.plot_traj <- function(tr) {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(3, 1), mar = c(3, 4, 2, 1))

  plot(tr$years, tr$iamse$catch, type = "n", lwd = 2,
       ylim = range(tr$iamse$catch, tr$ss3$catch, na.rm = TRUE),
       xlab = "Year", ylab = "Catch")
  lines(tr$years, tr$ss3$catch, col = "black", lwd = 2, lty = 1)
  lines(tr$years, tr$iamse$catch, col = "dodgerblue2", lwd = 2, lty = 2)
  legend("topright", legend = c("SS3", "iamse"),
         col = c("black", "dodgerblue2"), lty = c(1,2),
         lwd = 2,
         bg = "white")

  plot(tr$years, tr$iamse$F, type = "n", lwd = 2,
       ylim = range(tr$iamse$F, tr$ss3$F, na.rm = TRUE),
       xlab = "Year", ylab = "F", col = "black")
  lines(tr$years, tr$ss3$F, col = "black", lwd = 2, lty = 1)
  lines(tr$years, tr$iamse$F, col = "dodgerblue2", lwd = 2, lty = 2)

  plot(tr$years, tr$iamse$SSB, type = "n", lwd = 2,
       ylim = range(tr$iamse$SSB, tr$ss3$SSB, na.rm = TRUE),
       xlab = "Year", ylab = "SSB", col = "black")
  lines(tr$years, tr$ss3$SSB, col = "black", lwd = 2, lty = 1)
  lines(tr$years, tr$iamse$SSB, col = "dodgerblue2", lwd = 2, lty = 2)
}
