

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
    L_at_Amin <- min(sim_dat$LA[is.finite(sim_dat$LA) & sim_dat$LA > 0], na.rm = TRUE)
    L_at_Amax <- sim_dat$Linf
    mat_slope <- -log(19) / (sim_dat$Lm95 - sim_dat$Lm50)

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "NatM_p_1_Fem", M_value,
      lo = 0.001, hi = max(2, M_value * 3),
      phase = -1, pr_type = 3, pr_sd = 10
    )

    min_pop_len <- min(ss_dat$lbin_vector_pop)
    max_pop_len <- max(ss_dat$lbin_vector_pop)

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "L_at_Amin_Fem",
      0, ## L_at_Amin,
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
      ctl$MG_parms, "CV_young_Fem", 0.1,
      lo = 0.001, hi = 5,
      phase = -4, pr_sd = 0.5
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "CV_old_Fem", 0.1,
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
      ctl$MG_parms, "Mat50", sim_dat$Lm50,
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
      ctl$MG_parms, "Eggs_alpha", sim_dat$lwa,
      lo = 0,
      hi = 3,
      phase = -3, pr_sd = 0.8
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Eggs_beta", sim_dat$lwb,
      lo = 0,
      hi = 10,
      phase = -3, pr_sd = 0.8
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "NatM_p_1_Mal", M_value,
      lo = 0.001, hi = max(2, M_value * 3),
      phase = -1, pr_type = 3, pr_sd = 10
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "L_at_Amin_Mal",
      0, ## L_at_Amin,
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
      ctl$MG_parms, "CV_young_Mal", 0.1,
      lo = 0.001, hi = 5,
      phase = -4, pr_sd = 0.5
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "CV_old_Mal", 0.1,
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
      ctl$SR_parms[ii[1], "LO"] <- 0.0001
      ctl$SR_parms[ii[1], "HI"] <- 20
      ctl$SR_parms[ii[1], "INIT"] <- log(sim_dat$R0)
      ctl$SR_parms[ii[1], "PRIOR"] <- log(sim_dat$R0)
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


  ## ---------------- return ----------------

  ## TODO causes error in fitting (maybe fix later)
  ## ss_dat$Nsexes <- 1

  ## lengthbin 4-38

  inputs$dat <- ss_dat
  inputs$ctl <- ctl
  inputs$fore <- fore

  inputs
}


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
