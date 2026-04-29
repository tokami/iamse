

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

  ## TODO: 60 years of catch
  ## TODO: 10 years of index
  ## TODO: fixed uncertainty in ss3 for parameters
  ## TODO: catch starts at 0

  ## TODO: differences (dat)
  ## #_Nfleets 1

  ## #_fleetinfo
 ## [16] "#_type\tsurveytiming\tarea\tunits\tneed_catch_mult\tfleetname"
  ## [17] "1\t-1\t1\t1\t0\tFishery\t#_1"

  ## #_CPUE and surveyabundance_observations
 ##   [96] "#_fleet\tunits\terrtype\tSD_report"
  ## [97] "1\t1\t0\t0\t#_Fishery"

  ## lengthbin 4-38


## [140] "1 #_use_lencomp"
## [141] "#"
## [142] "#_len_info"
## [143] "#_mintailcomp\taddtocomp\tcombine_M_F\tCompressBins\tCompError\tParmSelect\tminsamplesize"
## [144] "-1\t0.001\t0\t0\t0\t0\t0.001\t#_Fishery"
## [145] "21 #_N_lbins"
## [146] "#_lbin_vector"
## [147] "12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 #_lbin_vector"
## [148] "#"
## [149] "#_lencomp"
## [150] "#_year\tmonth\tfleet\tsex\tpart\tNsamp\tf12\tf13\tf14\tf15\tf16\tf17\tf18\tf19\tf20\tf21\tf22\tf23\tf24\tf25\tf26\tf27\tf28\tf29\tf30\tf31\tf32\tm12\tm13\tm14\tm15\tm16\tm17\tm18\tm19\tm20\tm21\tm22\tm23\tm24\tm25\tm26\tm27\tm28\tm29\tm30\tm31\tm32"
## [151] " 1974\t1\t1\t1\t0\t64.232884\t 9\t100\t264\t690\t1297\t1752\t2439\t2585\t2113\t1759\t1441\t1122\t1174\t1282\t1016\t809\t514\t339\t167\t54\t47\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t#_2         "
## [152] " 1976\t1\t1\t1\t0\t23.269988\t39\t 44\t 35\t 56\t 194\t 443\t 943\t1013\t1068\t 959\t 761\t 560\t 494\t 365\t 272\t174\t 91\t 41\t 27\t14\t 5\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t#_3         "
## [153] " 1977\t1\t1\t1\t0\t54.989817\t 0\t  8\t 18\t 81\t 202\t 687\t2148\t3040\t2974\t2334\t1841\t1293\t1132\t 733\t 551\t497\t228\t122\t 41\t18\t 7\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t#_4         "
## [154] " 1980\t1\t1\t1\t0\t22.495138\t 3\t  3\t 13\t 76\t 248\t 629\t 810\t 901\t 869\t 755\t 777\t 677\t 520\t 403\t 273\t169\t104\t 59\t 38\t 7\t11\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t#_5         "
## [155] " 1981\t1\t1\t1\t0\t35.012174\t 0\t  5\t 35\t157\t 466\t 924\t1144\t1286\t1511\t1299\t1152\t 886\t 686\t 618\t 466\t335\t207\t121\t 76\t42\t16\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t#_6         "
  ## [156] "-9999\t0\t0\t0\t0\t        0\t 0\t  0\t  0\t  0\t   0\t   0\t   0\t   0\t   0\t   0\t   0\t   0\t   0\t   0\t   0\t  0\t  0\t  0\t  0\t 0\t 0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t#_terminator"

## [157] "32 #_N_agebins"
## [158] "#"
## [159] "#_agebin_vector"
## [160] "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 #_agebin_vector"
## [161] "#"
## [162] "#_ageing_error"
## [163] "1 #_N_ageerror_definitions"
## [164] "#_age0\tage1\tage2\tage3\tage4\tage5\tage6\tage7\tage8\tage9\tage10\tage11\tage12\tage13\tage14\tage15\tage16\tage17\tage18\tage19\tage20\tage21\tage22\tage23\tage24\tage25\tage26\tage27\tage28\tage29\tage30\tage31\tage32"
## [165] "   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t   -1\t#_1"
## [166] "0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t0.001\t#_2"
## [167] "#"
## [168] "#_age_info"
## [169] "#_mintailcomp\taddtocomp\tcombine_M_F\tCompressBins\tCompError\tParmSelect\tminsamplesize"
## [170] "-1\t0.001\t0\t0\t0\t0\t0.001\t#_Fishery"
## [171] "3 #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths"
## [172] " #_combine males into females at or below this bin number"
## [173] "#_X.9999\tX0\tX0.1\tX0.2\tX0.3\tX0.4\tX0.5\tX0.6\tX0.7\tX0.8\tX0.9\tX0.10\tX0.11\tX0.12\tX0.13\tX0.14\tX0.15\tX0.16\tX0.17\tX0.18\tX0.19\tX0.20\tX0.21\tX0.22\tX0.23\tX0.24\tX0.25\tX0.26\tX0.27\tX0.28\tX0.29\tX0.30\tX0.31\tX0.32\tX0.33\tX0.34\tX0.35\tX0.36\tX0.37\tX0.38\tX0.39\tX0.40\tX0.41\tX0.42\tX0.43\tX0.44\tX0.45\tX0.46\tX0.47\tX0.48\tX0.49\tX0.50\tX0.51\tX0.52\tX0.53\tX0.54\tX0.55\tX0.56\tX0.57\tX0.58\tX0.59\tX0.60\tX0.61\tX0.62\tX0.63\tX0.64\tX0.65\tX0.66\tX0.67\tX0.68\tX0.69\tX0.70\tX0.71"
## [174]

  ss_dat <- inputs$dat
  ctl <- inputs$ctl
  fore <- inputs$fore


  ## ---------------- helper ----------------
  set_par <- function(tab, pattern, init, lo = NULL, hi = NULL,
                      prior = init, phase = -99) {
    ii <- grep(pattern, rownames(tab))
    if (length(ii) == 0 || !is.finite(init)) return(tab)

    if (!is.null(lo)) tab[ii, "LO"] <- lo
    if (!is.null(hi)) tab[ii, "HI"] <- hi

    tab[ii, "INIT"] <- init
    tab[ii, "PRIOR"] <- prior
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
    row_ind <- sort(sample(1:nrow(obsL), ny_lencomp))

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

  ss_dat$N_agebins <- 0
  ss_dat$agebin_vector <- NULL
  ss_dat$N_ageerror_definitions <- 0
  ss_dat$ageerror <- NULL
  ss_dat$age_info <- NULL
  ss_dat$agecomp <- NULL
  ss_dat$use_MeanSize_at_Age_obs <- 0
  ss_dat$MeanSize_at_Age_obs <- NULL

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
      ctl$MG_parms, "NatM_uniform_Fem", M_value,
      lo = 0.001, hi = max(1, M_value * 3),
      phase = -99
    )

    min_pop_len <- min(ss_dat$lbin_vector_pop)
    max_pop_len <- max(ss_dat$lbin_vector_pop)

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "L_at_Amin_Fem",
      L_at_Amin,
      lo = min_pop_len,
      hi = max_pop_len,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "L_at_Amax_Fem",
      sim_dat$Linf,
      lo = min_pop_len,
      hi = max_pop_len,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "VonBert_K_Fem", sim_dat$K,
      lo = 0.001, hi = max(1, sim_dat$K * 5),
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Wtlen_1_Fem", sim_dat$lwa,
      lo = sim_dat$lwa / 100,
      hi = sim_dat$lwa * 100,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Wtlen_2_Fem", sim_dat$lwb,
      lo = 1,
      hi = 5,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Wtlen_1_Mal", sim_dat$lwa,
      lo = sim_dat$lwa / 100,
      hi = sim_dat$lwa * 100,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Wtlen_2_Mal", sim_dat$lwb,
      lo = 1,
      hi = 5,
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Mat50", sim_dat$Lm50,
      lo = 0.001,
      hi = max(100, sim_dat$Lm50 * 3),
      phase = -99
    )

    ctl$MG_parms <- set_par(
      ctl$MG_parms, "Mat_slope", mat_slope,
      lo = -10,
      hi = -0.001,
      phase = -99
    )
  }


  ctl$First_Mature_Age <- 0
  ctl$fecundity_option <- 2

 ## [48] "#_growth_parms"
 ## [49] "#_LO\tHI\tINIT\tPRIOR\tPR_SD\tPR_type\tPHASE\tenv_var&link\tdev_link\tdev_minyr\tdev_maxyr\tdev_PH\tBlock\tBlock_Fxn"
 ## [50] "0.001\t   2\t     0.17\t    -2.92\t0.22\t3\t -1\t0\t0\t0\t0\t0\t0\t0\t#_NatM_p_1_Fem_GP_1  "
 ## [51] "  -50\t 100\t        0\t        0\t  10\t0\t -3\t0\t0\t0\t0\t0\t0\t0\t#_L_at_Amin_Fem_GP_1 "
 ## [52] "    1\t 500\t    28.97\t    28.97\t  10\t0\t -2\t0\t0\t0\t0\t0\t0\t0\t#_L_at_Amax_Fem_GP_1 "
 ## [53] "0.001\t   2\t     0.11\t     0.11\t0.05\t0\t -3\t0\t0\t0\t0\t0\t0\t0\t#_VonBert_K_Fem_GP_1 "
 ## [54] "0.001\t   5\t      0.1\t      0.1\t 0.5\t0\t -4\t0\t0\t0\t0\t0\t0\t0\t#_CV_young_Fem_GP_1  "
 ## [55] "0.001\t   5\t      0.1\t      0.1\t 0.5\t0\t -4\t0\t0\t0\t0\t0\t0\t0\t#_CV_old_Fem_GP_1    "
 ## [56] "    0\t   3\t  4.6e-05\t  4.6e-05\t  99\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_Wtlen_1_Fem_GP_1   "
 ## [57] "    2\t   4\t  2.94493\t  2.94493\t  99\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_Wtlen_2_Fem_GP_1   "
 ## [58] "1e-04\t1000\t     22.1\t     22.1\t  99\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_Mat50%_Fem_GP_1    "
 ## [59] "   -2\t   4\t-0.640095\t-0.640095\t  99\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_Mat_slope_Fem_GP_1 "
 ## [60] "    0\t   3\t  4.6e-05\t  4.6e-05\t 0.8\t0\t -3\t0\t0\t0\t0\t0\t0\t0\t#_Eggs_alpha_Fem_GP_1"
 ## [61] "    0\t  10\t  2.94493\t  2.94493\t 0.8\t0\t -3\t0\t0\t0\t0\t0\t0\t0\t#_Eggs_beta_Fem_GP_1 "
 ## [62] "0.001\t   2\t     0.17\t    -2.92\t0.22\t3\t -1\t0\t0\t0\t0\t0\t0\t0\t#_NatM_p_1_Mal_GP_1  "
 ## [63] "  -50\t 100\t        0\t        0\t  10\t0\t -3\t0\t0\t0\t0\t0\t0\t0\t#_L_at_Amin_Mal_GP_1 "
 ## [64] "    1\t 500\t    28.97\t    28.97\t  10\t0\t -2\t0\t0\t0\t0\t0\t0\t0\t#_L_at_Amax_Mal_GP_1 "
 ## [65] "0.001\t   2\t     0.11\t     0.11\t0.05\t0\t -3\t0\t0\t0\t0\t0\t0\t0\t#_VonBert_K_Mal_GP_1 "
 ## [66] "0.001\t   5\t      0.1\t      0.1\t 0.5\t0\t -4\t0\t0\t0\t0\t0\t0\t0\t#_CV_young_Mal_GP_1  "
 ## [67] "0.001\t   5\t      0.1\t      0.1\t 0.5\t0\t -4\t0\t0\t0\t0\t0\t0\t0\t#_CV_old_Mal_GP_1    "
 ## [68] "    0\t   3\t  4.6e-05\t  4.6e-05\t  99\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_Wtlen_1_Mal_GP_1   "
 ## [69] "    2\t   4\t  2.94493\t  2.94493\t  99\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_Wtlen_2_Mal_GP_1   "
 ## [70] "  0.1\t  10\t        1\t        1\t   1\t0\t -1\t0\t0\t0\t0\t0\t0\t0\t#_CohortGrowDev      "
 ##  [71] " 0.01\t0.99\t      0.5\t      0.5\t 0.5\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_FracFemale_GP_1    "


 ## [82] "#_LO\tHI\tINIT\tPRIOR\tPR_SD\tPR_type\tPHASE\tenv-var\tuse_dev\tdev_mnyr\tdev_mxyr\tdev_PH\tBlock\tBlk_Fxn # parm_name"
 ## [83] "1e-04\t20\t  7\t  7\t  99\t0\t  1\t0\t0\t0\t0\t0\t0\t0\t#_SR_LN(R0)  "
 ## [84] "  0.2\t 1\t0.7\t0.7\t0.24\t3\t -1\t0\t0\t0\t0\t0\t0\t0\t#_SR_BH_steep"
 ## [85] "    0\t 2\t0.6\t0.6\t  99\t0\t -6\t0\t0\t0\t0\t0\t0\t0\t#_SR_sigmaR  "
 ## [86] "   -5\t 5\t  0\t  0\t  99\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_SR_regime  "
 ## [87] "    0\t 2\t  0\t  1\t  99\t0\t-99\t0\t0\t0\t0\t0\t0\t0\t#_SR_autocorr"


  ## ---------------- stock-recruitment ----------------

  ctl$SR_function <- 3
  ctl$MainRdevYrFirst <- ss_dat$styr
  ctl$MainRdevYrLast <- ss_dat$endyr
  ctl$F_ballpark_year <- ss_dat$styr


  if (!is.null(ctl$SR_parms)) {
    rn <- rownames(ctl$SR_parms)

    ii <- grep("R0|LnR0|log_R0|SR_LN\\(R0\\)", rn)
    if (length(ii) > 0 && !is.null(sim_dat$R0)) {
      ctl$SR_parms[ii[1], "INIT"] <- log(sim_dat$R0)
    }

    ii <- grep("steep|Steep", rn)
    if (length(ii) > 0 && !is.null(sim_dat$h)) {
      ctl$SR_parms[ii[1], "INIT"] <- sim_dat$h
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

## [124] "#_Q_setup for fleets with cpue or survey data"
## [125] "#_fleet\tlink\tlink_info\textra_se\tbiasadj\tfloat  #  fleetname"
## [126] "    1\t1\t0\t0\t0\t1\t#_1         "
## [127] "-9999\t0\t0\t0\t0\t0\t#_terminator"
## [128] "#_Q_parms(if_any);Qunits_are_ln(q)"
## [129] "#_LO\tHI\tINIT\tPRIOR\tPR_SD\tPR_type\tPHASE\tenv-var\tuse_dev\tdev_mnyr\tdev_mxyr\tdev_PH\tBlock\tBlk_Fxn  #  parm_name"
## [130] "-15\t15\t1\t0\t1\t0\t-1\t0\t0\t0\t0\t0\t0\t0\t#_LnQ_base_1(1)"
## [131] "#_no timevary Q parameters"


## [133] "#_size_selex_patterns"
## [134] "#_Pattern\tDiscard\tMale\tSpecial"
## [135] "24\t0\t0\t0\t#_1 Fishery"
## [136] "#"
## [137] "#_age_selex_patterns"
## [138] "#_Pattern\tDiscard\tMale\tSpecial"
## [139] "0\t0\t0\t0\t#_1 Fishery"
## [140] "#"
## [141] "#_SizeSelex"
## [142] "#_LO\tHI\tINIT\tPRIOR\tPR_SD\tPR_type\tPHASE\tenv-var\tuse_dev\tdev_mnyr\tdev_mxyr\tdev_PH\tBlock\tBlk_Fxn  #  parm_name"
## [143] "  12\t32\t    28\t    28\t99\t0\t 3\t0\t0\t0\t0\t0\t0\t0\t#_SizeSel_P_1_Fishery(1)"
## [144] " -15\t15\t    15\t    15\t99\t0\t-1\t0\t0\t0\t0\t0\t0\t0\t#_SizeSel_P_2_Fishery(1)"
## [145] "  -4\t12\t3.1391\t3.1391\t99\t0\t 3\t0\t0\t0\t0\t0\t0\t0\t#_SizeSel_P_3_Fishery(1)"
## [146] " -15\t 6\t   -15\t   -15\t99\t0\t-1\t0\t0\t0\t0\t0\t0\t0\t#_SizeSel_P_4_Fishery(1)"
## [147] "-999\t15\t   -15\t   -10\t99\t0\t-2\t0\t0\t0\t0\t0\t0\t0\t#_SizeSel_P_5_Fishery(1)"
## [148] " -15\t20\t    15\t    15\t99\t0\t-1\t0\t0\t0\t0\t0\t0\t0\t#_SizeSel_P_6_Fishery(1)"

  ctl$maxlambdaphase <- 15

## [165] "# read 2 changes to default Lambdas (default value is 1.0)"
## [166] "#_like_comp\tfleet\tphase\tvalue\tsizefreq_method"
## [167] "    8\t1\t1\t1\t1\t#_catch_Fishery1_Phz1        "
## [168] "    9\t1\t1\t0\t1\t#_init_equ_catch_Fishery_Phz1"
## [169] "-9999\t0\t0\t0\t0\t#_terminator                 "


  ctl$more_stddev_reporting <- 0



  ## ---------------- selectivity ----------------

  ## TODO
  ## if (isTRUE(fix_selectivity) && !is.null(ctl$size_selex_parms)) {

  ##   sel_width <- sim_dat$Ls95 - sim_dat$Ls50

  ##   ctl$size_selex_parms <- set_par(
  ##     ctl$size_selex_parms, "SizeSel_P_1_FISHERY",
  ##     sim_dat$Ls50,
  ##     lo = 0.001,
  ##     hi = max(100, sim_dat$Ls50 * 3),
  ##     phase = -99
  ##   )

  ##   ctl$size_selex_parms <- set_par(
  ##     ctl$size_selex_parms, "SizeSel_P_2_FISHERY",
  ##     sel_width,
  ##     lo = 0.001,
  ##     hi = max(100, sel_width * 5),
  ##     phase = -99
  ##   )

  ##   ctl$size_selex_parms <- set_par(
  ##     ctl$size_selex_parms, "SizeSel_P_1_SURVEY1",
  ##     sim_dat$Ls50,
  ##     lo = 0.001,
  ##     hi = max(100, sim_dat$Ls50 * 3),
  ##     phase = -99
  ##   )

  ##   ctl$size_selex_parms <- set_par(
  ##     ctl$size_selex_parms, "SizeSel_P_2_SURVEY1",
  ##     sel_width,
  ##     lo = 0.001,
  ##     hi = max(100, sel_width * 5),
  ##     phase = -99
  ##   )

  ## }

  ## Remove stale age selectivity parameters, if present
  if (!is.null(ctl$age_selex_parms)) {
    ctl$age_selex_parms <- NULL
  }

  ## ---------------- survey Q ----------------

  if (!is.null(ctl$Q_parms)) {

    iq <- grep("LnQ", rownames(ctl$Q_parms))

    ctl$Q_parms[iq, "INIT"] <- log(sim_dat$q)
    ctl$Q_parms[iq, "LO"] <- log(sim_dat$q / 100)
    ctl$Q_parms[iq, "HI"] <- log(sim_dat$q * 100)
    ctl$Q_parms[iq, "PHASE"] <- 1

  }

  ## ---------------- forecast ----------------

  fore$FirstYear_for_caps_and_allocations <- ss_dat$endyr + 1
  fore$Ydecl <- ss_dat$endyr + 1
  fore$Yinit <- ss_dat$endyr + 1


  ## ---------------- return ----------------

  ## TODO causes error in fitting (maybe fix later)
  ## ss_dat$Nsexes <- 1

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
