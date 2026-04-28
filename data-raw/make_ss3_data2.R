## data-raw/make_ss3_template.R

if (!requireNamespace("r4ss", quietly = TRUE)) {
  stop("Package 'r4ss' is required.")
}

pkg_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
template_dir <- file.path(pkg_root, "inst", "extdata", "ss3", "simple_small")

example_dir <- system.file("extdata", "simple_small", package = "r4ss")
if (example_dir == "") {
  stop("Could not find the example model in the r4ss package.")
}

if (dir.exists(template_dir)) {
  unlink(template_dir, recursive = TRUE, force = TRUE)
}
dir.create(template_dir, recursive = TRUE, showWarnings = FALSE)

r4ss::copy_SS_inputs(dir.old = example_dir, dir.new = template_dir)

inputs <- r4ss::SS_read(dir = template_dir, verbose = FALSE)

years <- 1:30





## ---------------- data ----------------
dat <- inputs$dat

dat$styr <- min(years)
dat$endyr <- max(years)
dat$Nfleets <- 1

dat$Nages <- 32

dat$fleetinfo <- data.frame(
  type = 1,
  surveytiming = -1,
  area = 1,
  units = 1,
  need_catch_mult = 0,
  fleetname = c("Fishery")
)
dat$fleetnames <- c("Fishery")
dat$surveytiming <- -1
dat$units_of_catch <- 1
dat$areas <- 1

dat$catch <- rbind(
  data.frame(year = -999, seas = 1, fleet = 1, catch = 1e-20, catch_se = 0.05),
  data.frame(year = years, seas = 1, fleet = 1, catch = rep(100, length(years)), catch_se = 0.05)
)

dat$CPUEinfo <- data.frame(
  fleet = 1,
  units = 1,
  errtype = 0,
  SD_report = 0,
  stringsAsFactors = FALSE
)
rownames(dat$CPUEinfo) <- "Fishery"

dat$CPUE <- data.frame(
  year = years,
  month = 6,
  index = 1,
  obs = rep(50, length(years)),
  se_log = 0.3
)


## keep length comps on, but keep structure from example
dat$use_lencomp <- 1
dat$len_info <- dat$len_info[1, , drop = FALSE]
rownames(dat$len_info) <- "Fishery"
dat$len_info$addtocomp <- 0.001
dat$len_info$minsamplesize <- 0.001
dat$N_lbins <- 21
dat$lbin_vector <- 12:32
p <- rep(0, dat$N_lbins)
p[10] <- 1

lencomp <- data.frame(
  year = min(years),
  month = 1,
  fleet = 1,
  sex = 1,
  part = 0,
  Nsamp = 10
)

lencomp <- cbind(
  lencomp,
  as.data.frame(t(c(p, rep(0, dat$N_lbins))))
)

colnames(lencomp) <- c(
  "year", "month", "fleet", "sex", "part", "Nsamp",
  paste0("f", dat$lbin_vector),
  paste0("m", dat$lbin_vector)
)

dat$lencomp <- lencomp

## ## alternative:
## dat$use_lencomp <- 0
## dat$len_info <- NULL
## dat$lencomp <- NULL


## turn off age-related data
dat$N_agebins <- 32
dat$agebin_vector <- 0:(dat$N_agebins-1)
dat$N_ageerror_definitions <- 1
dat$ageerror <- rbind(
  rep(-1, dat$N_agebins + 1),
  rep(0.001, dat$N_agebins + 1)
)
colnames(dat$ageerror) <- paste0("age", 0:dat$N_agebins)
dat$age_info <- dat$age_info[1, , drop = FALSE]
rownames(dat$age_info) <- "Fishery"
dat$age_info$addtocomp <- 0.001
dat$age_info$minsamplesize <- 0.001

dat$ageerror <- rbind(
  -1,
  rep(0.001, dat$N_agebins + 1)
)
colnames(dat$ageerror) <- paste0("age", 0:dat$N_agebins)

dat$age_info <- data.frame(
  mintailcomp = -1,
  addtocomp = 0.001,
  combine_M_F = 0,
  CompressBins = 0,
  CompError = 0,
  ParmSelect = 0,
  minsamplesize = 0.001
)
rownames(dat$age_info) <- "Fishery"

dat$lbin_method <- 3

lmax_pop <- ceiling(30 * 1.4)
dat$binwidth <- 1
dat$lbin_vector_pop <- seq(
  from = 0,
  to = lmax_pop,
  by = 1
)
dat$N_lbinspop <- length(dat$lbin_vector_pop)
dat$minimum_size <- 0
dat$maximum_size <- lmax_pop

## no age composition observations
dat$agecomp <- NULL
dat$use_MeanSize_at_Age_obs <- 0
dat$MeanSize_at_Age_obs <- NULL



## turn off tagging completely
dat$do_tags <- 0
dat$N_tag_groups <- 0
dat$N_recap_events <- 0

## remove tag-related tables if present
dat$tagging_data <- NULL
dat$tag_recaps <- NULL
dat$recap <- NULL
dat$TG_release <- NULL
dat$TG_recap_data <- NULL
dat$TG_mixperiod <- NULL
dat$TG_save_for_report <- NULL


inputs$dat <- dat






## ---------------- control ----------------
ctl <- inputs$ctl


ctl$Nfleets <- 1
ctl$fleetnames <- "Fishery"
ctl$MainRdevYrFirst <- min(years)
ctl$MainRdevYrLast <- max(years)
ctl$F_ballpark_year <- min(years)
ctl$TG_custom <- 0
ctl$TG_parms <- NULL
ctl$TG_loss_init <- NULL
ctl$TG_loss_chronic <- NULL
ctl$TG_overdisp <- NULL

ctl$TG_custom <- 0

## if these exist, drop them
for (nm in c("TG_parms", "TG_loss_init", "TG_loss_chronic",
             "TG_overdisp", "TG_Report_fleet", "TG_custom_parms")) {
  if (nm %in% names(ctl)) ctl[[nm]] <- NULL
}

## keep only first survey Q setup if needed
if (!is.null(ctl$Q_options) && nrow(ctl$Q_options) >= 1) {
  ctl$Q_options <- ctl$Q_options[1, , drop = FALSE]
  ctl$Q_options$fleet <- 1
  ctl$Q_options$extra_se <- 0
  ctl$Q_options$float <- 1
  rownames(ctl$Q_options) <- ""
}

if (!is.null(ctl$Q_parms) && nrow(ctl$Q_parms) >= 2) {
  ctl$Q_parms <- ctl$Q_parms[1, , drop = FALSE]
  ctl$Q_parms[1] <- -15
  ctl$Q_parms[2] <- 15
  ctl$Q_parms[3] <- 1
  ctl$Q_parms$PHASE <- -1
  rownames(ctl$Q_parms) <- "LnQ_base_1(1)"
}

## keep original selex blocks as much as possible, just trim to two fleets
if (!is.null(ctl$size_selex_types) && nrow(ctl$size_selex_types) >= 2) {
  ctl$size_selex_types <- ctl$size_selex_types[1, , drop = FALSE]
  ctl$size_selex_types[1] <- 24
  rownames(ctl$size_selex_types) <- "Fishery"
}
if (!is.null(ctl$age_selex_types) && nrow(ctl$age_selex_types) >= 2) {
  ctl$age_selex_types <- ctl$age_selex_types[1, , drop = FALSE]
  rownames(ctl$age_selex_types) <- "Fishery"
}
if (!is.null(ctl$size_selex_parms) && nrow(ctl$size_selex_parms) >= 4) {
  tmp <- ctl$size_selex_parms[1, , drop = FALSE]
  tmp2 <- rbind(tmp, tmp, tmp, tmp, tmp, tmp)
  rownames(tmp2) <- paste0("SizeSel_P_",1:6,"_Fishery(1)")
  tmp2$LO <- c(12,-15,-4,-15,-999,-15)
  tmp2$HI <- c(32,15,12,6,15,20)
  tmp2$INIT <- c(28,15,3.1391,-15,-15,15)
  tmp2$PRIOR <- c(28,15,3.1391,-15,-10,15)
  tmp2$PR_SD <- rep(99, 6)
  tmp2$PR_type <- rep(0, 6)
  tmp2$PHASE <- c(3,-1,3,-1,-2,-1)
  ctl$size_selex_parms <- tmp2
}

## no age selectivity for either fleet
ctl$age_selex_types <- data.frame(
  Pattern = 0,
  Discard = 0,
  Male = 0,
  Special = 0
)
rownames(ctl$age_selex_types) <- "Fishery"

## remove leftover age selex parameters if present
if ("age_selex_parms" %in% names(ctl)) ctl$age_selex_parms <- NULL
if ("age_selex_parms_reindexed" %in% names(ctl)) ctl$age_selex_parms_reindexed <- NULL


## do not touch lambdas for now

ctl$maxlambdaphase <- 15
ctl$sd_offset <- 1
ctl$N_lambdas <- 2

tmp <- ctl$lambdas[1,,drop = FALSE]
tmp2 <- rbind(tmp, tmp)
rownames(tmp2) <- c("catch_Fishery1_Phz1",
                    "init_equ_catch_Fishery_Phz1")

tmp2[,1] <- c(8,9)
tmp2[,2] <- c(1,1)
tmp2[,3] <- c(1,1)
tmp2[,4] <- c(1,0)
tmp2[,5] <- c(1,1)

ctl$lambdas <- tmp2

ctl$more_stddev_reporting <- 0
ctl$stddev_reporting_specs <- NULL
ctl$stddev_reporting_specs <- NULL
ctl$stddev_reporting_selex <- NULL
ctl$stddev_reporting_growth <- NULL
ctl$stddev_reporting_N_at_A <- NULL

inputs$ctl <- ctl






## ---------------- forecast ----------------
fore <- inputs$fore
if ("benchmarks" %in% names(fore)) fore$benchmarks <- 1
if ("MSY" %in% names(fore)) fore$MSY <- 2
if ("Forecast" %in% names(fore)) fore$Forecast <- 2
if ("Nforecastyrs" %in% names(fore)) fore$Nforecastyrs <- 1
if ("FirstYear_for_caps_and_allocations" %in% names(fore)) {
  fore$FirstYear_for_caps_and_allocations <- max(years) + 1
}
if ("Ydecl" %in% names(fore)) fore$Ydecl <- max(years) + 1
if ("Yinit" %in% names(fore)) fore$Yinit <- max(years) + 1

inputs$fore <- fore





## ---------------- write + verify ----------------
r4ss::SS_write(inputs, dir = template_dir, overwrite = TRUE)

## readLines(file.path(template_dir, "data.ss"))
## readLines(file.path(template_dir, "control.ss"))

## verify the template is readable
check_inputs <- r4ss::SS_read(dir = template_dir, verbose = FALSE)
rm(check_inputs)

## verify the template can be fitted
if ( FALSE ) {

  r4ss::run(
    dir = template_dir,
    exe = "~/ss3/bin/ss3",
    extras = "-nohess",
    skipfinished = FALSE,
    verbose = TRUE,
    show_in_console = TRUE
  )

  ## get results
  replist <- try(
    r4ss::SS_output(
      dir = template_dir,
      forecast = TRUE,
      verbose = FALSE,
      printstats = FALSE
    ),
    silent = TRUE
  )

  ## Diagnostics
  r4ss::SS_plots(
    replist = replist,
    plot = 1:25,
    print = TRUE,
    dir = template_dir
  )

  dir(template_dir)
  unlink(list.files(template_dir, full.names = TRUE), recursive = TRUE)

}


message("SS3 template written to: ", template_dir)
