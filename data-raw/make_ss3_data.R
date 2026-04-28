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
dat$Nfleets <- 2

dat$fleetinfo <- data.frame(
  type = c(1, 3),
  surveytiming = c(-1, 1),
  area = c(1, 1),
  units = c(1, 2),
  need_catch_mult = c(0, 0),
  fleetname = c("FISHERY", "SURVEY1")
)
dat$fleetnames <- c("FISHERY", "SURVEY1")
dat$surveytiming <- c(-1, 1)
dat$units_of_catch <- c(1, 2)
dat$areas <- c(1, 1)

dat$catch <- rbind(
  data.frame(year = -999, seas = 1, fleet = 1, catch = 1e-20, catch_se = 0.05),
  data.frame(year = years, seas = 1, fleet = 1, catch = rep(100, length(years)), catch_se = 0.05)
)

dat$CPUEinfo <- data.frame(
  fleet = 2,
  units = 1,
  errtype = 0,
  SD_report = 1,
  stringsAsFactors = FALSE
)
rownames(dat$CPUEinfo) <- "SURVEY1"

dat$CPUE <- data.frame(
  year = years,
  month = 6,
  index = 2,
  obs = rep(50, length(years)),
  se_log = 0.3
)


## keep length comps on, but keep structure from example
dat$use_lencomp <- 1
dat$len_info <- dat$len_info[1:2, , drop = FALSE]
dat$lencomp <- dat$lencomp[1, , drop = FALSE]
dat$lencomp$year <- min(years)
dat$lencomp$fleet <- 1
dat$lencomp$month <- 7
dat$lencomp$Nsamp <- 50

## ## alternative:
dat$use_lencomp <- 0
dat$len_info <- NULL
dat$lencomp <- NULL


## turn off age-related data
dat$N_agebins <- 0
dat$agebin_vector <- NULL
dat$N_ageerror_definitions <- 0
dat$ageerror <- NULL
dat$age_info <- NULL
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


ctl$Nfleets <- 2
ctl$fleetnames <- c("FISHERY", "SURVEY1")
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
  ctl$Q_options$fleet <- 2
}

if (!is.null(ctl$Q_parms) && nrow(ctl$Q_parms) >= 2) {
  ctl$Q_parms <- ctl$Q_parms[1:2, , drop = FALSE]
}

## keep original selex blocks as much as possible, just trim to two fleets
if (!is.null(ctl$size_selex_types) && nrow(ctl$size_selex_types) >= 2) {
  ctl$size_selex_types <- ctl$size_selex_types[1:2, , drop = FALSE]
}
if (!is.null(ctl$age_selex_types) && nrow(ctl$age_selex_types) >= 2) {
  ctl$age_selex_types <- ctl$age_selex_types[1:2, , drop = FALSE]
}
if (!is.null(ctl$size_selex_parms) && nrow(ctl$size_selex_parms) >= 4) {
  ctl$size_selex_parms <- ctl$size_selex_parms[1:4, , drop = FALSE]
}

grep("age", names(ctl), value = TRUE)

## no age selectivity for either fleet
ctl$age_selex_types <- data.frame(
  Pattern = c(0, 0),
  Discard = c(0, 0),
  Male = c(0, 0),
  Special = c(0, 0)
)
rownames(ctl$age_selex_types) <- c("FISHERY", "SURVEY1")

## remove leftover age selex parameters if present
if ("age_selex_parms" %in% names(ctl)) ctl$age_selex_parms <- NULL
if ("age_selex_parms_reindexed" %in% names(ctl)) ctl$age_selex_parms_reindexed <- NULL



## do not touch lambdas for now
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


## fix issue manually in data

datfile <- file.path(template_dir, "data.ss")
x <- readLines(datfile)

## find CPUE fleet setup header
ind <- grep("^#_fleet[[:space:]]+units[[:space:]]+errtype[[:space:]]+SD_report", x)

stopifnot(length(ind) == 1)

## find the CPUEinfo data row after that header
row_ind <- ind + 1
while (row_ind <= length(x) && grepl("^\\s*#", x[row_ind])) row_ind <- row_ind + 1

## insert terminator right after the CPUEinfo row if missing
next_ind <- row_ind + 1
while (next_ind <= length(x) && trimws(x[next_ind]) == "") next_ind <- next_ind + 1

if (!grepl("^-9999", trimws(x[next_ind]))) {
  x <- append(x, values = "-9999 0 0 0 # terminator", after = row_ind)
}

writeLines(x, datfile)


ctlfile <- file.path(template_dir, "control.ss")
x <- readLines(ctlfile)

ind <- grep("TG_custom", x, fixed = TRUE)
stopifnot(length(ind) == 1)

## keep TG_custom line
## remove any numeric placeholder lines immediately after it
keep <- rep(TRUE, length(x))

for (j in seq(ind + 1, min(length(x), ind + 5))) {
  if (grepl("^-6\\s+6\\s+1\\s+1\\s+2\\s+0\\.01\\s+-4", trimws(x[j]))) {
    keep[j] <- FALSE
  }
}

x <- x[keep]

## ensure the commented conditional placeholder is present right after TG_custom
j <- ind + 1
if (!grepl("^#_Cond\\s+-6\\s+6\\s+1\\s+1\\s+2\\s+0\\.01\\s+-4", trimws(x[j]))) {
  x <- append(
    x,
    values = "#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters",
    after = ind
  )
}

writeLines(x, ctlfile)



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

## ## save the template inputs object as RDS
## saveRDS(
##   inputs,
##   file = file.path(template_dir, "simple_small_inputs.rds")
## )

message("SS3 template written to: ", template_dir)
