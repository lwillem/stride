#!/usr/bin/env Rscript
#############################################################################
#  This file is part of the Stride software.
#  It is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or any
#  later version.
#  The software is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License,
#  along with the software. If not, see <http://www.gnu.org/licenses/>.
#  see http://www.gnu.org/licenses/.
#
#  Copyright 2026, Manansala R, Willem L
#############################################################################
#
# STEP-BY-STEP CALIBRATION TO THE SPARTANBURG (SC) 2025-26 MEASLES OUTBREAK
#
# Rerun this after the population, contact matrix or immunity files change.
# Documentation, targets and sources: doc/markdown/measles_sc_calibration.md
#
#   Step 1  population: renumber household ids 1..N ("hhdense")
#   Step 2  R0 -> transmission probability (new disease file)
#   Step 3  adult susceptibility factor  (target: adult share of cases)
#   Step 4  school-cluster rate x tracing efficiency "other"
#           (targets: cumulative cases by report date, share of cases traced
#            as contacts while infected = "developed illness in quarantine")
#   Step 5  check run with many seeds: model vs all targets
#
# Every step saves its result in `calib_dir` (calibration_state.rds) and the
# next step reads it, so the steps can be run in separate R sessions.
#
# Call this script from the main project folder (containing bin, config, data, ...)
# E.g.: path/to/stride $ ./bin/rStride_measles_calibrate_sc.R
# Set `steps` below first.
#############################################################################

rm(list = ls())

################################################ #
## 0. SETTINGS                               ####
################################################ #

steps                <- 1:5                    # step(s) to run now, e.g. 1, 2, c(3,4), 5
calib_dir            <- "sim_output/calibration_SC33"
num_parallel_workers <- 3

## ---- input files (the ones that may be updated) ----
pop_file_in      <- "data/population_6region/pop_usa_sc_spartanburg_c1000.csv"     # from social_contacts_usa2026.R
contact_file     <- "data/contacts_usa/contact_matrix_usa_sc_spartanburg_c1000_v2.xml"
holidays_file    <- "data/calendar_USA_2023_2026_measles.csv"
disease_template <- "data/disease_measles_usa.xml"                                  # natural history; b0/b1 are refitted in step 2
adult_raw_file   <- "data/immunity_6region/immunity_measles_adult_SC.xml"           # natural immunity, ages 18+
child_file       <- "data/immunity_6region/immunity_measles_child_SC.xml"           # vaccine coverage, ages 0-17

## ---- fixed inputs (sources: table 2 of the documentation) ----
fixed <- list(
  r0                           = 12,            # DPH: "12 to 18 non-immune contacts"; lower end
  start_date                   = "2025-09-15",  # assumption: ~1 week before the first outbreak cases
  num_days                     = 240,
  num_infected_seeds           = 1,
  seeding_age_min              = 1,
  seeding_age_max              = 17,
  contact_tracing_date         = "2025-10-02",  # outbreak declared (DPH)
  vaccine_hesitancy_rate       = 0.14,          # ASPE/CDC q9mh-h2tw, Spartanburg 0.1443
  vaccine_hesitancy_uptake     = 0.5,           # assumption: county child coverage ~88-90%
  vaccine_hesitancy_cluster_pooltype = "School",
  vaccine_hesitancy_num_clusters     = 33,       # 33 schools impacted with exposures from SC DPH
  vaccine_hesitancy_cluster_uptake   = 0,       # assumption
  detection_probability        = 0.9,           # assumption
  tracing_efficiency_household = 0.95,          # assumption
  delay_isolation_index        = 3,             # assumption
  delay_contact_tracing        = 1,             # assumption
  case_finding_capacity        = 10000,         # keep < 1e5: R writes 100000 as "1e+05", C++ reads 1
  test_false_negative          = 0,
  isolation_duration_index     = 7,             # model choice (CDC: isolate through 4 days after rash)
  isolation_duration_contact   = 21             # DPH/CDC: 21-day quarantine
)

## ---- start values of the calibrated parameters (October 2026 calibration) ----
start_values <- list(adult_factor = 0.35, cluster_rate = 0.45, tracing_other = 0.4)

## ---- targets (SC DPH; sources: table 1 of the documentation) ----
# cumulative outbreak cases by DPH REPORT date
obs <- data.frame(date = as.Date(c(
  "2025-10-02","2025-10-07","2025-10-09","2025-10-14","2025-10-17","2025-10-21","2025-10-24","2025-10-28",
  "2025-10-31","2025-11-07","2025-11-12","2025-11-14","2025-11-18","2025-11-21","2025-11-25","2025-11-28",
  "2025-12-02","2025-12-05","2025-12-09","2025-12-12","2025-12-16","2025-12-19","2025-12-23","2025-12-26",
  "2025-12-30","2026-01-02","2026-01-06","2026-01-09","2026-01-13","2026-01-16","2026-01-20","2026-01-23",
  "2026-01-27","2026-01-30","2026-02-03","2026-02-06","2026-02-10","2026-02-13","2026-02-17","2026-02-20",
  "2026-02-24","2026-02-27","2026-03-03","2026-03-06","2026-03-10","2026-03-13","2026-03-17")),
  cum_cases = c(5,7,8,13,16,20,22,30,34,35,43,44,49,52,55,62,76,84,111,126,135,144,153,156,
                176,185,211,310,434,558,646,700,789,847,876,920,933,950,962,973,979,985,990,991,
                993,996,997))
target_final_size          <- 997          # DPH, outbreak declared over 2026-04-27
target_adult_share         <- 87 / 997     # 87 adults of 997 cases with known age (DPH data via WSPA)
target_share_in_quarantine <- 496 / 989    # 496 of 2,219 quarantined fell ill in quarantine; 989 cases (DPH, 2026-03-02)
target_child_coverage      <- 0.90         # Spartanburg K-12 "just under 90%" (DPH school data, 2024-25)

## ---- how model output is compared with the targets ----
report_lag        <- 3      # days from symptom onset to DPH report (assumption)
min_outbreak_size <- 200    # runs below this did not take off and are excluded
min_p_outbreak    <- 0.5    # settings where fewer runs take off are penalised in the score (the observed outbreak should not be a rare event)
min_outbreak_runs <- 10     # settings with fewer outbreak runs are not scored (medians too noisy)
num_obs_cases     <- 997    # adult share is measured over the first 997 cases

## ---- scan grids and seeds ----
grid_r0             <- seq(0, 15, 3)
num_seeds_r0        <- 10
# grids widened 2026-11-18 after the 7- and 33-cluster runs (values at the grid edges)
grid_adult_factor   <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)          # chosen: 0.45 (7 clusters), 0.60 (33)
grid_cluster_rate   <- if (fixed$vaccine_hesitancy_num_clusters >= 20) {
                         c(0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40)          # 33 clusters: 0.30 was the top of the old grid
                       } else {
                         c(0.20, 0.25, 0.30, 0.35, 0.45)          # 7 clusters: 0.25 chosen; < 0.20 gave no outbreaks
                       }
# previous grids: adult c(0.4, 0.5, 0.6, 0.7, 0.8, 1.0); cluster c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
grid_tracing_other  <- c(0.15, 0.20, 0.25, 0.30, 0.40)          # 0.3 was the bottom of the old grid c(0.3, 0.4, 0.5)
num_seeds_scan      <- 40
num_seeds_check     <- 100

################################################ #
## FUNCTIONS                                 ####
################################################ #

state_file <- file.path(calib_dir, "calibration_state.rds")
load_state <- function() if (file.exists(state_file)) readRDS(state_file) else start_values
save_state <- function(st) { dir.create(calib_dir, showWarnings = FALSE, recursive = TRUE); saveRDS(st, state_file) }

## ---- population ----
# household ids 1..N in the order of the original ids (same results, ~200x faster)
make_hhdense <- function(pop_in) {
  pop_out <- sub("\\.csv$", "_hhdense.csv", pop_in)
  map_out <- sub("\\.csv$", "_hhdense_household_map.csv", pop_in)
  pop <- read.csv(pop_in, stringsAsFactors = FALSE)
  ids <- sort(unique(pop$household_id))
  pop$household_id <- match(pop$household_id, ids)
  write.csv(pop, pop_out, row.names = FALSE, quote = FALSE, na = "NA")
  write.csv(data.frame(household_id = seq_along(ids), original_household_id = ids),
            map_out, row.names = FALSE, quote = FALSE)
  cat(sprintf("[step 1] %d persons, %d households (max original id %d) -> %s\n",
              nrow(pop), length(ids), max(ids), pop_out))
  pop_out
}

## ---- immunity files ----
read_immunity_xml <- function(file) {
  txt  <- readLines(file, warn = FALSE)
  m    <- regmatches(txt, regexec("<age([0-9]+)>([-0-9.eE]+)</age[0-9]+>", txt))
  m    <- m[lengths(m) == 3]
  get_tag <- function(tag) {
    x <- regmatches(txt, regexec(sprintf("<%s>(.*)</%s>", tag, tag), txt)); x <- x[lengths(x) == 2]
    if (length(x) == 0) NA_character_ else x[[1]][2]
  }
  list(header = txt[grepl("^<!--", txt)][1], data_source = get_tag("data_source"),
       data_manipulation = get_tag("data_manipulation"),
       age = as.integer(vapply(m, `[`, "", 2)), immunity = as.numeric(vapply(m, `[`, "", 3)))
}
# immunity_new = 1 - f * (1 - immunity_raw) for ages >= 18
write_adult_factor_file <- function(raw_file, f, note = "") {
  out <- sub("\\.xml$", sprintf("_calib_f%03d.xml", as.integer(round(100 * f))), raw_file)
  x   <- read_immunity_xml(raw_file)
  imm <- ifelse(x$age >= 18, 1 - f * (1 - x$immunity), x$immunity)
  writeLines(c(if (!is.na(x$header)) x$header, "<immunity>",
               sprintf("  <data_source>%s</data_source>", x$data_source),
               sprintf("  <data_manipulation>%s; susceptibility (1 - immunity) at ages 18+ multiplied by %.2f%s</data_manipulation>",
                       x$data_manipulation, f, note),
               sprintf("  <age%d>%.4f</age%d>", x$age, imm, x$age), "</immunity>"), out)
  out
}

## ---- STRIDE parameters for the SC outbreak scenario ----
get_sc_param_list <- function(st, num_seeds) {
  p <- get_default_param()
  for (k in names(fixed)) p[[k]] <- fixed[[k]]
  p$num_rng_seeds              <- num_seeds
  p$disease_config_file        <- st$disease_file
  p$population_file            <- st$pop_file
  p$age_contact_matrix_file    <- contact_file
  p$holidays_file              <- holidays_file
  p$contact_probability_rule   <- "Mean"
  p$event_log_level            <- "Transmissions"     # [PRIM]/[TRAN]/[TRACE] lines for the metrics
  p$output_pop_snapshot        <- TRUE                # writes susceptibles_by_age.csv (coverage check)
  p$reference_hospital_data_file <- NA
  p$reference_serology_data_file <- NA
  p$immunity_profile           <- "AgeDependent"
  p$immunity_distribution_file <- st$adult_file
  p$immunity_link_probability  <- 0
  p$contact_survey_dates       <- p$start_date
  p$vaccine_profile            <- "AgeDependent"
  p$vaccine_distribution_file  <- child_file
  p$vaccine_min_age            <- 0
  p$vaccine_max_age            <- 17
  p$mass_immunize              <- TRUE
  p$vaccine_hesitancy_cluster_rate <- st$cluster_rate
  p$tracing_efficiency_other   <- st$tracing_other
  p$hosp_probability_factor    <- 1
  p$hospital_length_of_stay    <- 14
  p$hospital_probability_age   <- paste(0.18,0.06,0.12,sep=',')
  p$hospital_category_age      <- paste(0,5,20,sep=',')
  p$hospital_mean_delay_age    <- paste(2,2,5,sep=',')
  p
}
run_design <- function(p, postfix) {
  big <- vapply(p, function(v) is.numeric(v) && any(abs(v) >= 1e5, na.rm = TRUE), TRUE)
  if (any(big)) stop("numeric parameter(s) >= 1e5 are written as 1e+05: ", paste(names(p)[big], collapse = ", "))
  d <- .rstride$get_full_grid_exp_design(exp_param_list = p, num_rng_seeds = p$num_rng_seeds)
  run_rStride(d, postfix, remove_run_output = FALSE, num_parallel_workers = num_parallel_workers)
}

## ---- model output per run ----
# event log fields after the tag ([PRIM]/[TRAN]): 1 local_id, 2 infector_id, 3 part_age,
# 4 infector_age, 5 pool_type, 6 sim_day, 7 id_index_case, 8 start_infectiousness,
# 9 end_infectiousness, 10 start_symptoms, ... (8-10 are in days since infection)
# [TRACE]: 1 person_id, 2 age, 3 is_infected, 4 is_symptomatic, 5 pool ("Index" = detected case), 6 case_id, 7 case_age, 8 sim_day
read_run_log <- function(event_log) {
  txt <- readLines(event_log, warn = FALSE)
  tr  <- txt[startsWith(txt, "[PRIM]") | startsWith(txt, "[TRAN]")]
  p   <- strsplit(tr, " ", fixed = TRUE)
  cases <- data.frame(id = as.integer(vapply(p, `[`, "", 2)), age = as.integer(vapply(p, `[`, "", 4)),
                      sim_day = as.integer(vapply(p, `[`, "", 7)), onset = as.integer(vapply(p, `[`, "", 11)))
  tc  <- strsplit(txt[startsWith(txt, "[TRACE]")], " ", fixed = TRUE)
  traced_contacts <- unique(as.integer(vapply(tc, `[`, "", 2))[vapply(tc, `[`, "", 6) != "Index" &
                                                                vapply(tc, `[`, "", 4) == "true"])
  list(cases = cases, traced_contacts = traced_contacts)
}
child_coverage <- function(run_dir, ages = 5:17) {
  f <- file.path(run_dir, "susceptibles_by_age.csv")
  if (!file.exists(f)) return(NA_real_)
  s <- read.csv(f); s <- s[s$age %in% ages, ]
  sum(s$immune) / sum(s$immune + s$susceptible)
}
run_metrics <- function(project_dir, stride_dir = ".") {
  summ <- read.csv(file.path(stride_dir, project_dir, paste0(basename(project_dir), "_summary.csv")),
                   stringsAsFactors = FALSE)
  col <- function(k, i) if (k %in% names(summ)) summ[[k]][i] else NA
  out <- lapply(seq_len(nrow(summ)), function(i) {
    run_dir <- file.path(stride_dir, summ$output_prefix[i])
    if (!file.exists(file.path(run_dir, "event_log.txt"))) return(NULL)
    lg <- read_run_log(file.path(run_dir, "event_log.txt"))
    cs <- lg$cases[order(lg$cases$sim_day), ]
    first <- head(cs, num_obs_cases)
    report_date <- as.Date(summ$start_date[i]) + cs$sim_day + cs$onset + report_lag   # onset = days since infection
    cum <- vapply(obs$date, function(d) sum(report_date <= d), 0)
    data.frame(project = basename(project_dir), run = basename(summ$output_prefix[i]),
               adult_file = col("immunity_distribution_file", i),
               cluster_rate = col("vaccine_hesitancy_cluster_rate", i),
               tracing_other = col("tracing_efficiency_other", i),
               num_cases = nrow(cs),
               adult_share_first = mean(first$age >= 18),
               share_traced = mean(cs$id %in% lg$traced_contacts),
               share_traced_first = mean(first$id %in% lg$traced_contacts),
               coverage_5_17 = child_coverage(run_dir),
               t(setNames(cum, paste0("cum_", format(obs$date, "%Y%m%d")))),
               check.names = FALSE)
  })
  do.call(rbind, out)
}
# median over outbreak runs, per parameter combination
summarise_runs <- function(m, by = c("adult_file", "cluster_rate", "tracing_other")) {
  cum_cols <- grep("^cum_", names(m), value = TRUE)
  key <- do.call(paste, c(m[, by, drop = FALSE], sep = " | "))
  out <- do.call(rbind, lapply(split(m, key), function(d) {
    ok <- d[d$num_cases >= min_outbreak_size, ]
    r  <- d[1, by, drop = FALSE]
    r$runs <- nrow(d); r$n_outbreak <- nrow(ok); r$p_outbreak <- nrow(ok) / nrow(d)
    r$coverage_5_17 <- median(d$coverage_5_17, na.rm = TRUE)
    for (k in c("num_cases", "adult_share_first", "share_traced", "share_traced_first", cum_cols))
      r[[k]] <- if (nrow(ok) > 0) median(ok[[k]]) else NA
    r
  }))
  rownames(out) <- NULL
  out
}
# fit score (lower = better), sum of four terms:
#   score_curve    mean squared log error of the cumulative curve over the DPH report dates
#   score_traced   squared log error of the share of cases traced as contacts
#   score_final    squared log error of the final size (median of outbreak runs, day num_days) vs 997;
#                  penalises outbreaks that keep growing after the last DPH report
#   score_outbreak squared log of min_p_outbreak / p_outbreak if p_outbreak < min_p_outbreak, else 0
# settings with fewer than min_outbreak_runs outbreak runs get score NA and are listed last
score_fit <- function(s) {
  cum_cols <- grep("^cum_", names(s), value = TRUE)
  curve  <- apply(s[, cum_cols, drop = FALSE], 1, function(x) mean((log1p(x) - log1p(obs$cum_cases))^2))
  traced <- (log(s$share_traced_first) - log(target_share_in_quarantine))^2
  final  <- (log(s$num_cases) - log(target_final_size))^2
  outbr  <- ifelse(s$p_outbreak < min_p_outbreak, log(min_p_outbreak / pmax(s$p_outbreak, 0.01))^2, 0)
  s$score_curve <- curve; s$score_traced <- traced; s$score_final <- final; s$score_outbreak <- outbr
  s$score <- ifelse(s$n_outbreak >= min_outbreak_runs, curve + traced + final + outbr, NA)
  s[order(s$score, na.last = TRUE), ]
}
# warn when the best value of a scanned parameter is the lowest or highest grid value
warn_grid_edge <- function(value, grid, name) {
  if (isTRUE(value == min(grid)) || isTRUE(value == max(grid)))
    message(sprintf("%s = %s is at the edge of its grid (%s): widen the grid and rerun", name, value,
                    paste(grid, collapse = ", ")))
}
compare_with_targets <- function(s) {
  sel <- format(as.Date(c("2025-11-14","2025-12-16","2026-01-09","2026-01-27","2026-02-03","2026-03-03","2026-03-17")), "%Y%m%d")
  data.frame(target = c(paste("cumulative cases", sel), "final size (median of outbreak runs)",
                        "share of runs with an outbreak (>= min_outbreak_size)", "adult share (first 997)",
                        "share of cases traced as contacts (first 997)", "coverage ages 5-17"),
             observed = c(obs$cum_cases[match(sel, format(obs$date, "%Y%m%d"))], target_final_size,
                          min_p_outbreak, target_adult_share, target_share_in_quarantine, target_child_coverage),
             model = c(unlist(s[1, paste0("cum_", sel)]), s$num_cases[1], s$p_outbreak[1],
                       s$adult_share_first[1], s$share_traced_first[1], s$coverage_5_17[1]))
}

################################################ #
## STEP 1: POPULATION (hhdense)              ####
################################################ #
if (1 %in% steps) {
  st <- load_state()
  st$pop_file <- make_hhdense(pop_file_in)
  save_state(st)
}

################################################ #
## STEP 2: R0 -> TRANSMISSION PROBABILITY    ####
################################################ #
if (2 %in% steps) {
  source('./bin/rstride/rStride.R')
  st <- load_state()
  set.seed(num_seeds_r0)
  d <- expand.grid(r0 = grid_r0, num_days = 40, rng_seed = seq(num_seeds_r0), start_date = '2023-01-01',
                   num_infected_seeds = 20, seeding_age_min = 1, seeding_age_max = 99,
                   disease_config_file = disease_template, population_file = st$pop_file,
                   age_contact_matrix_file = contact_file, holidays_file = holidays_file,
                   contact_probability_rule = "Mean", stringsAsFactors = FALSE)
  d$rng_seed <- sample(nrow(d))
  project_dir <- run_rStride(exp_design = d, dir_postfix = '_r0', ignore_stdout = TRUE,
                             remove_run_output = TRUE, get_transmission_rdata = TRUE,
                             num_parallel_workers = num_parallel_workers)
  analyse_transmission_data_for_r0(project_dir)
  fit  <- list.files(project_dir, pattern = "_disease_.*\\.xml$", full.names = TRUE)[1]
  dest <- file.path("data", paste0("SC_", basename(fit)))
  file.copy(fit, dest, overwrite = FALSE)
  x  <- readLines(dest)
  b0 <- as.numeric(sub(".*<b0>(.*)</b0>.*", "\\1", grep("<b0>", x, value = TRUE)))
  b1 <- as.numeric(sub(".*<b1>(.*)</b1>.*", "\\1", grep("<b1>", x, value = TRUE)))
  cat(sprintf("[step 2] %s: R0 = %.2f + %.2f p -> p(R0 = %g) = %.4f\n", dest, b0, b1, fixed$r0, (fixed$r0 - b0) / b1))
  st$disease_file <- dest
  save_state(st)
}

################################################ #
## STEP 3: ADULT SUSCEPTIBILITY FACTOR       ####
################################################ #
if (3 %in% steps) {
  source('./bin/rstride/rStride.R')
  st <- load_state()
  files <- vapply(grid_adult_factor, function(f) write_adult_factor_file(adult_raw_file, f, " (calibration scan)"), "")
  p <- get_sc_param_list(st, num_seeds_scan)
  p$immunity_distribution_file <- files
  project_dir <- run_design(p, '_calib_SC_step3_adult')
  m <- run_metrics(file.path("sim_output", basename(project_dir)))
  s <- summarise_runs(m)
  s$factor <- grid_adult_factor[match(s$adult_file, files)]
  s <- s[order(s$factor), ]
  print(s[, c("factor", "runs", "n_outbreak", "p_outbreak", "num_cases", "adult_share_first")], row.names = FALSE)
  if (is.unsorted(s$adult_share_first)) warning("adult share not increasing in the factor: add seeds")
  f_hat <- approx(s$adult_share_first, s$factor, xout = target_adult_share, ties = mean)$y
  if (is.na(f_hat)) stop("target adult share outside the scanned range: extend grid_adult_factor")
  st$adult_factor <- round(f_hat / 0.05) * 0.05
  st$adult_file   <- write_adult_factor_file(adult_raw_file, st$adult_factor,
                       sprintf(", calibrated so that %.1f%% of SC outbreak cases are adults (%s)", 100 * target_adult_share, Sys.Date()))
  cat(sprintf("[step 3] interpolated factor %.3f -> %.2f: %s\n", f_hat, st$adult_factor, st$adult_file))
  write.csv(m, file.path(calib_dir, "step3_runs.csv"), row.names = FALSE)
  save_state(st)
}

################################################ #
## STEP 4: CLUSTER RATE x TRACING "OTHER"    ####
################################################ #
if (4 %in% steps) {
  source('./bin/rstride/rStride.R')
  st <- load_state()
  if (is.null(st$adult_file)) st$adult_file <- write_adult_factor_file(adult_raw_file, st$adult_factor)
  p <- get_sc_param_list(st, num_seeds_scan)
  p$vaccine_hesitancy_cluster_rate <- grid_cluster_rate
  p$tracing_efficiency_other       <- grid_tracing_other
  project_dir <- run_design(p, '_calib_SC_step4_cluster_tracing')
  m <- run_metrics(file.path("sim_output", basename(project_dir)))
  s <- score_fit(summarise_runs(m))
  print(s[, c("cluster_rate", "tracing_other", "n_outbreak", "p_outbreak", "num_cases", "adult_share_first",
              "share_traced_first", "score_curve", "score_traced", "score_final", "score_outbreak", "score")],
        row.names = FALSE, digits = 3)
  if (is.na(s$score[1])) stop("no setting has >= min_outbreak_runs outbreak runs: add seeds or change the grids")
  st$cluster_rate  <- s$cluster_rate[1]
  st$tracing_other <- s$tracing_other[1]
  cat(sprintf("[step 4] best: cluster_rate %.2f, tracing_other %.2f\n", st$cluster_rate, st$tracing_other))
  warn_grid_edge(st$cluster_rate, grid_cluster_rate, "cluster_rate")
  warn_grid_edge(st$tracing_other, grid_tracing_other, "tracing_other")
  write.csv(m, file.path(calib_dir, "step4_runs.csv"), row.names = FALSE)
  write.csv(s, file.path(calib_dir, "step4_scores.csv"), row.names = FALSE)
  save_state(st)
}

################################################ #
## STEP 5: CHECK RUN                         ####
################################################ #
if (5 %in% steps) {
  source('./bin/rstride/rStride.R')
  st <- load_state()
  project_dir <- run_design(get_sc_param_list(st, num_seeds_check), '_calib_SC_step5_check')
  m <- run_metrics(file.path("sim_output", basename(project_dir)))
  s <- score_fit(summarise_runs(m))
  cmp <- compare_with_targets(s)
  print(cmp, row.names = FALSE, digits = 3)
  cat(sprintf("[step 5] P(outbreak) %.2f, median final size %.0f (observed %d)\n",
              s$p_outbreak, s$num_cases, target_final_size))
  if (s$n_outbreak < min_outbreak_runs)
    message(sprintf("only %d outbreak runs: the medians above are unreliable", s$n_outbreak))
  if (abs(log(s$num_cases / target_final_size)) > log(1.5))
    message("median final size is off by more than a factor 1.5: check the cumulative curve after March")
  if (abs(s$adult_share_first - target_adult_share) > 0.02)
    message("adult share is off by > 2 percentage points: rerun step 3 with the step-4 values, then step 4 again")
  # cumulative curves: runs (grey), median (black), DPH (red)
  pdf(file.path(calib_dir, "step5_cumulative_cases.pdf"), 8, 5)
  cum_cols <- grep("^cum_", names(m), value = TRUE); ok <- m[m$num_cases >= min_outbreak_size, ]
  matplot(obs$date, t(ok[, cum_cols]), type = "l", lty = 1, col = "grey80", xlab = "report date",
          ylab = "cumulative cases", ylim = c(0, max(c(obs$cum_cases, unlist(ok[, cum_cols])))))
  lines(obs$date, unlist(s[1, cum_cols]), lwd = 2); points(obs$date, obs$cum_cases, col = "red", pch = 16)
  legend("topleft", c("outbreak runs", "model median", "DPH"), col = c("grey70", "black", "red"),
         lty = c(1, 1, NA), pch = c(NA, NA, 16), bty = "n")
  dev.off()
  write.csv(m, file.path(calib_dir, "step5_runs.csv"), row.names = FALSE)
  write.csv(cmp, file.path(calib_dir, "step5_model_vs_targets.csv"), row.names = FALSE)
  st$check_project <- basename(project_dir)
  save_state(st)
}

if (file.exists(state_file)) str(readRDS(state_file))
