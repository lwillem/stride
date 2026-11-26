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
# CALIBRATION TO THE SPARTANBURG (SC) 2025-26 MEASLES OUTBREAK, STEP BY STEP
#
# Calibration written out as one sequence of steps without helper functions,
# so the order is visible and the output of each step can be checked before
# the next one is run.
# Key mechanisms: vaccine hesitancy (households) and school clustering of
# children from hesitant households. Exemption policy acts through the uptake
# in hesitant households: children left unvaccinated by hesitancy = children
# with a non-medical exemption (CDC SchoolVaxView 2025). SC is calibrated in
# steps 2b-5; step 6 transfers the hesitancy values to the other regions.
# Targets and sources: doc/markdown/measles_sc_calibration.md
#
#   STEP 0  settings: input files, fixed parameters, targets, grids  (run every session)
#   STEP 1  population: household ids renumbered 1..N ("hhdense")
#   STEP 2  R0 -> transmission probability (new disease file)
#   STEP 2b hesitancy uptake and school-cluster rate (SC)
#           (number of cluster schools fixed as a share of schools)
#                                           targets: hesitancy gap = exemption rate (5.6%),
#                                           share of unvaccinated students in the
#                                           cluster schools (46.8%)
#   STEP 3  adult susceptibility factor     target: adult share of cases
#   STEP 4  tracing efficiency "other"      targets: case curve, final size,
#                                            share of runs that take off,
#                                            share of cases traced
#   STEP 5  check run with 100 seeds: model vs all targets, and the values to
#           use in the simulation scripts
#   STEP 6  transfer to CT, ME, MN, OH, TX: uptake from each region's exemption
#           rate with the same rule as SC; clusters from SC (CT: own school data)
#
# HOW TO USE
#   - Run from the main project folder (opt/stride), e.g. open this file in
#     RStudio with opt/stride as working directory.
#   - Run STEP 0, then ONE step at a time (RStudio: Ctrl+Alt+T runs the
#     current section). Read the "CHECK BEFORE CONTINUING" list of the step
#     before running the next one.
#   - Every step writes its results to `calib_dir` and the next step reads
#     them from there. You can close R between steps: run STEP 0 again first.
#   - County: set county_code in STEP 0. Steps 1 and 2 work for every county
#     (SC, CT, ME, MN, OH, TX); steps 2b-6 run with county_code "SC" (step 6
#     needs step 1 to have been run once for each other region).
#   - calib_dir contains the date (Sys.Date()): if you continue on another day,
#     set calib_dir to the folder of the earlier steps.
#   - When the population or contact matrix change: steps 1-6 (incl. 2b).
#     When the child immunity file, exemption rates or school targets change: steps 2b-6.
#     When only the adult immunity file changes: steps 3-5.
#   - Second round (step 5 says the adult share is off): run step 3 again
#     (it then uses the step-4 values), then steps 4 and 5.
#
# rStride functions used (the model interface, not defined here):
#   get_default_param(), .rstride$get_full_grid_exp_design(), run_rStride(),
#   analyse_transmission_data_for_r0()
#############################################################################


#############################################################################
## STEP 0: SETTINGS  (run at the start of every R session)               ####
#############################################################################
# PURPOSE  Define every input of the calibration in one place.
# INPUTS   the files below must exist in opt/stride/data
# OUTPUTS  objects in the R session only; creates calib_dir

rm(list = ls())
source('./bin/rstride/rStride.R')

county_code          <- "SC"                            # SC, CT, ME, MN, OH or TX
calib_dir            <- paste0("sim_output/", Sys.Date(), "_calibration_", county_code)   # all step results are written here; one folder per county
num_parallel_workers <- 3
dir.create(calib_dir, showWarnings = FALSE, recursive = TRUE)

## ---- input files and rates per county (the ones that may be updated) ----
# hesitancy_rate: ASPE county "estimated hesitant" (COVID-19, Household Pulse Survey; data.cdc.gov q9mh-h2tw)
# exemption_rate: non-medical exemptions among kindergarteners, state level, CDC SchoolVaxView 2025
#   (policy: CT, ME medical only; SC religious; MN personal; OH, TX personal and religious;
#    ME has no non-medical exemptions since 2021 -> 0)
county_inputs <- list(
  SC = list(pop = "pop_usa_sc_spartanburg_c1000.csv", contacts = "contact_matrix_usa_sc_spartanburg_c1000_v2.xml", hesitancy_rate = 0.14,  exemption_rate = 0.056),
  CT = list(pop = "pop_usa_ct_fairfield_c1000.csv",   contacts = "contact_matrix_usa_ct_fairfield_c1000_v2.xml",   hesitancy_rate = 0.05,  exemption_rate = 0.001),
  ME = list(pop = "pop_usa_me_piscataquis_c1000.csv", contacts = "contact_matrix_usa_me_piscataquis_c1000_v2.xml", hesitancy_rate = 0.105, exemption_rate = 0.000),
  MN = list(pop = "pop_usa_mn_hennepin_c1000.csv",    contacts = "contact_matrix_usa_mn_hennepin_c1000_v2.xml",    hesitancy_rate = 0.062, exemption_rate = 0.068),
  OH = list(pop = "pop_usa_oh_cuyahoga_c1000.csv",    contacts = "contact_matrix_usa_oh_cuyahoga_c1000_v2.xml",    hesitancy_rate = 0.122, exemption_rate = 0.055),
  TX = list(pop = "pop_usa_tx_gaines_c1000.csv",      contacts = "contact_matrix_usa_tx_gaines_c1000_v2.xml",      hesitancy_rate = 0.117, exemption_rate = 0.050))
exemption_ages <- c(5, 6)   # kindergarten: the ages at which the exemption rate is compared (steps 2b and 6)
pop_file_in      <- file.path("data/population_6region", county_inputs[[county_code]]$pop)  # from social_contacts_usa2026.R
contact_file     <- file.path("data/contacts_usa", county_inputs[[county_code]]$contacts)
holidays_file    <- "data/calendar_USA_2023_2026_measles.csv"
disease_template <- "data/disease_measles_usa.xml"                                  # natural history; b0/b1 are refitted in step 2
adult_raw_file   <- sprintf("data/immunity_6region/immunity_measles_adult_%s.xml", county_code)  # natural immunity, ages 18+ (uncalibrated)
child_file       <- sprintf("data/immunity_6region/immunity_measles_child_%s.xml", county_code)  # immunity by age, ages 0-17

## ---- fixed parameters (not calibrated; sources: table 2 of the documentation) ----
fixed <- list(
  r0                                 = 12,            # DPH: "12 to 18 non-immune contacts"; lower end
  start_date                         = "2025-09-15",  # assumption: ~1 week before the first outbreak cases
  num_days                           = 240,           # to 2026-05-13
  num_infected_seeds                 = 1,
  seeding_age_min                    = 1,
  seeding_age_max                    = 17,
  contact_tracing_date               = "2025-10-02",  # outbreak declared (DPH)
  vaccine_hesitancy_rate             = county_inputs[[county_code]]$hesitancy_rate,   # ASPE (SC: Spartanburg 0.1443)
  vaccine_hesitancy_cluster_pooltype = "School",
  vaccine_hesitancy_cluster_uptake   = 0,             # calibration baseline; a scenario parameter later
                                                      # (how do cases change when uptake rises in the clusters?)
  # used only if step 2b has not been run (step 2b fits these two):
  vaccine_hesitancy_uptake           = 0.5,
  vaccine_hesitancy_cluster_rate     = 0.20,
  # vaccine_hesitancy_num_clusters is not set here: step 1 computes it as
  # cluster_school_share x the number of schools in the population file
  detection_probability              = 0.9,           # assumption
  tracing_efficiency_household       = 0.95,          # assumption
  delay_isolation_index              = 3,             # assumption
  delay_contact_tracing              = 1,             # assumption
  case_finding_capacity              = 10000,         # keep < 1e5: R writes 100000 as "1e+05", C++ reads 1
  test_false_negative                = 0,
  isolation_duration_index           = 7,             # model choice (guidance: isolate through 4 days after rash)
  isolation_duration_contact         = 21             # DPH: 21-day quarantine
)

## ---- school coverage targets per county (2025-26 school year) ----
# One row per target: coverage (with the ages it refers to; a check only),
# cluster_school_share and unvaccinated_share_in_clusters (SC: step 2b target).
# Sources are in the file. Counties without school-level data use the SC share.
school_targets_file <- "data/immunity_6region/school_coverage_targets_2025-26.csv"
school_targets      <- read.csv(school_targets_file, stringsAsFactors = FALSE)
school_targets      <- school_targets[school_targets$county_code == county_code, ]
print(school_targets[, c("target", "age_min", "age_max", "value", "measure")], row.names = FALSE)

## ---- share of schools that are hesitant clusters ----
# Step 1 turns this into the number of cluster schools (vaccine_hesitancy_num_clusters)
# = round(share x schools in the population file); SC: 24 of 93 schools below 85%
# -> 25.8% of 90 schools = 23.
cluster_school_share <- school_targets$value[school_targets$target == "cluster_school_share"]

## ---- tracing "other" used in steps 2b and 3 before step 4 has been run ----
# (in a second round, step 3 reads the step-4 result from calib_dir instead)
start_tracing_other <- 0.20

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
target_adult_share         <- 87 / 997     # 87 adults among 997 cases (DPH data via WSPA)
target_share_in_quarantine <- 496 / 989    # 496 of 2,219 quarantined fell ill in quarantine; 989 cases (DPH, 2026-03-02)
target_child_coverage      <- 0.888        # Spartanburg PK-12, 57,280 students in 93 schools with complete immunizations (DPH 45-Day report 2025-26); check only

## ---- how model output is compared with the targets ----
report_lag        <- 3      # days from symptom onset to DPH report (assumption)
min_outbreak_size <- 200    # runs with fewer cases did not take off; excluded from the medians
min_p_outbreak    <- 0.5    # settings where fewer runs take off are penalised in the step-4 score
min_outbreak_runs <- 10     # settings with fewer outbreak runs are not scored
num_obs_cases     <- 997    # adult share and share traced are measured over the first 997 cases

## ---- scan grids and seeds ----
grid_r0            <- seq(0, 15, 3)
num_seeds_r0       <- 10
grid_adult_factor  <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
grid_tracing_other <- c(0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60)
num_seeds_scan     <- 40     # step 3
num_seeds_tracing  <- 100    # step 4 (one parameter only, so more seeds per value)
num_seeds_check    <- 100    # step 5

## ---- step 2b: grids and school data ----
grid_hesitancy_uptake    <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0)
grid_cluster_rate_school <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0)
num_seeds_school         <- 5       # each seed draws other cluster schools
min_school_size          <- 10      # model schools with fewer students are left out (DPH reports them as "<10")
# SC only: school-level data for the extra check plot in step 2b (not used in the score)
school_data_file         <- "data/immunity_6region/45_Day_Report_25-26_cleaned.csv"   # cleaned by clean_sc_45day_report.R
school_data_county       <- "Spartanburg"
step2b_delete_snapshots  <- TRUE    # population_snapshot.csv is ~14 MB per run; deleted after reading (steps 2b and 6)

## ---- step 6: transfer to the other regions ----
regions_transfer     <- c("CT", "ME", "MN", "OH", "TX")
grid_uptake_transfer <- seq(0, 1, by = 0.1)
num_seeds_transfer   <- 5

## ---- STRIDE parameters shared by steps 3-5 ----
# Steps 3-5 copy this list and add: population_file, disease_config_file,
# immunity_distribution_file, vaccine_hesitancy_cluster_rate,
# tracing_efficiency_other and num_rng_seeds.
sc_param_base <- get_default_param()
for (k in names(fixed)) sc_param_base[[k]] <- fixed[[k]]
sc_param_base$age_contact_matrix_file      <- contact_file
sc_param_base$holidays_file                <- holidays_file
sc_param_base$contact_probability_rule     <- "Mean"
sc_param_base$event_log_level              <- "Transmissions"   # [PRIM]/[TRAN]/[TRACE] lines are read in steps 3-5
sc_param_base$output_pop_snapshot          <- TRUE              # writes susceptibles_by_age.csv (coverage check)
sc_param_base$reference_hospital_data_file <- NA
sc_param_base$reference_serology_data_file <- NA
sc_param_base$immunity_profile             <- "AgeDependent"
sc_param_base$immunity_link_probability    <- 0
sc_param_base$contact_survey_dates         <- sc_param_base$start_date
sc_param_base$vaccine_profile              <- "AgeDependent"
sc_param_base$vaccine_distribution_file    <- child_file
sc_param_base$vaccine_min_age              <- 0
sc_param_base$vaccine_max_age              <- 17
sc_param_base$mass_immunize                <- TRUE
sc_param_base$hosp_probability_factor      <- 1
sc_param_base$hospital_length_of_stay      <- 14
sc_param_base$hospital_probability_age     <- paste(0.18, 0.06, 0.12, sep = ',')
sc_param_base$hospital_category_age        <- paste(0, 5, 20, sep = ',')
sc_param_base$hospital_mean_delay_age      <- paste(2, 2, 5, sep = ',')

# Event-log fields (after splitting a line on spaces; field 1 is the tag):
#   [PRIM]/[TRAN]: 2 local_id, 3 infector_id, 4 age, 5 infector_age, 6 pool_type,
#                  7 sim_day (infection day), 8 id_index_case, 9 start_infectiousness,
#                  10 end_infectiousness, 11 start_symptoms (days since infection), ...
#   [TRACE]:       2 person_id, 3 age, 4 is_infected, 5 is_symptomatic,
#                  6 pool ("Index" = detected case, otherwise a traced contact), 7 case_id, ...


#############################################################################
## STEP 1: POPULATION (household ids 1..N)                                ####
#############################################################################
# PURPOSE  Renumber household ids 1..N in the order of the original ids. Same
#          results as the original file, ~200x faster (STRIDE allocates one
#          pool per id up to the largest id).
# INPUTS   pop_file_in
#          Also counts the schools and sets the number of cluster schools
#          = round(cluster_school_share x schools); used by steps 2b-5.
# OUTPUTS  <pop>_hhdense.csv and <pop>_hhdense_household_map.csv next to the input;
#          calib_dir/step1_population.csv (file name, persons, households,
#          schools, cluster_school_share, cluster_schools)
#
# CHECK BEFORE CONTINUING
#   - persons and households are as expected for the county
#   - after renumbering, the largest household id equals the number of households
#   - schools and cluster_schools (SC: 90 schools -> 23 cluster schools)

pop <- read.csv(pop_file_in, stringsAsFactors = FALSE)
hh_ids_original  <- sort(unique(pop$household_id))
pop$household_id <- match(pop$household_id, hh_ids_original)

step1_pop_file <- sub("\\.csv$", "_hhdense.csv", pop_file_in)
step1_map_file <- sub("\\.csv$", "_hhdense_household_map.csv", pop_file_in)
write.csv(pop, step1_pop_file, row.names = FALSE, quote = FALSE, na = "NA")
write.csv(data.frame(household_id = seq_along(hh_ids_original), original_household_id = hh_ids_original),
          step1_map_file, row.names = FALSE, quote = FALSE)

school_ids   <- unique(pop$school_id[!is.na(pop$school_id) & pop$school_id != 0])
step1_result <- data.frame(pop_file = step1_pop_file, persons = nrow(pop),
                           households = length(hh_ids_original),
                           max_original_household_id = max(hh_ids_original),
                           max_new_household_id = max(pop$household_id),
                           schools = length(school_ids),
                           cluster_school_share = cluster_school_share,
                           cluster_schools = round(cluster_school_share * length(school_ids)))
write.csv(step1_result, file.path(calib_dir, "step1_population.csv"), row.names = FALSE)
print(step1_result)
rm(pop)


#############################################################################
## STEP 2: R0 -> TRANSMISSION PROBABILITY                                 ####
#############################################################################
# PURPOSE  Fit the link between R0 and the transmission probability for this
#          population and contact matrix: secondary cases of 20 index cases
#          (ages 1-99) in a fully susceptible population over 40 days, for a
#          grid of R0 values, then a linear fit R0 = b0 + b1 * p.
# INPUTS   step 1: calib_dir/step1_population.csv; contact_file; disease_template
# OUTPUTS  data/<county_code>_<run_tag>_disease_measles_usa.xml (new b0, b1);
#          calib_dir/step2_disease.csv (file, b0, b1, transmission probability at R0 = 12)
#
# CHECK BEFORE CONTINUING
#   - the fit plot <project>/<run_tag>_fit_r0.pdf: means close to the line, R^2 high
#   - p at R0 = 12 is plausible (16 Nov and 20 Nov 2026: 0.5346; Oct 2026: 0.5566).
#     If the population and contact matrix did not change, it should be the same.

step1_result <- read.csv(file.path(calib_dir, "step1_population.csv"), stringsAsFactors = FALSE)

set.seed(num_seeds_r0)
step2_design <- expand.grid(r0 = grid_r0, num_days = 40, rng_seed = seq(num_seeds_r0),
                            start_date = '2023-01-01', num_infected_seeds = 20,
                            seeding_age_min = 1, seeding_age_max = 99,
                            disease_config_file = disease_template,
                            population_file = step1_result$pop_file,
                            age_contact_matrix_file = contact_file,
                            holidays_file = holidays_file,
                            contact_probability_rule = "Mean",
                            stringsAsFactors = FALSE)
step2_design$rng_seed <- sample(nrow(step2_design))

step2_project <- run_rStride(exp_design = step2_design, dir_postfix = '_r0', ignore_stdout = TRUE,
                             remove_run_output = TRUE, get_transmission_rdata = TRUE,
                             num_parallel_workers = num_parallel_workers)
analyse_transmission_data_for_r0(step2_project)

step2_fit_file     <- list.files(step2_project, pattern = "_disease_.*\\.xml$", full.names = TRUE)[1]
step2_disease_file <- file.path("data", paste0(county_code, "_", basename(step2_fit_file)))
file.copy(step2_fit_file, step2_disease_file, overwrite = FALSE)

step2_xml <- readLines(step2_disease_file)
step2_b0  <- as.numeric(sub(".*<b0>(.*)</b0>.*", "\\1", grep("<b0>", step2_xml, value = TRUE)))
step2_b1  <- as.numeric(sub(".*<b1>(.*)</b1>.*", "\\1", grep("<b1>", step2_xml, value = TRUE)))
step2_result <- data.frame(disease_file = step2_disease_file, project = basename(step2_project),
                           b0 = step2_b0, b1 = step2_b1, r0 = fixed$r0,
                           transmission_probability = (fixed$r0 - step2_b0) / step2_b1)
write.csv(step2_result, file.path(calib_dir, "step2_disease.csv"), row.names = FALSE)
print(step2_result)


#############################################################################
## STEP 2b: VACCINE HESITANCY AND SCHOOL CLUSTERS (SC)                    ####
#############################################################################
# PURPOSE  Choose how vaccine hesitancy is spread over households and schools
#          in Spartanburg at the start of the simulation (t = 0):
#            vaccine_hesitancy_uptake        uptake in hesitant households outside the
#                                            cluster schools (the policy lever)      fitted
#            vaccine_hesitancy_cluster_rate  share of the households linked to the
#                                            cluster schools that is hesitant          fitted
#            vaccine_hesitancy_num_clusters  number of cluster schools, from step 1
#                                            (cluster_school_share x schools)          fixed
#            vaccine_hesitancy_cluster_uptake  uptake in clustered hesitant households:
#                                            0 here; a scenario parameter later         fixed
#          vaccine_hesitancy_rate stays at the ASPE value.
#
#          Targets:
#            hesitancy_gap   = immunity at ages 5-6 without hesitancy (from the child
#                              immunity file) minus immunity at ages 5-6 with hesitancy
#                            = the non-medical exemption rate (kindergarten, CDC
#                              SchoolVaxView 2025; SC 5.6%).
#                              Assumption: children unvaccinated because of hesitancy
#                              are the children with a non-medical exemption.
#            concentration   = share of all non-immune school attendees who are in the
#                              lowest-coverage cluster_school_share of the schools
#                              (SC 45-Day report: 46.8% in the lowest 24 of 93 schools)
#          Score (lower = better):
#            score_gap           ((hesitancy_gap - exemption rate) / 0.005)^2
#            score_concentration ((concentration - target) / 0.05)^2
#          so 0.5 percentage point on the gap weighs as much as 5 points on the
#          concentration. County coverage (DPH 88.8%, PK-12) is printed as a check only.
#          Only t = 0 is needed: every run lasts 1 day and writes
#          population_snapshot.csv (age, school and immunity of every person).
# INPUTS   step 1 (population, number of cluster schools), step 2 (disease file,
#          only needed to start STRIDE); exemption rate and school targets (STEP 0)
# OUTPUTS  calib_dir/step2b_schools.csv (coverage per school per run),
#          calib_dir/step2b_runs.csv (one row per run),
#          calib_dir/step2b_scores.csv (one row per combination, best first),
#          calib_dir/step2b_best.csv (chosen uptake and cluster rate; number of clusters),
#          calib_dir/step2b_school_coverage.pdf (sorted school coverage: DPH vs best combination)
#
# CHECK BEFORE CONTINUING
#   - best row: hesitancy_gap within ~0.5 point of the exemption rate (5.6%),
#     concentration within ~5 points of 46.8%
#   - check: coverage_3_18 close to the DPH coverage (88.8%)
#   - the chosen values are not at the edge of their grids
#   - compare the top rows: if their scores are close, seed noise can decide the
#     order; prefer the middle of that group or rerun with more seeds

if (county_code != "SC") stop("step 2b calibrates SC; other counties get their values in step 6")
step1_result <- read.csv(file.path(calib_dir, "step1_population.csv"), stringsAsFactors = FALSE)
step2_result <- read.csv(file.path(calib_dir, "step2_disease.csv"), stringsAsFactors = FALSE)

# 2b-1. targets
target_gap           <- county_inputs[[county_code]]$exemption_rate
target_concentration <- school_targets$value[school_targets$target == "unvaccinated_share_in_clusters"]
check_coverage       <- school_targets[school_targets$target == "coverage", ]
cov_cols             <- sprintf("coverage_%d_%d", check_coverage$age_min, check_coverage$age_max)

# immunity at ages 5-6 without hesitancy: child immunity file, weighted by the
# number of school attendees of each age in the population
child_txt   <- readLines(child_file, warn = FALSE)
child_m     <- regmatches(child_txt, regexec("<age([0-9]+)>([-0-9.eE]+)</age[0-9]+>", child_txt))
child_m     <- child_m[lengths(child_m) == 3]
child_imm   <- setNames(as.numeric(sapply(child_m, `[`, 3)), sapply(child_m, `[`, 2))
pop_ages    <- read.csv(step1_result$pop_file, colClasses = c("integer", "NULL", "integer", "NULL", "NULL", "NULL"))
pop_ages    <- pop_ages$age[!is.na(pop_ages$school_id) & pop_ages$school_id != 0 & pop_ages$age %in% exemption_ages]
base_immunity_5_6 <- mean(child_imm[as.character(pop_ages)])
cat(sprintf("[step 2b] SC targets: hesitancy gap at ages %s = %.3f (exemption rate); concentration %.3f\n",
            paste(exemption_ages, collapse = "-"), target_gap, target_concentration))
cat(sprintf("[step 2b] immunity at ages %s without hesitancy (child file): %.4f -> target immunity with hesitancy %.4f\n",
            paste(exemption_ages, collapse = "-"), base_immunity_5_6, base_immunity_5_6 - target_gap))

# 2b-2. number of cluster schools: from step 1 (cluster_school_share x schools)
if (is.null(step1_result$cluster_schools)) stop("step1_population.csv has no cluster_schools: rerun step 1")
step2b_num_clusters <- step1_result$cluster_schools
cat(sprintf("[step 2b] %d schools in the population; %.1f%% of them -> %d cluster schools\n",
            step1_result$schools, 100 * step1_result$cluster_school_share, step2b_num_clusters))

# 2b-3. run STRIDE for 1 day for every combination
step2b_param <- sc_param_base
step2b_param$num_rng_seeds                  <- num_seeds_school
step2b_param$num_days                       <- 1
step2b_param$population_file                <- step1_result$pop_file
step2b_param$disease_config_file            <- step2_result$disease_file
step2b_param$immunity_distribution_file     <- adult_raw_file      # adults do not matter here
step2b_param$tracing_efficiency_other       <- start_tracing_other
step2b_param$vaccine_hesitancy_uptake       <- grid_hesitancy_uptake       # scanned
step2b_param$vaccine_hesitancy_num_clusters <- step2b_num_clusters         # fixed (2b-2)
step2b_param$vaccine_hesitancy_cluster_rate <- grid_cluster_rate_school    # scanned
step2b_numeric <- unlist(step2b_param[sapply(step2b_param, is.numeric)])
if (any(abs(step2b_numeric) >= 1e5, na.rm = TRUE)) stop("a numeric parameter >= 1e5 would be written as 1e+05 (read as 1 by STRIDE)")
set.seed(2026)   # fixes which rng seeds go to which combination, so a rerun gives the same result
step2b_design  <- .rstride$get_full_grid_exp_design(exp_param_list = step2b_param, num_rng_seeds = step2b_param$num_rng_seeds)
step2b_project <- run_rStride(step2b_design, paste0('_calib_', county_code, '_step2b_school_coverage'),
                              remove_run_output = FALSE, num_parallel_workers = num_parallel_workers)

# 2b-4. per run: hesitancy gap at ages 5-6, concentration in the lowest schools, school coverage
step2b_summary <- read.csv(file.path(step2b_project, paste0(basename(step2b_project), "_summary.csv")), stringsAsFactors = FALSE)
step2b_runs    <- data.frame()
step2b_schools <- data.frame()
for (i in seq_len(nrow(step2b_summary))) {
  snap_file <- file.path(step2b_summary$output_prefix[i], "population_snapshot.csv")
  snap <- read.csv(snap_file, colClasses = c("integer", "NULL", "integer", "NULL", "NULL", "NULL", "character"))
  snap <- snap[snap$school_id != 0, ]                       # school attendees
  snap$immune <- snap$immunity_status == "immune"
  # schools with >= min_school_size students, sorted by coverage
  sch_students <- as.numeric(table(snap$school_id))
  sch_immune   <- as.numeric(tapply(snap$immune, snap$school_id, sum))
  sch_id       <- sort(unique(snap$school_id))
  keep         <- sch_students >= min_school_size
  sch_students <- sch_students[keep]; sch_immune <- sch_immune[keep]; sch_id <- sch_id[keep]
  sch_coverage <- sch_immune / sch_students
  is_low       <- seq_along(sch_id) %in% order(sch_coverage)[seq_len(round(cluster_school_share * length(sch_id)))]  # lowest schools
  # coverage at the ages of the check targets
  run_cov <- numeric(0)
  for (k in seq_len(nrow(check_coverage))) {
    at_age  <- snap$age >= check_coverage$age_min[k] & snap$age <= check_coverage$age_max[k] & snap$school_id %in% sch_id
    run_cov <- c(run_cov, mean(snap$immune[at_age]))
  }
  immunity_5_6 <- mean(snap$immune[snap$age %in% exemption_ages])
  step2b_runs <- rbind(step2b_runs, data.frame(
    run                 = basename(step2b_summary$output_prefix[i]),
    uptake              = step2b_summary$vaccine_hesitancy_uptake[i],
    num_clusters        = step2b_summary$vaccine_hesitancy_num_clusters[i],
    cluster_rate        = step2b_summary$vaccine_hesitancy_cluster_rate[i],
    immunity_5_6        = immunity_5_6,
    hesitancy_gap       = base_immunity_5_6 - immunity_5_6,
    concentration       = sum(sch_students[is_low] - sch_immune[is_low]) / sum(sch_students - sch_immune),
    coverage_all        = sum(sch_immune) / sum(sch_students),
    coverage_low        = sum(sch_immune[is_low]) / sum(sch_students[is_low]),
    coverage_rest       = sum(sch_immune[!is_low]) / sum(sch_students[!is_low]),
    below_90            = sum(sch_coverage < 0.90),
    below_85            = sum(sch_coverage < 0.85),
    lowest              = min(sch_coverage),
    t(setNames(run_cov, cov_cols))))
  step2b_schools <- rbind(step2b_schools, data.frame(
    run = basename(step2b_summary$output_prefix[i]), uptake = step2b_summary$vaccine_hesitancy_uptake[i],
    cluster_rate = step2b_summary$vaccine_hesitancy_cluster_rate[i],
    school_id = sch_id, students = sch_students, coverage = sch_coverage))
  if (step2b_delete_snapshots) file.remove(snap_file)
}
write.csv(step2b_runs, file.path(calib_dir, "step2b_runs.csv"), row.names = FALSE)
write.csv(step2b_schools, file.path(calib_dir, "step2b_schools.csv"), row.names = FALSE)

# 2b-5. medians per combination and score
step2b_scores <- aggregate(step2b_runs[, c("hesitancy_gap", "concentration", cov_cols, "coverage_low", "coverage_rest",
                                           "below_90", "below_85", "lowest")],
                           by = list(uptake = step2b_runs$uptake, cluster_rate = step2b_runs$cluster_rate),
                           FUN = median)
step2b_scores$score_gap           <- ((step2b_scores$hesitancy_gap - target_gap) / 0.005)^2
step2b_scores$score_concentration <- ((step2b_scores$concentration - target_concentration) / 0.05)^2
step2b_scores$score               <- step2b_scores$score_gap + step2b_scores$score_concentration
step2b_scores <- step2b_scores[order(step2b_scores$score), ]
write.csv(step2b_scores, file.path(calib_dir, "step2b_scores.csv"), row.names = FALSE)
cat(sprintf("targets: hesitancy gap %.3f, concentration %.3f; check: %s\n", target_gap, target_concentration,
            paste(sprintf("coverage ages %d-%d %.3f", check_coverage$age_min, check_coverage$age_max, check_coverage$value), collapse = ", ")))
print(head(step2b_scores, 10), row.names = FALSE, digits = 3)

# 2b-6. best combination
step2b_best <- data.frame(county_code                      = county_code,
                          vaccine_hesitancy_rate           = fixed$vaccine_hesitancy_rate,
                          vaccine_hesitancy_uptake         = step2b_scores$uptake[1],
                          vaccine_hesitancy_num_clusters   = step2b_num_clusters,
                          vaccine_hesitancy_cluster_rate   = step2b_scores$cluster_rate[1],
                          vaccine_hesitancy_cluster_uptake = fixed$vaccine_hesitancy_cluster_uptake,
                          cluster_school_share = step1_result$cluster_school_share, schools_in_population = step1_result$schools,
                          exemption_rate = target_gap, base_immunity_5_6 = base_immunity_5_6,
                          hesitancy_gap = step2b_scores$hesitancy_gap[1], concentration = step2b_scores$concentration[1],
                          score = step2b_scores$score[1], project = basename(step2b_project))
write.csv(step2b_best, file.path(calib_dir, "step2b_best.csv"), row.names = FALSE)
print(step2b_best)
if (step2b_best$vaccine_hesitancy_uptake %in% range(grid_hesitancy_uptake))          message("CHECK: uptake is at the edge of its grid")
if (step2b_best$vaccine_hesitancy_cluster_rate %in% range(grid_cluster_rate_school)) message("CHECK: cluster rate is at the edge of its grid")

# 2b-7. sorted school coverage, DPH (red) vs the runs of the best combination (grey)
school_data <- read.csv(school_data_file, stringsAsFactors = FALSE, check.names = FALSE)
school_data <- school_data[!is.na(school_data$County) & school_data$County == school_data_county, ]
school_data$students <- suppressWarnings(as.numeric(school_data$`Total Students`))   # "<10" becomes NA
school_data <- school_data[!is.na(school_data$students) & !is.na(school_data$Immunizations), ]
dph_coverage <- sort(as.numeric(school_data$Immunizations))
best_schools <- step2b_schools[step2b_schools$uptake == step2b_best$vaccine_hesitancy_uptake &
                               step2b_schools$cluster_rate == step2b_best$vaccine_hesitancy_cluster_rate, ]
pdf(file.path(calib_dir, "step2b_school_coverage.pdf"), 8, 5)
plot(seq_along(dph_coverage) / length(dph_coverage), dph_coverage, type = "n", ylim = c(0, 1),
     xlab = "schools, sorted by coverage (share of schools)", ylab = "coverage",
     main = sprintf("uptake %.2f, %d clusters, cluster rate %.2f", step2b_best$vaccine_hesitancy_uptake,
                    step2b_best$vaccine_hesitancy_num_clusters, step2b_best$vaccine_hesitancy_cluster_rate))
for (r in unique(best_schools$run)) {
  cv <- sort(best_schools$coverage[best_schools$run == r])
  lines(seq_along(cv) / length(cv), cv, col = "grey60")
}
lines(seq_along(dph_coverage) / length(dph_coverage), dph_coverage, col = "red", lwd = 2)
abline(h = 0.85, lty = 3); abline(v = cluster_school_share, lty = 3)
legend("bottomright", c("model runs, best combination", "DPH 45-Day report 2025-26"),
       col = c("grey60", "red"), lwd = c(1, 2), bty = "n")
dev.off()


#############################################################################
## STEP 3: ADULT SUSCEPTIBILITY FACTOR                                    ####
#############################################################################
# PURPOSE  Scale adult susceptibility so that adults are 8.7% of the first 997
#          cases: immunity_new = 1 - f * (1 - immunity_raw) for ages >= 18.
#          Simulates every factor in grid_adult_factor and interpolates the
#          factor at the target share.
# INPUTS   SC only. step 1 (population), step 2 (disease file), step 2b (hesitancy and
#          school clusters); adult_raw_file; tracing "other": start value (first
#          round) or calib_dir/step4_best.csv (second round)
# OUTPUTS  one adult file per factor (..._calib_fXXX.xml, scan only);
#          the calibrated adult file (factor rounded to 0.05);
#          calib_dir/step3_runs.csv (one row per run),
#          calib_dir/step3_by_factor.csv (medians per factor),
#          calib_dir/step3_adult_factor.csv (chosen factor and file)
#
# CHECK BEFORE CONTINUING
#   - step3_by_factor: adult share rises with the factor
#   - the target (8.7%) lies between the lowest and highest share; otherwise widen the grid
#   - each factor has >= 10 outbreak runs (n_outbreak); otherwise add seeds
#   - the chosen factor is not the lowest or highest grid value

if (county_code != "SC") stop("steps 3-5 calibrate to the SC outbreak data: for other counties run steps 1, 2 and 2b only")
step1_result <- read.csv(file.path(calib_dir, "step1_population.csv"), stringsAsFactors = FALSE)
step2_result <- read.csv(file.path(calib_dir, "step2_disease.csv"), stringsAsFactors = FALSE)

if (file.exists(file.path(calib_dir, "step2b_best.csv"))) {
  step2b_best <- read.csv(file.path(calib_dir, "step2b_best.csv"))
} else {                                                             # step 2b not run: fixed values
  step2b_best <- data.frame(vaccine_hesitancy_uptake       = fixed$vaccine_hesitancy_uptake,
                            vaccine_hesitancy_num_clusters = step1_result$cluster_schools,
                            vaccine_hesitancy_cluster_rate = fixed$vaccine_hesitancy_cluster_rate)
}
step3_tracing_other <- start_tracing_other
if (file.exists(file.path(calib_dir, "step4_best.csv"))) {          # second round
  step3_tracing_other <- read.csv(file.path(calib_dir, "step4_best.csv"))$tracing_other
}
step3_cluster_rate <- step2b_best$vaccine_hesitancy_cluster_rate
cat(sprintf("[step 3] uptake %.2f, %d clusters, cluster rate %.2f, tracing other %.2f\n",
            step2b_best$vaccine_hesitancy_uptake, step2b_best$vaccine_hesitancy_num_clusters,
            step3_cluster_rate, step3_tracing_other))

# 3a. one adult immunity file per factor
adult_raw_txt <- readLines(adult_raw_file, warn = FALSE)
adult_age_m   <- regmatches(adult_raw_txt, regexec("<age([0-9]+)>([-0-9.eE]+)</age[0-9]+>", adult_raw_txt))
adult_age_m   <- adult_age_m[lengths(adult_age_m) == 3]
adult_age     <- as.integer(sapply(adult_age_m, `[`, 2))
adult_imm_raw <- as.numeric(sapply(adult_age_m, `[`, 3))
adult_header  <- adult_raw_txt[grepl("^<!--", adult_raw_txt)][1]
adult_source  <- sub(".*<data_source>(.*)</data_source>.*", "\\1", grep("<data_source>", adult_raw_txt, value = TRUE))
adult_manip   <- sub(".*<data_manipulation>(.*)</data_manipulation>.*", "\\1", grep("<data_manipulation>", adult_raw_txt, value = TRUE))

step3_files <- character(0)
for (f in grid_adult_factor) {
  imm <- ifelse(adult_age >= 18, 1 - f * (1 - adult_imm_raw), adult_imm_raw)
  out <- sub("\\.xml$", sprintf("_calib_f%03d.xml", as.integer(round(100 * f))), adult_raw_file)
  writeLines(c(adult_header, "<immunity>",
               sprintf("  <data_source>%s</data_source>", adult_source),
               sprintf("  <data_manipulation>%s; susceptibility (1 - immunity) at ages 18+ multiplied by %.2f (calibration scan)</data_manipulation>",
                       adult_manip, f),
               sprintf("  <age%d>%.4f</age%d>", adult_age, imm, adult_age), "</immunity>"), out)
  step3_files <- c(step3_files, out)
}

# 3b. run STRIDE for every factor
step3_param <- sc_param_base
step3_param$num_rng_seeds                  <- num_seeds_scan
step3_param$population_file                <- step1_result$pop_file
step3_param$disease_config_file            <- step2_result$disease_file
step3_param$vaccine_hesitancy_uptake       <- step2b_best$vaccine_hesitancy_uptake
step3_param$vaccine_hesitancy_num_clusters <- step2b_best$vaccine_hesitancy_num_clusters
step3_param$vaccine_hesitancy_cluster_rate <- step3_cluster_rate
step3_param$tracing_efficiency_other       <- step3_tracing_other
step3_param$immunity_distribution_file     <- step3_files          # the scanned parameter
step3_numeric <- unlist(step3_param[sapply(step3_param, is.numeric)])
if (any(abs(step3_numeric) >= 1e5, na.rm = TRUE)) stop("a numeric parameter >= 1e5 would be written as 1e+05 (read as 1 by STRIDE)")
step3_design  <- .rstride$get_full_grid_exp_design(exp_param_list = step3_param, num_rng_seeds = step3_param$num_rng_seeds)
step3_project <- run_rStride(step3_design, '_calib_SC_step3_adult', remove_run_output = FALSE,
                             num_parallel_workers = num_parallel_workers)

# 3c. adult share among the first 997 cases of every run
step3_summary <- read.csv(file.path(step3_project, paste0(basename(step3_project), "_summary.csv")), stringsAsFactors = FALSE)
step3_runs <- data.frame()
for (i in seq_len(nrow(step3_summary))) {
  log_txt <- readLines(file.path(step3_summary$output_prefix[i], "event_log.txt"), warn = FALSE)
  cases   <- strsplit(log_txt[startsWith(log_txt, "[PRIM]") | startsWith(log_txt, "[TRAN]")], " ", fixed = TRUE)
  age     <- as.integer(sapply(cases, `[`, 4))
  inf_day <- as.integer(sapply(cases, `[`, 7))
  first   <- head(order(inf_day), num_obs_cases)                   # first 997 cases by infection day
  step3_runs <- rbind(step3_runs, data.frame(
    run               = basename(step3_summary$output_prefix[i]),
    adult_file        = step3_summary$immunity_distribution_file[i],
    factor            = grid_adult_factor[match(step3_summary$immunity_distribution_file[i], step3_files)],
    num_cases         = length(age),
    outbreak          = length(age) >= min_outbreak_size,
    adult_share_first = mean(age[first] >= 18)))
}
write.csv(step3_runs, file.path(calib_dir, "step3_runs.csv"), row.names = FALSE)

# 3d. medians per factor (outbreak runs only) and interpolation at the target
step3_by_factor      <- aggregate(cbind(n_outbreak = outbreak) ~ factor, data = step3_runs, FUN = sum)
step3_by_factor$runs <- aggregate(outbreak ~ factor, data = step3_runs, FUN = length)$outbreak
step3_by_factor <- merge(step3_by_factor,
  aggregate(cbind(num_cases, adult_share_first) ~ factor, data = step3_runs[step3_runs$outbreak, ], FUN = median),
  all.x = TRUE)
step3_by_factor <- step3_by_factor[order(step3_by_factor$factor), ]
write.csv(step3_by_factor, file.path(calib_dir, "step3_by_factor.csv"), row.names = FALSE)
print(step3_by_factor, row.names = FALSE, digits = 3)

step3_factor_interp <- approx(step3_by_factor$adult_share_first, step3_by_factor$factor,
                              xout = target_adult_share, ties = mean)$y
if (is.na(step3_factor_interp)) stop("target adult share outside the simulated range: widen grid_adult_factor")
step3_factor <- round(step3_factor_interp / 0.05) * 0.05

# 3e. write the calibrated adult file
imm <- ifelse(adult_age >= 18, 1 - step3_factor * (1 - adult_imm_raw), adult_imm_raw)
step3_adult_file <- sub("\\.xml$", sprintf("_calib_f%03d.xml", as.integer(round(100 * step3_factor))), adult_raw_file)
writeLines(c(adult_header, "<immunity>",
             sprintf("  <data_source>%s</data_source>", adult_source),
             sprintf("  <data_manipulation>%s; susceptibility (1 - immunity) at ages 18+ multiplied by %.2f, calibrated so that %.1f%% of SC outbreak cases are adults (%s)</data_manipulation>",
                     adult_manip, step3_factor, 100 * target_adult_share, Sys.Date()),
             sprintf("  <age%d>%.4f</age%d>", adult_age, imm, adult_age), "</immunity>"), step3_adult_file)

step3_result <- data.frame(factor_interpolated = step3_factor_interp, factor = step3_factor,
                           adult_file = step3_adult_file, cluster_rate_used = step3_cluster_rate,
                           tracing_other_used = step3_tracing_other, project = basename(step3_project))
write.csv(step3_result, file.path(calib_dir, "step3_adult_factor.csv"), row.names = FALSE)
print(step3_result)
if (is.unsorted(na.omit(step3_by_factor$adult_share_first))) message("CHECK: adult share does not rise with the factor (noise): add seeds")
if (any(step3_by_factor$n_outbreak < min_outbreak_runs)) message("CHECK: some factors have < 10 outbreak runs: add seeds")


#############################################################################
## STEP 4: TRACING EFFICIENCY "OTHER"                                     ####
#############################################################################
# PURPOSE  Choose the share of non-household contacts that is traced
#          (tracing_efficiency_other), by simulating every value and scoring it
#          against the outbreak data. The school clusters are fixed by step 2b.
#          (Before 2026-11-24 this step also scanned the cluster rate.)
#          Score (lower = better) = sum of
#            score_curve    mean squared log error of the cumulative cases at the 47 DPH report dates
#            score_traced   squared log error of the share of the first 997 cases traced as contacts
#            score_final    squared log error of the median final size (day num_days) vs 997
#            score_outbreak (log(0.5 / p_outbreak))^2 if fewer than half of the runs take off, else 0
#          Settings with < 10 outbreak runs get no score.
# INPUTS   steps 1-3 (population, disease file, hesitancy and clusters, calibrated adult file)
# OUTPUTS  calib_dir/step4_runs.csv (one row per run),
#          calib_dir/step4_scores.csv (one row per tracing value, best first),
#          calib_dir/step4_best.csv (chosen tracing; cluster rate copied from step 2b)
#
# CHECK BEFORE CONTINUING
#   - the best row has >= 10 outbreak runs and a final size near 997
#     (tracing alone may not reach 997: then the school clusters from step 2b,
#     not tracing, decide the outbreak size; report the gap rather than force it)
#   - the chosen value is not at the edge of its grid; otherwise widen the grid
#   - compare the top 3-5 rows: if their scores are close (difference < ~0.2),
#     seed noise can decide the order; prefer values in the middle of that group
#     or rerun the top rows with more seeds

if (county_code != "SC") stop("steps 3-5 calibrate to the SC outbreak data: for other counties run steps 1, 2 and 2b only")
step1_result <- read.csv(file.path(calib_dir, "step1_population.csv"), stringsAsFactors = FALSE)
step2_result <- read.csv(file.path(calib_dir, "step2_disease.csv"), stringsAsFactors = FALSE)
step3_result <- read.csv(file.path(calib_dir, "step3_adult_factor.csv"), stringsAsFactors = FALSE)
if (file.exists(file.path(calib_dir, "step2b_best.csv"))) {
  step2b_best <- read.csv(file.path(calib_dir, "step2b_best.csv"))
} else {
  step2b_best <- data.frame(vaccine_hesitancy_uptake       = fixed$vaccine_hesitancy_uptake,
                            vaccine_hesitancy_num_clusters = step1_result$cluster_schools,
                            vaccine_hesitancy_cluster_rate = fixed$vaccine_hesitancy_cluster_rate)
}

# 4a. run STRIDE for every tracing value
step4_param <- sc_param_base
step4_param$num_rng_seeds                  <- num_seeds_tracing
step4_param$population_file                <- step1_result$pop_file
step4_param$disease_config_file            <- step2_result$disease_file
step4_param$immunity_distribution_file     <- step3_result$adult_file
step4_param$vaccine_hesitancy_uptake       <- step2b_best$vaccine_hesitancy_uptake
step4_param$vaccine_hesitancy_num_clusters <- step2b_best$vaccine_hesitancy_num_clusters
step4_param$vaccine_hesitancy_cluster_rate <- step2b_best$vaccine_hesitancy_cluster_rate
step4_param$tracing_efficiency_other       <- grid_tracing_other   # scanned
step4_numeric <- unlist(step4_param[sapply(step4_param, is.numeric)])
if (any(abs(step4_numeric) >= 1e5, na.rm = TRUE)) stop("a numeric parameter >= 1e5 would be written as 1e+05 (read as 1 by STRIDE)")
step4_design  <- .rstride$get_full_grid_exp_design(exp_param_list = step4_param, num_rng_seeds = step4_param$num_rng_seeds)
step4_project <- run_rStride(step4_design, '_calib_SC_step4_cluster_tracing', remove_run_output = FALSE,
                             num_parallel_workers = num_parallel_workers)

# 4b. per run: size, adult share, share traced, cumulative cases by report date, coverage
step4_summary <- read.csv(file.path(step4_project, paste0(basename(step4_project), "_summary.csv")), stringsAsFactors = FALSE)
cum_cols   <- paste0("cum_", format(obs$date, "%Y%m%d"))
step4_runs <- data.frame()
for (i in seq_len(nrow(step4_summary))) {
  run_dir <- step4_summary$output_prefix[i]
  log_txt <- readLines(file.path(run_dir, "event_log.txt"), warn = FALSE)
  cases   <- strsplit(log_txt[startsWith(log_txt, "[PRIM]") | startsWith(log_txt, "[TRAN]")], " ", fixed = TRUE)
  id      <- as.integer(sapply(cases, `[`, 2))
  age     <- as.integer(sapply(cases, `[`, 4))
  inf_day <- as.integer(sapply(cases, `[`, 7))
  onset   <- as.integer(sapply(cases, `[`, 11))                    # days since infection
  trace   <- strsplit(log_txt[startsWith(log_txt, "[TRACE]")], " ", fixed = TRUE)
  traced  <- unique(as.integer(sapply(trace, `[`, 2))[sapply(trace, `[`, 6) != "Index" & sapply(trace, `[`, 4) == "true"])
  first   <- head(order(inf_day), num_obs_cases)
  report_date <- as.Date(step4_summary$start_date[i]) + inf_day + onset + report_lag
  cum     <- findInterval(obs$date, sort(report_date))              # cases reported on or before each DPH date
  sus     <- read.csv(file.path(run_dir, "susceptibles_by_age.csv"))
  sus     <- sus[sus$age %in% 5:17, ]
  step4_runs <- rbind(step4_runs, data.frame(
    run                = basename(run_dir),
    cluster_rate       = step4_summary$vaccine_hesitancy_cluster_rate[i],
    tracing_other      = step4_summary$tracing_efficiency_other[i],
    num_cases          = length(id),
    outbreak           = length(id) >= min_outbreak_size,
    adult_share_first  = mean(age[first] >= 18),
    share_traced_first = mean(id[first] %in% traced),
    coverage_5_17      = sum(sus$immune) / sum(sus$immune + sus$susceptible),
    t(setNames(cum, cum_cols)), check.names = FALSE))
}
write.csv(step4_runs, file.path(calib_dir, "step4_runs.csv"), row.names = FALSE)

# 4c. medians per combination (outbreak runs only) and score
step4_all      <- aggregate(cbind(n_outbreak = outbreak) ~ cluster_rate + tracing_other, data = step4_runs, FUN = sum)
step4_all$runs <- aggregate(outbreak ~ cluster_rate + tracing_other, data = step4_runs, FUN = length)$outbreak
step4_cov <- aggregate(coverage_5_17 ~ cluster_rate + tracing_other, data = step4_runs, FUN = median)
step4_ok  <- step4_runs[step4_runs$outbreak, ]
step4_med <- aggregate(step4_ok[, c("num_cases", "adult_share_first", "share_traced_first", cum_cols)],
                       by = list(cluster_rate = step4_ok$cluster_rate, tracing_other = step4_ok$tracing_other),
                       FUN = median)
step4_scores <- merge(merge(step4_all, step4_cov), step4_med, all.x = TRUE)
step4_scores$p_outbreak     <- step4_scores$n_outbreak / step4_scores$runs
step4_scores$score_curve    <- rowMeans((log1p(as.matrix(step4_scores[, cum_cols])) -
                                         matrix(log1p(obs$cum_cases), nrow(step4_scores), length(cum_cols), byrow = TRUE))^2)
step4_scores$score_traced   <- (log(step4_scores$share_traced_first) - log(target_share_in_quarantine))^2
step4_scores$score_final    <- (log(step4_scores$num_cases) - log(target_final_size))^2
step4_scores$score_outbreak <- ifelse(step4_scores$p_outbreak < min_p_outbreak,
                                      log(min_p_outbreak / pmax(step4_scores$p_outbreak, 0.01))^2, 0)
step4_scores$score          <- ifelse(step4_scores$n_outbreak >= min_outbreak_runs,
                                      step4_scores$score_curve + step4_scores$score_traced +
                                      step4_scores$score_final + step4_scores$score_outbreak, NA)
step4_scores <- step4_scores[order(step4_scores$score, na.last = TRUE), ]
write.csv(step4_scores, file.path(calib_dir, "step4_scores.csv"), row.names = FALSE)
print(head(step4_scores[, c("cluster_rate", "tracing_other", "n_outbreak", "p_outbreak", "num_cases",
                            "adult_share_first", "share_traced_first", "coverage_5_17", "score_curve",
                            "score_traced", "score_final", "score_outbreak", "score")], 10),
      row.names = FALSE, digits = 3)

# 4d. best combination
if (is.na(step4_scores$score[1])) stop("no tracing value has >= 10 outbreak runs: add seeds or change the grid")
step4_best <- data.frame(cluster_rate = step4_scores$cluster_rate[1], tracing_other = step4_scores$tracing_other[1],
                         score = step4_scores$score[1], project = basename(step4_project))
write.csv(step4_best, file.path(calib_dir, "step4_best.csv"), row.names = FALSE)
print(step4_best)
if (step4_best$tracing_other %in% range(grid_tracing_other)) message("CHECK: tracing other is at the edge of its grid")


#############################################################################
## STEP 5: CHECK RUN AND CALIBRATED VALUES                                ####
#############################################################################
# PURPOSE  Run the chosen values with 100 seeds and compare the model with every
#          target; write the values to use in the simulation scripts.
# INPUTS   steps 1-4
# OUTPUTS  calib_dir/step5_runs.csv, calib_dir/step5_model_vs_targets.csv,
#          calib_dir/step5_cumulative_cases.pdf,
#          calib_dir/calibrated_parameters.csv  <- copy these into
#          rStride_measles_explore.R / rStride_measles_baseline.R /
#          rStride_measles_scenarios_common.R
#
# CHECK BEFORE USING THE VALUES
#   - n_outbreak >= 10 (otherwise the medians are unreliable)
#   - final size within a factor 1.5 of 997
#   - adult share within 2 percentage points of 8.7%; otherwise go back to
#     step 3 (it uses the step-4 values), then steps 4 and 5
#   - share traced close to 50%
#   - the cumulative curve: the model is not expected to reproduce the January
#     surge (accepted limitation); the plot shows how far it is off
#   - coverage ages 5-17 vs 88.8% (DPH 2025-26; 82% before step 2b was added, Nov 2026)

if (county_code != "SC") stop("steps 3-5 calibrate to the SC outbreak data: for other counties run steps 1, 2 and 2b only")
step1_result <- read.csv(file.path(calib_dir, "step1_population.csv"), stringsAsFactors = FALSE)
step2_result <- read.csv(file.path(calib_dir, "step2_disease.csv"), stringsAsFactors = FALSE)
step3_result <- read.csv(file.path(calib_dir, "step3_adult_factor.csv"), stringsAsFactors = FALSE)
step4_best   <- read.csv(file.path(calib_dir, "step4_best.csv"), stringsAsFactors = FALSE)
if (file.exists(file.path(calib_dir, "step2b_best.csv"))) {
  step2b_best <- read.csv(file.path(calib_dir, "step2b_best.csv"))
  step2b_source <- "step 2b"
} else {
  step2b_best <- data.frame(vaccine_hesitancy_uptake       = fixed$vaccine_hesitancy_uptake,
                            vaccine_hesitancy_num_clusters = step1_result$cluster_schools,
                            vaccine_hesitancy_cluster_rate = fixed$vaccine_hesitancy_cluster_rate)
  step2b_source <- "fixed"
}

# 5a. run STRIDE with the calibrated values
step5_param <- sc_param_base
step5_param$num_rng_seeds                  <- num_seeds_check
step5_param$population_file                <- step1_result$pop_file
step5_param$disease_config_file            <- step2_result$disease_file
step5_param$immunity_distribution_file     <- step3_result$adult_file
step5_param$vaccine_hesitancy_uptake       <- step2b_best$vaccine_hesitancy_uptake
step5_param$vaccine_hesitancy_num_clusters <- step2b_best$vaccine_hesitancy_num_clusters
step5_param$vaccine_hesitancy_cluster_rate <- step2b_best$vaccine_hesitancy_cluster_rate
step5_param$tracing_efficiency_other       <- step4_best$tracing_other
step5_numeric <- unlist(step5_param[sapply(step5_param, is.numeric)])
if (any(abs(step5_numeric) >= 1e5, na.rm = TRUE)) stop("a numeric parameter >= 1e5 would be written as 1e+05 (read as 1 by STRIDE)")
step5_design  <- .rstride$get_full_grid_exp_design(exp_param_list = step5_param, num_rng_seeds = step5_param$num_rng_seeds)
step5_project <- run_rStride(step5_design, '_calib_SC_step5_check', remove_run_output = FALSE,
                             num_parallel_workers = num_parallel_workers)

# 5b. per run: same measures as in step 4
step5_summary <- read.csv(file.path(step5_project, paste0(basename(step5_project), "_summary.csv")), stringsAsFactors = FALSE)
cum_cols   <- paste0("cum_", format(obs$date, "%Y%m%d"))
step5_runs <- data.frame()
for (i in seq_len(nrow(step5_summary))) {
  run_dir <- step5_summary$output_prefix[i]
  log_txt <- readLines(file.path(run_dir, "event_log.txt"), warn = FALSE)
  cases   <- strsplit(log_txt[startsWith(log_txt, "[PRIM]") | startsWith(log_txt, "[TRAN]")], " ", fixed = TRUE)
  id      <- as.integer(sapply(cases, `[`, 2))
  age     <- as.integer(sapply(cases, `[`, 4))
  inf_day <- as.integer(sapply(cases, `[`, 7))
  onset   <- as.integer(sapply(cases, `[`, 11))
  trace   <- strsplit(log_txt[startsWith(log_txt, "[TRACE]")], " ", fixed = TRUE)
  traced  <- unique(as.integer(sapply(trace, `[`, 2))[sapply(trace, `[`, 6) != "Index" & sapply(trace, `[`, 4) == "true"])
  first   <- head(order(inf_day), num_obs_cases)
  report_date <- as.Date(step5_summary$start_date[i]) + inf_day + onset + report_lag
  cum     <- findInterval(obs$date, sort(report_date))
  sus     <- read.csv(file.path(run_dir, "susceptibles_by_age.csv"))
  sus     <- sus[sus$age %in% 5:17, ]
  step5_runs <- rbind(step5_runs, data.frame(
    run                = basename(run_dir),
    num_cases          = length(id),
    outbreak           = length(id) >= min_outbreak_size,
    last_infection_day = max(inf_day),
    adult_share_first  = mean(age[first] >= 18),
    share_traced_first = mean(id[first] %in% traced),
    coverage_5_17      = sum(sus$immune) / sum(sus$immune + sus$susceptible),
    t(setNames(cum, cum_cols)), check.names = FALSE))
}
write.csv(step5_runs, file.path(calib_dir, "step5_runs.csv"), row.names = FALSE)

# 5c. model vs targets (medians of the outbreak runs)
step5_ok  <- step5_runs[step5_runs$outbreak, ]
step5_med <- sapply(step5_ok[, c("num_cases", "adult_share_first", "share_traced_first", cum_cols)], median)
step5_key_dates <- c("20251114", "20251216", "20260109", "20260127", "20260203", "20260303", "20260317")
step5_compare <- data.frame(
  target   = c(paste("cumulative cases", step5_key_dates), "final size (median of outbreak runs)",
               "outbreak runs (share of runs >= 200 cases)", "adult share (first 997 cases)",
               "share traced as contacts (first 997 cases)", "coverage ages 5-17"),
  observed = c(obs$cum_cases[match(step5_key_dates, format(obs$date, "%Y%m%d"))], target_final_size,
               min_p_outbreak, target_adult_share, target_share_in_quarantine, target_child_coverage),
  model    = c(step5_med[paste0("cum_", step5_key_dates)], step5_med["num_cases"], mean(step5_runs$outbreak),
               step5_med["adult_share_first"], step5_med["share_traced_first"], median(step5_runs$coverage_5_17)))
write.csv(step5_compare, file.path(calib_dir, "step5_model_vs_targets.csv"), row.names = FALSE)
print(step5_compare, row.names = FALSE, digits = 3)
cat(sprintf("[step 5] %d outbreak runs of %d; final size p10-p90 %.0f-%.0f; median last infection %s\n",
            nrow(step5_ok), nrow(step5_runs), quantile(step5_ok$num_cases, 0.1), quantile(step5_ok$num_cases, 0.9),
            as.Date(fixed$start_date) + median(step5_ok$last_infection_day)))
if (nrow(step5_ok) < min_outbreak_runs) message("CHECK: fewer than 10 outbreak runs")
if (abs(log(step5_med["num_cases"] / target_final_size)) > log(1.5)) message("CHECK: final size off by more than a factor 1.5")
if (abs(step5_med["adult_share_first"] - target_adult_share) > 0.02) message("CHECK: adult share off by > 2 points: rerun steps 3-5")

# 5d. cumulative curves: outbreak runs (grey), median (black), DPH (red)
pdf(file.path(calib_dir, "step5_cumulative_cases.pdf"), 8, 5)
matplot(obs$date, t(step5_ok[, cum_cols]), type = "l", lty = 1, col = "grey80", xlab = "report date",
        ylab = "cumulative cases", ylim = c(0, max(c(obs$cum_cases, unlist(step5_ok[, cum_cols])))))
lines(obs$date, step5_med[cum_cols], lwd = 2)
points(obs$date, obs$cum_cases, col = "red", pch = 16)
legend("topleft", c("outbreak runs", "model median", "DPH"), col = c("grey70", "black", "red"),
       lty = c(1, 1, NA), pch = c(NA, NA, 16), bty = "n")
dev.off()

# 5e. calibrated values for the simulation scripts
calibrated_parameters <- data.frame(
  parameter = c("population_file", "age_contact_matrix_file", "disease_config_file", "transmission_probability (r0 = 12)",
                "immunity_distribution_file", "adult susceptibility factor", "vaccine_distribution_file",
                "non-medical exemption rate (CDC SchoolVaxView 2025)", "cluster_school_share",
                "vaccine_hesitancy_rate", "vaccine_hesitancy_uptake", "vaccine_hesitancy_num_clusters",
                "vaccine_hesitancy_cluster_rate", "vaccine_hesitancy_cluster_uptake", "tracing_efficiency_other",
                "tracing_efficiency_household", "detection_probability", "case_finding_capacity",
                "start_date", "num_days", "contact_tracing_date"),
  value = c(step1_result$pop_file, contact_file, step2_result$disease_file, round(step2_result$transmission_probability, 4),
            step3_result$adult_file, step3_result$factor, child_file,
            county_inputs[[county_code]]$exemption_rate, cluster_school_share,
            fixed$vaccine_hesitancy_rate, step2b_best$vaccine_hesitancy_uptake, step2b_best$vaccine_hesitancy_num_clusters,
            step2b_best$vaccine_hesitancy_cluster_rate, fixed$vaccine_hesitancy_cluster_uptake, step4_best$tracing_other,
            fixed$tracing_efficiency_household, fixed$detection_probability, fixed$case_finding_capacity,
            fixed$start_date, fixed$num_days, fixed$contact_tracing_date),
  source = c("step 1", "input", "step 2", "step 2", "step 3", "step 3", "input",
             "input", "input (school data)",
             "fixed", step2b_source, step2b_source, step2b_source, "fixed (scenario parameter)", "step 4", "fixed", "fixed", "fixed",
             "fixed", "fixed", "fixed"))
write.csv(calibrated_parameters, file.path(calib_dir, "calibrated_parameters.csv"), row.names = FALSE)
print(calibrated_parameters, row.names = FALSE)


#############################################################################
## STEP 6: TRANSFER TO THE OTHER REGIONS (vaccine hesitancy uptake)       ####
#############################################################################
# PURPOSE  Set the hesitancy parameters of CT, ME, MN, OH and TX from the SC
#          calibration, with exemption policy acting through the uptake:
#            vaccine_hesitancy_rate          the region's ASPE value          (data)
#            cluster_school_share            the region's school data (CT) or
#                                            the SC value (others)            (data / SC)
#            vaccine_hesitancy_cluster_rate  the SC value from step 2b        (SC)
#            vaccine_hesitancy_cluster_uptake 0, as in SC (scenario parameter) (SC)
#            vaccine_hesitancy_uptake        chosen so that the hesitancy gap at
#                                            ages 5-6 equals the region's
#                                            non-medical exemption rate       (same rule as SC)
#          The same rule as in step 2b: children left unvaccinated by hesitancy
#          = children with a non-medical exemption. A stricter policy (fewer
#          exemptions) therefore gives a higher uptake.
#          For every region: 1-day runs for every uptake in grid_uptake_transfer;
#          the uptake is interpolated where the median gap equals the exemption rate.
#          If the clusters alone already leave more children unvaccinated than the
#          exemption rate (gap at uptake 1 > exemption rate), the uptake is set to
#          1 and flagged: the clusters, not the policy, then set the gap.
#          The regions' coverage data are printed as a check, not used as targets.
# INPUTS   SC step 2b (calib_dir/step2b_best.csv) and step 2 (disease file, only
#          needed to start STRIDE); for every region: its _hhdense population
#          (run step 1 with that county_code once), contact matrix, child and adult
#          immunity files, ASPE rate and exemption rate (county_inputs, STEP 0)
#          and its school targets (school_coverage_targets_2025-26.csv)
# OUTPUTS  calib_dir/step6_runs.csv (one row per run),
#          calib_dir/step6_region_parameters.csv (one row per region, SC included):
#          the hesitancy values to use in the region and policy scenario scripts
#
# CHECK BEFORE USING THE VALUES
#   - uptake is ordered by policy: medical-only states (CT, ME) highest
#   - flagged regions (uptake set to 1): their gap is set by the clusters
#   - coverage at the data ages next to the regions' data: large differences
#     point to the ASPE rate, the child immunity file or the cluster pattern

if (county_code != "SC") stop("step 6 starts from the SC calibration: set county_code <- \"SC\" and rerun STEP 0")
step2_result <- read.csv(file.path(calib_dir, "step2_disease.csv"), stringsAsFactors = FALSE)
step2b_best  <- read.csv(file.path(calib_dir, "step2b_best.csv"), stringsAsFactors = FALSE)
all_targets  <- read.csv(school_targets_file, stringsAsFactors = FALSE)

step6_runs    <- data.frame()
step6_regions <- data.frame(region = "SC", hesitancy_rate = step2b_best$vaccine_hesitancy_rate,
                            exemption_rate = step2b_best$exemption_rate,
                            uptake = step2b_best$vaccine_hesitancy_uptake,
                            cluster_school_share = step2b_best$cluster_school_share,
                            num_clusters = step2b_best$vaccine_hesitancy_num_clusters,
                            cluster_rate = step2b_best$vaccine_hesitancy_cluster_rate,
                            cluster_uptake = step2b_best$vaccine_hesitancy_cluster_uptake,
                            gap_at_uptake_1 = NA, flag = "calibrated (step 2b)")

for (region in regions_transfer) {
  # 6a. inputs of this region
  r_pop      <- sub("\\.csv$", "_hhdense.csv", file.path("data/population_6region", county_inputs[[region]]$pop))
  if (!file.exists(r_pop)) stop(r_pop, " not found: run step 1 with county_code <- \"", region, "\" first")
  r_child    <- sprintf("data/immunity_6region/immunity_measles_child_%s.xml", region)
  r_adult    <- sprintf("data/immunity_6region/immunity_measles_adult_%s.xml", region)
  r_targets  <- all_targets[all_targets$county_code == region, ]
  r_share    <- r_targets$value[r_targets$target == "cluster_school_share"]
  r_coverage <- r_targets[r_targets$target == "coverage", ]
  r_exempt   <- county_inputs[[region]]$exemption_rate
  r_rate     <- county_inputs[[region]]$hesitancy_rate

  # schools and cluster schools in the region's population; immunity at ages 5-6 without hesitancy
  r_pop_ages <- read.csv(r_pop, colClasses = c("integer", "NULL", "integer", "NULL", "NULL", "NULL"))
  r_pop_ages <- r_pop_ages[!is.na(r_pop_ages$school_id) & r_pop_ages$school_id != 0, ]
  r_schools  <- length(unique(r_pop_ages$school_id))
  r_clusters <- round(r_share * r_schools)
  child_txt  <- readLines(r_child, warn = FALSE)
  child_m    <- regmatches(child_txt, regexec("<age([0-9]+)>([-0-9.eE]+)</age[0-9]+>", child_txt))
  child_m    <- child_m[lengths(child_m) == 3]
  child_imm  <- setNames(as.numeric(sapply(child_m, `[`, 3)), sapply(child_m, `[`, 2))
  r_base     <- mean(child_imm[as.character(r_pop_ages$age[r_pop_ages$age %in% exemption_ages])])
  cat(sprintf("[step 6] %s: ASPE %.3f, exemption rate %.3f, %d schools -> %d cluster schools (share %.3f), immunity at ages %s without hesitancy %.4f\n",
              region, r_rate, r_exempt, r_schools, r_clusters, r_share, paste(exemption_ages, collapse = "-"), r_base))

  # 6b. run STRIDE for 1 day for every uptake
  step6_param <- sc_param_base
  step6_param$num_rng_seeds                    <- num_seeds_transfer
  step6_param$num_days                         <- 1
  step6_param$population_file                  <- r_pop
  step6_param$age_contact_matrix_file          <- file.path("data/contacts_usa", county_inputs[[region]]$contacts)
  step6_param$disease_config_file              <- step2_result$disease_file   # only needed to start STRIDE
  step6_param$immunity_distribution_file       <- r_adult
  step6_param$vaccine_distribution_file        <- r_child
  step6_param$tracing_efficiency_other         <- start_tracing_other
  step6_param$vaccine_hesitancy_rate           <- r_rate
  step6_param$vaccine_hesitancy_num_clusters   <- r_clusters
  step6_param$vaccine_hesitancy_cluster_rate   <- step2b_best$vaccine_hesitancy_cluster_rate
  step6_param$vaccine_hesitancy_cluster_uptake <- step2b_best$vaccine_hesitancy_cluster_uptake
  step6_param$vaccine_hesitancy_uptake         <- grid_uptake_transfer         # scanned
  if (r_clusters == 0) step6_param$vaccine_hesitancy_cluster_pooltype <- "None"
  set.seed(2026)
  step6_design  <- .rstride$get_full_grid_exp_design(exp_param_list = step6_param, num_rng_seeds = step6_param$num_rng_seeds)
  step6_project <- run_rStride(step6_design, paste0('_calib_', region, '_step6_transfer'),
                               remove_run_output = FALSE, num_parallel_workers = num_parallel_workers)

  # 6c. per run: hesitancy gap at ages 5-6 and coverage at the ages of the region's data
  step6_summary <- read.csv(file.path(step6_project, paste0(basename(step6_project), "_summary.csv")), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(step6_summary))) {
    snap_file <- file.path(step6_summary$output_prefix[i], "population_snapshot.csv")
    snap <- read.csv(snap_file, colClasses = c("integer", "NULL", "integer", "NULL", "NULL", "NULL", "character"))
    snap <- snap[snap$school_id != 0, ]
    snap$immune <- snap$immunity_status == "immune"
    immunity_5_6 <- mean(snap$immune[snap$age %in% exemption_ages])
    run_cov <- numeric(0)
    for (k in seq_len(nrow(r_coverage))) {
      run_cov <- c(run_cov, mean(snap$immune[snap$age >= r_coverage$age_min[k] & snap$age <= r_coverage$age_max[k]]))
    }
    step6_runs <- rbind(step6_runs, data.frame(
      region = region, run = basename(step6_summary$output_prefix[i]),
      uptake = step6_summary$vaccine_hesitancy_uptake[i],
      immunity_5_6 = immunity_5_6, hesitancy_gap = r_base - immunity_5_6,
      coverage_check = paste(sprintf("ages %d-%d: model %.3f / data %.3f", r_coverage$age_min, r_coverage$age_max,
                                     run_cov, r_coverage$value), collapse = "; ")))
    if (step2b_delete_snapshots) file.remove(snap_file)
  }

  # 6d. uptake where the median gap equals the exemption rate
  r_runs <- step6_runs[step6_runs$region == region, ]
  r_gap  <- aggregate(hesitancy_gap ~ uptake, data = r_runs, FUN = median)
  print(r_gap, row.names = FALSE, digits = 3)
  gap_at_1 <- r_gap$hesitancy_gap[r_gap$uptake == max(r_gap$uptake)]
  gap_at_0 <- r_gap$hesitancy_gap[r_gap$uptake == min(r_gap$uptake)]
  if (r_exempt <= gap_at_1) {
    r_uptake <- 1;  r_flag <- "uptake 1: the clusters alone leave more children unvaccinated than the exemption rate"
  } else if (r_exempt >= gap_at_0) {
    r_uptake <- 0;  r_flag <- "uptake 0: even no uptake leaves fewer children unvaccinated than the exemption rate"
  } else {
    r_uptake <- approx(r_gap$hesitancy_gap, r_gap$uptake, xout = r_exempt, ties = mean)$y
    r_flag   <- "interpolated"
  }
  step6_regions <- rbind(step6_regions, data.frame(region = region, hesitancy_rate = r_rate, exemption_rate = r_exempt,
                                                   uptake = round(r_uptake, 3), cluster_school_share = r_share,
                                                   num_clusters = r_clusters,
                                                   cluster_rate = step2b_best$vaccine_hesitancy_cluster_rate,
                                                   cluster_uptake = step2b_best$vaccine_hesitancy_cluster_uptake,
                                                   gap_at_uptake_1 = round(gap_at_1, 4), flag = r_flag))
  cat(sprintf("[step 6] %s: uptake %.3f (%s)\n", region, r_uptake, r_flag))
  print(unique(r_runs$coverage_check[r_runs$uptake == r_runs$uptake[which.min(abs(r_runs$uptake - r_uptake))]])[1])
}

write.csv(step6_runs, file.path(calib_dir, "step6_runs.csv"), row.names = FALSE)
write.csv(step6_regions, file.path(calib_dir, "step6_region_parameters.csv"), row.names = FALSE)
print(step6_regions, row.names = FALSE)
