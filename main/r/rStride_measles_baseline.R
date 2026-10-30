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
#
#  Copyright 2026, Manansala R, Willem L
#############################################################################
#
# Measles scenario: Baseline: each county's own exemption policy, baseline clustering (C1).
# See the Claude Doc "Measles-Stride: baseline inputs and exemption-policy scenarios by county".
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right.
#
# E.g.: path/to/stride $ ./bin/rStride_measles_baseline.R
#
#############################################################################

# Clear work environment
rm(list=ls())

# Load rStride and the county helpers (input files per county, run loop)
source('./bin/rStride_measles_scenarios_common.R')

# scenario name (used in the output folder names)
scenario_name <- "baseline"

##################################
## DEFAULT PARAMETERS        ##
##################################

# get default parameters and values
exp_param_list <- get_default_param()

################################################ #
## GENERAL PARAMETERS (ALL COUNTIES)        ####
################################################ #

# simulation horizon and start: 1 year from the start of the school year
exp_param_list$num_days   <- 365
exp_param_list$start_date <- "2025-09-01"

# measles (transmission probability from the county-specific R0 fit, see disease_config_file per county)
exp_param_list$r0 <- 12

# hospital settings
exp_param_list$hosp_probability_factor  <- 1
exp_param_list$hospital_length_of_stay  <- 14  #TODO: set real value
exp_param_list$hospital_probability_age <- paste(0.18,0.06,0.12,sep=',')
exp_param_list$hospital_category_age    <- paste(0,5,20,sep=',')
exp_param_list$hospital_mean_delay_age  <- paste(2,2,5,sep=',')

# calendar
exp_param_list$holidays_file <- "data/calendar_USA_2023_2026_measles.csv"

# Combine the two age-specific contact probabilities of a candidate pair by their mean.
exp_param_list$contact_probability_rule <- "Mean"

# population-wide immunity snapshot at t = 0 (input for MeaslesClustering.R)
exp_param_list$output_pop_snapshot <- TRUE

# initial conditions
exp_param_list$num_infected_seeds <- 1
exp_param_list$seeding_age_min    <- 1   # infants <1y cannot be selected as index case
exp_param_list$seeding_age_max    <- 17

# log level(s). Options: Incidence, Transmissions, Contacts, None
exp_param_list$event_log_level <- "Transmissions"

# reference data file(s)
exp_param_list$reference_hospital_data_file <- NA
exp_param_list$reference_serology_data_file <- NA

# immunity (adults: data/immunity_6region/immunity_measles_adult_<ST>_calib_f035.xml, set per county)
exp_param_list$immunity_profile          <- "AgeDependent"
exp_param_list$immunity_link_probability <- 0

# virtual survey
exp_param_list$contact_survey_dates <- exp_param_list$start_date

# vaccination (children: data/immunity_6region/immunity_measles_child_<ST>.xml, set per county)
exp_param_list$vaccine_profile                    <- "AgeDependent"
exp_param_list$vaccine_min_age                    <- 0
exp_param_list$vaccine_max_age                    <- 17
exp_param_list$mass_immunize                      <- TRUE
exp_param_list$vaccine_hesitancy_cluster_pooltype <- "School"

# public health response (Spartanburg calibration)
exp_param_list$contact_tracing_date         <- as.character(as.Date(exp_param_list$start_date) + 17) # same delay after seeding as in SC (2025-09-15 -> 2025-10-02)
exp_param_list$detection_probability        <- 0.9
exp_param_list$tracing_efficiency_household <- 0.95
exp_param_list$tracing_efficiency_other     <- 0.4
exp_param_list$delay_isolation_index        <- 3
exp_param_list$delay_contact_tracing        <- 1
exp_param_list$case_finding_capacity        <- 10000  # keep < 1e5: R writes 100000 as "1e+05", which C++ reads as 1
exp_param_list$test_false_negative          <- 0
exp_param_list$isolation_duration_index     <- 7     # days
exp_param_list$isolation_duration_contact   <- 21    # days (measles quarantine)

################################################ #
## SCENARIO PARAMETERS BY COUNTY            ####
################################################ #
# one row per county; column names are rSTRIDE parameter names
# vaccine_hesitancy_rate: ASPE county "estimated hesitant" (COVID-19, Household Pulse Survey; data.cdc.gov q9mh-h2tw)
# own policy: CT, ME = P1 medical only (0.8 / 0.5); MN, TX, SC = P2 notarized (0.5 / 0); OH = P3 easy statement (0.3 / 0)
# C1 baseline clustering: PROVISIONAL num_clusters (SC: 7 of 90 schools scaled to the school count;
#    CT from school data), cluster rate 0.45; to be tuned to county K MMR coverage
# num_rng_seeds: 100 for all counties (same precision for P(outbreak) across counties)
# remove rows to run fewer counties
county_param <- rbind(
  data.frame(county = "CT", vaccine_hesitancy_rate = 0.05 , vaccine_hesitancy_uptake = 0.8, vaccine_hesitancy_cluster_uptake = 0.5, vaccine_hesitancy_num_clusters =  6, vaccine_hesitancy_cluster_rate = 0.45, num_rng_seeds = 100),  # own policy: P1 medical only
  data.frame(county = "ME", vaccine_hesitancy_rate = 0.105, vaccine_hesitancy_uptake = 0.8, vaccine_hesitancy_cluster_uptake = 0.5, vaccine_hesitancy_num_clusters =  1, vaccine_hesitancy_cluster_rate = 0.45, num_rng_seeds = 100),  # own policy: P1 medical only
  data.frame(county = "MN", vaccine_hesitancy_rate = 0.062, vaccine_hesitancy_uptake = 0.5, vaccine_hesitancy_cluster_uptake = 0, vaccine_hesitancy_num_clusters = 31, vaccine_hesitancy_cluster_rate = 0.45, num_rng_seeds = 100),  # own policy: P2 notarized exemption
  data.frame(county = "OH", vaccine_hesitancy_rate = 0.122, vaccine_hesitancy_uptake = 0.3, vaccine_hesitancy_cluster_uptake = 0, vaccine_hesitancy_num_clusters = 30, vaccine_hesitancy_cluster_rate = 0.45, num_rng_seeds = 100),  # own policy: P3 easy statement
  data.frame(county = "TX", vaccine_hesitancy_rate = 0.117, vaccine_hesitancy_uptake = 0.5, vaccine_hesitancy_cluster_uptake = 0, vaccine_hesitancy_num_clusters =  1, vaccine_hesitancy_cluster_rate = 0.45, num_rng_seeds = 100),  # own policy: P2 notarized exemption
  data.frame(county = "SC", vaccine_hesitancy_rate = 0.15 , vaccine_hesitancy_uptake = 0.5, vaccine_hesitancy_cluster_uptake = 0, vaccine_hesitancy_num_clusters =  7, vaccine_hesitancy_cluster_rate = 0.45, num_rng_seeds = 100)  # own policy: P2 notarized exemption
)

##################################
## RUN rSTRIDE                  ##
##################################
project_dirs <- run_measles_counties(exp_param_list,
                                     county_param,
                                     scenario_name,
                                     num_parallel_workers = 3)

## load existing project
# project_dir <- smd_file_path('sim_output','<run_tag>')
