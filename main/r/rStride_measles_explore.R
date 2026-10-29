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
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_explore.R 
#
#############################################################################

# Clear work environment
rm(list=ls())

# Load rStride
source('./bin/rstride/rStride.R')

# set directory postfix (optional)
dir_postfix <- '_expl_measles'

##################################
## DEFAULT PARAMETERS        ##
##################################

# get default parameters and values
exp_param_list <- get_default_param()

################################################ #
## UPDATE PARAMETERS VALUES         ####
################################################ #

# number of stochastic realisations
exp_param_list$num_rng_seeds <- 20

# simulation horizon and start
# calibrated to the Spartanburg 2025-26 outbreak (declared 2025-10-02, last case 2026-03-15)
exp_param_list$num_days   <- 240
exp_param_list$start_date <- "2025-09-15"

# case study for measles
exp_param_list$r0 <- 12
exp_param_list$disease_config_file <- "data/OH_20261022_232515_r0_disease_measles_usa.xml"

# hospital settings
# age-specific values taken from measles_usa_rm (2b927f4, "fixed age specific hosp prob")
exp_param_list$hosp_probability_factor  <- 1
exp_param_list$hospital_length_of_stay  <- 14  #TODO: set real value
exp_param_list$hospital_probability_age <- paste(0.18,0.06,0.12,sep=',')
exp_param_list$hospital_category_age    <- paste(0,5,20,sep=',')
exp_param_list$hospital_mean_delay_age  <- paste(2,2,5,sep=',')

# USA population
# household ids renumbered 1..N in the order of the original ids (same results, ~200x faster);
# original ids: data/population_6region/pop_usa_sc_spartanburg_c1000_hhdense_household_map.csv
exp_param_list$population_file         <- "data/population_6region/pop_usa_oh_cuyahoga_c1000_hhdense.csv"
exp_param_list$age_contact_matrix_file <- "data/contacts_usa/contact_matrix_usa_oh_cuyahoga_c1000.xml"
exp_param_list$holidays_file           <- "data/calendar_USA_2023_2026_measles.csv"

# Combine the two age-specific contact probabilities of a candidate pair by their mean.
exp_param_list$contact_probability_rule <- "Mean"

# Write the population-wide immunity/susceptibility snapshot at t = 0. This is the input
# MeaslesClustering.R reads; without it that analysis has nothing to work from. Off by
# default because it costs roughly 16% of run time and writes three files per experiment.
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

# immunity profile = starting condition (time consuming!)
exp_param_list$immunity_profile <- "AgeDependent"
# note: immunity_measles_WI.xml and the data/immunity_6region/ files used on
# measles_usa_rm are absent from both branches; only the dummy profile is in the tree.
# adult susceptibility x 0.35 (calibrated: ~10% of outbreak cases are adults)
exp_param_list$immunity_distribution_file <- "data/immunity_6region/immunity_measles_adult_OH.xml"
exp_param_list$immunity_link_probability <- 0 # immunity is distributed by household, this is the chance of continuing to immunize the next shuffled household member instead of jumping to a new random household
# NOTE: 0 selects the independent per-age sampler since 2026-10-06 (f57bfd4). Before that, 0 over-sampled
# people in small households, so earlier outputs of this script are not reproduced in distribution.

# virtual survey for immunity levels
exp_param_list$num_participants_survey
exp_param_list$contact_survey_dates <- exp_param_list$start_date # if log level is not "Contacts", only health data is obtained

# vaccine rate
# exp_param_list$vaccine_link_probability <- 0
# exp_param_list$vaccine_profile <- "Random"
# exp_param_list$vaccine_rate <- c(0.8)
# exp_param_list$vaccine_min_age <- 0
# exp_param_list$vaccine_max_age <- 17

# Vaccination with vaccine hesitancy (ImmunitySeeder). Note: vaccine_min_age/max_age and
# mass_immunize are not used with the AgeDependent profile (ages come from the file).
# exp_param_list$vaccine_link_probability   <- 0
exp_param_list$vaccine_profile            <- "AgeDependent" ## MAKE NEW VAC PROF HESITANCY??
exp_param_list$vaccine_distribution_file  <- "data/immunity_6region/immunity_measles_child_OH.xml"
exp_param_list$vaccine_min_age            <- 0
exp_param_list$vaccine_max_age            <- 17
exp_param_list$mass_immunize              <- TRUE
exp_param_list$vaccine_hesitancy_rate     <- 0.15   # share of ALL households that is hesitant (CDC/ASPE COVID-19 adult hesitancy estimate)
exp_param_list$vaccine_hesitancy_uptake   <- 0.5    # relative uptake in hesitant households (0 = refuse all)
# cluster part of the hesitant households around a few schools (low-coverage schools)
exp_param_list$vaccine_hesitancy_cluster_pooltype <- "School"  # "None" = spread at random
exp_param_list$vaccine_hesitancy_num_clusters     <- 7         # number of schools, drawn at random per run
exp_param_list$vaccine_hesitancy_cluster_rate     <- 0.45      # share of households linked to those schools that is hesitant
exp_param_list$vaccine_hesitancy_cluster_uptake   <- 0         # uptake in clustered hesitant households

# household clustering
# exp_param_list$household_clustering_date <- "2023-01-01"
# exp_param_list$household_clustering_ratio <- 4/7
# exp_param_list$household_clustering_delay <- 0
# exp_param_list$cnt_intensity_householdCluster <- c(0, 1)

# contact reduction when symptomatic: read from the DISEASE file
# (disease.sympt_cnt_reduction_*), default 1 = stay home from school/work/community.
# Setting them here has no effect.
# exp_param_list$sympt_cnt_reduction_workplace_school <- 1
# exp_param_list$sympt_cnt_reduction_community <- 1

# public health response: isolation of detected cases, tracing and quarantine of contacts
exp_param_list$contact_tracing_date         <- "2025-10-02" # outbreak declared
exp_param_list$detection_probability        <- 0.9   # share of cases detected (index cases)
exp_param_list$tracing_efficiency_household <- 0.95  # share of household contacts traced
exp_param_list$tracing_efficiency_other     <- 0.4   # share of other contacts traced (calibrated)
exp_param_list$delay_isolation_index        <- 3     # days after symptom onset (~rash) until isolation/tracing
exp_param_list$delay_contact_tracing        <- 1     # days until traced contacts are quarantined
exp_param_list$case_finding_capacity        <- 10000  # index cases traced per day; keep < 1e5: R writes 100000 as "1e+05", which C++ reads as 1 (0 would also stop after 1)
exp_param_list$test_false_negative          <- 0
exp_param_list$isolation_duration_index     <- 7     # days
exp_param_list$isolation_duration_contact   <- 21    # days (measles quarantine)

################################################ #
## GENERATE DESIGN OF EXPERIMENT GRID         ####
################################################ #

# get grid-based design of experiments
exp_design <- .rstride$get_full_grid_exp_design(exp_param_list = exp_param_list,
                                                num_rng_seeds  = exp_param_list$num_rng_seeds)
dim(exp_design)

# tmp fix:
if(any(is.na(exp_param_list$immunity_distribution_file))){
  exp_param_list$immunity_distribution_file <- NULL
}


##################################
## RUN rSTRIDE                  ##
##################################
project_dir <- run_rStride(exp_design,
                           dir_postfix,
                           remove_run_output = FALSE,   ## THIS NEEDS TO STAY 'FALSE' TO GET HOUSEHOLD AGG
                           num_parallel_workers = 3)

## load existing project
# project_dir <- smd_file_path('sim_output','20260817_172151_expl_measles')

#####################################
## EXPLORE INPUT-OUTPUT BEHAVIOR   ##
#####################################
inspect_summary(project_dir)


#####################################
## EXPLORE SURVEY PARTICIPANT DATA ##
#####################################
inspect_participant_data(project_dir)


##################################
## EXPLORE TRANSMISSION         ##
##################################
inspect_transmission_dynamics(project_dir)


############################# #
## INCIDENCE DATA          ####
############################# #
inspect_incidence_data(project_dir)


############################# #
## PREVALENCE              ####
############################# #
inspect_prevalence_data(project_dir)


# TMP: explore burden of disease
project_summary    <- .rstride$load_project_summary(project_dir)
data_incidence_all <- .rstride$load_aggregated_output(project_dir,'data_incidence')

