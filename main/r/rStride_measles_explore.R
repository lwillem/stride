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
exp_param_list$num_rng_seeds <- 2

# simulation horizon and start
exp_param_list$num_days   <- 40
exp_param_list$start_date <- "2023-01-01"

# case study for measles
exp_param_list$r0 <- 12
exp_param_list$disease_config_file <- "data/disease_measles_usa.xml"

# hospital settings
exp_param_list$hosp_probability_factor <- 0.1 #TODO: set real value
exp_param_list$hospital_length_of_stay <- 14  #TODO: set real value  
exp_param_list$hospital_probability_age       #TODO: set age-specific adjustment of 'hosp_probability_factor'
exp_param_list$hospital_category_age          #FYI: age categories for 'hosp_probability_factor' adjustments

# USA population
exp_param_list$age_contact_matrix_file <- "data/contact_matrix_usa_conditional.xml"
exp_param_list$holidays_file           <- "data/calendar_dane_2023_2026_measles.csv"
exp_param_list$population_file         <- "data/pop_usa_wisconsin_dane474k_c1000.csv"

# initial conditions
exp_param_list$num_infected_seeds <- 20
exp_param_list$seeding_age_min    <- 1   # infants <1y cannot be selected as index case
exp_param_list$seeding_age_max    <- 99

# log level(s). Options: Incidence, Transmissions, Contacts, None
exp_param_list$event_log_level <- "Transmissions"

# reference data file(s)
exp_param_list$reference_hospital_data_file <- NA
exp_param_list$reference_serology_data_file <- NA

# # immunity profile = starting condition (time consuming!)
 exp_param_list$immunity_profile <- "AgeDependent"
 exp_param_list$immunity_distribution_file <- "data/immunity_measles_dummy.xml" # "data/immunity_measles_belgium.xml" #c("data/immunity_measles_belgium.xml","data/immunity_measles_belgium_dummy.xml")
 exp_param_list$immunity_link_probability <- 0 # immunity is distributed by household, this is the chance of continuing to immunize the next shuffled household member instead of jumping to a new random household

# virtual survey for immunity levels
exp_param_list$num_participants_survey
exp_param_list$contact_survey_dates <- exp_param_list$start_date # if log level is not "Contacts", only health data is obtained

################################################ #
## GENERATE DESIGN OF EXPERIMENT GRID         ####
################################################ #

# add parameters and values to combine in a full-factorial grid
exp_design <- expand.grid(r0                            = seq(12,14,2),
                          num_days                      = c(40,50),
                          rng_seed                      = seq(num_seeds),
                          age_contact_matrix_file       = "sim_output/20260710_111118_WI-Dane_conditional_social_contacts/contact_matrix_usa_conditional.xml",
                          disease_config_file           = "data/disease_measles.xml",
                          holidays_file                 = "data/calendar_dane_2023_2026_measles.csv",
                          immunity_profiles             = "AgeDependent",
                          immunity_distribution_file    = "data/immunity_measles_WI.xml", 
                          # immunity_rate                 = 0.8,
                          # num_participants_survey       = 5000,
                          population_file               = "sim_output/20260710_111118_WI-Dane_conditional_social_contacts/20260710_111118_population_WI-Dane.csv",
                          seeding_age_max               = 99,
                          seeding_age_min               = 1,
                          seeding_rate                  = 0.00002,
                          num_infected_seeds            = 9,
                          start_date                    = "2023-02-01",
                          track_index_case              = 'false',
                          event_log_level               = "Transmissions",
                          stride_log_level              = "true",
                          adaptive_symptomatic_behavior = 'true',
                          stringsAsFactors = F)

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
                           remove_run_output = FALSE)


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

