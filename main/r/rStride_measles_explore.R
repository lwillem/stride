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
#  Copyright 2018, Willem L, Kuylen E & Broeckhove J
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
dir_postfix <- '_expl'

##################################
## DESIGN OF EXPERIMENTS        ##
##################################

# uncomment the following line to inspect the config xml tags
#names(xmlToList('./config/run_default.xml'))

# add default parameters and values
# exp_param_list <- get_measles_default_param()

# set the number of realisations per configuration set
num_seeds<- 5

################################################ #
## GENERATE DESIGN OF EXPERIMENT GRID         ####
################################################ #

# # add parameters and values to combine in a full-factorial grid
# exp_design <- expand.grid(r0                            = seq(12,14,2),
#                           num_days                      = c(40,50),
#                           rng_seed                      = seq(num_seeds),
#                           track_index_case              = 'false',
#                           contact_log_level             = "Transmissions",
#                           seeding_rate                  = 0.00002,
#                           disease_config_file           = "data/disease_measles_adaptive_behavior.xml",
#                           population_file               = "data/pop_flanders600.csv",
#                           age_contact_matrix_file       = "data/contact_matrix_flanders_subpop.xml",
#                           adaptive_symptomatic_behavior = 'true',
#                           stringsAsFactors = F)

# get grid-based design of experiments
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

# add a unique seed for each run
set.seed(125)
exp_design$rng_seed <- sample(1e4,nrow(exp_design))

##################################
## RUN rSTRIDE                  ##
##################################
project_dir <- run_rStride(exp_design,dir_postfix)


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

inspect_transmission_data(project_dir)
