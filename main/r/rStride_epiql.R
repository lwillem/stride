#!/usr/bin/env Rscript
############################################################################ #
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
#  Copyright 2025, Willem L.
############################################################################ #
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_epiql.R 
#
############################################################################ #

# Clear work environment
rm(list=ls())

# Load rStride
source('./bin/rstride/rStride.R')

# Load default parameter configurations
source('./bin/rStride_covid19_default_param.R')

# set directory postfix (optional)
dir_postfix <- '_epiql'

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# add default parameters and values
exp_param_list <- get_covid19_default_param()

# change population file
exp_param_list$population_file <- 'data/pop_belgium600k_c500_teachers_censushh.csv'
#exp_param_list$population_file <- 'data/pop_belgium1000k_c500_teachers_censushh.csv'

# change parameters and values to combine in a full-factorial grid
exp_param_list$num_days <- 240
exp_param_list$num_parallel_workers <- 8
exp_param_list$event_log_level <- c("Transmissions")
exp_param_list$num_seeds <- 2
 
names(exp_param_list)

exp_param_list$r0 <- 1.4
exp_param_list$hosp_probability_factor <- 0
exp_param_list$num_infected_seeds <- 500

exp_param_list$detection_probability <- c(0, 0.9)   # if symptomatic!
exp_param_list$tracing_efficiency_household <- 0.9
exp_param_list$tracing_efficiency_other <- 0.9
exp_param_list$test_false_negative   = 0 
exp_param_list$contact_tracing_date <- exp_param_list$start_date
exp_param_list$is_isolated_from_household <- 1

exp_param_list <- exp_param_list[!grepl('distancing',names(exp_param_list))]
exp_param_list$holidays_file <- 'data/holidays_none.csv'

# check period
range(as.Date(exp_param_list$start_date), as.Date(exp_param_list$start_date) + exp_param_list$num_days)

# vaccine

################################################ #
## GENERATE DESIGN OF EXPERIMENT GRID         ####
################################################ #

# get grid-based design of experiments
exp_design <- .rstride$get_full_grid_exp_design(exp_param_list = exp_param_list,
                                                num_seeds      = exp_param_list$num_seeds)
dim(exp_design)

################################## #
## RUN rSTRIDE                  ####
################################## #
project_dir <- run_rStride(exp_design               = exp_design,
                           dir_postfix              = dir_postfix,
                           get_tracing_rdata        = TRUE,
                           num_parallel_workers     = exp_param_list$num_parallel_workers,
                           remove_run_output        = FALSE)


############################# #
## INPUT-OUTPUT BEHAVIOR   ####
############################# #
inspect_summary(project_dir)


############################# #
## SURVEY PARTICIPANT DATA ####
############################# #
inspect_participant_data(project_dir)


########################################### #
## PARAMETER ESTIMATION (optional)       ####
########################################### #
#estimate_parameters(project_dir)


############################# #
## INCIDENCE DATA          ####
############################# #
inspect_incidence_data(project_dir)


############################# #
## PREVALENCE              ####
############################# #
inspect_prevalence_data(project_dir)


############################# #
## TRANSMISSION            ####
############################# #
inspect_transmission_dynamics(project_dir)
 

############################# #
## CONTACT TRACING         ####
############################# #
inspect_tracing_data(project_dir)

############################# #
## CONTACT SURVEY          ####
############################# #
#inspect_contact_data(project_dir)


 
