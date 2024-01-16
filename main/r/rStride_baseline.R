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
#  Copyright 2024, Willem L.
############################################################################ #
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_baseline.R 
#
############################################################################ #

# Clear work environment
rm(list=ls())

# Load rStride
source('./bin/rstride/rStride.R')

# Load default parameter configurations
source('./bin/rStride_intervention_baseline.R')

# set directory postfix (optional)
dir_postfix <- '_baseline'

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# add default parameters and values
exp_param_list <- get_covid19_default_param()

# change population file
exp_param_list$population_file <- 'pop_belgium600k_c500_teachers_censushh.csv'
#exp_param_list$population_file <- 'pop_belgium1000k_c500_teachers_censushh.csv'

# change parameters and values to combine in a full-factorial grid
exp_param_list$num_days <- 300
exp_param_list$num_parallel_workers <- 8
#exp_param_list$event_log_level <- c("Incidence")
exp_param_list$event_log_level <- c("Transmissions")
exp_param_list$num_seeds<- 2
 
exp_param_list$distancing_workplace_ratio    <- c_str(exp_param_list$distancing_workplace_ratio,0.2)
exp_param_list$distancing_workplace_date     <- c_str(exp_param_list$distancing_workplace_date,'2020-09-01')
exp_param_list$distancing_workplace_delay    <- c_str(exp_param_list$distancing_workplace_delay,7)

exp_param_list$distancing_school_ratio       <- c_str(exp_param_list$distancing_school_ratio,0) # end school closure
exp_param_list$distancing_school_date        <- c_str(exp_param_list$distancing_school_date,'2020-05-18')
exp_param_list$distancing_school_delay       <- c_str(exp_param_list$distancing_school_delay,0)

exp_param_list$distancing_community_ratio    <- c_str(exp_param_list$distancing_community_ratio,0.1)
exp_param_list$distancing_community_date     <- c_str(exp_param_list$distancing_community_date,'2020-09-01')
exp_param_list$distancing_community_delay    <- c_str(exp_param_list$distancing_community_delay,7)
 
exp_param_list$imported_cases_number         <- 10
exp_param_list$imported_cases_date           <- '2020-08-25'
exp_param_list$imported_cases_delay          <- 5

# check period
range(as.Date(exp_param_list$start_date), as.Date(exp_param_list$start_date)+ exp_param_list$num_days)

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
#inspect_tracing_data(project_dir)





 
