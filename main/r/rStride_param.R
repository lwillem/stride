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
#  Copyright 2026, Willem L
############################################################################ #
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_param.R 
#
############################################################################ #

# Clear work environment
rm(list=ls())

# Load rStride
source('./bin/rstride/rStride.R')

# Load default parameter configurations
source('./bin/rStride_covid19_default_param.R')

# set directory postfix (optional)
dir_postfix <- '_param'

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# set number of experiments
num_experiments  <- 16

# get default parameters and values to combine in a LHS design
exp_param_list <- get_default_param()

# # to explore model parameters
#View(get_default_param_info())

# option to change parameter values
exp_param_list$population_file <- 'data/pop_belgium600k_c500_teachers_censushh.csv'
exp_param_list$num_days <- 34      #34, 74, 134
exp_param_list$r0 <- c(0.5,4)
exp_param_list$num_infected_seeds <- c(0,560)
exp_param_list$num_rng_seeds <- 3

# # add maximum for log parsing (memory boundaries)
exp_param_list$logparsing_cases_upperlimit <- 1.5e6

# local cluster
exp_param_list$num_parallel_workers <- 8

# check period
range(as.Date(exp_param_list$start_date), as.Date(exp_param_list$start_date)+ exp_param_list$num_days)

################################################ #
## GENERATE DESIGN OF EXPERIMENT GRID         ####
################################################ #

# get LHS design of experiments
exp_design <- .rstride$get_lhs_exp_design(exp_param_list,
                                          num_experiments = num_experiments,
                                          num_rng_seeds    = exp_param_list$num_rng_seeds)

################################## #
## RUN rSTRIDE                  ####
################################## #
project_dir <- run_rStride(exp_design               = exp_design,
                           dir_postfix              = dir_postfix,
                           num_parallel_workers     = exp_param_list$num_parallel_workers)


############################# #
## INPUT-OUTPUT BEHAVIOR   ####
############################# #
inspect_summary(project_dir)


############################# #
## SURVEY PARTICIPANT DATA ####
############################# #
# inspect_participant_data(project_dir)



################################# #
## PARAMETER ESTIMATION        ####
################################# #
estimate_parameters(project_dir)


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




 
