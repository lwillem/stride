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
#  Copyright 2020, Willem L, Libin P
############################################################################ #
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_explore.R 
#
############################################################################ #

# Clear work environment
rm(list=ls())

# Load rStride
source('./bin/rstride/rStride.R')

# Load default parameter configurations
source('./bin/rStride_intervention_baseline.R')

# set directory postfix (optional)
dir_postfix <- '_int'

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# add default parameters and values to combine in a full-factorial grid
exp_param_list <- get_exp_param_default(bool_age_specific_param = T,bool_min_restrictive = T)

# change parameters and values to combine in a full-factorial grid
exp_param_list$num_days <- 300 #300
exp_param_list$num_parallel_workers <- 8
#exp_param_list$event_log_level <- c("Incidence")
exp_param_list$num_seeds<- 2

exp_param_list$temporal_distancing_workplace    <- c_str(seq(0.85,0.75,length=8))
exp_param_list$dates_distancing_workplace       <- NA

# exp_param_list$temporal_distancing_community    <- paste(seq(0.85,0.80,length=8),collapse=',')
# exp_param_list$dates_distancing_community       <- NA
exp_param_list$temporal_distancing_community      <- c_str(0.70,0.4,0.8)
exp_param_list$dates_distancing_community         <- c_str('2020-05-01','2020-08-15','2020-10-01','2020-11-01')

exp_param_list$temporal_distancing_collectivity <- c_str(seq(0.70,0.60,length=8))
exp_param_list$dates_distancing_collectivity    <- NA
#exp_param_list$temporal_distancing_school       <- paste(seq(0.70,0.60,length=5),collapse=',')

exp_param_list$temporal_imported_cases          <- c_str(0,1,1,0)
exp_param_list$dates_imported_cases             <- c_str('2020-06-01','2020-07-31','2020-08-01','2020-08-31','2020-09-01')
exp_param_list$num_infected_seeds               <- 50

exp_param_list$cnt_reduction_school_exit           <- 0.8
exp_param_list$cnt_reduction_school_exit_secondary <- 0.1
exp_param_list$cnt_reduction_school_exit_tertiary  <- 0.6

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
                           num_parallel_workers     = exp_param_list$num_parallel_workers)


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





 
