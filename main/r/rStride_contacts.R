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
#  Copyright 2024, Willem et al.
#############################################################################
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_contacts.R 
#
#############################################################################

# Clear work environment
rm(list=ls())

# load rStride
source('./bin/rstride/rStride.R')

# Load default parameter configurations
source('./bin/rStride_covid19_default_param.R')

# set directory postfix (optional)
dir_postfix <- '_cnt'

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# get default parameters and values to combine in a full-factorial grid
exp_param_list <- get_covid19_default_param()

# contact parameters
exp_param_list$event_log_level           <- "Participants"
exp_param_list$num_days                  <- 1
exp_param_list$num_infected_seeds        <- 1
exp_param_list$num_participants_survey   <- 4999
exp_param_list$num_seeds                 <- 2

# # specific start dates for COVID-19 in BEL
# exp_param_list$start_date <- c('2020-02-17',#'2020-02-22', # weekday, weekend
#                                '2020-02-24',              # holiday
#                                "2020-03-24",#"2020-03-28", # lockdown: weekday, weekend
#                                "2020-05-05",              # first deconfinement 
#                                "2020-06-10",#"2020-06-11", # hh bubble start (weekend, weekday)
#                                "2020-06-02",#"2020-06-07", # exit, with school
#                                "2020-07-07"##"2020-06-12") # exit, without school

exp_param_list$contact_survey_dates          <- c_str(exp_param_list$start_date)
exp_param_list$contact_survey_ages           <- c_str(seq(0,90,10))
exp_param_list$contact_survey_resample       <- 1

# change parameters for development (BEL)
exp_param_list$population_file <- 'data/pop_belgium600k_c500_teachers_censushh.csv'

# change parameters for USA
exp_param_list$start_date              <- c('2020-02-10', # Monday
                                            '2020-02-16') # Sunday
exp_param_list$population_file         <- 'data/pop_usa_wisconsin_dane474k_c1000.csv'
exp_param_list$age_contact_matrix_file <- 'data/contact_matrix_usa_conditional.xml'
exp_param_list$holidays_file           <- 'data/holidays_none.csv'
exp_param_list$contact_survey_dates    <- c_str(exp_param_list$start_date)

# check period
range(as.Date(exp_param_list$start_date), as.Date(exp_param_list$start_date)+ exp_param_list$num_days)

################################################ #
## GENERATE DESIGN OF EXPERIMENT GRID         ####
################################################ #

# get grid-based design of experiments
exp_design <- .rstride$get_full_grid_exp_design(exp_param_list = exp_param_list,
                                                num_seeds      = exp_param_list$num_seeds)
dim(exp_design)

##################################
## RUN rSTRIDE                  ##
##################################
project_dir <- run_rStride(exp_design  = exp_design,
                           dir_postfix = dir_postfix,
                           ignore_stdout            = TRUE,
                           remove_run_output        = FALSE)


#####################################################
## EXPLORE SOCIAL CONTACT PATTERNS                 ##
#####################################################
inspect_contact_data(project_dir)


