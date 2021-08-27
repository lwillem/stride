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
#  Copyright 2020, Willem et al.
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
source('./bin/rStride_intervention_baseline.R')

# set directory postfix (optional)
dir_postfix <- '_cnt_long'

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# add default parameters and values to combine in a full-factorial grid
exp_param_list <- get_exp_param_default(bool_revised_model_param = T,
                                        bool_min_restrictive = F)

# contact parameters
exp_param_list$event_log_level           <- "All"
#exp_param_list$num_days                  <- 7
exp_param_list$num_infected_seeds        <- 1
exp_param_list$num_participants_survey   <- 4999
exp_param_list$num_seeds                 <- 1

exp_param_list$cnt_reduction_workplace_exit <- exp_param_list$cnt_reduction_workplace_exit[-2]
exp_param_list$cnt_reduction_other_exit <- exp_param_list$cnt_reduction_other_exit[-2]


# exp_param_list$compliance_delay_other     <- 0
# exp_param_list$compliance_delay_workplace <- 0
# 
# exp_param_list$start_date <- c('2020-02-17',#'2020-02-22', # weekday, weekend
#                                '2020-02-24',              # holiday
#                                "2020-03-24",#"2020-03-28", # lockdown: weekday, weekend
#                                "2020-05-05",              # first deconfinement 
#                                "2020-06-10",#"2020-06-11", # hh bubble start (weekend, weekday)
#                                "2020-06-02",#"2020-06-07", # exit, with school
#                                "2020-07-07"##"2020-06-12"
#                                ) # exit, without school


# change parameters
#exp_param_list$population_file <- 'pop_belgium600k_c500_teachers_censushh.csv'
exp_param_list$population_file                <- "pop_belgium600k_c500_teachers_censushh_extended3_size2.csv"
exp_param_list$cnt_intensity_householdCluster <- c(0,2/7)


# check period
range(as.Date(exp_param_list$start_date), as.Date(exp_param_list$start_date)+ exp_param_list$num_days)

################################################ #
## GENERATE DESIGN OF EXPERIMENT GRID         ####
################################################ #

# get grid-based design of experiments
exp_design <- .rstride$get_full_grid_exp_design(exp_param_list = exp_param_list,
                                                num_seeds      = exp_param_list$num_seeds)
dim(exp_design)



# ##################################
# ## DESIGN OF EXPERIMENTS        ##
# ##################################
# 
# # uncomment the following line to inspect the config xml tags
# #names(xmlToList('./config/run_default.xml'))
# 
# # set directory postfix (optional)
# dir_postfix <- '_cnt'
# 
# # set the number of realisations per configuration set
# num_seeds  <- 1
# 
# # add parameters and values to combine in a full-factorial grid
# exp_design <- expand.grid(event_log_level           = "All",
#                           num_days                  = 1,
#                           num_infected_seeds        = 1,
#                           num_participants_survey   = 4999,
#                           start_date                = c('2020-02-16','2020-02-17','2020-02-25',"2020-04-10","2020-04-11","2020-05-03","2020-05-04"),
#                           #start_date                = c('2020-03-1'),
#                           rng_seed                  = 1:num_seeds,
#                           disease_config_file       = "disease_covid19_age.xml", 
#                           population_file               = "pop_belgium600k_c500_teachers_censushh.csv",
#                           # population_file               = c("pop_belgium600k_c500_teachers_censushh_extended3.csv",
#                           #                                   "pop_belgium600k_c500_teachers_censushh.csv"),
#                           age_contact_matrix_file   = "contact_matrix_flanders_conditional_teachers.xml",
#                           #holidays_file             = "calendar_belgium_2020_covid19_may_workplace.json",
#                           holidays_file             = "calendar_belgium_2020_covid19_exit_school_adjusted.json",
#                           # holidays_file             = c("calendar_belgium_2020_covid19_may_school.json",
#                           #                               "calendar_belgium_2020_covid19_may_workplace.json"),
#                           school_system_adjusted        = 1,
#                           cnt_reduction_workplace       = c(0.8),
#                           cnt_reduction_other           = c(0.85),
#                           compliance_delay_workplace    = c(0),
#                           compliance_delay_other        = c(0),
#                           num_daily_imported_cases      = c(0),
#                           cnt_intensity_householdCluster = c(0,0.99),
#                           stringsAsFactors = F)
# 
# # add a unique seed for each run
# set.seed(125)
# exp_design$rng_seed <- sample(nrow(exp_design))

##################################
## RUN rSTRIDE                  ##
##################################
project_dir <- run_rStride(exp_design  = exp_design,
                           dir_postfix = dir_postfix,
                           ignore_stdout            = TRUE,
                           remove_run_output        = TRUE)


#####################################################
## EXPLORE SOCIAL CONTACT PATTERNS                 ##
#####################################################
inspect_contact_data(project_dir)


