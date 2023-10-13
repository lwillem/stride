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
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
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

# set the number of realisations per configuration set
num_seeds  <- 4

# add parameters and values to combine in a full-factorial grid
exp_design <- expand.grid(r0                            = seq(2.7,2.8,0.1),
                          num_days                      = 20,
                          rng_seed                      = seq(num_seeds),
                          num_participants_survey       = 3000,
                          num_infected_seeds            = 180, 
                          disease_config_file           = "disease_covid19_lognorm.xml",
                          population_file               = c("pop_belgium600k_c500_teachers_censushh.csv"),
                          age_contact_matrix_file       = "contact_matrix_flanders_conditional_teachers.xml",
                          start_date                    = c('2020-02-22'),
                          holidays_file                 = "holidays_belgium_2019_2021.csv",
                          num_daily_imported_cases      = 0,
                          stringsAsFactors = F)

# add a unique seed for each run
set.seed(125)
exp_design$rng_seed <- sample(nrow(exp_design))
dim(exp_design)

##################################
## RUN rSTRIDE                  ##
##################################
project_dir <- run_rStride(exp_design               = exp_design,
                           dir_postfix              = dir_postfix)


#####################################
## EXPLORE INPUT-OUTPUT BEHAVIOR   ##
#####################################
inspect_summary(project_dir)


#####################################
## EXPLORE SURVEY PARTICIPANT DATA ##
#####################################
inspect_participant_data(project_dir)


##################################
## EXPLORE INCIDENCE DATA       ##
##################################
inspect_incidence_data(project_dir)


##################################
## EXPLORE TRANSMISSION         ##
##################################
inspect_transmission_dynamics(project_dir)


##################################
## EXPLORE PREVALENCE           ##
##################################
inspect_prevalence_data(project_dir)




