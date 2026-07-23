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
#  Copyright 2024, Willem L
#############################################################################
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_r0.R 
#
#############################################################################

# Clear work environment
rm(list=ls())

# Load rStride
source('./bin/rstride/rStride.R')

##################################
## DESIGN OF EXPERIMENTS        ##
##################################

# uncomment the following line to inspect the config xml tags
#names(xmlToList('./config/run_default.xml'))

# set directory postfix (optional)
dir_postfix <- '_r0'

# set the number of realisations per configuration set
num_seeds  <- 5

# add parameters and values to combine in a full-factorial grid
# exp_design <- expand.grid(r0                            = sort(c(seq(1,5,length=5))),
#                           num_days                      = c(20),
#                           rng_seed                      = seq(num_seeds),
#                           start_date                    = c('2020-02-01','2020-02-02','2020-02-03','2020-02-04','2020-02-05','2020-02-06','2020-02-07'),
#                           num_infected_seeds            = 20,
#                           seeding_age_min               = 1,
#                           seeding_age_max               = 99,
#                           disease_config_file           = "disease_covid19_lognorm.xml",
#                           population_file               = "pop_belgium600k_c500_teachers_censushh.csv",
#                           age_contact_matrix_file       = "contact_matrix_flanders_conditional_teachers.xml",
#                           holidays_file                 = "holidays_none.csv",
#                           stringsAsFactors = F)

exp_design <- expand.grid(r0                            = seq(0,24,3),
                          num_days                      = c(40),
                          rng_seed                      = seq(num_seeds),
                          start_date                    = c('2023-01-01'),
                          num_infected_seeds            = 20,
                          seeding_age_min               = 1,
                          seeding_age_max               = 99,
                          disease_config_file           = "data/disease_measles.xml",
                          population_file               = "sim_output/20260710_111118_WI-Dane_conditional_social_contacts/20260710_111118_population_WI-Dane.csv",
                          age_contact_matrix_file       = "sim_output/20260710_111118_WI-Dane_conditional_social_contacts/contact_matrix_usa_conditional.xml",
                          holidays_file                 = "data/calendar_dane_2023_2026_measles.csv",
                          stringsAsFactors = F)

# add a unique seed for each run
set.seed(num_seeds)
exp_design$rng_seed <- sample(nrow(exp_design))
dim(exp_design)

##################################
## RUN rSTRIDE                  ##
##################################
project_dir <- run_rStride(exp_design = exp_design,
                           dir_postfix = dir_postfix,
                           ignore_stdout = T,
                           remove_run_output = T,
                           get_transmission_rdata = T)


############################# #
## INPUT-OUTPUT BEHAVIOR   ####
############################# #
# inspect_summary(project_dir)
inspect_participant_data(project_dir)
# inspect_incidence_data(project_dir)
# inspect_prevalence_data(project_dir)
inspect_transmission_dynamics(project_dir)

##################################
## REPRODUCTION NUMBER          ##
##################################
analyse_transmission_data_for_r0(project_dir)


################################### #
## HOSPITAL ADMISSIONS BY AGE    ####
#####################################

# covid-19 specific!!
analyse_transmission_data_for_hospital_admissions(project_dir)



