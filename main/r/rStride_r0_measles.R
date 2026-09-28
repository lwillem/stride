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
#  Copyright 2024, Manansala R, Willem L
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
exp_design <- expand.grid(r0                            = seq(0,16,3),
                          num_days                      = c(40),
                          rng_seed                      = seq(num_seeds),
                          start_date                    = c('2023-01-01'),
                          num_infected_seeds            = 20,
                          seeding_age_min               = 1,
                          seeding_age_max               = 99,
                          disease_config_file           = "data/disease_measles_usa.xml",
                          population_file               = "data/pop_usa_tx_gaines_c1000.csv",
                          age_contact_matrix_file       = "data/contact_matrix_usa_tx_gaines_c1000.xml",
                          holidays_file                 = "data/calendar_USA_2023_2026_measles.csv",
                          stringsAsFactors = F)


# combine the two age-specific contact probabilities of a candidate pair by
# their mean rather than their minimum (the kernel default).
# NOTE: the resulting fit is valid for this rule only. The b0/b1/b2 currently
# in disease_measles_usa.xml were produced under "Min" and do not apply here.
exp_design$contact_probability_rule <- "Mean"

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



