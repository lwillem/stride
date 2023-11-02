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
source('./bin/rStride_intervention_baseline.R')


args = commandArgs(trailingOnly=TRUE)
exp_id = "test1"

exp_param_list <- get_exp_param_default(bool_child_param = TRUE, bool_min_restrictive = TRUE, bool_revised_model_param=TRUE)

exp_param_list$start_date <- c('2021-01-01')
exp_param_list$num_days <- 181 #365  ##60
exp_param_list$num_threads <- 16
exp_param_list$rng_seed <- 0

exp_param_list$output_prefix <- "config/vsc_1/"
dir <- exp_param_list$output_prefix
smd_print("MDP start date", exp_param_list$start_date)

# immunity
exp_param_list$immunity_link_probability = 0
exp_param_list$immunity_profile = "AgeDependent"
exp_param_list$immunity_distribution_file = "../FullPop/immunity_covid_belgium.xml"


# contact tracing
exp_param_list$detection_probability         = .7

exp_param_list$tracing_efficiency_household  = .9
exp_param_list$tracing_efficiency_other      = .5

exp_param_list$case_finding_capacity         = 10000

exp_param_list$delay_contact_tracing         = 2
exp_param_list$test_false_negative           = .1

# TODO:
# exp_param_list$detection_probability         = .0
#
# exp_param_list$tracing_efficiency_household  = .0
# exp_param_list$tracing_efficiency_other      = .0
#
# exp_param_list$case_finding_capacity         = 0
#
# exp_param_list$delay_contact_tracing         = 0
# exp_param_list$test_false_negative           = .0



run_tag <- exp_id
# generate grid
exp_design <- expand.grid(exp_param_list, stringsAsFactors = F)
exp_design$id <- exp_id

# check period
range(as.Date(exp_param_list$start_date), as.Date(exp_param_list$start_date)+ exp_param_list$num_days)


exp_dir <- paste0(dir,"/")
xml_fn <- smd_file_path(exp_dir,"config.xml")
exp_row <- match(TRUE, exp_design$id == exp_id)

config_default_filename <- './config/run_default.xml'
config_default <- create_default_config(config_default_filename, run_tag)

config_exp <- create_config_exp(config_default, exp_dir, exp_design, exp_row)

config_exp$cnt_reduction_school <- 0
config_exp$cnt_reduction_school_secondary <- 0.5
config_exp$cnt_reduction_school_tertiary <- 1.0

config_exp$cnt_reduction_workplace <- 0.7  # 0.5
config_exp$cnt_reduction_other <- 0.7

# config_exp$cnt_reduction_workplace     <- 0.0
# config_exp$cnt_reduction_other         <- 0.0
# config_exp$cnt_baseline_collectivity   <- 0.0

file_name <- smd_file_path(config_exp$output_prefix,'calendar.csv')
end_date <- "2022-12-31"
school_holidays <- T  # TODO
# smd_print("end date:", end_date)
# smd_print("school holidays:", school_holidays)
create_new_cnt_calendar_file(file_name, config_exp, end_date, school_holidays)  # TODO: added

file_name <- paste0("../", file_name)
config_exp$holidays_file <- file_name
config_exp$output_prefix <- "runs/vcs_1/"

save_config_xml(config_exp,xml_fn)
