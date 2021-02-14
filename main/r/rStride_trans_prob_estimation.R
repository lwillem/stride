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
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
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

run_simulations <- function(scenario_name, num_runs,
                            tp_distribution = "Constant", tp_overdispersion = "10") {
  
  exp_design <- expand.grid(
    track_index_case                              = "true",
    transmission_probability                      = c(0.00, 0.025, 0.05, 0.075, 0.10),
    num_days                                      = 30,
    rng_seed                                      = seq(num_runs),
    num_participants_survey                       = 10,
    num_infected_seeds                            = 1,
    disease_config_file                           = "disease_covid19_lognorm_nocntreduction.xml",
    population_file                               = "pop_belgium3000k_c500_teachers_censushh.csv",
    age_contact_matrix_file                       = "contact_matrix_flanders_conditional_teachers.xml",
    start_date                                    = "2020-11-01",
    holidays_file                                 = "holidays_none.json",
    school_system_adjusted                        = 0,
    cnt_reduction_workplace                       = 0,
    cnt_reduction_other                           = 0,
    compliance_delay_workplace                    = 0,
    compliance_delay_other                        = 0,
    num_daily_imported_cases                      = 0,
    cnt_reduction_workplace_exit                  = 0,
    cnt_reduction_other_exit                      = 0,
    cnt_reduction_school_exit                     = 0,
    cnt_intensity_householdCluster                = 0,
    detection_probability                         = 0, 
    tracing_efficiency_household                  = 0,
    tracing_efficiency_other                      = 0,
    case_finding_capacity                         = 0,
    delay_isolation_index                         = 0,
    delay_contact_tracing                         = 0,
    test_false_negative                           = 0,
    cnt_other_exit_delay                          = 0,
    event_log_level                               = "Transmissions",
    hosp_probability_factor                       = 1,
    
    
    transmission_probability_distribution                 = tp_distribution,
    transmission_probability_distribution_overdispersion  = tp_overdispersion,
    
    stringsAsFactors                              = F
    
  )
  
  # add a unique seed for each run
  #set.seed(125)
  #exp_design$rng_seed <- sample(nrow(exp_design))
  exp_design$rng_seed <- sample(0:100000000, nrow(exp_design))
  
  # run rSTRIDE
  project_dir <- run_rStride(exp_design          = exp_design,
                             dir_postfix         = scenario_name,
                             remove_run_output   = FALSE,
                             parse_log_data      = FALSE,
                             use_date_prefix     = FALSE, 
                             stdout_fn                = 'out.txt' )
  
}

# Load rStride
source('./bin/rstride/rStride.R')


run_simulations("index_case_only_baseline", 200)

run_simulations("index_case_only_superspreading_30", 200, tp_distribution = "Gamma", tp_overdispersion = 0.3)
run_simulations("index_case_only_superspreading_40", 200, tp_distribution = "Gamma", tp_overdispersion = 0.4)
run_simulations("index_case_only_superspreading_50", 200, tp_distribution = "Gamma", tp_overdispersion = 0.5)

run_simulations("index_case_only_superspreading_1000", 200, tp_distribution = "Gamma", tp_overdispersion = 10)


