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

# TODO cnt_reduction_school?
# TODO contact tracing
# TODO contact rate distribution 
run_simulations <- function(scenario_name, 
                            tp_distribution, tp_overdispersion, 
                            cnt_reduction_workplace,
                            cnt_reduction_other,
                            holidays_file,
                            num_infected_seeds,
                            num_runs) {
  exp_design <- expand.grid(
      transmission_probability                      = 0.08,
      num_days                                      = 120,
      rng_seed                                      = seq(num_runs),
      num_infected_seeds                            = num_infected_seeds, 
      disease_config_file                           = "disease_covid19_lognorm.xml",
      population_file                               = "pop_belgium3000k_c500_teachers_censushh.csv",
      age_contact_matrix_file                       = "contact_matrix_flanders_conditional_teachers.xml",
      start_date                                    = "2020-11-01",
      holidays_file                                 = holidays_file,
      event_log_level                               = "Transmissions",
      num_participants_survey                       = 10,
      school_system_adjusted                        = 0,
      cnt_reduction_workplace                       = cnt_reduction_workplace,
      cnt_reduction_other                           = cnt_reduction_other,
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
      test_false_negative                                   = 0,
      cnt_other_exit_delay                                  = 0,
      hosp_probability_factor                               = 1,
      
      transmission_probability_distribution                 = tp_distribution,
      transmission_probability_distribution_overdispersion  = tp_overdispersion,
      
      stringsAsFactors                                      = F
      
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
                           use_date_prefix     = FALSE)
    
}
 
# Load rStride
source('./bin/rstride/rStride.R')

# Simulations without interventions 
run_simulations(scenario_name = "baseline", 
                 tp_distribution = "None", 
                 tp_overdispersion =  0, 
                 cnt_reduction_workplace = 0,
                 cnt_reduction_other = 0,
                 holidays_file = "holidays_none.json",
                 num_infected_seeds = 1,
                 num_runs = 200)
run_simulations(scenario_name = "superspreading_30", 
                tp_distribution = "Gamma", 
                tp_overdispersion =  0.3, 
                cnt_reduction_workplace = 0,
                cnt_reduction_other = 0,
                holidays_file = "holidays_none.json",
                num_infected_seeds = 1,
                num_runs = 200)
run_simulations(scenario_name = "superspreading_40", 
                 tp_distribution = "Gamma", 
                 tp_overdispersion =  0.4, 
                 cnt_reduction_workplace = 0,
                 cnt_reduction_other = 0,
                 holidays_file = "holidays_none.json",
                num_infected_seeds = 1,
                 num_runs = 200)
run_simulations(scenario_name = "superspreading_50", 
                tp_distribution = "Gamma", 
                tp_overdispersion =  0.5, 
                cnt_reduction_workplace = 0,
                cnt_reduction_other = 0,
                holidays_file = "holidays_none.json",
                num_infected_seeds = 1,
                num_runs = 200)
run_simulations(scenario_name = "superspreading_1000", 
                 tp_distribution = "Gamma", 
                 tp_overdispersion =  10, 
                 cnt_reduction_workplace = 0,
                 cnt_reduction_other = 0,
                 holidays_file = "holidays_none.json",
                num_infected_seeds = 1,
                 num_runs = 200)

# TODO test with different number of infected seeds
# TODO test with some pre-existing immunity in the population
# TODO effect on contact tracing 
# TODO heterogeneity in contact behavior