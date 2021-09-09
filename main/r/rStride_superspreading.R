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

# Load rStride
source('./bin/rstride/rStride.R')

# Function to run simulations for 1 scenario
run_simulations <- function(scenario_name,                        # Label for output
                            track_index_case,                     # Track only index case?
                            
                            tp_distribution,                      # Distribution applied to individual transmission probability
                            tp_mean,                              # Mean transmission probability (of non-truncated distribution)
                            tp_overdispersion,                    # Overdispersion of distribution for individual transmission probability
                            
                            contact_distribution,                 # Distribution for individual contact rates
                            contact_distribution_overdispersion,  # Overdispersion of distribution for individual contact rates
                            
                            disease_config_file,                  # File with disease parameters
                            holidays_file,                        # File with holidays + periods of social distancing
                            num_days,                             # Number of days to run simulation
                            num_infected_seeds,                   # Number of infected to seed at beginning of simulation
                            num_runs                              # Number of simulations to run
) {
  
  exp_design <- expand.grid(
    
    # General parameters 
    age_contact_matrix_file                        = "contact_matrix_flanders_conditional_teachers.xml",
    cnt_intensity_householdCluster                 = 0,
    disease_config_file                            = disease_config_file,
    event_log_level                                = "Transmissions",
    hosp_probability_factor                        = 1,
    population_file                                = "pop_belgium3000k_c500_teachers_censushh.csv",
    num_daily_imported_cases                       = 0,
    num_days                                       = num_days,
    num_infected_seeds                             = num_infected_seeds,
    num_participants_survey                        = 0,
    output_cases                                   = "false",
    rng_seed                                       = seq(num_runs),
    start_date                                     = "2020-11-01",
    track_index_case                               = track_index_case,
    
    # Parameters relating to social distancing measures
    holidays_file                                 = holidays_file,
    compliance_delay_workplace                    = 7,
    compliance_delay_other                        = 7,
    cnt_other_exit_delay                          = 0,
    
    # Parameters relating to individual transmission probability distribution 
    transmission_probability_distribution         = tp_distribution,
    transmission_probability                      = tp_mean,
    transmission_probability_distribution_overdispersion = tp_overdispersion,
    
    # Parameters relating to individual contact rate distribution 
    contact_distribution = contact_distribution,
    contact_distribution_overdispersion = contact_distribution_overdispersion,
    
    # Parameters relating to contact tracing
    detection_probability                         = 0, 
    tracing_efficiency_household                  = 0,
    tracing_efficiency_other                      = 0,
    case_finding_capacity                         = 0,
    delay_isolation_index                         = 0,
    delay_contact_tracing                         = 0,
    test_false_negative                           = 0,
    
    # threshold for log parsing (default is NA == no threshold)
    logparsing_cases_upperlimit                   = NA,
    
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
                             use_date_prefix     = FALSE)
  
}


# 
# ################################## #
# ## DESIGN OF EXPERIMENTS        ####
# ################################## #
# 
# # get default parameter values to combine in a full-factorial grid
# get_exp_param_default <- function(){
#   
#   out <- list(
#               holidays_file                 = 'calendar_belgium_2020_covid19_exit_school_adjusted.csv',
#               compliance_delay_workplace    = 6,
#               compliance_delay_other        = 6,
#               compliance_delay_collectivity = 1, # dummy, since not used in original setting
#   )


###############################################################################################
# Run simulations without interventions.                                                      #
###############################################################################################

num_days <- 200
num_infected_seeds <- 1
num_runs <- 8
 
run_simulations(
  scenario_name = "baseline_infectiousness_baseline_contacts", track_index_case = "false",
  tp_distribution = "Constant", tp_mean = 0.08, tp_overdispersion = 0,
  contact_distribution = "Constant", contact_distribution_overdispersion = 0,
  disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.csv",
  num_days = num_days, num_infected_seeds = num_infected_seeds, num_runs = num_runs)

run_simulations(
  scenario_name = "superspreading_40_infectiousness_baseline_contacts", track_index_case = "false",
  tp_distribution = "Gamma", tp_mean = 0.081264, tp_overdispersion = 0.4,
  contact_distribution = "Constant", contact_distribution_overdispersion = 0,
  disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.csv",
  num_days = num_days, num_infected_seeds = num_infected_seeds, num_runs = num_runs)

run_simulations(
  scenario_name = "baseline_infectiousness_superspreading_40_contacts", track_index_case = "false",
  tp_distribution = "Constant", tp_mean = 0.08, tp_overdispersion = 0,
  contact_distribution = "Gamma", contact_distribution_overdispersion = 0.4,
  disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.csv",
  num_days = num_days, num_infected_seeds = num_infected_seeds, num_runs = num_runs)

run_simulations(
  scenario_name = "superspreading_40_infectiousness_superspreading_40_contacts", track_index_case = "false",
  tp_distribution = "Gamma", tp_mean = 0.081264, tp_overdispersion = 0.4,
  contact_distribution = "Gamma", contact_distribution_overdispersion = 0.4,
  disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.csv",
  num_days = num_days, num_infected_seeds = num_infected_seeds, num_runs = num_runs)


# run_simulations(
#   scenario_name = "superspreading_1000", track_index_case = "false",
#   tp_distribution = "Gamma", tp_mean = 0.08, tp_overdispersion = 10,
#   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
#   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,  
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.json",
#   num_days = 180, num_infected_seeds = 1, num_runs = 200)
# 
# run_simulations(
#   scenario_name = "superspreading_60", track_index_case = "false",
#   tp_distribution = "Gamma", tp_mean = 0.080166, tp_overdispersion =  0.6, 
#   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
#   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.json",
#   num_days = 180, num_infected_seeds = 1, num_runs = 200)
# 
# run_simulations(
#   scenario_name = "superspreading_20", track_index_case = "false",
#   tp_distribution = "Gamma", tp_mean = 0.094622, tp_overdispersion =  0.2, 
#   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
#   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.json",
#   num_days = 180, num_infected_seeds = 1, num_runs = 200)
# 
# ###############################################################################################
# # Run simulations for 1100 days with a period of social distancing.                           #
# ###############################################################################################
# 
# # Holidays file:
# #     - Run for 30 days without interventions
# #     - 'Lockdown' for 70 days (schools closed, workplace cnt reduction = 86%, community cnt reduction = 85%)
# #     - Partial relaxation after day 100 (school cnt reduction = 50%, workplace cnt reduction = 75%, community cnt reduction = 85%)
# 
# run_simulations(
#   scenario_name = "social_distancing_baseline", track_index_case = "false",
#   tp_distribution = "Constant", tp_mean = 0.08, tp_overdispersion =  0, 
#   cnt_reduction_workplace = 0.86, cnt_reduction_other = 0.85, 
#   cnt_reduction_workplace_exit = 0.75, cnt_reduction_other_exit = 0.85, cnt_reduction_school_exit = 0.5,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "calendar_social_distancing_no_holidays.csv",
#   num_days = 1100, num_infected_seeds = 1, num_runs = 200)
# 
# run_simulations(
#   scenario_name = "social_distancing_superspreading_1000", track_index_case = "false",
#   tp_distribution = "Gamma", tp_mean = 0.08, tp_overdispersion = 10, 
#   cnt_reduction_workplace = 0.86, cnt_reduction_other = 0.85, 
#   cnt_reduction_workplace_exit = 0.75, cnt_reduction_other_exit = 0.85, cnt_reduction_school_exit = 0.5,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "calendar_social_distancing_no_holidays.csv",
#   num_days = 1100, num_infected_seeds = 1, num_runs = 200)
# 
# run_simulations(
#   scenario_name = "social_distancing_superspreading_60", track_index_case = "false",
#   tp_distribution = "Gamma", tp_mean = 0.080166, tp_overdispersion = 0.6, 
#   cnt_reduction_workplace = 0.86, cnt_reduction_other = 0.85, 
#   cnt_reduction_workplace_exit = 0.75, cnt_reduction_other_exit = 0.85, cnt_reduction_school_exit = 0.5,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "calendar_social_distancing_no_holidays.csv",
#   num_days = 1100, num_infected_seeds = 1, num_runs = 200)
# 
# run_simulations(
#   scenario_name = "social_distancing_superspreading_40", track_index_case = "false",
#   tp_distribution = "Gamma", tp_mean = 0.081264, tp_overdispersion = 0.4, 
#   cnt_reduction_workplace = 0.86, cnt_reduction_other = 0.85, 
#   cnt_reduction_workplace_exit = 0.75, cnt_reduction_other_exit = 0.85, cnt_reduction_school_exit = 0.5,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "calendar_social_distancing_no_holidays.csv",
#   num_days = 1100, num_infected_seeds = 1, num_runs = 200)
# 
# run_simulations(
#   scenario_name = "social_distancing_superspreading_20", track_index_case = "false",
#   tp_distribution = "Gamma", tp_mean = 0.094622, tp_overdispersion = 0.2, 
#   cnt_reduction_workplace = 0.86, cnt_reduction_other = 0.85, 
#   cnt_reduction_workplace_exit = 0.75, cnt_reduction_other_exit = 0.85, cnt_reduction_school_exit = 0.5,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "calendar_social_distancing_no_holidays.csv",
#   num_days = 1100, num_infected_seeds = 1, num_runs = 200)
# 
# ###############################################################################################
# # Run simulations, tracking only the index case, for 40 days.                                 #
# ###############################################################################################
# 
# run_simulations(
#   scenario_name = "index_case_only_baseline", track_index_case = "true",
#   tp_distribution = "Constant", tp_mean = c(0.0, 0.025, 0.05, 0.075, 0.1), tp_overdispersion = 0, 
#   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
#   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.json",
#   num_days = 40, num_infected_seeds = 1, num_runs = 1000)
# 
# run_simulations(
#   scenario_name = "index_case_only_superspreading_1000", track_index_case = "true",
#   tp_distribution = "Gamma", tp_mean = c(0.0, 0.025, 0.05, 0.075, 0.1), tp_overdispersion = 10, 
#   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
#   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.json",
#   num_days = 40, num_infected_seeds = 1, num_runs = 1000)
# 
# run_simulations(
#   scenario_name = "index_case_only_superspreading_60", track_index_case = "true",
#   tp_distribution = "Gamma", tp_mean = c(0.0, 0.025002, 0.050004, 0.075102, 0.10086), tp_overdispersion = 0.6, 
#   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
#   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.json",
#   num_days = 40, num_infected_seeds = 1, num_runs = 1000)
# 
# run_simulations(
#   scenario_name = "index_case_only_superspreading_40", track_index_case = "true",
#   tp_distribution = "Gamma", tp_mean = c(0.0, 0.025, 0.050044, 0.075852, 0.10438), tp_overdispersion = 0.4, 
#   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
#   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.json",
#   num_days = 40, num_infected_seeds = 1, num_runs = 1000)
# 
# run_simulations(
#   scenario_name = "index_case_only_superspreading_20", track_index_case = "true",
#   tp_distribution = "Gamma", tp_mean = c(0.0, 0.025014, 0.051518, 0.085884, 0.141108), tp_overdispersion = 0.2, 
#   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
#   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
#   disease_config_file = "disease_covid19_lognorm.xml", holidays_file = "holidays_none.json",
#   num_days = 40, num_infected_seeds = 1, num_runs = 1000)
# 
# ###############################################################################################
# # Run simplified simulations, tracking only the index case, for 40 days.                      #
# # Results for these are used to compare to estimates obtained with theoretical description.   #
# ###############################################################################################
# 
# # Some 'hacks' needed in C++ code to run these:
# #     - In Person::Update(): everyone is always present in 
# #       household / school / workplace / primary and secondary community, 
# #       and always absent in HouseholdCluster
# #     - In DiseaseSeeder::ImportInfectedCases(): fill in id of person being looked at on line 98,
# #       so that the same individual is chosen each time as the index case.
# #       I.e..: Person& p = pop->at(static_cast<size_t>(person_id));
# #       Persons I ran simulations for:
# #           * Adult (person_id = 4)
# #           * Child (person_id = 3)
# #           * Elderly (person_id = 5)
# #
# # To run simulations for child / elderly, only the scenario_name needs to be updated in the R script below. 
# 
# # run_simulations(
# #   scenario_name = "index_case_only_simple_adult_baseline", track_index_case = "true",
# #   tp_distribution = "Constant", tp_mean = c(0.0, 0.025, 0.05, 0.075, 0.1), tp_overdispersion = 0, 
# #   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
# #   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
# #   disease_config_file = "disease_covid19_lognorm_nocntreduction.xml", holidays_file = "holidays_none.json",
# #   num_days = 40, num_infected_seeds = 1, num_runs = 200)
# # 
# # run_simulations(
# #   scenario_name = "index_case_only_simple_adult_superspreading_1000", track_index_case = "true",
# #   tp_distribution = "Gamma", tp_mean = c(0.0, 0.025, 0.05, 0.075, 0.1), tp_overdispersion = 10, 
# #   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
# #   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
# #   disease_config_file = "disease_covid19_lognorm_nocntreduction.xml", holidays_file = "holidays_none.json",
# #   num_days = 40, num_infected_seeds = 1, num_runs = 200)
# # 
# # run_simulations(
# #   scenario_name = "index_case_only_simple_adult_superspreading_60", track_index_case = "true",
# #   tp_distribution = "Gamma", tp_mean = c(0.0, 0.025002, 0.050004, 0.075102, 0.10086), tp_overdispersion = 0.6, 
# #   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
# #   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
# #   disease_config_file = "disease_covid19_lognorm_nocntreduction.xml", holidays_file = "holidays_none.json",
# #   num_days = 40, num_infected_seeds = 1, num_runs = 200)
# # 
# # run_simulations(
# #   scenario_name = "index_case_only_simple_adult_superspreading_40", track_index_case = "true",
# #   tp_distribution = "Gamma", tp_mean = c(0.0, 0.025, 0.050044, 0.075852, 0.10438), tp_overdispersion = 0.4, 
# #   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
# #   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
# #   disease_config_file = "disease_covid19_lognorm_nocntreduction.xml", holidays_file = "holidays_none.json",
# #   num_days = 40, num_infected_seeds = 1, num_runs = 200)
# # 
# # run_simulations(
# #   scenario_name = "index_case_only_simple_adult_superspreading_20", track_index_case = "true",
# #   tp_distribution = "Gamma", tp_mean = c(0.0, 0.025014, 0.051518, 0.085884, 0.141108), tp_overdispersion = 0.2, 
# #   cnt_reduction_workplace = 0, cnt_reduction_other = 0, 
# #   cnt_reduction_workplace_exit = 0, cnt_reduction_other_exit = 0, cnt_reduction_school_exit = 0,
# #   disease_config_file = "disease_covid19_lognorm_nocntreduction.xml", holidays_file = "holidays_none.json",
# #   num_days = 40, num_infected_seeds = 1, num_runs = 200)