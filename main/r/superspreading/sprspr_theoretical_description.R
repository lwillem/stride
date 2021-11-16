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
#  Copyright 2021, Kuylen E
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
 
# Load functions to run simulations
source('./bin/sprspr_util.R')

sprspr_theoretical_description <- function(output_prefix, 
                                           infected_seed_id, 
                                           num_runs) 
{
  event_log_level <- "Transmissions"
  run_simplified <- "true"
  track_index_case <- "true"
  
  # No contact reduction for symptomatic individuals, 7 day infectious period for all infected 
  disease_config_file <- "disease_covid19_lognorm_nocntreduction.xml" 
  # No holidays / social distancing 
  holidays_file <- "holidays_none.csv"
  population_file <- "pop_belgium3000k_c500_teachers_censushh.csv"
  
  num_days <- 40
  start_date <- "2020-02-17"
  
  mean_transmission_probabilities <- c(0.0, 0.025, 0.05, 0.075, 0.1)
  
  # Baseline
  run_simulations(scenario_name = paste0(output_prefix, "theoretical_baseline"), event_log_level=event_log_level, 
                  run_simplified=run_simplified, track_index_case=track_index_case,
                  tp_distribution = "Constant", tp_mean = mean_transmission_probabilities, tp_overdispersion = 0,
                  contact_distribution = "Constant", contact_distribution_overdispersion = 0,
                  disease_config_file = disease_config_file, holidays_file = holidays_file, 
                  num_days=num_days, population_file = population_file, start_date = start_date,
                  infected_seed_id = infected_seed_id, num_runs = num_runs)
  
  # Vary dispersion paramter of infectiousness distribution 
  scenario_names <- c("theoretical_infectiousness_overdispersion_1000", 
                      "theoretical_infectiousness_overdispersion_100", 
                      "theoretical_infectiousness_overdispersion_60", 
                      "theoretical_infectiousness_overdispersion_40", 
                      "theoretical_infectiousness_overdispersion_20")
  tp_overdispersions <- c(10, 1, 0.6, 0.4, 0.2) 
  
  for (i in seq_along(scenario_names)) {
    tp_means <- sapply(mean_transmission_probabilities, get_mean_non_truncated_gamma, shape=tp_overdispersions[i])
    
    run_simulations(scenario_name = paste0(output_prefix, scenario_names[i]), event_log_level=event_log_level, 
                    run_simplified=run_simplified, track_index_case=track_index_case,
                    tp_distribution = "Gamma", tp_mean = tp_means, tp_overdispersion = tp_overdispersions[i],
                    contact_distribution = "Constant", contact_distribution_overdispersion = 0,
                    disease_config_file = disease_config_file, holidays_file = holidays_file, 
                    num_days=num_days, population_file = population_file, start_date = start_date,
                    infected_seed_id = infected_seed_id, num_runs = num_runs)
  }
  
  # Vary dispersion parameter of contacts distribution 
  scenario_names <- c("theoretical_contacts_overdispersion_1000", 
                      "theoretical_contacts_overdispersion_100", 
                      "theoretical_contacts_overdispersion_60", 
                      "theoretical_contacts_overdispersion_40", 
                      "theoretical_contacts_overdispersion_20")
  contact_overdispersions <- c(10, 1, 0.6, 0.4, 0.2) 
  
  for (i in seq_along(scenario_names)) {
    run_simulations(scenario_name = paste0(output_prefix, scenario_names[i]), event_log_level=event_log_level, 
                    run_simplified=run_simplified, track_index_case=track_index_case,
                    tp_distribution = "Constant", tp_mean = mean_transmission_probabilities, tp_overdispersion = 0,
                    contact_distribution = "Gamma", contact_distribution_overdispersion = contact_overdispersions[i],
                    disease_config_file = disease_config_file, holidays_file = holidays_file, 
                    num_days=num_days, population_file = population_file, start_date = start_date,
                    infected_seed_id = infected_seed_id, num_runs = num_runs)
  }
}

###############################################################################################
# Run simplified simulations, tracking only the index case, for 40 days.                      #
# Results for these are used to compare to estimates obtained with theoretical description.   #
###############################################################################################

# Run simulations for adult
sprspr_theoretical_description(output_prefix = "adult_", 
                               infected_seed_id = 4, 
                               num_runs = 8) 

# Run simulations for child 
sprspr_theoretical_description(output_prefix = "child_", 
                               infected_seed_id = 3, 
                               num_runs = 8) 

# Run simulations for elderly 
sprspr_theoretical_description(output_prefix = "elderly_", 
                               infected_seed_id = 5, 
                               num_runs = 8) 

# TODO combinations of infectiousness + contacts overdispersion? 
