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

###################################################################################################
# Run simulations for 40 days, only following the infections of the index case.                   #
###################################################################################################

sprspr_index_case_only <- function(num_infected_seeds, num_runs) 
{
  event_log_level <- "Transmissions"
  track_index_case <- "true"
  
  disease_config_file <- "disease_covid19_lognorm.xml"
  holidays_file <- "holidays_none.csv"
  population_file <- "pop_belgium3000k_c500_teachers_censushh.csv"
  
  num_days <- 40
  start_date <- "2020-02-17"
  
  mean_transmission_probabilities <- c(0.0, 0.025, 0.05, 0.075, 0.1)
  
  # Baseline 
  run_simulations(scenario_name="ico_baseline", track_index_case=track_index_case, event_log_level=event_log_level, 
                  tp_distribution="Constant", tp_mean=mean_transmission_probabilities, tp_overdispersion=0, 
                  contact_distribution="Constant", contact_distribution_overdispersion=0,
                  disease_config_file=disease_config_file, holidays_file=holidays_file, num_days=num_days, 
                  num_infected_seeds=num_infected_seeds, population_file=population_file, start_date=start_date, 
                  num_runs=num_runs)
  
  # Vary dispersion parameter of infectiousness distribution
  scenario_names <- c("ico_infectiousness_overdispersion_1000", 
                      "ico_infectiousness_overdispersion_100", 
                      "ico_infectiousness_overdispersion_60", 
                      "ico_infectiousness_overdispersion_40", 
                      "ico_infectiousness_overdispersion_20")
  
  tp_overdispersions <- c(10, 1, 0.6, 0.4, 0.2) 
  
  for (i in seq_along(scenario_names)) {
    tp_means <- sapply(mean_transmission_probabilities, get_mean_non_truncated_gamma, shape=tp_overdispersions[i])
    
    run_simulations(scenario_name=scenario_names[i], track_index_case = track_index_case, event_log_level=event_log_level,
                    tp_distribution="Gamma", tp_mean=tp_means, tp_overdispersion = tp_overdispersions[i], 
                    contact_distribution = "Constant", contact_distribution_overdispersion = 0,
                    disease_config_file = disease_config_file, holidays_file = holidays_file, num_days=num_days,
                    num_infected_seeds = num_infected_seeds, population_file = population_file, start_date = start_date,
                    num_runs=num_runs)
    
  }
  
  # Vary dispersion parameter of contacts distribution 
  scenario_names <- c("ico_contacts_overdispersion_1000", 
                      "ico_contacts_overdispersion_100", 
                      "ico_contacts_overdispersion_60", 
                      "ico_contacts_overdispersion_40", 
                      "ico_contacts_overdispersion_20")
  
  contact_overdispersions <- c(10, 1, 0.6, 0.4, 0.2)
  
  for (i in seq_along(scenario_names)) {
    run_simulations(scenario_name=scenario_names[i], track_index_case=track_index_case, event_log_level=event_log_level,
                    tp_distribution="Constant", tp_mean=mean_transmission_probabilities, tp_overdispersion=0,
                    contact_distribution="Gamma", contact_distribution_overdispersion=contact_overdispersions[i],
                    disease_config_file=disease_config_file, holidays_file=holidays_file, num_days=num_days,
                    num_infected_seeds=num_infected_seeds, population_file=population_file, start_date=start_date,
                    num_runs=num_runs)
  }
}

sprspr_index_case_only (num_infected_seeds = 1, num_runs = 8) 

# TODO sensitivity analysis (vary number of infected seeds)?