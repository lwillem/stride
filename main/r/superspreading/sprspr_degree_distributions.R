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
 
################################################################################################
# Run simulations for 7 days to get degree distribution.                                       #
################################################################################################

sprspr_degree_distributions <- function(mean_transmission_probability, 
                                        num_infected_seeds, 
                                        num_runs)
{
  event_log_level <- "All"
  track_index_case <- "false"
  
  disease_config_file <- "disease_covid19_lognorm.xml"
  holidays_file <- "holidays_none.csv"
  
  num_days <- 7 # To get degree distribution over week + weekend days 
  start_date <- "2020-02-17"
  
  scenario_names = c("degree_dist_baseline_infectiousness_baseline_contacts_part", 
                     "degree_dist_baseline_infectiousness_superspreading_20_contacts",
                     "degree_dist_baseline_infectiousness_superspreading_40_contacts_part",
                     "degree_dist_baseline_infectiousness_superspreading_60_contacts",
                     "degree_dist_baseline_infectiousness_superspreading_100_contacts",
                     "degree_dist_baseline_infectiousness_superspreading_1000_contacts")
  
  contact_distributions <- c("Constant", "Gamma", "Gamma", "Gamma")
  contact_overdispersions <- c(0, 0.2, 0.4, 0.6, 1, 10)

  for (i in seq_along(scenario_names)) {
    run_simulations(scenario_name = scenario_names[i], track_index_case = track_index_case, event_log_level=event_log_level,
                    tp_distribution = "Constant", tp_mean = mean_transmission_probability, tp_overdispersion=0,
                    contact_distribution = contact_distributions[i], contact_distribution_overdispersion = contact_overdispersions[i],
                    disease_config_file = disease_config_file, holidays_file = holidays_file, num_days = num_days,
                    num_infected_seeds = num_infected_seeds, population_file = population_file, start_date = start_date, 
                    num_runs = num_runs )
  }
  
}

sprspr_degree_distributions (mean_transmission_probability = 0.08, 
                             num_infected_seeds = 1, # FIXME or 0 seeds for this experiment? 
                             num_runs = 1) 

# TODO degree distribution when social distancing is practiced 
# TODO sensitivity analysis