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

# parse command line interface (CLI) arguments
cli_args = commandArgs(trailingOnly=TRUE)

# clear workspace, but leave the command line arguments
rm(list=ls()[ls()!='cli_args'])

# Load rStride
source('./bin/rstride/rStride.R')

# Load functions to run simulations
source('./bin/sprspr_util.R')

sprspr_full_runs <- function(output_prefix, 
                             holidays_file,
                             mean_transmission_probability, 
                             num_days,
                             num_infected_seeds,
                             num_runs) 
{
  if (missing(output_prefix)) {
    output_prefix <- ""
  }
  
  # Parameters 
  event_log_level <- "Transmissions"
  track_index_case <- "false"
  
  start_date <- "2020-02-17"
  
  disease_config_file <- "disease_covid19_lognorm.xml"
  population_file <- "pop_belgium11M_c500_teachers_censushh.csv"
  
  # Baseline 
  run_simulations(scenario_name=paste0(output_prefix, "baseline"), event_log_level=event_log_level, track_index_case=track_index_case, 
                  tp_distribution="Constant", tp_mean=mean_transmission_probability, tp_overdispersion=0, 
                  contact_distribution="Constant", contact_distribution_overdispersion=0,
                  disease_config_file=disease_config_file, holidays_file=holidays_file, num_days=num_days, 
                  num_infected_seeds=num_infected_seeds, population_file=population_file, start_date=start_date, 
                  num_runs=num_runs)
  
  # Vary dispersion parameter of infectiousness distribution 
  scenario_names <- c("infectiousness_overdispersion_1000", 
                      "infectiousness_overdispersion_100", 
                      "infectiousness_overdispersion_60", 
                      "infectiousness_overdispersion_40", 
                      "infectiousness_overdispersion_20")
  tp_overdispersions <- c(10, 1, 0.6, 0.4, 0.2) 
  
  for (i in seq_along(scenario_names)) {
    tp_mean <- get_mean_non_truncated_gamma(mean_transmission_probability, tp_overdispersions[i])
    
    run_simulations(scenario_name=paste0(output_prefix, scenario_names[i]), track_index_case=track_index_case, event_log_level=event_log_level, 
                    tp_distribution="Gamma", tp_mean=tp_mean, tp_overdispersion=tp_overdispersions[i], 
                    contact_distribution="Constant", contact_distribution_overdispersion=0,
                    disease_config_file=disease_config_file, holidays_file=holidays_file, num_days=num_days, 
                    num_infected_seeds=num_infected_seeds, population_file=population_file, start_date=start_date, 
                    num_runs=num_runs)
  }
  
  # Vary dispersion parameter of contacts distribution 
  scenario_names <- c("contacts_overdispersion_1000", 
                      "contacts_overdispersion_100", 
                      "contacts_overdispersion_60", 
                      "contacts_overdispersion_40", 
                      "contacts_overdispersion_20")
  
  contact_overdispersions <- c(10, 1, 0.6, 0.4, 0.2) 
  
  for (i in seq_along(scenario_names)) {
    run_simulations(scenario_name=paste0(output_prefix, scenario_names[i]), track_index_case=track_index_case, event_log_level=event_log_level, 
                    tp_distribution="Constant", tp_mean=mean_transmission_probability, tp_overdispersion=0, 
                    contact_distribution="Gamma", contact_distribution_overdispersion=contact_overdispersions[i],
                    disease_config_file=disease_config_file, holidays_file=holidays_file, num_days=num_days, 
                    num_infected_seeds=num_infected_seeds, population_file=population_file, start_date=start_date, 
                    num_runs=num_runs)
  }
}

use_interventions <- FALSE

if(length(cli_args) >= 1) {
  use_interventions <- as.logical(cli_args[[1]])
}

if (use_interventions) {
  # Run simulations with social distancing intervention 
  sprspr_full_runs(output_prefix = "sd_",
                   holidays_file = "calendar_social_distancing_comm_85_65_work_85_65.csv",
                   mean_transmission_probability = 0.08,
                   num_days = 200, 
                   num_infected_seeds = 1, 
                   num_runs = 200) 
} else {
  # Run simulations without interventions
  sprspr_full_runs(holidays_file = "holidays_belgium_2019_2021.csv",
                   mean_transmission_probability = 0.08, 
                   num_days = 200,
                   num_infected_seeds = 1, 
                   num_runs = 200)   
}

# TODO combinations of infectiousness + contacts overdispersion? 
# TODO sensitivity analysis (number of infected seeds + mean transmission probability)
