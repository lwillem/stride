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

# Function to run simulations for 1 scenario
run_simulations <- function(
                              scenario_name,                        # Label for output
                              
                              event_log_level,                      # Log level
                              run_simplified,                       # Run simplified simulations for comparison to theor. description
                              track_index_case,                     # Track only index case?
                              
                              
                              tp_distribution,                      # Distribution applied to individual transmission probability
                              tp_mean,                              # Mean transmission probability (of non-truncated distribution)
                              tp_overdispersion,                    # Overdispersion of distribution for individual transmission probability
                              
                              contact_distribution,                 # Distribution for individual contact rates
                              contact_distribution_overdispersion,  # Overdispersion of distribution for individual contact rates
                              
                              disease_config_file,                  # File with disease parameters
                              holidays_file,                        # File with holidays + periods of social distancing
                              num_days,                             # Number of days to run simulation
                              population_file,                      # File with full population 
                              start_date,                           # Start date for the simulation
                              
                              num_infected_seeds,                   # Number of infected to seed at beginning of simulation
                              infected_seed_id,                     # Or alternatively, the id of the person to infect at beginning of simulation
                              
                              num_runs                             # Number of simulations to run
) {
  
  if (missing(run_simplified)) {
    run_simplified <- "false"
  }
  
  exp_design <- expand.grid(
        #######################
        # Variable parameters #
        #######################
        disease_config_file                                  = disease_config_file,
        event_log_level                                      = event_log_level,
        holidays_file                                        = holidays_file, 
        num_days                                             = num_days, 
        population_file                                      = population_file, 
        rng_seed                                             = seq(num_runs),
        run_simplified                                       = run_simplified,
        start_date                                           = start_date,
        track_index_case                                     = track_index_case,
        
        # Parameters relating to individual transmission probability distribution 
        transmission_probability_distribution                = tp_distribution,
        transmission_probability                             = tp_mean,
        transmission_probability_distribution_overdispersion = tp_overdispersion,
         
        # Parameters relating to individual contact rate distribution 
        contact_distribution                                 = contact_distribution,
        contact_distribution_overdispersion                  = contact_distribution_overdispersion,
    
        ####################
        # Fixed parameters #
        ####################
        age_contact_matrix_file                              = "contact_matrix_flanders_conditional_teachers.xml",
        cnt_intensity_householdCluster                       = 0,
        hosp_probability_factor                              = 1,
        logparsing_cases_upperlimit                          = NA,
        num_daily_imported_cases                             = 0,
        num_participants_survey                              = 0,
        output_cases                                         = "false",
        
        # Parameters relating to contact tracing (fixed)
        detection_probability                                = 0, 
        tracing_efficiency_household                         = 0,
        tracing_efficiency_other                             = 0,
        case_finding_capacity                                = 0,
        delay_isolation_index                                = 0,
        delay_contact_tracing                                = 0,
        test_false_negative                                  = 0,
    
        stringsAsFactors                                     = F
  )
  
  if (!(missing(num_infected_seeds))) {
    exp_design$num_infected_seeds <- rep(num_infected_seeds, nrow(exp_design))
  }
  
  if (!(missing(infected_seed_id))) {
    exp_design$infected_seed_id <- rep(infected_seed_id, nrow(exp_design))
  }

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

get_mean_non_truncated_gamma <- function(target_mean, shape) 
{
  tolerance <- 1.49e-4
  scale_est <- target_mean / shape 
  scale_params <- seq(scale_est / 2, scale_est * 2, by=0.00001)
  
  best_scale <- NaN
  best_mean <- 0 # FIXME is this ok to start? NaN not working 
  
  for (scale in scale_params) {
    cdf1 <- pgamma(0, shape=shape, scale=scale)
    cdf2 <- pgamma(1, shape=shape, scale=scale)
    
    f<-function(t){return(t*dgamma(t,shape = shape, scale=scale))}
    
    mean.tr<-(integrate(f,lower = 0,upper = 1)$value)/(cdf2-cdf1)
    
    if (abs(mean.tr - target_mean) < tolerance) {
      if (abs(mean.tr - target_mean) < abs(best_mean - target_mean)) {
        best_mean <- mean.tr
        best_scale <- scale 
      }
      
    }
  }
  
  return(best_scale * shape)
}

