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
#  Copyright 2020, Willem L, Libin P
############################################################################ #
#
# Baseline settings for rStride intervention scenarios
#
############################################################################ #

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# get default parameter values to combine in a full-factorial grid
get_exp_param_default <- function(bool_child_param = FALSE, 
                                  bool_min_restrictive = FALSE,
                                  bool_revised_model_param = FALSE){
   
   # create calendar files
   create_calendar_files()
   
   out <- list(r0                            = seq(3.4,3.4,0.1),
                num_days                      = 196,
                num_seeds                     = 10,
                num_participants_survey       = 30,
                num_infected_seeds            = 510,
                disease_config_file           = "disease_covid19_age.xml",
                population_file               = c("pop_belgium11M_c500_teachers_censushh.csv"),
                age_contact_matrix_file       = "contact_matrix_flanders_conditional_teachers.xml",
                start_date                    = c('2020-02-17'),
                holidays_file                 = 'calendar_belgium_2020_covid19_exit_school_adjusted.csv',
                telework_probability          = 0,
                cnt_reduction_workplace       = 0.8,
                cnt_reduction_other           = 0.85,
                compliance_delay_workplace    = 6,
                compliance_delay_other        = 6,
                num_daily_imported_cases      = 0,
                cnt_reduction_workplace_exit  = c(0.50,0.75),  
                cnt_reduction_other_exit      = c(0.70,0.85),
                cnt_reduction_school_exit     = 0.5,
                cnt_reduction_collectivity    = 0.0, # no reduction, during lockdown
                cnt_baseline_collectivity     = 0.0, # no reduction, by default
                cnt_intensity_householdCluster = 0,
                detection_probability          = 0,
                tracing_efficiency_household   = 0.9, 
                tracing_efficiency_other       = 0.7,
                case_finding_capacity          = 10000, # no limit at this stage
                delay_isolation_index          = 1,
                delay_contact_tracing          = 1, 
                test_false_negative            = 0.1,
               
                # log level
                event_log_level                 = "Transmissions",

                # factor for parameter estimation and fitting
                hosp_probability_factor        = 1,
               
                # unversal testing
                unitest_pool_allocation       = c("data/pop_belgium11M_c500_pool_allocation_$unitest_pool_size.csv"),
                unitest_fnr                   = c(0.01),
                unitest_n_tests_per_day       = 0, #c(25000),
                unitest_pool_size             = c(32),
                unitest_test_compliance       = c(0.9),
                unitest_isolation_compliance  = c(0.8),
               
               # hospital admissions
                hospital_category_age         = paste(0,19,60,80,sep=','),
                hospital_probability_age      = paste(0.049,0.03024,0.1197,0.5922,sep=','),
                hospital_mean_delay_age       = paste(3,7,7,6,sep=','),
                
               # threshold for log parsing (default is NA == no threshold)
               logparsing_cases_upperlimit    = NA
               
          )
   
    if(bool_revised_model_param){
      # relative proportions
      # reference: hospital survey data by age (faes et al) 
      # update on 19/10 : hospital admissions in week 11-13 / simulated sympt cases by age in R0 callibration 2020-09-17
      out$hospital_category_age         = paste(c(seq(0,80,10)),collapse=',')
      out$hospital_probability_age      = paste(c(0.091,0.009,0.044,0.033,0.057,0.075,0.143,0.373,1.000 ),collapse=',') # still requires rescaling
      out$hospital_mean_delay_age       = paste(3,3,7,7,7,7,6,6,1,sep=',')
      
      # disease history: literature based distributions
      out$disease_config_file <- 'disease_covid19_lognorm.xml'
      
      ## parameters from a20201031_132800_param4_d73_05k_n10_parameter_pareto_incidence_single_hosp
      out$r0 <- 3.42
      out$num_infected_seeds <- 263 
      out$hosp_probability_factor <- 0.40
      out$cnt_reduction_workplace <- 0.86 
      out$cnt_reduction_other     <- 0.85 
      out$compliance_delay_workplace <- 7 
      out$compliance_delay_other  <- 7

      #out$num_seeds <- NA
      
   }
   
   # change parameters if childrens infectiousness is 1/2 compared to adults
   if(bool_child_param){ 
      
      # update on 19/10 : hospital admissions in week 11-13 / simulated sympt cases by age in R0_child callibration 2020-09-17
      out$hospital_category_age         = paste(c(seq(0,80,10)),collapse=',')
      out$hospital_probability_age      = paste(c(0.221,0.018,0.041,0.033,0.056,0.071,0.139,0.367,1.000),collapse=',') # still requires rescaling
      out$hospital_mean_delay_age       = paste(3,3,7,7,7,7,6,6,1,sep=',')
      
      # parameters from a20201031_132751_param4_child_d73_05k_n10_parameter_pareto_incidence_single_hosp
      out$disease_config_file <- "disease_covid19_lognorm_child.xml"
      out$r0 <- 3.37
      out$num_infected_seeds <- 255
      out$hosp_probability_factor <- 0.35
      out$cnt_reduction_workplace <- 0.76
      out$cnt_reduction_other     <- 0.86
      out$compliance_delay_workplace <- 6
      out$compliance_delay_other  <- 7
   }
   
   # select least stringent social mixing assumptions
   if(bool_min_restrictive){
      out$cnt_reduction_workplace_exit <- min(out$cnt_reduction_workplace_exit)
      out$cnt_reduction_other_exit <- min(out$cnt_reduction_other_exit)
   }  
   
   
   # number of paralel workers
   out$num_parallel_workers <- 50
   
   # return parameters
   return(out)
}

