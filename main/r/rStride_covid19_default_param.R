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
#  Copyright 2024
############################################################################ #
#
# Baseline settings for rStride COVID-19 intervention scenarios
#
############################################################################ #

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# get default parameter values to combine in a full-factorial grid
get_covid19_default_param <- function(){
   
   ## parameters from a20201031_132800_param4_d73_05k_n10_parameter_pareto_incidence_single_hosp
   out <- list(r0                             = 3.42,
                num_days                      = 196,
                num_seeds                     = 10,
                num_participants_survey       = 30,
                num_infected_seeds            = 263,
                disease_config_file           = "disease_covid19_lognorm.xml",
                population_file               = "pop_belgium11M_c500_teachers_censushh.csv",
                age_contact_matrix_file       = "contact_matrix_flanders_conditional_teachers.xml",
               
                # update 2024-01-12: shifted from 17/2 to 16/2 because infected cases are now introduced on day 0 instead of -1 
                start_date                    = '2020-02-16', 
                # holidays_file                 = 'calendar_belgium_2020_covid19_exit_school_adjusted.csv',
                holidays_file                 = 'calendar_belgium_2019_2021.csv',
                cnt_intensity_householdCluster = 0,
               
                # tracing 
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
                hosp_probability_factor        = 0.40,
               
               # hospital admissions (relative proportions)
               # reference: hospital survey data by age (Faes et al) 
               # update on 19/10 : hospital admissions in week 11-13 / simulated sympt cases by age in R0 calibration 2020-09-17
                hospital_category_age         = paste(c(seq(0,80,10)),collapse=','),
                hospital_probability_age      = paste(c(0.091,0.009,0.044,0.033,0.057,0.075,0.143,0.373,1.000 ),collapse=','),
                hospital_mean_delay_age       = paste(3,3,7,7,7,7,6,6,1,sep=','),
               
               # stochastic compartment model ('OCTO')  
               hospital_length_of_stay        = 12,
                
               # threshold for log parsing (default is NA == no threshold)
               logparsing_cases_upperlimit    = NA
          )
   
   
   # 2020 lock down parameters
   date_t0                    <- as.Date('2020-03-13')
   cnt_reduction_workplace    <- 0.86
   cnt_reduction_school       <- 1
   cnt_reduction_other        <- 0.85
   compliance_delay_workplace <- 7
   compliance_delay_school    <- 0
   compliance_delay_other     <- 7
   
   # include 2020 lock-down parameters
   out$distancing_workplace_ratio         <- c_str(cnt_reduction_workplace)
   out$distancing_workplace_date          <- c_str(paste(date_t0))
   out$distancing_workplace_delay         <- c_str(7)
   
   out$distancing_community_ratio        <- c_str(cnt_reduction_other)
   out$distancing_community_date         <- c_str(paste(date_t0))
   out$distancing_community_delay        <- c_str(7)
   
   out$distancing_school_ratio        <- 1
   out$distancing_school_date         <- paste(date_t0) 
   out$distancing_school_delay        <- 1 # school closure did not start on Friday 13/3.
   
   # number of parallel workers (on UA cluster)
   out$num_parallel_workers <- 50
   
   # # household clustering?
   # out$household_clustering_ratio_ratio <- c_str(4/7)
   # out$household_clustering_ratio_date  <- c_str(out$start_date)
   # out$household_clustering_ratio_delay <- c_str(0)
   
   # # start importing external infections ? (the population remains constant)
   # out$imported_cases_number   <- 10
   # out$imported_cases_date     <- '2020-08-25 # until the end of the simulation
   # out$imported_cases_delay    <- 5
   
   # # conduct social contact survey?
   # out$event_log_level         <- 'Participants'
   # out$num_participants_survey <- 5000        # number of participants
   # out$contact_survey_dates    <- c_str('2020-03-16','2020-03-17) # single days
   # out$contact_survey_ages     <- c_str(0,18,110)
   
   # return parameters
   return(out)
}

get_linear_approx <- function(level_start,level_end,delay){
   return(approx(y = c(level_start,level_end),x = c(0,delay),xout = 0:delay)$y)
}

