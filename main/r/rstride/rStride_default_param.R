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
#  Copyright 2026
############################################################################ #
#
# rStride parameter information and baseline COVID-19 values
#
############################################################################ #

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# get default parameter values (flat list) to combine in a full-factorial grid.
# kept as a plain name -> value list for backward compatibility: callers
# override individual elements directly (e.g. exp_param_list$num_days <- 2)
# and pass the whole list into expand.grid()-based design generation.
get_default_param <- function(){
   return(get_param_values(get_default_param_definitions()))
}

# get a data.frame with the default value and description of each parameter,
# for documentation purposes (e.g. to render a parameter table in a report).
get_default_param_info <- function(){
   return(get_param_info(get_default_param_definitions()))
}

# master list of default parameters: each parameter is a list with a
# 'value' and a 'description', so the two can never drift out of sync.
# note that c_str() is used to concatenate strings to keep them together in the design of experiments
get_default_param_definitions <- function(){

   ## parameters from a20201031_132800_param4_d73_05k_n10_parameter_pareto_incidence_single_hosp
   out <- list(
      r0 = list(value = 3.42,
                description = "Basic reproduction number (R0)"),
      num_days = list(value = 196,
                description = "Number of simulated days"),
      num_rng_seeds = list(value = 10,
                description = "Number of stochastic realizations (RNG seeds) per parameter combination"),
      num_infected_seeds = list(value = 263,
                description = "Number of initially infected/seed cases at day 0"),
      disease_config_file = list(value = "data/disease_covid19_lognorm.xml",
                description = "Path to the XML file with disease natural-history parameters"),
      population_file = list(value = "data/pop_belgium11M_c500_teachers_censushh.csv",
                description = "Path to the CSV file with the synthetic population"),
      age_contact_matrix_file = list(value = "data/contact_matrix_flanders_conditional_teachers.xml",
                description = "Path to the XML file with age-specific social contact rates"),
      start_date = list(value = '2020-02-16',
                description = "Simulation start date"),
      holidays_file = list(value = 'data/calendar_belgium_2019_2021.csv',
                description = "Path to the CSV file with the school holiday/calendar schedule"),
      cnt_intensity_householdCluster = list(value = 0,
                description = "Relative contact intensity for household cluster contacts (0 = disabled)"),

      # tracing
      detection_probability = list(value = 0,
                description = "Probability that a symptomatic case is detected/tested (0 = no detection/testing)"),
      tracing_efficiency_household = list(value = 0,
                description = "Probability that a household contact of a detected case is successfully traced"),
      tracing_efficiency_other = list(value = 0,
                description = "Probability that a non-household contact of a detected case is successfully traced"),
      case_finding_capacity = list(value = 0,
                description = "Maximum number of cases that can be processed by contact tracing per day (NA: no limit)"),
      delay_isolation_index = list(value = NA,
                description = "Delay (in days) between symptom onset/detection and isolation of the index case"),
      delay_contact_tracing = list(value = NA,
                description = "Delay (in days) between index case detection and notification of traced contacts"),
      test_false_negative = list(value = 0,
                description = "False-negative rate of the diagnostic test used to confirm cases"),

      # log level
      event_log_level = list(value = "Transmissions",
                description = "Level of detail for the STRIDE event log: None, Incidence, Transmissions, Participants"),

      # factor for parameter estimation and fitting
      hosp_probability_factor = list(value = 1,
                description = "Scaling factor applied to the age-specific hospitalization probabilities, used for calibration/fitting"),

      # hospital admissions (relative proportions)
      hospital_category_age = list(value = paste(c(seq(0,80,10)),collapse=','),
                description = "Age breakpoints (comma-separated) defining the age categories used for hospitalization probability/delay"),
      hospital_probability_age = list(value = paste(c(0.091,0.009,0.044,0.033,0.057,0.075,0.143,0.373,1.000),collapse=','),
                description = "Age-specific relative hospitalization probabilities (comma-separated), based on hospital survey data by age (Faes et al); updated 2020-10-19 using hospital admissions in week 11-13 vs. simulated symptomatic cases by age from the 2020-09-17 R0 calibration"),
      hospital_mean_delay_age = list(value = paste(3,3,7,7,7,7,6,6,1,sep=','),
                description = "Age-specific mean delay (days) between symptom onset and hospital admission (comma-separated)"),
      # stochastic compartment model 
      hospital_length_of_stay = list(value = NA,
                description = "Mean hospital length of stay (days)"),

      # threshold for log parsing (default is NA == no threshold)
      logparsing_cases_upperlimit = list(value = NA,
                description = "Upper limit on case counts used when parsing simulation logs (NA = no threshold)")
   )

   # set social contact restriction parameters
   date_t0                    <- NA # date to start for example lockdown
   cnt_reduction_workplace    <- NA
   cnt_reduction_school       <- NA
   cnt_reduction_other        <- NA
   compliance_delay_workplace <- NA
   compliance_delay_school    <- NA
   compliance_delay_other     <- NA

   # include reduced contact parameters: school
   out$distancing_school_date  <- list(value = c_str(paste(date_t0)),
                                       description = "Relative reduction in school contacts: start date(s)")
   out$distancing_school_ratio <- list(value = c_str(cnt_reduction_school),
                                       description = "Relative reduction in school contacts: value(s)")
   out$distancing_school_delay <- list(value = c_str(compliance_delay_school),
                                       description = "Compliance delay (in days) for reduced school contacts")
   
   # include reduced contact parameters: workplace
   out$distancing_workplace_date  <- list(value = c_str(paste(date_t0)),
                                          description = "Relative reduction in workplace contacts: start date(s)")
   out$distancing_workplace_ratio <- list(value = c_str(cnt_reduction_workplace),
                                          description = "Relative reduction in workplace contacts: value(s)")
   out$distancing_workplace_delay <- list(value = c_str(compliance_delay_workplace),
                                          description = "Compliance delay (in days) for reduced workplace contacts")
   
   # include reduced contact parameters: community
   out$distancing_community_date  <- list(value = c_str(paste(date_t0)),
                                          description = "Relative reduction in community contacts: start date(s)")
   out$distancing_community_ratio <- list(value = c_str(cnt_reduction_workplace),
                                          description = "Relative reduction in community contacts: value(s)")
   out$distancing_community_delay <- list(value = c_str(compliance_delay_other),
                                          description = "Compliance delay (in days) for reduced community contacts")
   
   # include temporal household clustering
   out$household_clustering_date  <- list(value = c_str(out$start_date$value),
                                          description = "Start date of household clustering")
   out$household_clustering_ratio <- list(value = NA,
                                          description = "Probability for social contacts within household cluster, for example 4/7 to suggest presence 4 out of 7 days")
   out$household_clustering_delay <- list(value = c_str(0),
                                          description = "Compliance delay (in days) for household clustering")
   
   # import of external infections (the population size remains constant)
   out$imported_cases_date   <- list(value = NA, 
                                     description = "Start date of imported infections, until the end of the simulation, or a new value is provided")
   out$imported_cases_number <- list(value = NA,
                                     description = "Number of imported infections per introduction event")
   out$imported_cases_delay  <- list(value = NA,
             description = "Delay (in days) untill the specified introductions are reached (linear increase)")
   
   # social contact survey parameters 
   out$num_participants_survey <- list(value = 5000, 
                                       description = "Number of participants in the social contact survey, only active if event_log_level = Participants")
   out$contact_survey_dates    <- list(value = NA, 
                                       description = "Single days on which the contact survey is held, if NA: all simulation days. Note that multiple dates should be oncatenated to keep them togheter in the experimental design. For example: c_str('2020-03-16','2020-03-17')")
   out$contact_survey_ages     <- list(value = NA, 
                                       description = "Age breakpoints for the contact survey. For example c_str(0,18,110)")
   out$contact_survey_resample <- list(value = 0, 
                                       description = "Sample new set of participants each survey day? (0: no, so keep the same participatns or 1: yes)")
   # survey quota [0-1] for symptomatic participants
   out$contact_survey_quota_symptomatic <- list(value = 0,
                                                description = "Survey quota in terms of fraction symptomatic participants (0 = no quota, fully at random, 1 = only symptomatic participants)")
   
   # number of parallel workers (on UA cluster)
   out$num_parallel_workers <- list(value = 8,
             description = "Number of parallel worker processes to use when running the design of experiments")

   # reference data, used in HealthAgencyData.R
   out$reference_hospital_data_file <- list(value = 'data/covid19_hospital_age_2020_full.csv',
                                            description = "File name for observed hospital admissions by age over time")
   out$reference_serology_data_file <- list(value = 'data/covid19_serology_BE_reference.csv',
                                            description = "File name for observed serology levels by age over time")
   
   # return parameter definitions (value + description)
   return(out)
}

# strip a parameter-definitions list (name -> list(value, description)) down
# to a plain name -> value list, e.g. for use as exp_param_list.
get_param_values <- function(param_def_list){
   return(setNames(lapply(param_def_list, function(x) x$value), names(param_def_list)))
}

# turn a parameter-definitions list (name -> list(value, description)) into
# a data.frame with one row per parameter, e.g. for documentation.
get_param_info <- function(param_def_list){
   return(data.frame(parameter   = names(param_def_list),
                      value       = sapply(param_def_list, function(x) paste(x$value, collapse = ',')),
                      description = sapply(param_def_list, function(x) x$description),
                      stringsAsFactors = FALSE,
                      row.names = NULL))
}

