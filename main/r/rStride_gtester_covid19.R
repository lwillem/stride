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
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_gtester_covid19.R 
#
############################################################################ #

# Clear work environment
rm(list=ls())

# Load rStride
source('./bin/rstride/rStride.R')
source('./bin/rStride_covid19_default_param.R')

# set directory postfix (optional)
dir_postfix <- '_gtester'

################################## #
## DESIGN OF EXPERIMENTS        ####
################################## #

# uncomment the following line to inspect the config xml tags
#names(xmlToList('./config/run_default.xml'))

# set the number of realisations per configuration set
num_rng_seeds  <- 5

# add parameters and values to combine in a full-factorial grid
exp_design_base <- expand.grid(r0                       = 2.5,
                          num_days                      = 31,
                          rng_seed                      = seq(num_rng_seeds),
                          num_participants_survey       = 10,   
                          num_infected_seeds            = 540,
                          disease_config_file           = 'data/disease_covid19_age.xml',
                          population_file               = 'data/pop_belgium600k_c500_teachers_censushh.csv',
                          age_contact_matrix_file       = 'data/contact_matrix_flanders_conditional_teachers.xml',
                          start_date                    = '2020-03-04',
                          holidays_file                 = 'data/holidays_belgium_2019_2021.csv',

                          detection_probability          = 0,
                          tracing_efficiency_household   = 0,
                          tracing_efficiency_other       = 0,
                          case_finding_capacity          = 0,
                          test_false_negative            = 0,
                          
                          gtester_label                  = 'covid_base',
                          event_log_level                = 'Transmissions',
                          
                          # hospital_category_age         = paste(0,19,60,80,sep=','),
                          # hospital_probability_age      = paste(0.049,0.03024,0.1197,0.5922,sep=','),
                          # hospital_mean_delay_age       = paste(3,7,7,6,sep=','),
                          hospital_category_age         = paste(0,sep=','),
                          hospital_probability_age      = paste(0,sep=','),
                          hospital_mean_delay_age       = paste(0,sep=','),
                          hospital_length_of_stay       = 0,
                          
                          disease_susceptibility_age      = NA,
                          disease_susceptibility_agecat   = NA,
                          transmission_probability_distribution = NA,
                          transmission_probability_distribution_overdispersion = NA,

                          transmission_probability = NA,
                          
                          distancing_workplace_ratio = NA,
                          distancing_workplace_date  = NA,
                          distancing_workplace_delay = NA,
                          distancing_school_ratio    = NA,
                          distancing_school_date     = NA,
                          distancing_school_delay    = NA,
                          distancing_community_ratio = NA,
                          distancing_community_date  = NA,
                          distancing_community_delay = NA,
                          
                          imported_cases_number = NA,
                          imported_cases_date   = NA,
                          imported_cases_delay  = NA,
                          
                          household_clustering_ratio = NA,
                          household_clustering_delay = NA,
                          household_clustering_date = NA,
                      
                          contact_survey_dates = NA,
                          contact_survey_ages = c_str(0,18,110),
                          
                          subpools_community_file = NA,
                          subpools_community_used = NA,
                          pool_characteristics_file = NA,
                          airborne_transmission = NA,
                          stringsAsFactors = F)

# Contacts: virtual survey ---- 
exp_design_all <- exp_design_base
exp_design_all$event_log_level            <- 'Participants'
exp_design_all$contact_survey_dates       <- c_str(paste(as.Date(exp_design_base$start_date[1])+(1:exp_design_base$num_days[1])-1))
exp_design_all$gtester_label              <- 'covid_logParticipants'

# no logging ----
exp_design_none <- exp_design_base
exp_design_none$event_log_level            <- 'None'
exp_design_none$gtester_label              <- 'covid_none'

# hospital admission ----
exp_design_hosp <- exp_design_base
exp_design_hosp$hospital_category_age         <- paste(0,19,60,80,sep=',')
exp_design_hosp$hospital_probability_age      <- paste(0.049,0.03024,0.1197,0.5922,sep=',')
exp_design_hosp$hospital_mean_delay_age       <- paste(3,7,7,6,sep=',')
exp_design_hosp$hospital_length_of_stay       <- 12
exp_design_hosp$gtester_label                 <- 'covid_hosp'

# daily seeding ----
exp_design_daily <- exp_design_base
exp_design_daily$imported_cases_number    <- c_str(10)
exp_design_daily$imported_cases_date      <- c_str('2020-03-02')
exp_design_daily$imported_cases_delay     <- c_str(0)
exp_design_daily$gtester_label            <- 'covid_daily'


# distancing ----
exp_design_dist <- exp_design_base
exp_design_dist$holidays_file              <- 'data/calendar_belgium_2020_covid19_exit_school_adjusted.csv'
exp_design_dist$distancing_workplace_ratio <- 0.3;
exp_design_dist$distancing_workplace_date  <- '2020-03-13';
exp_design_dist$distancing_workplace_delay <- 3  
exp_design_dist$distancing_community_ratio <- 0.4;
exp_design_dist$distancing_community_date  <- '2020-03-13';
exp_design_dist$distancing_community_delay <- 4 
exp_design_dist$gtester_label              <- 'covid_distancing'


# age_15min ----
exp_design_15min <- exp_design_base
exp_design_15min$disease_config_file     <- 'data/disease_covid19_age_15min.xml'
exp_design_15min$age_contact_matrix_file <- 'data/contact_matrix_flanders_conditional_teachers_15min.xml'
exp_design_15min$gtester_label           <- 'covid_15min'

# householdCluster ----
exp_design_hhcl <- exp_design_base
exp_design_hhcl$population_file       <- 'data/pop_belgium600k_c500_teachers_censushh_extended3_size2.csv'
exp_design_hhcl$holidays_file         <- 'data/calendar_belgium_2020_covid19_exit_schoolcategory_adjusted.csv'
# exp_design_hhcl$holidays_file         <- 'calendar_belgium_2020_covid19_exit_schoolcategory_adjusted_hhclustering.csv' # hardcoded ratio of 4/7
exp_design_hhcl$start_date            <- '2020-05-31'
exp_design_hhcl$household_clustering_ratio <- 4/7
exp_design_hhcl$household_clustering_delay <- 0
exp_design_hhcl$household_clustering_date <- '2020-06-01'
exp_design_hhcl$gtester_label         <- 'covid_hhcl'


# contact tracing ----
exp_design_cts <- exp_design_base
exp_design_cts$detection_probability        <- 0.5
exp_design_cts$holidays_file                <- 'data/calendar_belgium_2020_covid19_exit_schoolcategory_adjusted.csv'
exp_design_cts$start_date                   <- '2020-05-31'
exp_design_cts$tracing_efficiency_household <- 1.0
exp_design_cts$tracing_efficiency_other     <- 0.7
exp_design_cts$test_false_negative          <- 0.1
exp_design_cts$case_finding_capacity        <- 1000
exp_design_cts$event_log_level              <- 'Transmissions'
exp_design_cts$gtester_label                <- 'covid_tracing'

# contact tracing all
exp_design_cts_all <- exp_design_cts
exp_design_cts_all$event_log_level          <- 'ContactTracing'
exp_design_cts_all$gtester_label            <- 'covid_tracing_all'

# age-specific susceptibility: baseline ----
# note: this should provide exact the same results as 'covid_base'
exp_design_susceptible <- exp_design_base
exp_design_susceptible$gtester_label            <- 'covid_suscept'
tmp_susceptible  <- rep(1,100)
exp_design_susceptible$disease_susceptibility_agecat <- paste(0:99,collapse=',')
exp_design_susceptible$disease_susceptibility_age <- paste(tmp_susceptible,collapse=',')

# age-specific susceptibility: adapted
exp_design_susceptible_adapt <- exp_design_base
exp_design_susceptible_adapt$gtester_label            <- 'covid_suscept_adapt'
tmp_susceptible[-seq(1,91,9)] <- 0.90
exp_design_susceptible_adapt$disease_susceptibility_agecat <- paste(0:99,collapse=',')
exp_design_susceptible_adapt$disease_susceptibility_age <- paste(tmp_susceptible,collapse=',')

# individual-based transmission: baseline ----
exp_design_transm <- exp_design_base
exp_design_transm$gtester_label            <- 'covid_transm'
exp_design_transm$transmission_probability_distribution   <- 'Constant'
exp_design_transm$transmission_probability_distribution_overdispersion   <- 0

# individual-based transmission: adapted
exp_design_transm_adapt <- exp_design_base
exp_design_transm_adapt$gtester_label            <- 'covid_transm_gamma'
exp_design_transm_adapt$transmission_probability_distribution   <- 'Gamma'
exp_design_transm_adapt$transmission_probability_distribution_overdispersion   <- 0.8


# fitting transmission and susceptibility: baseline ----
exp_design_fitting <- exp_design_base
exp_design_fitting$gtester_label            <- 'covid_fitting'
# b0 <- 0.124492138353664; b1 <- 39.6458896077442            # from: disease_covid19_lognormal 
b0 <- 0.14743616688954;  b1 <- 43.9598287259418              # from: disease_covid19_age  
tmp_transmission <- rep(unique(exp_design_fitting$r0 - b0) / b1,100)
exp_design_fitting$disease_susceptibility_age <- paste(tmp_transmission,collapse=',')
exp_design_fitting$disease_susceptibility_agecat <- paste(0:99,collapse=',')
# exp_design_fitting$r0 = (1-b0)/b1
exp_design_fitting$transmission_probability = 1

# fitting transmission and susceptibility: modified
exp_design_fitting_adapt <- exp_design_fitting
exp_design_fitting_adapt$gtester_label            <- 'covid_fitting_adapt'
tmp_transmission[seq(0,18,1)]  <- 0.02
tmp_transmission[seq(60,70,1)] <- 0.08
exp_design_fitting_adapt$disease_susceptibility_age <- paste(tmp_transmission,collapse=',')

# fitting transmission and susceptibility: modified
exp_design_fitting_agegroup <- exp_design_fitting
exp_design_fitting_agegroup$gtester_label            <- 'covid_fitting_agegroup'
t_base <- unique(exp_design_fitting$r0 - b0) / b1
exp_design_fitting_agegroup$disease_susceptibility_age <- paste(c(0.02,t_base,0.08,t_base),collapse=',')
exp_design_fitting_agegroup$disease_susceptibility_agecat <- c('0,18,59,70')

# collectivity ----
exp_design_collectivity <- exp_design_base
exp_design_collectivity$population_file              <- 'data/pop_belgium600k_c500_teachers_censushh_collectivity.csv'
exp_design_collectivity$age_contact_matrix_file      <- 'data/contact_matrix_flanders_conditional_teachers_collectivity20.xml'
exp_design_collectivity$gtester_label                <- 'covid_collectivity'

# collectivity population, but in strict isolation
exp_design_collectivity_isolation <- exp_design_base
exp_design_collectivity_isolation$population_file    <- 'data/pop_belgium600k_c500_teachers_censushh_collectivity.csv'
exp_design_collectivity_isolation$gtester_label      <- 'covid_collectivity_isolation'

# collectivity mixing, default population
exp_design_collectivity_mixing <- exp_design_base
exp_design_collectivity_mixing$age_contact_matrix_file  <- 'data/contact_matrix_flanders_conditional_teachers_collectivity20.xml'
exp_design_collectivity_mixing$gtester_label            <- 'covid_collectivity_mixing'

# default param ----
exp_design_default_param <- exp_design_base
exp_design_default_param[,names(get_covid19_default_param())] <- get_covid19_default_param()
exp_design_default_param[,!names(exp_design_default_param) %in% names(exp_design_base)] <- NULL
exp_design_default_param$population_file              <- 'data/pop_belgium600k_c500_teachers_censushh.csv'
exp_design_default_param$num_days                     <- 61
exp_design_default_param$gtester_label                <- 'covid_default_param'
names(exp_design_base) %in% names(exp_design_default_param)
names(exp_design_default_param) %in% names(exp_design_base)

# subpools ----
exp_design_subpools <- exp_design_base
exp_design_subpools$population_file           <- 'data/pop_belgium10k_c500_teachers_censushh.csv'
exp_design_subpools$subpools_community_file   <- 'data/pop_belgium10k_c500_teachers_censushh_subpools_community.csv'
exp_design_subpools$subpools_community_used   <- 'true'
exp_design_subpools$gtester_label             <- 'covid_subpools'

# airborne transmission ----
exp_design_airborne <- exp_design_subpools
exp_design_airborne$pool_characteristics_file <- 'data/pool_characteristics.xml'
exp_design_airborne$airborne_transmission     <- 'true'
exp_design_airborne$gtester_label             <- 'covid_airborne'


# rbind all designs
exp_design <- rbind(exp_design_base, exp_design_all,
                    exp_design_cts_all, exp_design_cts,
                    exp_design_daily, exp_design_dist,
                    exp_design_15min, exp_design_hhcl,
                    exp_design_susceptible,exp_design_susceptible_adapt,
                    exp_design_transm,exp_design_transm_adapt,
                    exp_design_fitting,exp_design_fitting_adapt,
                    exp_design_fitting_agegroup,
                    exp_design_collectivity,exp_design_collectivity_isolation,
                    exp_design_collectivity_mixing,
                    exp_design_none, exp_design_hosp,
                    exp_design_default_param,
                    exp_design_subpools,exp_design_airborne)


# add a unique seed for each run
# note: the rng seeds don't change (anymore) with additional tests
exp_design$rng_seed <- 1:nrow(exp_design)
dim(exp_design)

# align rng seeds for "base" with "none", "suscept", "transm" and "fitting" tests
exp_design$rng_seed[grepl('covid_none',exp_design$gtester_label)]    <- exp_design$rng_seed[exp_design$gtester_label %in% c('covid_base')]
exp_design$rng_seed[grepl('covid_suscept',exp_design$gtester_label)] <- exp_design$rng_seed[exp_design$gtester_label %in% c('covid_base')]
exp_design$rng_seed[grepl('covid_transm',exp_design$gtester_label)]  <- exp_design$rng_seed[exp_design$gtester_label %in% c('covid_base')]
exp_design$rng_seed[grepl('covid_fitting',exp_design$gtester_label)] <- exp_design$rng_seed[exp_design$gtester_label %in% c('covid_base')]


# # selection? ----
#exp_design <- exp_design[exp_design$gtester_label %in% c('covid_base','covid_hosp','covid_subpools','covid_airborne'),]
#exp_design <- exp_design[exp_design$gtester_label %in% c('covid_base','covid_collectivity','covid_collectivity_isolation','covid_collectivity_mixing'),]
#exp_design <- exp_design[exp_design$gtester_label %in% c('covid_base','covid_fitting_base','covid_fitting_adapt'),]
#exp_design <- exp_design[exp_design$gtester_label %in% c('covid_base','covid_transm','covid_transm_gamma'),]
#exp_design <- exp_design[exp_design$gtester_label %in% c('covid_base','covid_default_param','covid_distancing'),]
 # exp_design <- exp_design[grepl('_base',exp_design$gtester_label) |
 #                            grepl('_collectivity',exp_design$gtester_label) |
 #                            grepl('_fitting',exp_design$gtester_label),]


table(exp_design$gtester_label)
################################## #
## RUN rSTRIDE                  ####
################################## #
project_dir <- run_rStride(exp_design               = exp_design,
                           dir_postfix              = dir_postfix,
                           ignore_stdout            = TRUE,
                           remove_run_output        = FALSE )


##################################### #
## RUN ABC METHODS ----
##################################### #
smd_print("START ABC FUNCTION TEST")
# get one parameter config, set workdir and save parameter RDS file
model_param_abc <- exp_design[exp_design$gtester_label %in% c('covid_base'),]
model_param_abc <- exp_design[1,]
model_param_abc$event_log_level <- "Incidence"
if(!exists('project_dir')){project_dir <- smd_file_path('sim_output/abc_test') }
setwd(project_dir)
saveRDS(model_param_abc,'model_param_update.rds')

# run rStride_abc
rstride_out_abc <- run_rStride_abc(c(rng_seed = 100, 
                                     r0 = 3, 
                                     num_infected_seeds= 400, 
                                     hosp_probability_factor=0.4,
                                     distancing_workplace_ratio=0.85,
                                     distancing_workplace_delay=7.4,
                                     distancing_community_ratio=0.85,
                                     distancing_community_delay=4.51
                                     )
                                   )

# restore workdir
setwd('../..')


##################################### #
## EXPLORE INPUT-OUTPUT BEHAVIOR   ####
##################################### #
inspect_summary(project_dir)
inspect_participant_data(project_dir)
inspect_incidence_data(project_dir)
inspect_prevalence_data(project_dir)
inspect_transmission_dynamics(project_dir)
inspect_tracing_data(project_dir)
#inspect_contact_data(project_dir)



##################################### #
## CHECK INPUT-OUTPUT              ####
##################################### #

# terminal message
smd_print('START REGRESSION TEST')

## Load project summary 
project_summary <- .rstride$load_project_summary(project_dir)

# CHECK summary: plot number of cases
plot_final_sizes <- function(project_summary){
  y_lim     <- range(pretty(c(project_summary$num_cases*0.9,project_summary$num_cases*1.1)))
  bplt_mean <- aggregate(num_cases ~ gtester_label,data=project_summary,mean)
  bplt_mean$num_cases <- round(bplt_mean$num_cases)
  par(mar=c(10,4,4,2))
  bplt <- boxplot(num_cases ~ gtester_label,data=project_summary,las=2,ylim=y_lim,xlab='')
  x_ticks_mean <- (1:ncol(bplt$stats))+0.2
  points(x = x_ticks_mean,
         y = bplt_mean$num_cases,
         pch = 8,
         col = 4)
  arrows(x0 = x_ticks_mean,
         y0 = bplt_mean$num_cases * 0.9,
         y1 = bplt_mean$num_cases * 1.1,
         col = 4, lwd = 2,length = 0
  )
  text(x = 1:ncol(bplt$stats),
       y = bplt_mean$num_cases*1.1,
       labels = bplt_mean$num_cases,
       cex=0.8,
       pos = 3,
       col=4)
  legend('bottom',
         c('mean',
           'mean ± 10%'),
         pch=c('*','I'),
         col=4,
         ncol=2)
  grid()
}
par(mfrow=c(1,1))
plot_final_sizes(project_summary)


str2id <- function(str){
  
  if(length(dim(str))==2){
    str_num <- matrix(0,nrow(str),ncol(str))
    for(i in 1:nrow(str)){
      for(j in 1:ncol(str)){
        str_num[i,j] <- str2id_base(str[i,j])
      }
    }
    return(str_num)
  } else{
    return(unlist(lapply(str,str2id_base)))
  }
}

str2id_base <- function(str){
  if(is.numeric(str)){
    return(str)
  }
  str <- tolower(as.character(str))
  return(sum(as.numeric(factor(unlist(strsplit(str, "")), levels = letters)),na.rm=T))
}

mean_by_exp_id <- function(project_output){
  
  # safety check
  if(all(is.na(project_output))){
    return(project_output)
  }
  
  # make sure all columns are numeric
  is_character <- grepl('Length',summary(project_output)[1,])
  if(any(is_character)){
    project_output[,is_character]   <- str2id(project_output[,is_character])
  }
 
  # replace 'NA' by '0' 
  project_output[is.na(project_output)] <- 0 
 
  # aggregate by calculating the mean
  project_output   <- aggregate(. ~ exp_id, data = project_output, mean)
  
  # return
  return(project_output)
}

compare_output <- function(project_dir,output_type){
  
  # load reference file names
  reference_file_names   <- dir('./tests',full.names = TRUE,pattern = '.rds')
  
  # make sure the provided output_type is valid
  if(!any(grepl(output_type,reference_file_names))){
    smd_print("ERROR in 'compare_output()' with invallid output_type:",output_type,WARNING = T)
    return(NULL)
  }
  
  ## Load project summary 
  project_summary                <- .rstride$load_project_summary(project_dir)

  # load new results
  project_output     <- .rstride$load_aggregated_output(project_dir,output_type)
  
  # define a boolean for the summary comparison
  bool_summary       <- output_type == 'summary'
  if(bool_summary){
    project_output <- project_summary
  }

  if(all(is.na(project_output)) || nrow(project_output) == 0){
    smd_print("NO OUTPUT TO COMPARE FOR:",output_type)
    return(NULL)
  }
      
  # load previous results  
  reference_output       <- readRDS(file=reference_file_names[grepl(output_type,reference_file_names)])
  reference_summary      <- readRDS(file=reference_file_names[grepl('summary.rds',reference_file_names)])
  
  # Do we have to select reference scenarios?
  if(nrow(project_summary) != nrow(reference_summary) && nrow(reference_output)>0){
    reference_summary <- reference_summary[reference_summary$gtester_label %in% unique(project_summary$gtester_label),]
    reference_output  <- reference_output[reference_output$exp_id %in% unique(reference_summary$exp_id),]
    if(bool_summary) smd_print("REGRESSION TEST DOES NOT CONTAIN ALL SCENARIOS",WARNING = T)
  }
  
  # Do we have to exclude new scenarios?
  if(nrow(project_summary) != nrow(reference_summary) && nrow(project_output)>0){
    project_summary <- project_summary[project_summary$gtester_label %in% unique(reference_summary$gtester_label),]
    project_output  <- project_output[project_output$exp_id %in% unique(project_summary$exp_id),]
    if(bool_summary)  smd_print("REGRESSION TEST HAS NEW SCENARIOS",WARNING = T)
  }
  

  # compare length and names, and adjust if possible
  if(length(project_output) != length(reference_output) || 
     length(setdiff(names(reference_output),names(project_output)))>0){
    smd_print(c('!! Model output has different columns or column names for type = ',output_type),WARNING = TRUE)
    
    smd_print(paste('!! NEW: ',paste(names(project_output)[!names(project_output) %in% names(reference_output)],collapse = ', ')), WARNING = TRUE)
    smd_print(paste('!! PREV: ',paste(names(reference_output)[!names(reference_output) %in% names(project_output)],collapse = ', ')), WARNING = TRUE)
  }
  
  # make sure that columns with identical names are compared   
  common_names     <- intersect(names(reference_output),names(project_output))
  project_output   <- project_output[,common_names]
  reference_output <- reference_output[,common_names]
  
  # make sure all columns are numeric
  is_character <- grepl('Length',summary(project_output)[1,])
  if(any(is_character)){
    project_output[,is_character]   <- str2id(project_output[,is_character])
    reference_output[,is_character] <- str2id(reference_output[,is_character])
  }
  
  # option to aggregate data
  if(output_type %in% c('contacts','participants') &&
     !any(is.na(project_output))){
    project_output <- mean_by_exp_id(project_output)
    reference_output <- mean_by_exp_id(reference_output)
  }
  
  # compare rows and adjust if possible
  if(nrow(project_output) != nrow(reference_output)){
    smd_print(paste0('!! Model output has different number of rows for type = ',output_type),
              paste0('[NEW: ', nrow(project_output),' -- PREV: ', nrow(reference_output),']'),
              WARNING = TRUE)
    smd_print('CONTINUE WITH AGGREGATED STATISTICS', WARNING = TRUE)
    
    project_output   <- mean_by_exp_id(project_output)
    reference_output <- mean_by_exp_id(reference_output)
  }
  
  # compare again length and names, but first select non-id columns
  col_select       <- which(!(grepl('_id',names(project_output))  | 
                              grepl('time',names(project_output)) |
                              grepl('run',names(project_output)) |
                              grepl('tag',names(project_output))))
  project_exp_id   <- project_output$exp_id
  project_output   <- project_output[,col_select]
  reference_output <- reference_output[,col_select]
  
  # compare output, both using compare.list() and '=='
  diff_list      <- !compare.list(project_output,reference_output)
  diff_operator  <- colSums(project_output != reference_output,na.rm=T) != 0
  elements_are_different  <- diff_list & diff_operator

  if(all(!elements_are_different)){
    smd_print(paste0("Model output '",output_type,"' did not change."))
  } else{
    
    #check for textual changes
    is_character[col_select]
    if(any(elements_are_different & is_character[col_select])){
      diff_character <- names(elements_are_different)[elements_are_different & is_character[col_select]]
      smd_print(paste0("Model output '",output_type,"' did change for element(s): ",diff_character), WARNING = TRUE)
      
      # make sure we can use rowSums (which requires 2 dimensions)
      diff_row           <- rowSums(as.matrix(project_output[,diff_character]) != as.matrix(reference_output[,diff_character])) > 0
      diff_gtester_label <- project_summary$gtester_label[project_summary$exp_id %in% project_exp_id[diff_row]]
      smd_print(paste0("Model output '",output_type,"' did change for gtester(s): "),paste(unique(diff_gtester_label),collapse =', '), WARNING = TRUE)
      
      # narrow down for numerical results
      elements_are_different[is_character[col_select]] <- FALSE
    }
 
    # if there are still changes, check for floating point issues
    if(any(elements_are_different)){
    digits_cutoff   <- 6
    elements_are_different_fp <- elements_are_different
    if(any(elements_are_different_fp)){
      for(i_elem in 1:length(project_output)){
        if(elements_are_different_fp[i_elem]){
          if(all(dim(project_output[[i_elem]]) == dim(reference_output[[i_elem]]))){
            order <- -(log10(abs(range(unlist(project_output[[i_elem]]) - unlist(reference_output[[i_elem]]),na.rm=T))))
            if(all(order > digits_cutoff)){
              elements_are_different_fp[i_elem] <- FALSE
            } 
          } 
        }
      }  
    }
    # report outcome of comparison
    if(!any(elements_are_different_fp)){
      smd_print(paste0("Model output '",output_type,"' did not differ more than 1e-",digits_cutoff))
    } else{
      
      # make sure we can use rowSums (which requires 2 dimensions)
      diff_row           <- rowSums(as.matrix(project_output[,elements_are_different_fp]) != as.matrix(reference_output[,elements_are_different_fp])) >0
      diff_gtester_label <- project_summary$gtester_label[project_summary$exp_id %in% project_exp_id[diff_row]]
      smd_print(paste0("Model output '",output_type,"' did substantially change for gtester(s): "),paste(unique(diff_gtester_label),collapse = ', '), WARNING = TRUE)
      
      if(sum(elements_are_different_fp)<10){
        smd_print(paste(c('with different results for:', names(project_output)[elements_are_different_fp]),collapse=' '), WARNING = TRUE
        )
      } else{
        smd_print("with at least more than 10 columns changed", WARNING = TRUE)
      }
      
      return(unique(diff_gtester_label))
    }
    }
  }
  return(NULL)
} # end function


# COMPARE SUMMARY ----
diff_gtester <- compare_output(project_dir,'summary')

# plot potential changes in number of cases
if(length(diff_gtester)>0){
  ref_project_summary  <- readRDS(file='tests/regression_rstride_summary.rds')

  summary_new <- project_summary[project_summary$gtester_label %in% diff_gtester,]
  summary_ref <- ref_project_summary[ref_project_summary$gtester_label %in% diff_gtester,]
  
  par(mar=c(8,4,4,2))
  y_lim <- range(pretty(c(summary_new$num_cases,summary_ref$num_cases)))
  bplt_ref <- boxplot(num_cases ~ gtester_label,
                data=summary_ref,main='CHANGES',ylim=y_lim, las=2,xlab='');grid()
  if(length(bplt_ref$name)==1) axis(1,at=1,labels=bplt_ref$names)
  bplt_new <- boxplot(num_cases ~ gtester_label,
                      data=summary_new,add=T,
                      border=2,
                      col=alpha(2,0.4),main='',ylim=y_lim,las=2,xlab='')  ;
  bool_different <- colSums(bplt_new$stats != bplt_ref$stats) >0
  legend('topleft',c('reference','new','changed'),col=c(1,alpha(2,0.4),4),pch=c('I','I','*'),cex=0.8)
  points(which(bool_different)+0.5,bplt_new$stats[3,bool_different],col=4,pch='*',cex=3)

  par(mfrow=c(1,1),mar=c(8,4,4,2))
}


## COMPARE OTHER MODEL OUTPUT ----
compare_output(project_dir,"incidence")
compare_output(project_dir,"prevalence")
compare_output(project_dir,"contacts")
compare_output(project_dir,"participants")

## COMPARE ABC ----
ref_rstride_out_abc <- readRDS(file='tests/regression_rstride_out_abc.rds')
if(setequal(rstride_out_abc,ref_rstride_out_abc)){
  smd_print("rSTRIDE ABC OK")
} else{
  smd_print("rSTRIDE ABC CHANGED!",WARNING = T)
  stride_diff <- setdiff(rstride_out_abc,ref_rstride_out_abc)
  smd_print(names(stride_diff),WARNING = T)
}



# COMPARE PERFORMANCE ----
ref_project_summary  <- readRDS(file='tests/regression_rstride_summary.rds')
common_gtester_label <- intersect(project_summary$gtester_label,ref_project_summary$gtester_label)
current_run_times    <- aggregate(run_time ~ gtester_label, data= project_summary[project_summary$gtester_label %in% common_gtester_label,],mean)
previous_run_times   <- aggregate(run_time ~ gtester_label, data= ref_project_summary[ref_project_summary$gtester_label %in% common_gtester_label,],mean)
run_time_diff        <- current_run_times$run_time - previous_run_times$run_time
smd_print('Total run time and abs. difference (s):',
          round(sum(current_run_times$run_time/1e3),1), '::',
          round(sum(run_time_diff/1e3),1)
)
smd_print('Average run time and abs. difference  (s):',
          round(mean(current_run_times$run_time/1e3),1), '::',
          round(mean(run_time_diff/1e3),1)
)
smd_print('Test with highest time differenct:',
          current_run_times$gtester_label[order(run_time_diff)[1]])


# terminal message
smd_print('REGRESSION TEST COMPLETE')

# short call for "reset reference values"
rrv <- function(stride_repo_dir = 'tests'){
  
  saveRDS(.rstride$load_project_summary(project_dir),
          file=file.path(stride_repo_dir,'regression_rstride_summary.rds'))
  saveRDS(.rstride$load_aggregated_output(project_dir,'data_incidence'), 
          file=file.path(stride_repo_dir,'regression_rstride_incidence.rds'))
  saveRDS(.rstride$load_aggregated_output(project_dir,'data_prevalence'),
          file=file.path(stride_repo_dir,'regression_rstride_prevalence.rds'))
  saveRDS(rstride_out_abc,
          file=file.path(stride_repo_dir,'regression_rstride_out_abc.rds'))
  
  # store aggregated social contact survey data
  saveRDS(mean_by_exp_id(.rstride$load_aggregated_output(project_dir,'data_contacts')),
          file=file.path(stride_repo_dir,'regression_rstride_contacts.rds'))
  saveRDS(mean_by_exp_id(.rstride$load_aggregated_output(project_dir,'data_participants')),
          file=file.path(stride_repo_dir,'regression_rstride_participants.rds'))
  
  pdf(file=file.path(stride_repo_dir,'regression_rstride_cases.pdf'),14,7)
  plot_final_sizes(project_summary)
  dev.off()
  
  smd_print('NEW REFERENCE VALES STORED IN FOLDER:',stride_repo_dir)
}

# update the repository and local rStride reference values (note: local function for LW)
rrv_repo <- function(){
  stride_repo_dir <- '~/Documents/university/research/stride/repo/stride_2023/main/resources/rstride_test'
  rrv(stride_repo_dir = stride_repo_dir)
  rrv()
}



