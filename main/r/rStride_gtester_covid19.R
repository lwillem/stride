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
num_seeds  <- 5

# add parameters and values to combine in a full-factorial grid
exp_design_base <- expand.grid(r0                       = 2.5,
                          num_days                      = 31,
                          rng_seed                      = seq(num_seeds),
                          num_participants_survey       = 10,   
                          num_infected_seeds            = 540,
                          disease_config_file           = 'disease_covid19_age.xml',
                          population_file               = 'pop_belgium600k_c500_teachers_censushh.csv',
                          age_contact_matrix_file       = 'contact_matrix_flanders_conditional_teachers.xml',
                          start_date                    = '2020-03-04',
                          holidays_file                 = 'holidays_belgium_2019_2021.csv',

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
                      
                          stringsAsFactors = F)

# Contacts: virtual survey ---- 
exp_design_all <- exp_design_base
exp_design_all$event_log_level            <- 'Participants'
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
exp_design_dist$holidays_file              <- 'calendar_belgium_2020_covid19_exit_school_adjusted.csv'
exp_design_dist$distancing_workplace_ratio <- 0.3;
exp_design_dist$distancing_workplace_date  <- '2020-03-13';
exp_design_dist$distancing_workplace_delay <- 3  
exp_design_dist$distancing_community_ratio <- 0.4;
exp_design_dist$distancing_community_date  <- '2020-03-13';
exp_design_dist$distancing_community_delay <- 4 
exp_design_dist$gtester_label              <- 'covid_distancing'


# age_15min ----
exp_design_15min <- exp_design_base
exp_design_15min$disease_config_file     <- 'disease_covid19_age_15min.xml'
exp_design_15min$age_contact_matrix_file <- 'contact_matrix_flanders_conditional_teachers_15min.xml'
exp_design_15min$gtester_label           <- 'covid_15min'

# householdCluster ----
exp_design_hhcl <- exp_design_base
exp_design_hhcl$population_file       <- 'pop_belgium600k_c500_teachers_censushh_extended3_size2.csv'
exp_design_hhcl$holidays_file         <- 'calendar_belgium_2020_covid19_exit_schoolcategory_adjusted.csv'
# exp_design_hhcl$holidays_file         <- 'calendar_belgium_2020_covid19_exit_schoolcategory_adjusted_hhclustering.csv' # hardcoded ratio of 4/7
exp_design_hhcl$start_date            <- '2020-05-31'
exp_design_hhcl$household_clustering_ratio <- 4/7
exp_design_hhcl$household_clustering_delay <- 0
exp_design_hhcl$household_clustering_date <- '2020-06-01'
exp_design_hhcl$gtester_label         <- 'covid_hhcl'


# contact tracing ----
exp_design_cts <- exp_design_base
exp_design_cts$detection_probability        <- 0.5
exp_design_cts$holidays_file                <- 'calendar_belgium_2020_covid19_exit_schoolcategory_adjusted.csv'
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
exp_design_collectivity$population_file              <- 'pop_belgium600k_c500_teachers_censushh_collectivity.csv'
exp_design_collectivity$age_contact_matrix_file      <- 'contact_matrix_flanders_conditional_teachers_collectivity20.xml'
exp_design_collectivity$gtester_label                <- 'covid_collectivity'

# collectivity population, but in strict isolation
exp_design_collectivity_isolation <- exp_design_base
exp_design_collectivity_isolation$population_file    <- 'pop_belgium600k_c500_teachers_censushh_collectivity.csv'
exp_design_collectivity_isolation$gtester_label      <- 'covid_collectivity_isolation'

# collectivity mixing, default population
exp_design_collectivity_mixing <- exp_design_base
exp_design_collectivity_mixing$age_contact_matrix_file  <- 'contact_matrix_flanders_conditional_teachers_collectivity20.xml'
exp_design_collectivity_mixing$gtester_label            <- 'covid_collectivity_mixing'

# default param ----
exp_design_default_param <- exp_design_base
exp_design_default_param[,names(get_covid19_default_param())] <- get_covid19_default_param()
exp_design_default_param[,!names(exp_design_default_param) %in% names(exp_design_base)] <- NULL
exp_design_default_param$population_file              <- 'pop_belgium600k_c500_teachers_censushh.csv'
exp_design_default_param$num_days                     <- 61
exp_design_default_param$gtester_label                <- 'covid_default_param'
names(exp_design_base) %in% names(exp_design_default_param)
names(exp_design_default_param) %in% names(exp_design_base)

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
                    exp_design_default_param)


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
# exp_design <- exp_design[exp_design$gtester_label %in% c('covid_base'),]
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
project_summary$run_time_id    <- project_summary$run_time
project_summary$output_prefix  <- NULL
project_summary$run_tag        <- NULL
project_summary$run_time       <- NULL
project_summary$total_time     <- NULL

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
# plot_final_sizes(ref_project_summary)

# load incidence output
data_incidence     <- .rstride$load_aggregated_output(project_dir,'data_incidence')
dim(data_incidence)

# get prevalence output
data_prevalence <- .rstride$load_aggregated_output(project_dir,'data_prevalence')
dim(data_prevalence)

## Load reference data
ref_project_summary  <- readRDS(file='tests/regression_rstride_summary.rds')
ref_data_incidence   <- readRDS(file='tests/regression_rstride_incidence.rds')
ref_data_prevalence  <- readRDS(file='tests/regression_rstride_prevalence.rds')

# Do we have to select reference scenarios?
if(nrow(project_summary) != nrow(ref_project_summary)){
  ref_project_summary <- ref_project_summary[ref_project_summary$gtester_label %in% unique(project_summary$gtester_label),]
  ref_data_incidence  <- ref_data_incidence[ref_data_incidence$exp_id %in% unique(ref_project_summary$exp_id),]
  ref_data_prevalence <- ref_data_prevalence[ref_data_prevalence$exp_id %in% unique(ref_project_summary$exp_id),]

  # adjust exp_id (i.e. this is based on the number of experiments, but make sure the other parameters are similar)
  ref_project_summary$exp_id <- 1:nrow(ref_project_summary)

  smd_print("REGRESSION TEST DOES NOT CONTAIN ALL SCENARIOS",WARNING = T)
}

# fix: use base name of file names (and get rid of project-specific file paths)
# note: this is not needed any more (2024-01-10)
ref_project_summary$holidays_file <- basename(ref_project_summary$holidays_file)
project_summary$holidays_file     <- basename(project_summary$holidays_file)



## COMPARE SUMMARY ----

if(!setequal(project_summary[,!grepl('_id',names(project_summary))],
             ref_project_summary[,!grepl('_id',names(ref_project_summary))])){ 
  
  smd_print("SUMMARY CHANGED",WARNING = T)
  
  # remove new rows
  select_project_summary <- project_summary[project_summary$gtester_label %in% unique(ref_project_summary$gtester_label),]
  
  # check columns
  if(all(dim(select_project_summary) == dim(ref_project_summary))){
    col_changed <- which(colSums(select_project_summary != ref_project_summary) > 0)
    smd_print('column(s) with changes:', paste(names(col_changed),collapse = ','),WARNING = T)
  } else{
    smd_print(paste(c('Summary dimensions changed:',setdiff(names(select_project_summary),names(ref_project_summary))),collapse='\n\t\t'),WARNING = T)
  }
  
  # remove new column names
  select_project_summary <- project_summary[,names(project_summary) %in% names(ref_project_summary)]
  
  # remove redundant exp
  ref_project_summary <- ref_project_summary[,names(ref_project_summary) %in% names(select_project_summary)]
  
  # get difference (excluding _id columns)
  diff_summary    <- setdiff(select_project_summary[,!grepl('_id',names(select_project_summary))],
                             ref_project_summary[,!grepl('_id',names(select_project_summary))])
  if(length(diff_summary)>0 && all(dim(select_project_summary) == dim(ref_project_summary))){
    smd_print('CHANGES: ',paste(names(diff_summary),collapse = ', '),WARNING = T)
    flag <- (select_project_summary[,names(diff_summary)] != ref_project_summary[,names(diff_summary)])
   if(length(diff_summary)>1) {
     flag <- rowSums(flag)>0
   } 
    smd_print('EXP_ID with changes:', paste(unique(select_project_summary$gtester_label[flag]),collapse = ','))
    select_project_summary[flag,c('gtester_label',names(diff_summary))] ==
    ref_project_summary[flag,c('gtester_label',names(diff_summary))]
    
    if("num_cases" %in% names(diff_summary)){
      #par(mfrow=c(1,2),mar=c(8,4,4,2))
      par(mar=c(8,4,4,2))
      y_lim <- range(pretty(c(ref_project_summary$num_cases,select_project_summary$num_cases)))
      bplt_ref <- boxplot(num_cases ~ gtester_label,
                          data=ref_project_summary,main='REFERENCE',ylim=y_lim, las=2,xlab='');grid()
      boxplot(num_cases ~ gtester_label,
              data=ref_project_summary,main='BOTH',ylim=y_lim, las=2,xlab='');grid()
      bplt_new <- boxplot(num_cases ~ gtester_label,
                          data=select_project_summary,add=T,
                          border=2,
                          col=alpha(2,0.4),main='',ylim=y_lim,las=2,xlab='')  ;
      bool_different <- colSums(bplt_new$stats != bplt_ref$stats) >0
      legend('topleft',c('reference','new','changed'),col=c(1,alpha(2,0.4),4),pch=c('I','I','*'),cex=0.8)
      points(which(bool_different)+0.5,bplt_new$stats[3,bool_different],col=4,pch='*',cex=3)
      
      par(mfrow=c(1,1),mar=c(8,4,4,2))
    }
  }
  #print(head(diff_summary))
} else{
  smd_print("SUMMARY OK")
}

## COMPARE INCIDENCE ----
if(setequal(data_incidence[,names(data_incidence) != 'exp_id'], 
            ref_data_incidence[,names(ref_data_incidence) != 'exp_id'])){
  smd_print("INCIDENCE OK")
} else{
  
  missing_colnames_new <- names(data_incidence)[!names(data_incidence) %in% names(ref_data_incidence)]
  if(length(missing_colnames_new)>0){
    smd_print('INCIDENCE columns added:', paste(missing_colnames_new,collapse = ','),WARNING = T)
  }
  
  missing_colnames_ref <- names(ref_data_incidence)[!names(ref_data_incidence) %in% names(data_incidence)]
  if(length(missing_colnames_ref)>0){
    smd_print('INCIDENCE columns missing:', paste(missing_colnames_ref,collapse = ','),WARNING = T)
  }
  
  compare_col <- names(ref_data_incidence)[!names(ref_data_incidence) %in% c('exp_id',missing_colnames_ref,missing_colnames_new)]
  diff_incidence  <- setdiff(data_incidence[,compare_col],ref_data_incidence[,compare_col])
  if(length(diff_incidence)>0){ 
    smd_print("INCIDENCE CHANGED",WARNING = T)
    
    diff_incidence_colnames <- names(diff_incidence)[names(diff_incidence) %in% names(ref_data_incidence)]
    diff_incidence_colnames <- unique(gsub('_age.*','',diff_incidence_colnames))
    smd_print(diff_incidence_colnames,WARNING = T)
    
    # include "exp_id" column, to make sure there are at least 2 columns for the rowSums
    diff_incidence$exp_id <- data_incidence$exp_id
    
    if(all(dim(data_incidence) == dim(ref_data_incidence))){
      flag <- rowSums(data_incidence[,names(diff_incidence)] != ref_data_incidence[,names(diff_incidence)],na.rm=T)>0
      smd_print('EXP_ID with changes:', paste(unique(data_incidence$exp_id[flag]),collapse = ','))
      # data_incidence[flag,names(diff_incidence)]
      smd_print('gtester_label with changes:', paste(unique(project_summary$gtester_label[project_summary$exp_id %in% data_incidence$exp_id[flag]]),collapse = ','))
      # ref_data_incidence[flag,names(diff_incidence)]    
      
      # db_changes <- data.frame(exp_id = data_incidence$exp_id,
      #                          num_changes = rowSums(data_incidence[,names(diff_incidence)] != ref_data_incidence[,names(diff_incidence)],na.rm=T))
      # db_changes <- merge(db_changes,project_summary,by='exp_id')
      # aggregate(num_changes ~ gtester_label,data=db_changes,sum)
      
      # head(data_incidence[,bool_colnames])
      # head(ref_data_incidence[,bool_colnames])
      #head(data_incidence[names(diff_incidence)])
      
    } else { # dimensions changed!!
      smd_print('INCIDENCE ISSUE: dimensions changed',WARNING = TRUE)

    }
  } else {
    smd_print("OVERLAPPING INCIDENCE OK",WARNING = FALSE)
  }
}

## COMPARE PREVALENCE ----
sel_col <- names(data_prevalence)[names(data_prevalence) %in% names(ref_data_prevalence) & names(data_prevalence) != 'exp_id'] # make sure the same columns are compared
if(length(sel_col)>0 && setequal(data_prevalence[,sel_col],
            ref_data_prevalence[,sel_col])){ 
  smd_print("PREVALENCE OK")
} else{
  smd_print("PREVALENCE CHANGED",WARNING = T)
  #diff_prevalence <- setdiff(data_prevalence,ref_data_prevalence)
  #print(head(diff_prevalence))
}


## COMPARE ABC ----
ref_rstride_out_abc <- readRDS(file='tests/regression_rstride_out_abc.rds')
if(setequal(rstride_out_abc,ref_rstride_out_abc)){
  smd_print("rSTRIDE ABC OK")
} else{
  smd_print("rSTRIDE ABC CHANGED!",WARNING = T)
  stride_diff <- setdiff(rstride_out_abc,ref_rstride_out_abc)
  smd_print(names(stride_diff),WARNING = T)
}

# terminal message
smd_print('REGRESSION TEST COMPLETE')

# short call for "reset reference values"
rrv <- function(stride_repo_dir = 'tests'){
  
  saveRDS(project_summary,
          file=file.path(stride_repo_dir,'regression_rstride_summary.rds'))
  saveRDS(.rstride$load_aggregated_output(project_dir,'data_incidence'), 
          file=file.path(stride_repo_dir,'regression_rstride_incidence.rds'))
  saveRDS(.rstride$load_aggregated_output(project_dir,'data_prevalence'),
          file=file.path(stride_repo_dir,'regression_rstride_prevalence.rds'))
  saveRDS(rstride_out_abc,
          file=file.path(stride_repo_dir,'regression_rstride_out_abc.rds'))
  
  pdf(file=file.path(stride_repo_dir,'regression_rstride_summary.pdf'),14,7)
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

# print total run-time
smd_print('total and mean run time (s):', 
          round(mean(project_summary$run_time_id)/1e3,1),
          ':',
          round(sum(project_summary$run_time_id)/1e3,1)
          )
