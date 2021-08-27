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
#  Copyright 2021, Willem L
############################################################################ #
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_abc.R 
#
############################################################################ #

# lw's developing option
if(any(grepl('lwillem',dir('~',full.names=T))) && exists('.rstride')){
  .rstride$set_wd()
}

# Clear work environment
rm(list=ls())

# Load rStride
source('./bin/rstride/rStride.R')

# Load default parameter configurations
source('./bin/rStride_intervention_baseline.R')

# load ABC package
library(EasyABC)

# set samples and cluster size
n_cluster = 8
n_sample = n_cluster * 4 #24

# set acceptance level
pacc=0.7

# weight of hospital data as reference
rel_importance_hosp_data <- 5

# command line arguments?
args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  job_id <- paste0(format(Sys.time(), format="%Y%m%d_"),args[1])
} else{
  job_id <-format(Sys.time(), format="%Y%m%d_%H%M_abc_age")
}
if(length(args)>=3){
  n_cluster = as.numeric(args[2])
  n_sample = n_cluster * as.numeric(args[3])
  pacc= as.numeric(args[4])
  rel_importance_hosp_data= as.numeric(args[5])
} 

# set directory postfix (== run tag)
dir_postfix <- paste0(job_id,'_c',n_cluster,'_n',n_sample,'_p',formatC(pacc*100,flag = 0,digits = 2),'_h',rel_importance_hosp_data)
dir_postfix

# set project directory
project_dir <- smd_file_path('./sim_output',dir_postfix)

################################################ #
## DESIGN OF EXPERIMENTS  ----
################################################ #

# add default parameters and values to combine in a full-factorial grid
model_param_update <- get_exp_param_default(bool_age_specific_param = T,
                                            bool_min_restrictive = T)

# additional model parameters
#model_param_update$population_file <- "pop_belgium600k_c500_teachers_censushh_collectivity.csv"
model_param_update$num_days        <- 319

ref_period <- seq(as.Date('2020-03-15'),
                  as.Date(model_param_update$start_date) + model_param_update$num_days-1,
                  1)
range(ref_period)
################################## #
## ABC PARAMETERS  ----
################################## #
model_param_update
# set priors
stride_prior <- list( num_daily_imported_cases     = c("unif",0,50),
                      
                      temporal_distancing_workplace_opt1 = c("unif",0.5,0.80),  # 0.8360346
                      temporal_distancing_workplace_opt2 = c("unif",0.5,0.80),
                      temporal_distancing_workplace_opt3 = c("unif",0.5,0.80),
                      
                      temporal_distancing_community_opt1 = c("unif",0.7,0.9),  # 0.8967415
                      temporal_distancing_community_opt2 = c("unif",0.7,0.9),
                      temporal_distancing_community_opt3 = c("unif",0.7,0.9),
                      
                      temporal_distancing_collectivity_opt1 = c("unif",0.5,0.9),  # 0.6798729
                      temporal_distancing_collectivity_opt2 = c("unif",0.5,0.9),
                      temporal_distancing_collectivity_opt3 = c("unif",0.5,0.9),
                      
                      cnt_reduction_school_exit             = c("unif",0.0,0.50),
                      cnt_reduction_school_exit_secondary   = c("unif",0.0,0.50),
                      cnt_reduction_school_exit_tertiary    = c("unif",0.0,0.50)
                     )  


temporal_parameters_timepoints <- c_str('2020-05-01','2020-08-15','2020-10-01','2020-11-01')
# exp_param_list$temporal_distancing_workplace    <- c_str(seq(0.85,0.75,length=8))
model_param_update$dates_distancing_workplace       <- temporal_parameters_timepoints

# exp_param_list$temporal_distancing_community      <- c_str(0.70,0.4,0.8)
model_param_update$dates_distancing_community         <- temporal_parameters_timepoints

# exp_param_list$temporal_distancing_collectivity <- c_str(seq(0.70,0.60,length=8))
model_param_update$dates_distancing_collectivity    <- temporal_parameters_timepoints

#exp_param_list$num_infected_seeds                  <- 50
model_param_update$temporal_imported_cases          <- c_str(0,1,1,0)
model_param_update$dates_imported_cases             <- c_str('2020-06-01','2020-07-31','2020-08-01','2020-08-31','2020-09-01')

# exp_param_list$cnt_reduction_school_exit           <- 0.8
# exp_param_list$cnt_reduction_school_exit_secondary <- 0.1
# exp_param_list$cnt_reduction_school_exit_tertiary  <- 0.6


length(stride_prior)

################################################ #
## REFERENCE DATA  ----
################################################ #

sum_stat_obs <- get_abc_reference_data(ref_period               = ref_period,
                                       bool_age                 = TRUE,
                                       bool_doubling_time       = FALSE,
                                       bool_hospital            = TRUE,
                                       bool_serology            = FALSE,
                                       rel_importance_hosp_data = rel_importance_hosp_data,
                                       age_cat_hosp_str         = model_param_update$hospital_category_age,
                                       bool_add_pop_stat        = TRUE,
                                       bool_truncate_serology   = FALSE)
table(sum_stat_obs$category)
table(is.na(sum_stat_obs$value))

################################################ #
## RUN ABC   ----
################################################ #

# create output folder and set workdir
run_file_path <- dirname(smd_file_path(project_dir,'test'))
setwd(run_file_path)
saveRDS(model_param_update,'model_param_update.rds')
saveRDS(sum_stat_obs,'sum_stat_obs.rds')
saveRDS(stride_prior,'stride_prior.rds')

# run_param  <- sample_param_from_prior(stride_prior)
# stride_out <- run_rStride_abc(run_param,remove_run_output = F)
# length(stride_out)
# dim(sum_stat_obs)

# ABC_stride<-ABC_rejection(model     = run_rStride_abc,
#                            prior    = stride_prior,
#                            nb_simul = n_sample,
#                            summary_stat_target=sum_stat_obs$value,
#                            tol=pacc,
#                            verbose = F,
#                            n_cluster=n_cluster,
#                            use_seed=TRUE,
#                            progress_bar=T)
# 
ABC_stride<-ABC_sequential(model=run_rStride_abc,
                           prior=stride_prior,
                           nb_simul=n_sample,
                           summary_stat_target=sum_stat_obs$value,
                           method = "Lenormand",
                           p_acc_min=pacc,
                           verbose = T,
                           n_cluster=n_cluster,
                           use_seed=TRUE,
                           progress_bar=T)
# 
# set back workdir
setwd('../..')


# par(mfrow=c(3,2))
saveRDS(ABC_stride,file=smd_file_path(project_dir,'ABC_stride.rds'))
save(list=ls(),file=smd_file_path(project_dir,'ABC_stride_all.RData'))

############################# #
## EXPLORE RESULTS         ####
############################# #

# # load results
# ABC_stride <- readRDS(smd_file_path(project_dir,'ABC_stride.rds'))
# load(file=smd_file_path(project_dir,'ABC_stride_all.RData'))

# or use intermediate output
#ABC_stride <- load_partial_results_abc(project_dir)

# # re-load rStride
# source('./bin/rstride/rStride.R')

print(ABC_stride$computime/3600)
print(ABC_stride$nsim)
print(length(ABC_stride$intermediary))

# plot (final) results
plot_abc_results(ABC_stride,project_dir)

# plot posterior distribution per iteration
plot_abc_posterior(ABC_stride,project_dir)

# plot parameter correlation
plot_abc_correlation(ABC_stride,project_dir)

# plot single best 
plot_abc_singleton(project_dir)

# # intermediate results
# plot_abc_intermediate(ABC_stride,project_dir)


