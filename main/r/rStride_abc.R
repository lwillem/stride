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
#  Copyright 2020, Willem L
############################################################################ #
#
# Call this script from the main project folder (containing bin, config, lib, ...)
# to get all relative data links right. 
#
# E.g.: path/to/stride $ ./bin/rStride_abc.R 
#
############################################################################ #

args = commandArgs(trailingOnly=TRUE)
if(length(args)>0){
  job_id <- paste0(format(Sys.time(), format="%Y%m%d_"),args[1])
} else{
  job_id <-''
}

# lw's developing option
if(any(grepl('lwillem',dir('~',full.names=T)))){
  .rstride$set_wd()
}

# Clear work environment
rm(list=ls()[ls() != 'job_id'])

# Load rStride
source('./bin/rstride/rStride.R')

# Load default parameter configurations
source('./bin/rStride_intervention_baseline.R')

# load ABC package
library(EasyABC)

# set samples and cluster size
n_cluster = 8
n_sample = n_cluster*2

# set acceptance level
pacc=0.5

# set directory postfix (optional)
dir_postfix <- paste0('_abc_n',n_sample,'_c',n_cluster,'_p',formatC(pacc*100,flag = 0,digits = 2))
dir_postfix

# create run tag using the current time if use_date_prefix == TRUE
# use_date_prefix <- TRUE
use_date_prefix <- job_id ==''
run_tag <- ifelse(use_date_prefix,format(Sys.time(), format="%Y%m%d_%H%M%S"),job_id)
run_tag_data <- paste0(run_tag,dir_postfix)

# set project directory
project_dir <- smd_file_path('./sim_output',run_tag_data)

################################################ #
## DESIGN OF EXPERIMENTS  ----
################################################ #

# add default parameters and values to combine in a full-factorial grid
model_param_update <- get_exp_param_default(bool_revised_model_param = T,
                                            bool_min_restrictive = T)

# TEMP
model_param_update$population_file <- "pop_belgium600k_c500_teachers_censushh.csv"
model_param_update$num_days        <- 74
#model_param_update$logparsing_cases_upperlimit <- 2.5e6


ref_period <- seq(as.Date('2020-03-15'),
                  as.Date(model_param_update$start_date) + model_param_update$num_days-1,
                  1)

# ################################################ #
# ## REFERENCE DATA  ----
# ################################################ #


sum_stat_obs <- get_abc_reference_data(ref_period,
                                       bool_age = FALSE,
                                       bool_doubling_time = TRUE)


################################## #
## RUN ABC  ----
################################## #

# set priors
stride_prior <- list(r0                         = c("unif",1.0,5.0),   
                     num_infected_seeds         = c("unif",200,600),
                     hosp_probability_factor    = c("unif",0.05,0.95),
                     cnt_reduction_workplace    = c("unif",0.60,0.95),
                     compliance_delay_workplace = c("unif",4.51,7.49),  # rounded: 5-7
                     cnt_reduction_other        = c("unif",0.60,0.95),
                     compliance_delay_other     = c("unif",4.51,7.49))  # rounded: 5-7



## for debugging
#.rstride$set_wd()                    

# create output folder and set workdir
run_file_path <- dirname(smd_file_path('./sim_output',run_tag_data,'test'))
setwd(run_file_path)
saveRDS(model_param_update,'model_param_update.rds')
saveRDS(sum_stat_obs,'sum_stat_obs.rds')
saveRDS(stride_prior,'stride_prior.rds')

# run_param <- c(rng_seed=20,4,400,0.4,0.85,7.4,0.85,4.51)
# names(run_param)[-1] <- names(stride_prior)
# stride_out <- run_rStride_abc(run_param)
# length(stride_out)
# dim(sum_stat_obs)

# ABC_stride<-ABC_rejection(model     = run_rStride_abc,
#                            prior    = stride_prior,
#                            nb_simul = n_sample,
#                            summary_stat_target=sum_stat_obs$value,
#                            tol=pacc,
#                            verbose = T,
#                            n_cluster=n_cluster,
#                            use_seed=TRUE,
#                            progress_bar=T)


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

# # intermediate results
# plot_abc_intermediate(ABC_stride,project_dir)



