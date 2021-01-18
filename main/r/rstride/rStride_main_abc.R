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
# ABC controller for the Stride model
#
############################################################################ #

#' Main rStride function for ABC
# abc_function_param <- c(15,3.4,256,0.4,0.85,7.4,0.85,4.51)
#abc_function_param <- c(41,4,400,0.4,0.85,7.4,0.85,4.51)
# abc_function_param <- run_param; remove_run_output <- FALSE

################################################ #
## RUN  ----
################################################ #

run_rStride_abc <- function(abc_function_param,
                            remove_run_output     = TRUE)
{
  # define rng seed
  rng_seed <- abc_function_param[1]

  # add process-id-specific delay 
  # note: to prevent multiple project-dirs when starting many parallel processes
  # note: "+1" to prevent "-Inf"
  Sys.sleep(unlist(log(rng_seed+1,10)))

  wd_start <- getwd()
  run_tag <- basename(getwd())
  setwd('../..')

  # load functions within parallel worker  
  source('./bin/rstride/rStride.R')

  ################################## #
  ## GENERAL OPTIONS              ####
  ################################## #
  stride_bin              <- './bin/stride'
  config_opt              <- '-c'
  config_default_filename <- './config/run_default.xml'
  output_dir              <- 'sim_output'
  
  ################################## #
  ## OUTPUT DIRECTORY             ####
  ################################## #

  # create project directory
  project_dir <- smd_file_path(output_dir,run_tag)

  
  ################################## #
  ## CONFIGURATION                 ####
  ################################## #
  
  # get default config
  config_exp <- create_default_config(config_default_filename, run_tag)
  
  # set event_log_level to "Incidence"
  config_exp$event_log_level <- 'Incidence'
  
  # incorportate experiment-specific parameter values
  model_param_update <- readRDS(file.path('./sim_output',run_tag,'model_param_update.rds'))
  config_exp[names(model_param_update)] <- model_param_update
 
  # check/adjust parameter names
  if(is.null(names(abc_function_param))){
     names(abc_function_param) <- c('rng_seed',names(readRDS(file.path('./sim_output',run_tag,'stride_prior.rds'))))
  }

  # aggregate age-specific parameters
  abc_param_aggr <- collapse_age_param(abc_function_param)

  # copy parameter values
   for(i_param in names(abc_param_aggr)){
      config_exp[i_param]  <- abc_param_aggr[i_param]
   }
  
  # make sure some input parameters are coded as integer value
  config_exp$num_infected_seeds         <- round(config_exp$num_infected_seeds )
  config_exp$compliance_delay_workplace <- round(config_exp$compliance_delay_workplace)
  config_exp$compliance_delay_other     <- round(config_exp$compliance_delay_other)
  

  ################################## #
  ## RUN                          ####
  ################################## #
  
  # create experiment tag
  i_exp   <- rng_seed
  exp_tag <- .rstride$create_exp_tag(i_exp)
  
   #save the config as XML file
   output_prefix       = smd_file_path(project_dir,exp_tag,.verbose=FALSE)
   config_exp$output_prefix <- output_prefix 
   
   config_exp_filename = paste0(output_prefix,".xml")
   save_config_xml(config_exp, config_exp_filename)

   # run stride (using the C++ Controller)
   cmd = paste(stride_bin,config_opt, paste0("../", config_exp_filename))
   system(cmd,ignore.stdout = TRUE)

   # load output summary
   summary_filename <- file.path(output_prefix,'summary.csv')
   run_summary      <- read.table(summary_filename,header=T,sep=',')
   
   # merge output summary with input param
   # note: do not use "merge" to prevent issues with decimal numbers
   config_df   <- as.data.frame(config_exp)
   config_df   <- config_df[,!names(config_df) %in% names(run_summary)]
   run_summary <- cbind(run_summary,config_df)
   
   # save summary
   write.table(run_summary,file=file.path(project_dir=output_prefix,
                                          paste0(run_tag,'_summary.csv')),sep=',',row.names=F)
   
   # parse event_log and process model output
   parse_log_file(config_exp, 
                  i_exp, 
                  # get_burden_rdata, 
                  get_transmission_rdata = FALSE, 
                  get_tracing_rdata = FALSE, 
                  project_dir_exp = config_exp$output_prefix)
   
   # get transmission output
   parsed_logfile <- dir(output_prefix,pattern = 'rds',full.names = T)
   data_incidence_all <- readRDS(parsed_logfile)$data_incidence
   
   names(data_incidence_all)
  
   # get ref data (temp)
   sim_stat_filename <- file.path('./sim_output',run_tag,'sum_stat_obs.rds')
   if(file.exists(file.path('./sim_output',run_tag,'sum_stat_obs.rds'))){
      sum_stat_obs <- readRDS(file.path('./sim_output',run_tag,'sum_stat_obs.rds'))
   } else{
      ref_period   <- unique(data_incidence_all$sim_date)
      sum_stat_obs <- get_abc_reference_data(ref_period)
   }
   
   
   # if doubling time is part of the reference output ==> add summary statistic for given period
   if(any(grepl('doubling_time',sum_stat_obs$category))){
      rstride_category <- as.character(unique(sum_stat_obs$category[grepl('doubling_time',sum_stat_obs$category)]))
      if(length(rstride_category)>1){
         smd_print('rSTRIDE CANNOT HANDLE MORE THAN ONE REFERENCE DOUBLING TIME (yet) => RETURN NA!',WARNING = T)
         data_incidence_all[,eval(rstride_category):=NA]
         #data_incidence_all[,get(rstride_category)]
      } else{
         
         dtime_period <- unique(sum_stat_obs$date[grepl('doubling_time',sum_stat_obs$category)])
         flag_exp_doubling   <- data_incidence_all$sim_date %in% dtime_period
         if(sum(flag_exp_doubling)>0){
            doubling_time_model     <- get_doubling_time(data_incidence_all$new_infections[flag_exp_doubling])
         } else{
            doubling_time_model <- NA
         }
         data_incidence_all[,eval(rstride_category):=doubling_time_model]
      }
      names(data_incidence_all)
   }
   
   # create model-based summary stats
   abc_out <- vector(length=nrow(sum_stat_obs))
   
   i_ref <- 6
   sum_stat_obs[i_ref,]
   for(i_ref in 1:nrow(sum_stat_obs)){
      
      # category and age ==>> column
      flag_col <- names(data_incidence_all) == sum_stat_obs$category[i_ref]
      table(flag_col)
      name_col <- names(data_incidence_all)[flag_col]
      
      # date => row
      data_incidence_all[,get(name_col)]
      data_incidence_all[sim_date == sum_stat_obs$date[i_ref],]
      
      # save value
      abc_out[i_ref] <- data_incidence_all[sim_date == sum_stat_obs$date[i_ref],get(name_col)]
   }
   head(abc_out)

   # tmp debug code
   #write.table((abc_out),paste0(output_prefix,'.csv'),sep=',',row.names=F)
   
   # remove experiment output and config
   if(remove_run_output){
     unlink(config_exp$output_prefix,recursive=TRUE)
     unlink(config_exp_filename,recursive = TRUE)
     #unlink(paste0(output_prefix,'.csv'),recursive = T) # tmp debug code
   }
  
   # reset wd
   setwd(wd_start)

  return(abc_out)
  
} # end run_rStride_abc function

################################################ #
## EXPLORE RESULTS  ----
################################################ #

# function to plot the ABC results: parameters & summary statistics (over time)
plot_abc_results <- function(ABC_out,project_dir,bool_pdf=TRUE){
   
   # model parameters
   stride_prior       <- readRDS(dir(project_dir,pattern = 'stride_prior',full.names = T))
   sum_stat_obs       <- readRDS(dir(project_dir,pattern = 'sum_stat_obs',full.names = T))

   # open pdf stream
   if(bool_pdf) .rstride$create_pdf(project_dir = project_dir,file_name = 'results_ABC')
   
   ## model parameters ----
   par(mfrow=c(4,3))
   # parameters
   for(i in 1:ncol(ABC_out$param)){
      
      hist(ABC_out$param[,i],20,
           xlim = as.numeric(stride_prior[[i]][-1]),
           xlab = names(stride_prior)[i],
           main = nameLabels(stride_prior)[i])
      legend('topright',
             title='mean',
             paste(round(mean(ABC_out$param[,i]),digits=2)),
             cex=0.5)
   }
   
   
   for(i in grep('compliance_delay',names(stride_prior))){
      hist(round(ABC_out$param[,i]),
           xlab = names(stride_prior)[i],
           main = paste(nameLabels(stride_prior)[i],'\n[DISCRETE]'))
   }
   
   
   ## output statistics ----
   par(mfrow=c(4,3))
   
   ## if  over time
   
   # get categories
   output_cat <- as.character(unique(sum_stat_obs$category))
   
   # iterate over categories
   for(i_cat in output_cat){
      flag_out <- sum_stat_obs$category == i_cat & sum_stat_obs$bool_orig == TRUE
      sum(flag_out)
      dim(sum_stat_obs)
      dim(ABC_out$stats)
      
      # if over time
      if(sum(flag_out)>1){
         y_lim <- range(pretty(c(sum_stat_obs$value[flag_out],
                                 sum_stat_obs$value_low[flag_out],
                                 sum_stat_obs$value_high[flag_out],
                                 ABC_out$stats[,flag_out])),
                        na.rm=T)
         
         plot(sum_stat_obs$date[flag_out],
              sum_stat_obs$value[flag_out],
              ylim = y_lim,
              main = string2label(i_cat),
              ylab = i_cat)
         
         if(all(!is.na(sum_stat_obs$value_low[flag_out]))){
            add_interval(x = sum_stat_obs$date[flag_out],
                         y1 = sum_stat_obs$value_low[flag_out],
                         y2 = sum_stat_obs$value_high[flag_out])
         }
         
         for(i_out in 1:nrow(ABC_out$stats)){
            lines(sum_stat_obs$date[flag_out],
                  ABC_out$stats[i_out,flag_out],
                  col=alpha(4,0.8))
         }
      } else{ # else: hist + reference
         
         # initial doubling time
         x_lim <- range(0,pretty(sum_stat_obs$value_high[flag_out]*1.1))
         hist(ABC_out$stats[,flag_out],40,
              xlim=x_lim,
              xlab=i_cat,
              main=i_cat,
              col=alpha(4,0.2),
              border = alpha(4,0.2))
      #   abline(v=mean(sum_stat_obs$value[flag_out]),col=1)
         add_interval_hor(x1  = sum_stat_obs$value_low[flag_out],
                          x2  = sum_stat_obs$value_high[flag_out],
                          x_mean = sum_stat_obs$value[flag_out],
                          col = 1)
      }
      
      
      
   }
   
   # close pdf stream
   if(bool_pdf) dev.off()
}


plot_abc_intermediate <- function(ABC_out,project_dir,bool_pdf=TRUE){
   
   # if intermediate results present ==>> plot
   if( 'intermediary' %in% names(ABC_out)){
   
      stride_prior       <- readRDS(dir(project_dir,pattern = 'stride_prior',full.names = T))

      # open pdf stream
      if(bool_pdf){ .rstride$create_pdf(project_dir = project_dir,file_name = 'results_ABC_intermediate') }
      
      for(i_seq in 1:length(ABC_out$intermediary)){
         ABC_out_temp <- ABC_out
         num_param <- length(stride_prior)
         ABC_out_temp$param <- ABC_out$intermediary[[i_seq]]$posterior[,2:(num_param+1)]
         ABC_out_temp$stats <- ABC_out$intermediary[[i_seq]]$posterior[,-(0:(num_param+1))]
         plot_abc_results(ABC_out_temp,project_dir,bool_pdf = FALSE)
      }
      
      if(bool_pdf) dev.off()
   }
}

plot_abc_posterior <- function(ABC_out,project_dir,bool_pdf=TRUE){
   
   # if intermediate results present ==>> plot
   if( 'intermediary' %in% names(ABC_out)){
      
      stride_prior       <- readRDS(dir(project_dir,pattern = 'stride_prior',full.names = T))

      ## open pdf steam
      if(bool_pdf) .rstride$create_pdf(project_dir = project_dir,file_name = 'results_ABC_posterior')
      
      get_stat <- function(x){
         c(quantile(x,0.025),
           quantile(x,0.25),
           mean(x),
           quantile(x,0.75),
           quantile(x,0.975))
      }
      
      param_num   <- length(stride_prior) 
      param_names <- names(stride_prior)
      
      foreach(i_iter = 1:length(ABC_out$intermediary),
              .combine = 'rbind') %do% {
                 
                 foreach(i_param = 1:length(param_names),
                         .combine = 'cbind') %do% {
                            get_stat(ABC_out$intermediary[[i_iter]]$posterior[,i_param+1])
                  } -> p_out_matrix
                 colnames(p_out_matrix) <- param_names
                 p_out <- as.double(p_out_matrix)
                 names(p_out) <- paste0(rep(param_names,each=nrow(p_out_matrix)),'_',rep(1:nrow(p_out_matrix),param_num))
                 
                 c(iter = i_iter,
                   n_simul_tot = ABC_out$intermediary[[i_iter]]$n_simul_tot,
                   tol_step = ABC_out$intermediary[[i_iter]]$tol_step,
                  p_out)
              } -> db_abc
      db_abc <- data.frame(db_abc)
      
      if(ncol(db_abc)==1){
         db_abc <- data.frame(t(db_abc)) 
      }
      
      par(mfrow=c(3,3))
      plot(db_abc$iter,(db_abc$n_simul_tot),main='num simulations')
      plot(db_abc$iter,log(db_abc$tol_step),main='log(tolerance)')
      
      
      for(i in order(names(stride_prior))){
         tmp_out <- db_abc[,paste0(param_names[i],'_',1:5)]
         
         plot(db_abc$iter,tmp_out[,3],type='b',
              ylim=range(as.numeric(stride_prior[[i]][-1])),
              main=nameLabels(stride_prior)[i],
              ylab=names(stride_prior)[i])
         lines(db_abc$iter,tmp_out[,1],lty=2)
         lines(db_abc$iter,tmp_out[,5],lty=2)
         
         polygon(x=c(db_abc$iter,rev(db_abc$iter)),
                 y=c(tmp_out[,2],rev(tmp_out[,4])),
                 col=alpha(1,0.1),
                 border = NA)
         
         text(db_abc$iter[2],
              tmp_out[2,5],
              '95%',
              pos=1,
              cex=0.5)
         
         text(db_abc$iter[2],
              tmp_out[2,4],
              '75%',
              pos=1,
              cex=0.5)
      }
      
      # close pdf stream
      if(bool_pdf) dev.off()
   }
}

plot_abc_correlation <- function(ABC_out,project_dir){
   
   stride_prior       <- readRDS(dir(project_dir,pattern = 'stride_prior',full.names = T))
   
   .rstride$create_pdf(project_dir = project_dir,file_name = 'results_ABC_correlation')
   posterior_param <- ABC_out$param
   colnames(posterior_param) <- nameLabels(stride_prior)
   corrplot(cor(posterior_param))
   dev.off()
}


################################################ #
## REFERENCE DATA  ----
################################################ #

get_abc_reference_data <- function(ref_period ,
                                   bool_age  = FALSE,
                                   bool_doubling_time = TRUE,
                                   bool_hospital = TRUE,
                                   rel_importance_hosp_data = NA){
   
   # set contribution hospital data vs other data sources
   # if 1: number of hospital admission data points == number of (e.g.) seroprevalence data points
   # if 2: importance of hospital addmission data increases...
   rel_factor_hosp_data <- ifelse(!is.na(rel_importance_hosp_data),1/rel_importance_hosp_data,1)
   
   ## hospital reference data ----
   # use (local version of) most recent SCIENSANO data (or backup version)
   hosp_ref_data          <- get_hospital_incidence_age(config_exp$hospital_category_age)
   hosp_ref_data$sim_date <- as.Date(hosp_ref_data$sim_date)
   dim(hosp_ref_data)
   hosp_ref_data          <- hosp_ref_data[hosp_ref_data$sim_date %in% ref_period,]
   
   if(bool_hospital){
      if(bool_age){
         abc_age_cat       <- seq(0,80,10) #TODO: make generic
         abc_hosp_stat     <- data.table(value      = unlist(hosp_ref_data[,grepl('hospital_admissions_',names(hosp_ref_data))]),
                                         value_low  = NA,
                                         value_high = NA,
                                         date       = rep(hosp_ref_data$sim_date,length(abc_age_cat)),
                                         age_min    = rep(abc_age_cat,each=nrow(hosp_ref_data)),
                                         category   = rep(paste0('new_hospital_admissions_age',1:length(abc_age_cat)),each=length(hosp_ref_data$sim_date)),
                                         bool_orig  = TRUE)
      } else{
         abc_hosp_stat     <- data.table(value      = hosp_ref_data$hospital_admissions,
                                         value_low  = NA,
                                         value_high = NA,
                                         date       = hosp_ref_data$sim_date,
                                         age_min    = NA,
                                         category   = 'new_hospital_admissions',
                                         bool_orig  = TRUE)
      }  
   } else{
      abc_hosp_stat <- NULL
   }
    
   dim(abc_hosp_stat)
   
   ## seroprevalence data ----
   prevalence_ref <- load_observed_seroprevalence_data(ref_period = ref_period,
                                                       analysis = ifelse(bool_age,'age','overall'))
   # temporary fix for 80-90 year olds
   prevalence_ref <- prevalence_ref[prevalence_ref$age_min!=90,]
   
   # set "level" as factor
   prevalence_ref$level <- as.factor(prevalence_ref$level)
   
   abc_sero_stat     <- data.table(value      = prevalence_ref$point_incidence_mean,
                                   value_low  = prevalence_ref$point_incidence_low,
                                   value_high = prevalence_ref$point_incidence_high,
                                   date       = prevalence_ref$seroprevalence_date,
                                   age_min    = prevalence_ref$age_min,
                                   category   = 'cumulative_infections',
                                   bool_orig  = TRUE)
   if(bool_age){
      abc_sero_stat[,category := paste0('cumulative_infections_age', as.numeric(prevalence_ref$level))]
   }
   
   ## doubling time 3.1 (2.4-4.4) \cite{pellis2020challenges} ----
   ref_doubling_time <- data.frame(dates    = seq(as.Date('2020-02-24'),as.Date('2020-03-08'),1),
                                   mean     = 3.1,
                                   CI_low   = 2.4,
                                   CI_upper = 4.4
   )
   ref_doubling_time <- ref_doubling_time[ref_doubling_time$dates %in% ref_period,]
   
   if(bool_doubling_time && nrow(ref_doubling_time) > 0){
      abc_dtime_stat    <- data.table(value      = ref_doubling_time$mean,
                                      value_low  = ref_doubling_time$CI_low,
                                      value_high = ref_doubling_time$CI_upper,
                                      date       = ref_doubling_time$dates,
                                      age_min    = NA,
                                      category   = 'doubling_time_march',
                                      bool_orig  = c(TRUE,rep(FALSE,length(ref_doubling_time$mean)-1)))
   } else{
      abc_dtime_stat <- NULL
   }
      
   # combine ----
   sum_stat_obs      <- rbind(abc_hosp_stat,
                              abc_sero_stat,
                              abc_dtime_stat)
   dim(sum_stat_obs)

   # duplicate serology data to give it the same weight in the reference data
   if(!is.null(abc_hosp_stat)){
      sero_rep_factor   <- floor(nrow(abc_hosp_stat) / nrow(abc_sero_stat) * rel_factor_hosp_data)
      for(i in 2:sero_rep_factor){
         sum_stat_obs <- rbind(sum_stat_obs,abc_sero_stat[, bool_orig := FALSE])
      }
      dim(sum_stat_obs); table(sum_stat_obs$bool_orig)
   }
   
   
   if(!is.null(abc_hosp_stat) && !is.null(abc_dtime_stat)){
      sero_dtime_factor   <- floor(nrow(abc_hosp_stat) / nrow(abc_dtime_stat) * rel_factor_hosp_data )
      for(i in 2:sero_dtime_factor){
         sum_stat_obs <- rbind(sum_stat_obs,abc_dtime_stat[,bool_orig := FALSE])
      }
      dim(sum_stat_obs); table(sum_stat_obs$bool_orig)
   }
   
   return(sum_stat_obs)
   
   # nrow(abc_hosp_stat) / nrow(abc_sero_stat)
   # 
   # sero_rep_factor   <- floor(nrow(abc_hosp_stat) / nrow(abc_sero_stat))
   # sum_stat_obs      <- rbind(abc_hosp_stat,abc_sero_stat)
   # sum_stat_obs['bool_orig'] <- TRUE
   # for(i in 2:sero_rep_factor){
   #    sum_stat_obs <- rbind(sum_stat_obs,
   #                          cbind(abc_sero_stat,
   #                                bool_orig = FALSE))
   # }
   
}


# get a sample from the given ABC parameter prior list
# development function to run the ABC procedure
sample_param_from_prior <- function(stride_prior){
   
   param_names_all <- names(stride_prior)
   
   # start with general parameters (not age-specific)
   param_value_out <- data.frame(matrix(0,ncol=length(param_names_all)+1))
   names(param_value_out) <- c('rng_seed',param_names_all)
   
   param_value_out$rng_seed <- 100
   
   # add mean for other parameters
   for(i_param in param_names_all){
      x_range <- as.numeric(stride_prior[[i_param]][-1])
      x_min <- x_range[1]
      x_max <- x_range[2]
      param_value_out[i_param] <- round(runif(1,x_min,x_max),digits=2)
   }
   
   # return result
   return(param_value_out)
}


# development function to run the ABC procedure
# param_list <- abc_function_param
collapse_age_param <- function(param_list){

   param_names_all <- names(param_list)

   # start with general parameters (not age-specific)
   param_general  <- param_names_all[!grepl('_opt',param_names_all)]
   param_list_out <- param_list[param_general]
   param_list_out <- as.list(param_list_out)
   cat('****',file = 'function_out.txt',append = T,fill = T)
   cat(typeof(param_list),file = 'function_out.txt',append = T,fill = T)
   cat(typeof(param_list_out),file = 'function_out.txt',append = T,fill = T)
   cat(length(param_list_out),file = 'function_out.txt',append = T,fill = T)
   cat('****',file = 'function_out.txt',append = T,fill = T)
   
   
   # age specific parameters
   param_age <- param_names_all[!param_names_all %in% param_general]
   param_age_main  <- gsub('_opt.','',param_age)
   param_age_level <- gsub('.*_opt','',param_age)

   for(i_param in unique(param_age_main)){
      age_values <- NULL
      for(i_param_age in param_age[param_age_main == i_param]){
         age_values <- c(age_values,param_list[i_param_age])
      }
      param_list_out[i_param] <- paste(age_values,collapse=',')
   }

   return(param_list_out)
}

################################################ #
## LOAD PARTIAL RESULTS  ----
################################################ #

load_partial_results_abc <- function(project_dir){
   
   # model parameters
   stride_prior       <- readRDS(dir(project_dir,pattern = 'stride_prior',full.names = T))
   sum_stat_obs       <- readRDS(dir(project_dir,pattern = 'sum_stat_obs',full.names = T))
   model_param_update <- readRDS(dir(project_dir,pattern = 'model_param_update',full.names = T))

   length(stride_prior)
   dim(sum_stat_obs)
   length(model_param_update)
   
   # set column names
   col_names    <- c('step_id','tab_weight',names(stride_prior),sum_stat_obs$category)
   param_names  <- names(stride_prior)
   output_names <- sum_stat_obs$category
   
   # model step
   model_step_files <- dir(project_dir,pattern = 'model_step',full.names = T)
   foreach(i_file = 1:length(model_step_files),
           .combine = 'rbind') %do%{
              step_id <- unlist(strsplit(model_step_files[i_file],'model_step'))[2]
              data.frame(step_id,read.table(model_step_files[i_file],sep=' '))
   } -> model_step
   dim(model_step)  
   names(model_step) <- col_names
   
   # output_step
   output_step_files <- dir(project_dir,pattern = 'output_step',full.names = T)
   foreach(i_file = 1:length(output_step_files),
           .combine = 'rbind') %do%{
              step_id <- unlist(strsplit(output_step_files[i_file],'output_step'))[2]
              data.frame(step_id,read.table(output_step_files[i_file],sep=' '))
   } -> output_step
   dim(output_step) 
   names(output_step) <- col_names
   
   # n_simul_tot_step
   n_simul_tot_step_files <- dir(project_dir,pattern = 'n_simul_tot_step',full.names = T)
   foreach(i_file = 1:length(n_simul_tot_step_files),
           .combine = 'rbind') %do%{
              step_id <- unlist(strsplit(n_simul_tot_step_files[i_file],'n_simul_tot_step'))[2]
              data.frame(step_id,read.table(n_simul_tot_step_files[i_file],sep=' '))
           } -> n_simul_tot_step
   dim(n_simul_tot_step) 
   names(n_simul_tot_step) <- c('step_id','n_simul_tot_step')
   
   # tolerance_step
   tolerance_step_files <- dir(project_dir,pattern = 'tolerance_step',full.names = T)
   foreach(i_file = 1:length(n_simul_tot_step_files),
           .combine = 'rbind') %do%{
              step_id <- unlist(strsplit(tolerance_step_files[i_file],'tolerance_step'))[2]
              data.frame(step_id,read.table(tolerance_step_files[i_file],sep=' '))
           } -> tolerance_step
   dim(tolerance_step) 
   names(tolerance_step) <- c('step_id','tolerance_step')
   
   # p_acc
   p_acc_step_files <- dir(project_dir,pattern = 'p_acc_step',full.names = T)
   if(length(p_acc_step_files)>0){
      foreach(i_file = 1:length(p_acc_step_files),
              .combine = 'rbind') %do%{
                 step_id <- unlist(strsplit(p_acc_step_files[i_file],'p_acc_step'))[2]
                 data.frame(step_id,read.table(p_acc_step_files[i_file],sep=' '))
              } -> p_acc_step
      dim(p_acc_step) 
      names(p_acc_step) <- c('step_id','p_acc_step')
   } else{
      p_acc_step <- data.frame(step_id=1,
                               p_acc_step = NA)
   }
   
   # compute time?
   all_file_info <- file.info(dir(project_dir,full.names = T))
   computime <- round(as.numeric(difftime(max(all_file_info$mtime),min(all_file_info$mtime),units='secs')))
   
   
   # reformat: ABC_out$param & ABC_out$stats
   model_step      <- list2double(model_step)
   flag_final_step <- model_step[,'step_id'] == max(model_step[,'step_id'])
   flag_output_col <- colnames(model_step) %in% sum_stat_obs$category

   ABC_out    <- list(param = model_step[flag_final_step,param_names],
                      stats = model_step[flag_final_step,flag_output_col],
                      computime = computime,
                      nsim = sum(n_simul_tot_step$n_simul_tot_step))
   #plot_abc_results(ABC_out,project_dir,bool_pdf = F)
   
   # reformat: ABC_out$intermediate
   ABC_out$intermediary <- list(1:nrow(tolerance_step))
   output_step$step_id <- as.numeric(output_step$step_id)
   i_step <- 1
   for(i_step in 1:max(output_step$step_id)){
      posterior <- list2double(output_step[output_step$step_id == i_step,-1])
      ABC_out$intermediary[[i_step]] <- list(n_simul_tot = as.numeric(n_simul_tot_step$n_simul_tot_step[n_simul_tot_step$step_id == i_step]),
                                             tol_step = as.numeric(tolerance_step$tolerance_step[tolerance_step$step_id == i_step]),
                                             p_acc_step = as.numeric(p_acc_step$p_acc_step[p_acc_step$step_id == i_step]),
                                             posterior=posterior)
   }
   #plot_abc_intermediate(ABC_out,project_dir,bool_pdf = F)
   
   # return result
   return(ABC_out)
   
}


list2double <- function(x_list){
   x_matrix <- as.matrix(x_list)
   x_double <- matrix(as.double(as.matrix(x_matrix)),nrow=nrow(x_matrix),ncol=ncol(x_matrix))
   colnames(x_double) <- names(x_list)
   return(x_double)
}

#prior_names <- names(stride_prior)
nameLabels <- function(data_list){

   name_list <- names(data_list)
   return(string2label(name_list))
}

string2label <- function(name_list){
   
   name_list <- gsub('cumulative','cum',name_list)
   name_list <- gsub('hospital','hosp',name_list)
   name_list <- gsub('admissions','adm',name_list)
   name_list <- gsub('probability','prob',name_list)
   name_list <- gsub('reduction','reduct',name_list)
   name_list <- gsub('susceptibility','susc',name_list)
   name_list <- gsub('disease','dis',name_list)
   name_list <- gsub('infections','infect',name_list)
   
   return(name_list)
}


process_partial_results <- function(output_folders = NA){
   
   if(any(is.na(output_folders))){
      output_folders <- dir('sim_output/tmp',pattern = 'abc',full.names = T)
      output_folders <- output_folders[!grepl('.csv',output_folders)]
   } else {
      
   }
   
   for(project_dir in output_folders){
      
      print(project_dir)
      
      # load partial results
      ABC_stride <- load_partial_results_abc(project_dir)
      
      # plot (final) results
      plot_abc_results(ABC_stride,project_dir)
      
      # plot posterior distribution per iteration
      plot_abc_posterior(ABC_stride,project_dir)
      
      # plot parameter correlation
      plot_abc_correlation(ABC_stride,project_dir)
   }
   
}


# project_dir <- 'sim_output/tmp/20210117_11068_abc_age1_n120_c40_p010/'
get_final_best <- function(project_dir){
   
   # model output
   ABC_out <- load_partial_results_abc(project_dir)
   
   # model parameters and reference
   stride_prior       <- readRDS(dir(project_dir,pattern = 'stride_prior',full.names = T))
   sum_stat_obs       <- readRDS(dir(project_dir,pattern = 'sum_stat_obs',full.names = T))
   
   # names(ABC_out)
   # dim(ABC_out$param)
   # dim(ABC_out$stats)
   # colnames(ABC_out$intermediary[[1]]$posterior)
   # 
   
   dim(ABC_out$stats)
   dim(sum_stat_obs)
   get_binom_poisson(sum_stat_obs$value,ABC_out$stats[1,])
   
   foreach(i_run = 1:nrow(ABC_out$stats),
           .combine='rbind') %do% {
              get_binom_poisson(sum_stat_obs$value,ABC_out$stats[i_run,])
           } -> ABC_binom_poisson
   
   hist(ABC_binom_poisson)
   sel_poisson <- which(ABC_binom_poisson == min(ABC_binom_poisson))
   
   # copy/paste...
   # get categories
   output_cat <- as.character(unique(sum_stat_obs$category))
   par(mfrow=c(3,3))
   # iterate over categories
   for(i_cat in output_cat){
      flag_out <- sum_stat_obs$category == i_cat & sum_stat_obs$bool_orig == TRUE
      sum(flag_out)
      dim(sum_stat_obs)
      dim(ABC_out$stats[sel_poisson,])
      
      # if over time
      if(sum(flag_out)>1){
         y_lim <- range(pretty(c(sum_stat_obs$value[flag_out],
                                 sum_stat_obs$value_low[flag_out],
                                 sum_stat_obs$value_high[flag_out],
                                 ABC_out$stats[sel_poisson,flag_out])),
                        na.rm=T)
         
         plot(sum_stat_obs$date[flag_out],
              sum_stat_obs$value[flag_out],
              ylim = y_lim,
              main = string2label(i_cat),
              ylab = i_cat)
         
         if(all(!is.na(sum_stat_obs$value_low[flag_out]))){
            add_interval(x = sum_stat_obs$date[flag_out],
                         y1 = sum_stat_obs$value_low[flag_out],
                         y2 = sum_stat_obs$value_high[flag_out])
         }
         
         # for(i_out in 1:nrow(ABC_out$stats)){
            i_out <- sel_poisson
            lines(sum_stat_obs$date[flag_out],
                  ABC_out$stats[i_out,flag_out],
                  col=alpha(4,0.8))
         # }
       } else{ # else: hist + reference
   #       
   #       # initial doubling time
   #       x_lim <- range(0,pretty(sum_stat_obs$value_high[flag_out]*1.1))
   #       hist(ABC_out$stats[,flag_out],40,
   #            xlim=x_lim,
   #            xlab=i_cat,
   #            main=i_cat,
   #            col=alpha(4,0.2),
   #            border = alpha(4,0.2))
   #       #   abline(v=mean(sum_stat_obs$value[flag_out]),col=1)
   #       add_interval_hor(x1  = sum_stat_obs$value_low[flag_out],
   #                        x2  = sum_stat_obs$value_high[flag_out],
   #                        x_mean = sum_stat_obs$value[flag_out],
   #                        col = 1)
   #    }
   #    
       }
   }
   
   
   
   
}



