############################################################################# #
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
############################################################################# #
#
# MODEL PARAMETER ESTIMATION
#
############################################################################# #
library(mgcv)


if(0==1){
  
  # Load rStride
  source('./bin/rstride/rStride.R')
  
  output_dir <- '/Users/lwillem/Documents/university/research/stride/results_covid19_revision/20200918_fitting/param_estim/'
  output_folders <- dir(output_dir,full.names = T)
  output_folders <- output_folders[!grepl('pdf',output_folders)]
  i_folder <- output_folders[6]
  for(i_folder in output_folders){
    print(basename(i_folder))
    estimate_parameters(i_folder)
    
  }
}

#' @param project_dir   name of the project folder
estimate_parameters <- function(project_dir)
{
  # command line message
  smd_print('START PARAMETER ESTIMATION...')
  
  # load project summary
  project_summary    <- .rstride$load_project_summary(project_dir)
  
  # remove "initially_infected" (copy of num_infected_seeds)
  project_summary$initially_infected <- NULL
  
  # retrieve all variable model parameters
  input_opt_design   <- .rstride$get_variable_model_param(project_summary)
  
  # get all transmission output
  data_incidence_all <- .rstride$load_aggregated_output(project_dir,'data_incidence')
  
  if(length(data_incidence_all) == 1 && is.na(data_incidence_all)){
    smd_print('NO INCIDENCE DATA AVAILABLE.')
    return(.rstride$no_return_value())
  }
  
  # add config_id (this can be removed in the future)
  project_summary$config_id  <- .rstride$get_config_id(project_summary)
  input_opt_design$config_id <- .rstride$get_config_id(input_opt_design)
  
  # remove string based input parameters (e.g. different calendars)
  flag_character  <- unlist(lapply(input_opt_design,is.character)) & names(input_opt_design) != "config_id"
  input_opt_design <- input_opt_design[,!flag_character]
  
  # add config_id to incidence data
  data_incidence_all            <- merge(data_incidence_all,project_summary[,c('exp_id','config_id')] )
  
  # remove rows missing sim_date
  data_incidence_all <- data_incidence_all[!is.na(data_incidence_all$sim_date),]
  
  # check for NA's, and replace by 0
  data_incidence_all[is.na(data_incidence_all)] <- 0

  ## HOSPITAL REFERENCE DATA COVID-19 ----
  # use (local version of) most recent SCIENSANO data (or backup version)
  ref_data          <- get_observed_incidence_data()
  ref_data$sim_date <- as.Date(ref_data$sim_date)
  dim(ref_data)
  
  # reformat reference data
  hosp_adm_data <- data.frame(date = ref_data$sim_date,
                              num_adm = ref_data$hospital_admissions,
                              cum_adm = ref_data$cumulative_hospital_admissions)
  
  # remove reference data if simulation period is shorter
  flag_compare  <- hosp_adm_data$date %in% data_incidence_all$sim_date
  #flag_compare  <- flag_compare & hosp_adm_data$date < as.Date("2020-05-15")
  hosp_adm_data <- hosp_adm_data[flag_compare,]
  
  ## SERO-PREVALENCE DATA ----
  prevalence_ref <- load_observed_seroprevalence_data()
  
  # select simulation period
  sel_ref_dates  <- prevalence_ref$seroprevalence_date %in% data_incidence_all$sim_date
  prevalence_ref <- prevalence_ref[sel_ref_dates,]
  
  ## DOUBLING TIME 3.1 (2.4-4.4) \cite{pellis2020challenges} ----
  sim_start_date <- as.Date(unique(project_summary$start_date))
  sim_ref_dates  <- seq(sim_start_date+10,sim_start_date+26,1)
  ref_doubling_time <- data.frame(dates = seq(as.Date('2020-02-24'),as.Date('2020-03-08'),1),
                                  mean = 3.1,
                                  CI_low = 2.4,
                                  CI_upper = 4.4
                                  )
  
  ## GENERATION INTERVAL 5.20 (3.78-6.78) days \cite{ganyani2020estimating} ----
  ref_generation_interval  <- data.frame(dates = sim_ref_dates,
                                         mean = 5.20,
                                         CI_low = 3.78,
                                         CI_upper =6.78
  )
  
  # create matrix to keep track of the scores
  df_loglike <- data.frame(config_id = rep(NA,nrow(project_summary)),
                           hospital_pois = NA,
                           incidence_pois = NA,
                           doubling_pois = NA,
                           generation_pois = NA,
                           doubling_time_baseline = NA
  )
  
  # reset doubling time
  data_incidence_all$doubling_time <- NA
  
  i_exp <- 1  
  # loop over each experiment ----
  for(i_exp in unique(data_incidence_all$exp_id)){

    # select one run
    flag_exp       <- data_incidence_all$exp_id == i_exp
    
    # 0. save config_id
    df_loglike$config_id[i_exp] <- unique(data_incidence_all$config_id[flag_exp])
    
    # 1. Hospital admissions
    flag_hosp_data       <- flag_exp & data_incidence_all$sim_date %in% hosp_adm_data$date
    flag_hosp_ref        <- hosp_adm_data$date %in% data_incidence_all$sim_date[flag_exp]
    
    hosp_adm_observed                <- hosp_adm_data$num_adm[flag_hosp_ref]
    hosp_adm_model                   <- data_incidence_all$new_hospital_admissions[flag_hosp_data]
    df_loglike$hospital_pois[i_exp]  <- get_binom_poisson(hosp_adm_observed,hosp_adm_model)
    
    # 2. sero-prevalence
    ref_dates                  <- prevalence_ref$seroprevalence_date
    total_incidence_observed   <- prevalence_ref$point_incidence_mean
    total_incidence_model      <- data_incidence_all$cumulative_infections[ flag_exp & data_incidence_all$sim_date %in% ref_dates]
    df_loglike$incidence_pois[i_exp]   <- get_binom_poisson(total_incidence_observed,total_incidence_model)
    
    # 3. doubling time
    ref_dates                  <- ref_doubling_time$dates
    doubling_time_observed     <- mean(ref_doubling_time$mean)
    
    flag_exp_doubling   <- flag_exp & data_incidence_all$sim_date %in% ref_dates
    if(sum(flag_exp_doubling)>0){
      doubling_time_model                      <- get_doubling_time(data_incidence_all$new_infections[flag_exp_doubling])
      df_loglike$doubling_pois[i_exp]          <- get_binom_poisson(doubling_time_observed,doubling_time_model)
      df_loglike$doubling_time_baseline[i_exp] <- as.numeric(doubling_time_model)
      data_incidence_all$doubling_time[flag_exp_doubling] <- as.numeric(doubling_time_model)

    } else{
      df_loglike$doubling_pois[i_exp]          <- 0
      df_loglike$doubling_time_baseline[i_exp] <- 0
      data_incidence_all$doubling_time[flag_exp_doubling] <- 0

    }
     
    # 4. generation interval
    ref_dates                          <- ref_generation_interval$dates
    generation_interval_observed       <- ref_generation_interval$mean
    generation_interval_model          <- data_incidence_all$gen_interval[flag_exp & data_incidence_all$sim_date %in% ref_dates]
    df_loglike$generation_pois[i_exp]  <- get_binom_poisson(generation_interval_observed,generation_interval_model)
    
  }
  
  # remove rows with NA's => incidence > threshold
  df_loglike <- df_loglike[!is.na(df_loglike$config_id),]
  
  # save as rds file
  df_loglike_filename <- smd_file_path(project_dir,paste0(basename(project_dir),'_loglike.RData'))
  saveRDS(df_loglike,df_loglike_filename)
  #df_loglike <- readRDS(df_loglike_filename)
  
  # GET RESULTS ####
  df_loglike_orig <- df_loglike
  # all results (including stochastic error)
  if(length(unique(df_loglike_orig$config_id))>1)
  select_ensemble_and_plot(df_loglike_orig,input_opt_design,hosp_adm_data,
                           data_incidence_all,project_summary,prevalence_ref)
  
  # based on average per parameter sets
  select_ensemble_and_plot(df_loglike_orig,input_opt_design,hosp_adm_data,
                           data_incidence_all,project_summary,prevalence_ref,
                           bool_average = T)
  
}

select_ensemble_and_plot <- function(df_loglike,input_opt_design,hosp_adm_data,
                                     data_incidence_all,project_summary,prevalence_ref,
                                     bool_average = F) {
  
  # if multiple realisations AND bool = TRUE  ==>> get average  
  if(bool_average) {
    df_loglike <- aggregate( .  ~ config_id,data=df_loglike,mean,na.action=na.pass)
  }
  dim(df_loglike)
  
  # merge with parameters
  df_loglike <- merge(df_loglike,input_opt_design,by.x='config_id',by.y='config_id')
  
  df_loglike$pareto_front <- FALSE
  q_value <- 0.13
  q_increase <- 0.02
  while((q_value+q_increase < 1) && 
         length(unique(df_loglike$config_id[df_loglike$pareto_front]))<10){
    q_value <- q_value+q_increase
    
    x_pq <- quantile(df_loglike$hospital_pois,q_value,na.rm=T)
    y_pq <- quantile(df_loglike$incidence_pois,q_value,na.rm=T)
    z_pq <- quantile(df_loglike$doubling_pois,q_value,na.rm=T)
    
    bool_pareto_front <- df_loglike$hospital_pois <= x_pq &
                          df_loglike$incidence_pois <= y_pq &
                          df_loglike$doubling_pois <= z_pq
    bool_pareto_front[is.na(bool_pareto_front)] <- FALSE
    df_loglike$pareto_front <- bool_pareto_front
    #print(length(unique(df_loglike$config_id[df_loglike$pareto_front])))
  }
  q_value <- round(q_value,digits=2)
  table(df_loglike$pareto_front )
  table(table(df_loglike$config_id[df_loglike$pareto_front]))
  
  
  # get names (without ids)
  param_names <- names(input_opt_design)
  param_names <- param_names[!grepl('id',param_names)]
  
  # create plot with scores
  .rstride$create_pdf(project_dir,paste0('parameter_inspection',ifelse(bool_average,'_avg','')),width = 6, height = 7)
  
  i_param <- param_names[1]
  par(mfrow=c(3,4))
  for(i_param in param_names){
  
    if(length(unique(df_loglike[,i_param]))>10){
      plot(formula(paste('hospital_pois ~ ',i_param)), data = df_loglike, xlab=i_param,ylab='Pois neg log like',main='Hospital Admissions')
      plot(formula(paste('incidence_pois ~ ',i_param)), data = df_loglike,xlab=i_param,ylab='Pois neg log like',main='Total Incidence')
      plot(formula(paste('doubling_pois ~ ',i_param)), data = df_loglike, xlab=i_param,ylab='Pois neg log like',main='Doubling Time')
    } else{
      boxplot(formula(paste('hospital_pois ~ ',i_param)), data = df_loglike,xlab=i_param,ylab='Pois neg log like',main='Hospital Admissions')
      boxplot(formula(paste('incidence_pois ~ ',i_param)), data = df_loglike,xlab=i_param,ylab='Pois neg log like',main='Total Incidence')
      boxplot(formula(paste('doubling_pois ~ ',i_param)), data = df_loglike,xlab=i_param,ylab='Pois neg log like',main='Doubling Time')
      
    }
     
    if(is.numeric(df_loglike[,i_param])){
      boxplot(formula(paste(i_param,' ~ pareto_front')), data = df_loglike,xlab=i_param,horizontal=T,ylab='Pareto front?',main=paste0('Pareto front\n(Q',q_value,')'))
      title(sub=)
    } else{
      boxplot(0,col=0,border=0)
    }
  }
  
  # Show pareto front (using 2 scales)
  plot_range_q <- 0.4
  for(plot_range_q in c(1,0.4)){
    par(mfrow=c(2,2))
    plot(df_loglike$hospital_pois,
         df_loglike$doubling_pois,
         xlab='Pois negloglike (hospital admissions)',
         ylab='Pois negloglike (doubling time)',
         xlim=quantile(df_loglike$hospital_pois,c(0,plot_range_q),na.rm=T),
         ylim=quantile(df_loglike$doubling_pois,c(0,plot_range_q),na.rm=T))
    points(df_loglike$hospital_pois[df_loglike$pareto_front],
           df_loglike$doubling_pois[df_loglike$pareto_front],
           col = 4,
           pch=19)
    legend('topright',
           paste(c('N_total','N_output','N_pareto'),c(nrow(project_summary),nrow(df_loglike),sum(df_loglike$pareto_front)),sep=': '),
           cex=0.8,
           bg='white')
    plot(df_loglike$incidence_pois,
         df_loglike$doubling_pois,
         xlab='Pois negloglike (total_incidence)',
         ylab='Pois negloglike (doubling time)',
         xlim=quantile(df_loglike$incidence_pois,c(0,plot_range_q),na.rm=T),
         ylim=quantile(df_loglike$doubling_pois,c(0,plot_range_q),na.rm=T))
    points(df_loglike$incidence_pois[df_loglike$pareto_front],
           df_loglike$doubling_pois[df_loglike$pareto_front],
           col = 4,
           pch=19)
    plot(df_loglike$hospital_pois,
         df_loglike$incidence_pois,
         xlab='Pois negloglike (hospital admissions)',
         ylab='Pois negloglike (total incidence)',
         xlim=quantile(df_loglike$hospital_pois,c(0,plot_range_q),na.rm=T),
         ylim=quantile(df_loglike$incidence_pois,c(0,plot_range_q),na.rm=T))
    points(df_loglike$hospital_pois[df_loglike$pareto_front],
           df_loglike$incidence_pois[df_loglike$pareto_front],
           col = 4,
           pch=19)
  }
  dev.off()
  
  # select incidence data
  config_tag_sel          <- unique(df_loglike$config_id[df_loglike$pareto_front])
  flag_plot               <- data_incidence_all$config_id %in% config_tag_sel
  data_incidence_ensemble <- data_incidence_all[flag_plot,]
  
  ## PARETO PLOTS ####
  .rstride$create_pdf(project_dir,paste0('parameter_pareto_incidence',ifelse(bool_average,'_avg','')),width = 6, height = 7)
  par(mfrow=c(4,1))
  
  bool_add_param <- TRUE
  
  # plot temporal patterns (ensemble)
  plot_incidence_data(data_incidence_ensemble,
                      project_summary,
                      hosp_adm_data,
                      input_opt_design,
                      prevalence_ref,
                      bool_add_param = bool_add_param,
                      bool_add_doubling_time = TRUE)
  
  ## ENSEMBLE (SINGLE RUNS) ####
  opt_config_id <- config_tag_sel
  i_config <- opt_config_id[1]
  for(i_config in opt_config_id){
    
    # select subset
    data_incidence_sel <- data_incidence_all[data_incidence_all$config_id == i_config,]
    data_incidence_sel$cumulative_hospital_cases <- data_incidence_sel$cumulative_hospital_cases - data_incidence_sel$cumulative_hospital_cases[1]
    # plot
    plot_incidence_data(data_incidence_sel,project_summary,
                        hosp_adm_data,input_opt_design,prevalence_ref,
                        bool_add_param = bool_add_param,
                        bool_add_doubling_time = TRUE)
  }
  
  # close pdf
  dev.off()
  
  
  ## BEST CONFIG FOR 2 TARGETS (mean) ----
  #table(df_loglike$pareto_front,df_loglike$pareto_num)
  df_loglike_mean <- aggregate( .  ~ config_id,data=df_loglike[df_loglike$pareto_front>0,!flag_character],mean,na.action=na.pass)
  #df_loglike_mean <- df_loglike
  order_table     <- data.frame(hosp   = order(df_loglike_mean$hospital_pois),
                                inc    = order(df_loglike_mean$incidence_pois),
                                double = order(df_loglike_mean$doubling_pois))
  
  i_row <- 1
  tbl <- table(unlist(order_table[1:i_row,-3]))
  while(!any(tbl>1)){
    i_row <- i_row+1
    tbl <- table(unlist(order_table[1:i_row,-3]))
    print(i_row)
  }
  
  # select subset
  i_config_single    <- df_loglike_mean$config_id[as.numeric(names(which(tbl>1)))[1]]
  data_incidence_sel <- data_incidence_all[data_incidence_all$config_id == i_config_single,]
  data_incidence_sel$cumulative_hospital_cases <- data_incidence_sel$cumulative_hospital_cases - data_incidence_sel$cumulative_hospital_cases[1]
  
  .rstride$create_pdf(project_dir,paste0('parameter_pareto_incidence_single',ifelse(bool_average,'_avg','')),width = 6, height = 7)
  par(mfrow=c(4,1))
  # plot
  plot_incidence_data(data_incidence_sel,project_summary,
                      hosp_adm_data,input_opt_design,prevalence_ref,
                      bool_add_param = bool_add_param,
                      bool_add_doubling_time = TRUE)
  dev.off()
  
  
  # select subset, best hospital score
  i_config_single    <- df_loglike_mean$config_id[order_table$hosp[1]]
  data_incidence_sel <- data_incidence_all[data_incidence_all$config_id == i_config_single,]
  data_incidence_sel$cumulative_hospital_cases <- data_incidence_sel$cumulative_hospital_cases - data_incidence_sel$cumulative_hospital_cases[1]
  
  .rstride$create_pdf(project_dir,paste0('parameter_pareto_incidence_single_hosp',ifelse(bool_average,'_avg','')),width = 6, height = 7)
  par(mfrow=c(4,1))
  # plot
  plot_incidence_data(data_incidence_sel,project_summary,
                      hosp_adm_data,input_opt_design,prevalence_ref,
                      bool_add_param = bool_add_param,
                      bool_add_doubling_time = TRUE)
  dev.off()
  
  
  ## save scores ----
  df_loglike_summary <- merge(project_summary,df_loglike)
  saveRDS(df_loglike_summary,smd_file_path(project_dir,paste0(basename(project_dir),paste0('_poison_neg_loglikelihood_scores',ifelse(bool_average,'_avg',''),'.RData'))))
  
  ## get parameter configuration of the pareto ensemble
  pareto_param <- .rstride$get_variable_model_param(df_loglike_summary[df_loglike_summary$pareto_front,])
  dim(pareto_param)
  
  # remove duplicates
  if('pareto_num' %in% names(pareto_param)){
    pareto_param <- aggregate(. ~ pareto_num, data=pareto_param,mean)
  }
  
  # add pareto number and project_dir as run info
  # note: columns with "num" will be handled as integers in a LHS
  pareto_param$pareto_num  <- 1:nrow(pareto_param)
  pareto_param$q_selection <- q_value
  pareto_param$num_df_loglike <- nrow(df_loglike)
  pareto_param$run_info    <- basename(project_dir)
  
  # # add boolean for "single best"
  pareto_param$config_id_single <- i_config_single
  
  # save as .RData and csv
  saveRDS(pareto_param,smd_file_path(project_dir,paste0(basename(project_dir),'_config_pareto_ensemble',ifelse(bool_average,'_avg',''),'.RData')))
  write.table(pareto_param,smd_file_path(project_dir,paste0(basename(project_dir),'_config_pareto_ensemble',ifelse(bool_average,'_avg',''),'.csv')),sep=',',row.names=F)
  
  # PARETO FRONT PARAMETERS
  param_design <- names(input_opt_design)
  param_design <- param_design[param_design != 'config_id']
  param_design
  
  ## PARETO PLOTS ####
  if(length(param_design)>2){
    .rstride$create_pdf(project_dir,paste0('parameter_pareto_config',ifelse(bool_average,'_avg','')),width = 6, height = 7)
    par(mfrow=c(3,3))
    for(i_param in 1:length(param_design))
      for(j_param in i_param:length(param_design)){
        if(i_param != j_param)
          plot(df_loglike_summary[df_loglike_summary$pareto_front,param_design[i_param]],
               df_loglike_summary[df_loglike_summary$pareto_front,param_design[j_param]])
      }
    
    df_loglike_cor <- cor(df_loglike_summary[df_loglike_summary$pareto_front,param_design])
    
    # parameter correlation?
    par(mfrow=c(1,1))
    corrplot(df_loglike_cor)  # using function from 'corrplot' library
    
    # add heatmap
    if(!any(is.na(df_loglike_cor))){
      palette = colorRampPalette(c("green", "white", "red")) (20)
      heatmap(x = df_loglike_cor, col = palette, symm = TRUE)
    }
    
    dev.off()
  }
  
  
  # command line message
  smd_print('PARAMETER ESTIMATION COMPLETE')
  
} # end function


# incidence_observed <- hosp_adm_observed;incidence_predicted <- hosp_adm_model
get_lsq_score <- function(incidence_observed,incidence_predicted){
  
  # # from the model estimate of the incidence and the data, we can 
  # # estimate the fraction of cases that were confirmed
  # frac_confirmed = sum(incidence_observed,na.rm=T)/sum(incidence_predicted,na.rm=T)
  # 
  # # normalize the model prediction so area under curve equals the sum of the data incidence
  # incidence_predicted = incidence_predicted*frac_confirmed 
  # 
  if (length(incidence_predicted)==length(incidence_observed)
      &!is.na(sum(incidence_predicted))){
    
    # calculate the least squares statistic
    leastsquares = sum((incidence_predicted-incidence_observed)^2,na.rm=T)
    
    return(leastsquares)
    
  } else{
    return(NA)
  }
}
#incidence_observed <- doubling_time_observed;incidence_predicted <- doubling_time_model
#incidence_observed <- hosp_adm_observed; incidence_predicted <- hosp_adm_model    
get_binom_poisson <- function(incidence_observed,incidence_predicted){

  # # from the model estimate of the incidence and the data, we can 
  # # estimate the fraction of cases that were confirmed
  #frac_confirmed = sum(incidence_observed,na.rm = T)/sum(incidence_predicted,na.rm = T)
  # 
  # # normalize the model prediction so area under curve
  # # equals the sum of the data incidence
  #incidence_predicted = incidence_predicted*frac_confirmed 
  
  # adjust
  incidence_predicted[incidence_predicted==0] <- 1
  
  # now calculate the Poisson neglog likelihood
  # statistic that compares the data to this model calculated
  # under a particular hypothesis of R0 and t0
  
  if (length(incidence_predicted)==length(incidence_observed)
      &&!is.na(sum(incidence_predicted))
      && all(incidence_predicted > 0) ){
    
    # calculate the Poisson negative log-likelihood statistic
    # see http://sherrytowers.com/2014/07/10/poisson-likelihood/
    pois_negloglike = sum(-incidence_observed*log(incidence_predicted)+incidence_predicted,na.rm = T)
    return(pois_negloglike)
    
  } else{
    return(NA)
  }
}


#vpois_negloglike <- df_loglike$hospital_pois; num_obs <- nrow(df_loglike)
#vpois_negloglike <- df_loglike$doubling_pois; num_obs <- nrow(df_loglike)
get_normal_negloglike <- function(num_obs,vpois_negloglike){
  
  ################################################################################# #
  ################################################################################# #
  # work out the values of the likelihood that represent the one
  # std dev uncertainty fmin + 1/2
  # 95% CI fmin + 1.96^2/2
  # with fmin equal to the minimum value of the likelihood
  # (the Normal negative log likelihood in this case)
  ##################################################################################
  lmin = which.min(vpois_negloglike)
  a = min(vpois_negloglike,na.rm=T)
  b = min(vpois_negloglike,na.rm=T)+0.5
  c = min(vpois_negloglike,na.rm=T)+1.96^2/2
  ib = which(vpois_negloglike<=b)
  ic = which(vpois_negloglike<=c)
  
  R0_best  <- df_loglike$r0[vpois_negloglike == a]
  vR0      <- df_loglike$r0
  t0_best  <- df_loglike$num_infected_seeds[which(vpois_negloglike == a)]
  vt0      <- df_loglike$num_infected_seeds
  
  plot(vR0,vpois_negloglike)
  
  around = 3
  eR0_best_1stddev = round((range(vR0[ib])[2]-range(vR0[ib])[1])/2,around)
  eR0_best_95CI_lo = round(range(vR0[ic])[1],around)
  eR0_best_95CI_hi = round(range(vR0[ic])[2],around)
  
  et0_best_1stddev = (range(vt0[ib])[2]-range(vt0[ib])[1])/2
  et0_best_95CI_lo = range(vt0[ic])[1]
  et0_best_95CI_hi = range(vt0[ic])[2]
  
  ################################################################################# #
  ################################################################################# #
  # print out the best-fit values and the one std dev uncertainty and
  # 95% CI
  ################################################################################# #
  cat("\n")
  cat(paste("The best fit R0 is ",round(R0_best,around)," with one std dev uncertainty ",eR0_best_1stddev," and 95% CI [",eR0_best_95CI_lo,",",eR0_best_95CI_hi,"]\n",sep=""))
  cat(paste("The best fit t0 is ",round(t0_best,around)," with one std dev uncertainty ",et0_best_1stddev," and 95% CI [",et0_best_95CI_lo,",",et0_best_95CI_hi,"]\n",sep=""))
  
  
}

# based on https://github.com/thomasallanhouse/covid19-growth
# dat <- summary_out$new_infections;npts=200; plt=TRUE; figtitle="db";bool_confidence=FALSE
get_doubling_time <- function(dat, npts=200,bool_confidence=FALSE){
  
  # if no new infections => doubling time is infinity, so return NA
  if(length(unique(dat))==1 || sum(!is.na(dat))==1){
    return(NA)
  }
  
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))#,expdt=rep(0,3))
  Tv <- seq(1,length(dat),1)
  
  MGLM <- glm(dat~(Tv), family=quasipoisson)
  doubling<-data.frame(mean=MGLM$coefficients[2])  # instantanious growth rate
  
  if(bool_confidence){
    doubling$ci_025= confint(MGLM)[2,1]
    doubling$ci_950 =confint(MGLM)[2,2]
  }
  
  # go from growth rate to doubling time
  doubling<-log(2)/doubling # doubling time
  
  # make numeric
  doubling <- as.numeric(doubling)
  
  return(doubling)
}

# dat <- summary_infections$new_infections
get_longitudinal_doubling_time <- function(dat){
  
  # if no new infections => doubling time is infinity, so return NA
  if(length(unique(dat))==1 || sum(!is.na(dat))==1){
    return(data.frame(mean=NA,sd=NA))
  }
  
  npts <- length(dat)
  n_knots <- min(10,length(dat))
  
  # fix if number of new infections is almost zero
  if(n_knots <= 2){
    return(data.frame(mean=NA,sd=NA))
  }
  
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))#,expdt=rep(0,3))
  Tv <- seq(1,length(dat),1)
  MGAM <- gam(dat~s(Tv,k=n_knots), family=quasipoisson)
  
  xv<-seq(1,length(dat), length=npts)
  newd <- data.frame(Tv=xv)
  X0 <- predict(MGAM, newd,type="lpmatrix")
  
  eps <- 1e-7 ## finite difference interval
  xv <- seq(1,length(dat), length=npts) +eps
  newd <- data.frame(Tv=xv)
  X1 <- predict(MGAM, newd,type="lpmatrix")
  
  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  
  Xi <- Xp*0 
  Xi[,-1] <- Xp[,-1] ## Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(MGAM)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  # plot(xv,log(2)/df,type="l", ylab='Days', xlab='Time', main='Doubling time', ylim=c(0,10)) #,ylim=range(log(2)/c(df+2*df.sd,df-2*df.sd)))
  # lines(xv,log(2)/(df+2*df.sd),lty=2);
  # lines(xv,log(2)/ifelse(df<2*df.sd,0.00001, df-2*df.sd),lty=2)
  
  doubling_longitudinal <- data.frame(mean=log(2)/df,
                                      sd = log(2)/df.sd)
  return(doubling_longitudinal)
}

#dat <- summary_infections$new_infections;npts=200; plt=TRUE; figtitle=""
analyse_doubling_time <- function(dat, npts=200, plt=FALSE, figtitle=""){
  
  res<- data.frame(sdt=rep(0,npts),sdtup=rep(0,npts),sdtlow=rep(0,npts))#,expdt=rep(0,3))
  Tv <- seq(1,length(dat),1)
  MGAM <- gam(dat~s(Tv), family=quasipoisson)
  
  xv<-seq(1,length(dat), length=npts)
  newd <- data.frame(Tv=xv)
  X0 <- predict(MGAM, newd,type="lpmatrix")
  
  eps <- 1e-7 ## finite difference interval
  xv <- seq(1,length(dat), length=npts) +eps
  newd <- data.frame(Tv=xv)
  X1 <- predict(MGAM, newd,type="lpmatrix")
  
  Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
  
  Xi <- Xp*0 
  Xi[,1:9+1] <- Xp[,1:9+1] ## Xi%*%coef(MGAM) = smooth deriv i
  df <- Xi%*%coef(MGAM)              ## ith smooth derivative 
  df.sd <- rowSums(Xi%*%MGAM$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
  
  res$sdt <- df
  res$sdtup <- df+2*df.sd
  res$sdtlow <- df-2*df.sd
  res$sdt <- ifelse(res$sdt < 0, 0.0001, res$sdt)
  res$sdt <- log(2)/res$sdt
  res$sdtup <- ifelse(res$sdtup < 0, 0.0001, res$sdtup)
  res$sdtup <- log(2)/res$sdtup
  res$sdtlow <- ifelse(res$sdtlow < 0, 0.0001, res$sdtlow)
  res$sdtlow <- log(2)/res$sdtlow
  
  MGLM <- glm(dat~(Tv), family=quasipoisson)
  Doubling<-c(MGLM$coefficients[2],confint(MGLM)[2,1],confint(MGLM)[2,2])
  #  res$expdt <- Doubling
  
  if(plt==TRUE){
    par(mfrow=c(1,3))
    plot(Tv, dat, main='Fit compared with model', ylab='Cases', xlab='Time')
    lines(Tv, fitted(MGAM))
    lines(Tv, fitted(MGLM), col=2)
    
    p <- predict(MGAM, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGAM$family$linkinv(upr)
    lwr <- MGAM$family$linkinv(lwr)
    lines(Tv, upr, col=1, lty=2)
    lines(Tv, lwr, col=1, lty=2)
    
    p <- predict(MGLM, type = "link", se.fit = TRUE)
    upr <- p$fit + (2 * p$se.fit)
    lwr <- p$fit - (2 * p$se.fit)
    upr <- MGLM$family$linkinv(upr)
    lwr <- MGLM$family$linkinv(lwr)
    lines(Tv, upr, col=2, lty=2)
    lines(Tv, lwr, col=2, lty=2)
    
    plot(xv,df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)), ylab='Instantaneous growth rate', xlab='Time', main=figtitle)
    lines(xv,df+2*df.sd,lty=2);
    lines(xv,df-2*df.sd,lty=2)
    abline(h=Doubling, lty=c(1,2,2), col=2)
    
    plot(xv,log(2)/df,type="l", ylab='Days', xlab='Time', main='Doubling time', ylim=c(0,10)) #,ylim=range(log(2)/c(df+2*df.sd,df-2*df.sd)))
    lines(xv,log(2)/(df+2*df.sd),lty=2);
    lines(xv,log(2)/ifelse(df<2*df.sd,0.00001, df-2*df.sd),lty=2)
    Doubling<-log(2)/Doubling
    abline(h=Doubling, lty=c(1,2,2), col=2)
    
    #    mtext(figtitle, outer=TRUE,  cex=1, line=0)
  }
  
  #res
  return(log(2)/Doubling)
}




