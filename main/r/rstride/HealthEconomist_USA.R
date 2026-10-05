#############################################################################
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
#  Copyright 2020, Willem L, Kuylen E & Broeckhove J
#############################################################################

#############################################################################
# ESTIMATE COST-EFFECTIVENESS                                              ##
#############################################################################

get_cea_param <- function(num_samples = 1)
{
  
  # # REFERENCE: Zwanziger, Szilagi & PadmaKaul, Health Service Research, 2001 Oct; 36(5): 885–909.
  # # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1089266/pdf/hsresearch00006-0069.pdf
  # data_cea <- data.frame(saving_per_averted_case    = c(2089,(2089/2)),              # dollar2001
  #                        QALYgain_per_averted_case  = c(0.086,(0.086/2)),
  #                        price_dose         = c(6.3,0),                            # Beutels & Gay 2003 
  #                        price_admin_dose   = c(7.3,0))                            # Beutels & Gay 2003 
  # rownames(data_cea) <- c('mean','se')
  # 
  # # source: Carabin et al. (2003) The cost of measles in industrialised countries
  # data_cea <- list( reporting_prob     = 0.75,      # Carabin et al 2003 (Upper Income Country)
  #                   severe_case_prob   = 0.8,       # Carabin et al 2003 (Upper Income Country)
  #                   hospital_stay_days = 4,         # Carabin et al 2003 (Upper Income Country)
  #                   QALY_lost          = 6.8/365,   # Thorrington et al 2014 (UK)
  #                         
  #                   age_cat            = c(0  , 5,   10, 15, 20, 100),   # Beutels & Gay 2003 
  #                   cfr_per100000      = c(18 , 9.6, 32, 87.2, 92.6),    # Beutels & Gay 2003 
  #                   hosp_rate          = c(5.8, 1.7, 1.7, 1.6, 5.3)/100, # Beutels & Gay 2003 
  #                   direct_cost_case   = c(123, 112, 110, 111,148),      # Beutels & Gay 2003 
  #                   life_expectancy    = c(73, 68, 63, 58, 27),          # Beutels & Gay 2003 
  #                   price_dose         = 6.3,                            # Beutels & Gay 2003 
  #                   price_admin_dose   = 7.3)                            # Beutels & Gay 2003 
  # 
  # 
  # burden <- list(age_cat  = c(0  , 5,   10, 20, 30, 100),    #Perry & Halsey 2004
  #                total    = c(28730,6492,18580,9161,4068),   #Perry & Halsey 2004
  #                any      = c(11883,1173,2369,2656,1399),    #Perry & Halsey 2004
  #                hospital = c(7470,612,1612,2075,1107),      #Perry & Halsey 2004
  #                death    = c(97,9,18,26,27))
  #
  # # convert case fatality rate
  # data_cea$cfr_per_case <- data_cea$cfr_per100000/1e5
  # 
  # # unit cost: 2001USD
  # data_cost2001 <- data.frame(cost_hospital_day = 253,      # Carabin et al 2003 (BE)
  #                             cost_gp_visit     = 16,       # Carabin et al 2003 (BE)
  #                             cost_mmr_delivery = 4.3,      # Carabin et al 2003 (BE)
  #                             cost_mmr_dose     = 5)        # Carabin et al 2003 (BE)
  # 
  # cpi_2001USD <- 81.20256846   #2001USDA
  # cpi_2017USD <- 112.4115573   #2017USD
  # conversion_2001USDA_2017EURO <- 0.80          #2017EURO == 2017USD
  # data_cost <- data_cost2001 / cpi_2001USD * cpi_2017USD * conversion_2001USDA_2017EURO
     
  # initialise CEA data
  data_cea_sample <- data.frame(cea_id = 1:num_samples)
  
  # outpatient cost
  cost_GP_visit        = 408.08 # $324.54 (Pike et al., 2022: outpatient cost $22,718 ÷ 70 cases, 2019 Q1$) - adjusted for inflation 2025
  outpatient_cost_mean = cost_GP_visit
  outpatient_cost_se   = 0
  data_cea_sample$outpatient_unit_cost <- rnorm(num_samples,outpatient_cost_mean,outpatient_cost_se)
  
  # inpatient cost
  # ref: KFF hospital expense per inpatient day
  ## How to make age specific?
  inpatient_cost_mean = 1700 # euro 
  inpatient_cost_se   = 0
  data_cea_sample$inpatient_unit_cost <- rnorm(num_samples,inpatient_cost_mean,inpatient_cost_se)
  
  
  # hospital probability
  # ref: CDC (https://www.cdc.gov/measles/data-research/index.html, accessed 31 July 2026)
  study_hosp_cases  = 106
  study_total_cases = 584
  data_cea_sample$hospital_probability <- rbeta(num_samples,study_hosp_cases,study_total_cases-study_hosp_cases)
    
  # program cost
  # ref: RIZIV + BCFI
  num_vfc      <- 0.522 # approximately 52.2% of U.S. children were VFC eligible. https://www.cdc.gov/mmwr/volumes/73/wr/mm7333e1.htm
  price_admin  <- cost_GP_visit #TODO?
  price_dose   <- 22.14
  data_cea_sample$price_dose_admin <- price_dose + price_admin
  
  # QALY
  # ref: Thorrington et al, Plos ONE, 2012
  qaly_loss_mean = 0.019
  qaly_loss_se   = 2 * sqrt((0.022-0.016)/203)
  # normal distribution => include ZERO
  #data_cea_sample$qaly_loss_case <- rnorm(num_samples,qaly_loss_mean,qaly_loss_se)
 
  # so, use beta distribution
  rbeta_alpha = (qaly_loss_mean^2*(1-qaly_loss_mean)/qaly_loss_se^2)-qaly_loss_mean 
  rbeta_beta  = rbeta_alpha * (1-qaly_loss_mean) / qaly_loss_mean
  data_cea_sample$qaly_loss_case <- rbeta(num_samples,rbeta_alpha,rbeta_beta)
  
  return(data_cea_sample)
}


calculate_cost_effectiveness <- function(project_dir){
  
  # command line message
  smd_print('START COST-EFFECTIVENESS ANALYSIS...')
  
  # load project output summary
  project_summary    <- .rstride$load_project_summary(project_dir)
  
  # get average burden   
  average_burden_all <- get_average_burden_averted(project_dir)  
  head(average_burden_all)
  
  # save average burden
  write.table(average_burden_all,file=file.path(project_dir,paste0(project_summary$run_tag[1],'_cea_burden_average.csv')),row.names=F,sep=',')
  
  # add burden id 
  average_burden_all$burden_id <- average_burden_all$reference_id

  # open pdf stream
  pdf(file=file.path(project_dir,paste0(project_summary$run_tag[1],'_uncertainty.pdf')),10,5)
  #par(mfrow=1:2)
  legend_cex <- 0.6
  
  # for each burden type
  i_exp_burden_tag <- levels(average_burden_all$burden_exp_tag)[1]
  for(i_exp_burden_tag in levels(average_burden_all$burden_exp_tag)){
    average_burden <- average_burden_all[average_burden_all$burden_exp_tag == i_exp_burden_tag,]

    # get CEA parameters
    num_cea_samples  <- unique(project_summary$num_cea_samples)
    if(length(num_cea_samples)>1){
      smd_print('MULTIPLE CEA UNCERTAINTY DIMENTIONS ==>> ERROR',WARNING = TRUE ,FORCED = TRUE)
    }
    num_samples      <- num_cea_samples
    cea_param_sample <- get_cea_param(num_samples = num_samples)
  
    # combine cea param with stride results
    cea_rstride <- cbind(cea_id    = cea_param_sample$cea_id,
                         burden_id = unique(average_burden$burden_id))
    # merge
    cea_rstride <- merge(cea_rstride,average_burden)
    cea_rstride <- merge(cea_rstride,cea_param_sample)
  
    # cases averted: outpatient / inpatient
    cea_rstride$cases_averted_inpatient  <- round(cea_rstride$incr_mean_num_cases_averted * cea_rstride$hospital_probability)
    cea_rstride$cases_averted_outpatient <- cea_rstride$incr_mean_num_cases_averted - cea_rstride$cases_averted_inpatient
    
    # costs averted: outpatient / inpatient / total
    cea_rstride$cost_averted_outpatient <- cea_rstride$cases_averted_outpatient * cea_rstride$outpatient_unit_cost
    cea_rstride$cost_averted_inpatient  <- cea_rstride$cases_averted_inpatient  * cea_rstride$inpatient_unit_cost
    cea_rstride$cost_averted_total      <- cea_rstride$cost_averted_outpatient + cea_rstride$cost_averted_inpatient 
    
    # program costs
    cea_rstride$cost_program <- cea_rstride$incr_mean_num_vaccines_admin * cea_rstride$price_dose
    
    # total cost
    cea_rstride$incr_cost    <- cea_rstride$cost_program - cea_rstride$cost_averted_total
    
    # burden averted
    cea_rstride$qaly_gain    <- cea_rstride$incr_mean_num_cases_averted * cea_rstride$qaly_loss_case
  
    # incremental cost effectives ratio (ICER)
    cea_rstride$icer <- cea_rstride$incr_cost / cea_rstride$qaly_gain
    cea_rstride$icer[cea_rstride$qaly_gain == 0] <- NA
    
    # plotting...
    cea_legend <- get_factor_legend(cea_rstride$intervention_tag)
    names(cea_legend) <- c('intervention_tag','intervention_color')
    cea_rstride <- merge(cea_rstride,cea_legend)
    
    cea_legend_pch <- data.frame(name = levels(cea_rstride$burden_exp_tag),
                                 pch  = 1:nlevels(cea_rstride$burden_exp_tag))
    
    # CEA PLANE
    plot(cea_rstride$qaly_gain,
         cea_rstride$incr_cost,
         col=alpha(cea_rstride$intervention_color,0.1),
         main=i_exp_burden_tag,
         xlab='QALY gain',
         ylab='Incremental cost',
         #pch=16,
         pch=as.numeric(cea_rstride$burden_exp_tag),
         cex=1)
    abline(h=0,v=0,lty=2)
    legend('topright',
           paste(cea_legend$intervention_tag),
           fill = cea_legend$intervention_color,
           cex=legend_cex,
           title='Scenario',
           bg='white')
    # legend('right',
    #        paste(cea_legend_pch$name),
    #        pch = cea_legend_pch$pch,
    #        col=1,
    #        cex=legend_cex,
    #        title='Scenario',
    #        bg='white')
    
    # ## ICER
    # par(mar=c(10,5,1,1))
    # boxplot(icer ~ intervention_tag + burden_exp_tag,
    #         data = cea_rstride,
    #         xlab='intervention',
    #         ylab='icer',
    #         las=2)
    # abline(h=0,lty=2)
    
    # create experiment id
    input_param_names         <- c(names(cea_param_sample),'burden_id','r0')
    num_scenario              <- nlevels(cea_rstride$intervention_tag)
    num_exp                   <- nlevels(as.factor(cea_rstride$cea_id))
    
    # set WTP levels
    wtp_opt <- seq(0,5e5,length=101)
    num_wtp <- length(wtp_opt)
    
    # initiate output parameters
    prob_high_net_benefit      <- matrix(NA,num_scenario,num_wtp)  # CEAC
    prob_high_mean_net_benefit <- prob_high_net_benefit            # CEAF
    net_benefit_all            <- array(NA,dim=c(num_wtp,num_scenario,num_exp))
    net_benefit_legend         <- get_factor_legend(cea_rstride$intervention_tag)
    
    i_wtp <- 3
    for(i_wtp in 1:num_wtp){
      cea_rstride_nb             <- cea_rstride
      cea_rstride_nb$net_benefit <- get_net_benefit(cea_rstride$qaly_gain,cea_rstride$incr_cost,wtp_opt[i_wtp])  
      
      net_benefit_matrix        <- data.frame(matrix(NA,ncol=num_scenario,nrow=num_exp))
      names(net_benefit_matrix) <- levels(cea_rstride_nb$intervention_tag)
      
      i_scen <- 2
      for(i_scen in 1:num_scenario){
        flag              <- as.numeric(cea_rstride_nb$intervention_tag) == i_scen
        
        net_benefit_value <- cea_rstride_nb$net_benefit[flag]
        net_benefit_index <- as.numeric(cea_rstride_nb$cea_id[flag])
        
        net_benefit_matrix[net_benefit_index,i_scen] <- net_benefit_value
        
      }
      head(net_benefit_matrix)
      dim(net_benefit_matrix)
      
      input_param_scen <- unique(cea_rstride_nb[,input_param_names])
      
      # get highest net benefit per simulation
      high_net_benefit <- (apply(X = net_benefit_matrix, MARGIN = 1,FUN=function(X){X==max(X,na.rm=T)}))
      # aggregate to get a probability
      prob <- rowSums(high_net_benefit,na.rm=T) / num_exp
      
      # store probability of highest net benefit => CEAC
      prob_high_net_benefit[,i_wtp] <- prob
      
      # if not the highest "mean net benefit", set NA => CEAF
      mean_nb   <- colMeans(net_benefit_matrix,na.rm = T)
      prob_mean <- prob
      prob_mean[mean_nb != max(mean_nb,na.rm = T)] <- NA
      
      # store
      prob_high_mean_net_benefit[,i_wtp] <- prob_mean
      
      # store net_benefit
      dim(t(net_benefit_matrix))
      dim(net_benefit_all)
      net_benefit_all[i_wtp,,] <- t(net_benefit_matrix)
    }
    
    plot(range(wtp_opt),c(0,1.2),col=0,
         xlab = 'Willingness to pay for a QALY (1000 EURO)',
         ylab = 'Probability highest net benefit',
         main=i_exp_burden_tag,
         xaxt='n')
    axis(1,pretty(wtp_opt),pretty(wtp_opt)/1e3)
    grid(col=alpha(1,0.5))
    #abline(v=seq(0,max(wtp_opt),1000),lty=3,col=alpha(1,0.5))
    # CEAC
    for(i in 1:dim(prob_high_net_benefit)[1])
    {
      lines(wtp_opt,prob_high_net_benefit[i,],col=alpha(cea_legend$intervention_color[i],0.8),lwd=4,pch=20)
    }
    # CEAF
    for(i in 1:dim(prob_high_mean_net_benefit)[1])
    {
      lines( wtp_opt,prob_high_mean_net_benefit[i,],col=cea_legend$intervention_color[i],lwd=4) # add line
      points(wtp_opt,prob_high_mean_net_benefit[i,],col=1,lwd=2,pch=1)
    }
    legend_ncol <- ceiling((length(net_benefit_legend$color)+1)/2)
    legend('topleft',
           c(net_benefit_legend$name,'CEAF'),
           pch=c(rep(NA,length(net_benefit_legend$color)),1),
           lty=c(rep(1,length(net_benefit_legend$color)),0),
           lwd=2,
           col=c(net_benefit_legend$color,1),
           bg='white',
           ncol=legend_ncol,
           cex=0.9)
    
    ## EVPI
    ## Based on code from Joke for Typhoid and McMarcel
  
    param_opt <- names(.rstride$get_variable_model_param(input_param_scen))
    param_opt <- param_opt[!grepl('id',param_opt)]
    num_param <- length(param_opt)
    
    flag_col         <- names(input_param_scen) %in% param_opt 
    country_param    <- input_param_scen[,flag_col]
    evppi <- matrix(NA,num_wtp,num_param)  
    j <- 1; i <- 1
    for(j in 1:num_wtp){
      NB_tmp <- t(net_benefit_all[j,,])
      for(i in 1:num_param){
        evppi[j,i] <- evpi_gam(NB_tmp,country_param[,i])
      }
    }  
    colnames(evppi) <- param_opt
    
    # plot EVPPI
    plot_evppi(evppi,wtp_opt)
    grid(col=alpha(1,0.5))

  } # end for-loop: burden type
  
  # close pdf stream
  dev.off()
  
  # command line message
  smd_print('COST-EFFECTIVENESS ANALYSIS COMPLETE')

}


# intervention_ref_name = 'vaccine_rate'; intervention_ref_value = 0.7
get_average_burden_averted <- function(project_dir,
                                       intervention_ref_name  = 'vaccine_rate',
                                       intervention_ref_value = 0){
  
  # load project summary
  project_summary    <- .rstride$load_project_summary(project_dir)
  
  # retrieve variable model parameters
  exp_design_opt         <- .rstride$get_variable_model_param(project_summary)
  exp_design_col_names   <- names(exp_design_opt)
  vaccine_col_names      <- exp_design_col_names[grepl('vaccine',exp_design_col_names)]
  exp_disease_col_names  <- exp_design_col_names[!exp_design_col_names %in% vaccine_col_names]
  
  # bugfix, if no disease parameters are varied
  if(length(exp_disease_col_names)==0){
    exp_disease_col_names <- 'r0'
  }
  
  # create burden tag 
  if(length(exp_disease_col_names)==1){
    values_format <- formatC(project_summary[,exp_disease_col_names],width = 2, flag ='0')
    project_summary$burden_exp_tag    <- paste(exp_disease_col_names,values_format,sep='_')
  } else{
    project_summary$burden_exp_tag    <- apply(project_summary[,exp_disease_col_names],1,paste,collapse='_') 
  }
  project_summary$burden_exp_tag      <- as.factor(project_summary$burden_exp_tag)
  
  # create experiment tag
  project_summary$experiment_tag      <- apply(project_summary[,exp_design_col_names],1,paste,collapse='_')
  project_summary$experiment_tag      <- as.factor(project_summary$experiment_tag)
  outbreak_total_runs <- as.data.frame(table(experiment_tag = project_summary$experiment_tag),
                                       responseName = 'num_runs')
  
  # add number of runs per experiment tag
  project_summary <- merge(project_summary,outbreak_total_runs)
  
  # create intervention tag
  project_summary[1,vaccine_col_names]
  project_summary$intervention_value  <- project_summary[,intervention_ref_name]
  project_summary$intervention_tag    <- paste0(project_summary$vaccine_profile,
                                                ': cov. ',project_summary$vaccine_rate*100,
                                                '% [',project_summary$vaccine_min_age,
                                                '-',project_summary$vaccine_max_age,']y')
  #project_summary$intervention_tag    <- apply(project_summary[,vaccine_col_names],1,paste,collapse='_') 
  
  # add 'current' strategy and create intervention tags
  flag_c                                   <- project_summary$intervention_value == intervention_ref_value
  project_summary$intervention_tag[flag_c] <- 'current'
  
  
  # get unique levels, and sort (with "current" first)
  intervention_levels                      <- unique(project_summary$intervention_tag)
  intervention_levels                      <- c('current',sort(intervention_levels[intervention_levels != 'current']))
  
  # create factor variable
  project_summary$intervention_tag         <- factor(project_summary$intervention_tag, levels = intervention_levels)
  
    # add intervention data => use "vaccination data"
  datafile_intervention <- file.path(project_dir,paste0(project_summary$run_tag[1],'_data_vaccination.RData'))
  if(file.exists(datafile_intervention)){
    data_intervention <- load(datafile_intervention)
    data_intervention <- get(data_intervention)
  } else{
    data_intervention <- data.frame(exp_id=0)
  }
  
  # get vaccines per scenario
  # outbreak_total_vaccines        <- data.frame(table(data_intervention$exp_id))
  # names(outbreak_total_vaccines) <- c('exp_id','num_vaccines_mean')
  data_intervention$num_vaccines_mean <- 1
  outbreak_total_vaccines <- aggregate(num_vaccines_mean ~ exp_id + pool_type, data = data_intervention, sum )

  # merge project and vaccine data
  project_summary <- merge(project_summary,outbreak_total_vaccines, all.x=T)
  project_summary$num_vaccines_mean[is.na(project_summary$num_vaccines_mean)] <- 0
  
  # select CE related output parameters
  ce_burden_param    <- c(exp_design_col_names,'intervention_tag','num_infected_seeds','burden_exp_tag')
  ce_burden_formula  <- as.formula(paste('num_cases ~ ',paste(ce_burden_param,collapse=' + ')))
  ce_burden_vaccines <- as.formula(paste('num_vaccines_mean ~ ',paste(ce_burden_param,collapse=' + ')))
 
  # calculate average results
  burden_mean        <- aggregate(ce_burden_formula, data=project_summary, mean, na.rm=TRUE)
  burden_lowerlimit  <- aggregate(ce_burden_formula, data=project_summary, quantile, 0.025, na.rm=TRUE)
  burden_upperlimit  <- aggregate(ce_burden_formula, data=project_summary, quantile, 0.975, na.rm=TRUE)
  burden_mean$num_cases_mean  <- burden_mean$num_cases
  burden_mean$num_cases_q2p5  <- burden_lowerlimit$num_cases
  burden_mean$num_cases_q97p5 <- burden_upperlimit$num_cases
  burden_mean$num_cases       <- NULL
  
  # add number of vaccine
  vaccines_mean <- aggregate(ce_burden_vaccines,data=project_summary, mean, na.rm=TRUE)
  burden_mean   <- merge(burden_mean,vaccines_mean)
  
  # fix factors
  burden_mean$intervention_tag  <- factor(burden_mean$intervention_tag,levels=levels(project_summary$intervention_tag))
  burden_mean$burden_exp_tag    <- factor(burden_mean$burden_exp_tag,levels=levels(project_summary$burden_exp_tag))
  
  # set disease-specific reference
  burden_mean$row_id       <- 1:nrow(burden_mean)
  burden_mean$reference_id <- NA
  # exp_disease_opt <- unique(burden_mean[,exp_disease_col_names])
  exp_disease_opt <- unique(burden_mean$burden_exp_tag)
  flag_reference  <- burden_mean$intervention_tag == 'current'
  
  for(i_opt in 1:length(exp_disease_opt)){
    # flag_exp_disease_opt <- burden_mean[,exp_disease_col_names] == exp_disease_opt[i_opt] 
    flag_exp_disease_opt <- burden_mean$burden_exp_tag == exp_disease_opt[i_opt] 
    reference_id         <- which(flag_exp_disease_opt & flag_reference)
    burden_mean$reference_id[flag_exp_disease_opt] <- reference_id
  }
  
  # averted burden
  burden_mean$incr_mean_num_cases_averted  <- burden_mean$num_cases_mean[burden_mean$reference_id] - burden_mean$num_cases_mean
  burden_mean$incr_mean_num_vaccines_admin <- burden_mean$num_vaccines_mean - burden_mean$num_vaccines_mean[burden_mean$reference_id] 
  
  # return burden estimates
  return(burden_mean)
}

##########################
# HELP FUNCTIONS         #
##########################

# standard error
std_err <- function(x) sqrt(var(x)/length(x))

# calculate the net benefit, given a price per DALY 'p'
get_net_benefit <- function(scen_daly_averted,scen_incr_cost_disc,p)
{
  NB <-   p*scen_daly_averted - scen_incr_cost_disc
  return(NB)
}


# NB <- NB_tmp; input_parameter_values <- country_param[,i]
evpi_gam <- function(NB,input_parameter_values){
  
  # note: current should have NB == 0
  D_opt  <- which(!colSums(NB) == 0)
  D <- ncol(NB)
  N <- nrow(NB)
  
  g.hat_new <- matrix(0,nrow=N,ncol=D)
  for(d in D_opt)
  {
    model <- gam(NB[,d] ~ input_parameter_values)
    g.hat_new[,d] <- model$fitted
  }   
  
  perfect.info  <- mean(apply(g.hat_new,1,max))
  baseline      <- max(colSums(g.hat_new)/N)
  
  partial.evpi  <- round(perfect.info - baseline, digits=4) ## estimate EVPPI 
  
  return(partial.evpi)
}

get_factor_legend <- function(factor_codes){
  
  factor_legend  <- data.frame(name=levels(factor_codes),
                               color=seq(1,nlevels(factor_codes))+1,
                               stringsAsFactors = F)
  
  # # manual corrections (to make plots uniform)
  # factor_legend$color[factor_legend$name == 'maternal'] <- 3
  # factor_legend$color[factor_legend$name == 'current']  <- 4
  
  return(factor_legend)
}


plot_evppi <- function(evppi,wtp_opt,plot_main=''){
  
  # modify evpi colnames for figure legend 
  evpi_legend <- colnames(evppi)
  # evpi_legend <- gsub('_',' ',evpi_legend)
  # evpi_legend <- gsub('hosp','hosp.',evpi_legend)
  # evpi_legend <- gsub('prob','probability',evpi_legend)
  # evpi_legend <- gsub('rsv','RSV',evpi_legend)
  # evpi_legend <- gsub('comm','comm.',evpi_legend)
  evpi_legend
  
  # rescale or transform the EVPPI
  #evppi_factor <- 1
  #evppi_plot   <- log(evppi+1) / evppi_factor
  evppi_plot   <- evppi
  
  ylim_plot <- range(c(0,evppi*1.15),na.rm=T)
  #ylim_plot <- range(c(0,8e6))/evppi_factor
  #ylim_plot <- c(0,8e6)
  
  plot(range(wtp_opt),range(evppi_plot),col=0,type='l',lwd=2,ylim=ylim_plot,
       #xlab='Willingness to pay for a DALY averted (USD)',ylab=paste0('log EVPPI (USD)'),
       xlab='Willingness to pay for a QALY (EURO)',ylab=paste0('EVPPI (EURO)'),
       main=plot_main,
       xaxt='n')
  axis(1,pretty(wtp_opt),pretty(wtp_opt)/1e3)
  
  for(i in 1:ncol(evppi_plot)){
    lines(wtp_opt,evppi_plot[,i],col=i,lwd=2)
  }
  
  legend('top',evpi_legend,
         col=(1:ncol(evppi_plot)),
         lty=1,
         lwd=2,
         pch=NA,
         ncol=3,
         cex=0.9)
}


