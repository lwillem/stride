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
# MODEL INCIDENCE EXPLORATION
#
############################################################################# #


#' @param project_dir   name of the project folder
inspect_incidence_data <- function(project_dir, num_selection = 4, bool_add_param=TRUE)
{
  # command line message
  smd_print('INSPECT INCIDENCE DATA...')
  
  #debug
  if(!exists('bool_add_param')) {bool_add_param = TRUE}
  
  # load project summary
  project_summary    <- .rstride$load_project_summary(project_dir)
  
  # reduce population filenames
  project_summary$population_file <- gsub('_belgium','',project_summary$population_file)
  project_summary$population_file <- gsub('_c500_teachers_censushh','',project_summary$population_file)
  project_summary$population_file <- gsub('.csv','',project_summary$population_file)

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
  
  # add config_id 
  project_summary$config_id <- .rstride$get_config_id(project_summary)
  
  # add contact id
  project_summary$contact_id <- .rstride$get_contact_id(project_summary)
  
  # add config_id and tracing_id to incidence data
  data_incidence_all         <- merge(data_incidence_all,project_summary[,c('exp_id','config_id','contact_id')] )

  # remove rows missing sim_date
  data_incidence_all <- data_incidence_all[!is.na(data_incidence_all$sim_date),]
  
  # check for NA's, and replace by 0
  data_incidence_all[is.na(data_incidence_all)] <- 0
  
  ## REFERENCE DATA COVID-19: new hospital admissions ----
  # use (local version of) most recent SCIENSANO data (or backup version)
  #ref_data          <- get_observed_incidence_data()
  ref_data          <- get_hospital_incidence_age()
  ref_data$sim_date <- as.Date(ref_data$sim_date)
  dim(ref_data)
  
  # reformat reference data
  hosp_adm_data <- data.frame(date = ref_data$sim_date,
                              num_adm = ref_data$hospital_admissions,
                              cum_adm = ref_data$cumulative_hospital_admissions,
                              ref_data[,grepl('admissions_',names(ref_data))])
  
  # remove reference data if simulation period is shorter
  flag_compare  <- hosp_adm_data$date %in% data_incidence_all$sim_date
  hosp_adm_data <- hosp_adm_data[flag_compare,]
  
  ## SEROPREVALENCE DATA
  prevalence_ref <- load_observed_seroprevalence_data()
  
  # select simulation period
  sel_ref_dates <- prevalence_ref$seroprevalence_date %in% data_incidence_all$sim_date
  prevalence_ref <- prevalence_ref[sel_ref_dates,]
  
  ## ALL PLOTS  ----
  .rstride$create_pdf(project_dir,'incidence_inspection',width = 6, height = 7)
  par(mfrow=c(4,1))
  
  opt_config_id <- unique(data_incidence_all$config_id)
  i_config <- opt_config_id[1]
  for(i_config in opt_config_id){
    
    # select subset
    data_incidence_sel <- data_incidence_all[data_incidence_all$config_id == i_config,]
    
    # plot
    plot_incidence_data(data_incidence_sel,project_summary,
                        hosp_adm_data,input_opt_design,prevalence_ref,
                        bool_add_param)
  }
  
  # close pdf
  dev.off()
  #--------------------------#
  
  
  ## R0     ####
  # add R0 to input opt design if not present
  if(any(is.null(input_opt_design$r0))){ 
    input_opt_design$r0 <- unique(project_summary$r0)
    
  }
  
  ## PER R0: plot temporal patterns
  input_opt_design$r0 <- round(input_opt_design$r0,digits=1)
  opt_r0 <- unique(input_opt_design$r0)
  if(length(opt_r0)>0){
    .rstride$create_pdf(project_dir,'incidence_R0',width = 6, height = 7)
    par(mfrow=c(4,1))
    
    
    i_r0 <- opt_r0[1]
    for(i_r0 in opt_r0){
      
      # select config_id
      opt_config_id <- unique(input_opt_design$config_id[input_opt_design$r0 ==  i_r0])

      # select subset
      if(is.null(opt_config_id)){
        data_incidence_sel <- data_incidence_all
      } else{
        data_incidence_sel <- data_incidence_all[data_incidence_all$config_id %in% opt_config_id,]
      }
      dim(data_incidence_sel)
      
      # check selection
      if(nrow(data_incidence_sel)>0){
        # plot
        plot_incidence_data(data_incidence_sel,project_summary,
                            hosp_adm_data,input_opt_design,prevalence_ref,
                            bool_add_param)
      }
    }
    
    # close pdf
    dev.off()
  }
  
  ## ALL TOGETHER (PDF) ####
  .rstride$create_pdf(project_dir,'incidence_all',width = 6, height = 2.5)
  par(mar=c(3,5,1,3))
  plot_incidence_data(data_incidence_all,project_summary,
                      hosp_adm_data,input_opt_design,prevalence_ref,
                      bool_add_param,bool_only_hospital_adm = TRUE) 
  dev.off()
  
  # all => polygon
  .rstride$create_pdf(project_dir,'incidence_reproduction',width = 5, height = 5)
  plot_incidence_reproduction(data_incidence = data_incidence_all,
                              hosp_adm_data = hosp_adm_data,
                              scen_color = 1)
  dev.off()
  
  ## ALL TOGETHER (JPEG) ####
  .rstride$create_jpg(project_dir,'incidence_all',width = 6, height = 2.5)
  par(mar=c(3,5,1,5))
  plot_incidence_data(data_incidence_all,project_summary,
                      hosp_adm_data,input_opt_design,prevalence_ref,
                      bool_add_param,bool_only_hospital_adm = TRUE) 
  dev.off()
  
  # all => polygon
  .rstride$create_jpg(project_dir,'incidence_reproduction',width = 5, height = 4)
  plot_incidence_reproduction(data_incidence = data_incidence_all,
                              hosp_adm_data = hosp_adm_data,
                              scen_color = 1)
  dev.off()
  
  ## ALL TOGETHER: NO PARAM ####
  .rstride$create_pdf(project_dir,'incidence_no_param',width = 6, height = 2.5)
  par(mar=c(3,5,1,3))

  # plot
  plot_incidence_data(data_incidence_all,project_summary,
                      hosp_adm_data,input_opt_design,prevalence_ref,
                      bool_add_param = FALSE,
                      bool_only_hospital_adm = FALSE) 
  
  # close pdf
  dev.off()
  
  ## PARETO ENSEMBLE
  filename_summary_score <- file.path(project_dir,paste0(basename(project_dir),'_poison_neg_loglikelihood_scores.RData'))
  if(file.exists(filename_summary_score)){
    
    
    summary_score      <- readRDS(filename_summary_score)
    config_selection   <- summary_score$config_id[summary_score$pareto_front]
    data_incidence_sel <- data_incidence_all[data_incidence_all$config_id %in% config_selection,]
    
    ## PARETO (PDF) ####
    
    # all
    .rstride$create_pdf(project_dir,'incidence_pareto_all',width = 6, height = 7)
    par(mfrow=c(4,1))
    plot_incidence_data(data_incidence_sel,project_summary,
                        hosp_adm_data,input_opt_design,prevalence_ref,
                        bool_add_param,bool_only_hospital_adm = FALSE) 
    dev.off()
    
    # hospital admissions
    .rstride$create_pdf(project_dir,'incidence_pareto_hosp',width = 6, height = 2.5)
    par(mar=c(3,5,1,3))
    plot_incidence_data(data_incidence_sel,project_summary,
                        hosp_adm_data,input_opt_design,prevalence_ref,
                        bool_add_param,bool_only_hospital_adm = TRUE) 
    dev.off()
    
    # polygon
    .rstride$create_pdf(project_dir,'incidence_pareto_reproduction',width = 5, height = 5)
    plot_incidence_reproduction(data_incidence = data_incidence_sel,
                                hosp_adm_data = hosp_adm_data,
                                scen_color = 1)
    dev.off()
    
    
  }
  
  ## AGE-SPECIFIC PLOTS
  
  par(mfrow=c(3,3))
  i_age <- 4
  names(data_incidence_sel)
  col_ind_hosp_age <- which(grepl('hospital_admissions',names(hosp_adm_data)))
  
  for(i_age in 1:9){
    
    data_incidence_sel <- data_incidence_all
    for(i_col in c('new_hospital_admissions','new_infections')){
      data_incidence_sel[,i_col] <- data_incidence_sel[,grepl(paste0(i_col,'_age',i_age),names(data_incidence_sel))]
    }
    
    hosp_adm_data_sel <- hosp_adm_data
    hosp_adm_data_sel$num_adm <- hosp_adm_data_sel[,col_ind_hosp_age[i_age]]
    
    data_incidence_sel$new_hospital_admissions
    
    plot_incidence_data(data_incidence_sel,
                        project_summary,
                        hosp_adm_data_sel,
                        input_opt_design,
                        prevalence_ref,
                        bool_add_param,
                        bool_only_hospital_adm = TRUE) 
  }
  
  
  
  
   # command line message
  smd_print('INSPECTION OF INCIDENCE DATA COMPLETE')
  
} # end function

plot_incidence_data <- function(data_incidence_sel,project_summary,
                                hosp_adm_data,input_opt_design,prevalence_ref,
                                bool_add_param,
                                bool_add_axis4 = TRUE,
                                bool_only_hospital_adm = FALSE,
                                bool_seroprev_limited = FALSE,
                                bool_add_doubling_time = FALSE){

  if(nrow(data_incidence_sel)==0){
    return(NULL)
  }
  
  # change figure margins
  if(!bool_only_hospital_adm){
    par(mar=c(3,5,1,5))
  }
  
  # set color definitions and other layout definitions
  pcolor <- data.frame(E = "black",  # exposed (or total infections)
                       I = "darkgoldenrod3",  # infectious
                       S = "red",  # symptomatic
                       H = 'blue',  # hospitalized
                       D = 'black',
                       alpha = 0.1,
                       lwd = 3,
                       pch = 20 , # if points are used
                       stringsAsFactors = F)  # data
  
  # change transparancey for low number of rng-runs
  if(max(table(project_summary$config_id)) < 10){
    pcolor$alpha <- 0.2
  }
  
  ## FIX FOR PLOTTING: Set all values for the last sim_day to NA
  data_incidence_sel[data_incidence_sel$sim_date %in% max(data_incidence_sel$sim_date,na.rm=T),] <- NA
  
  # set y-lim
  y_lim <- range(0,pretty(max(hosp_adm_data$num_adm,na.rm=T)*1.1),max(data_incidence_sel$new_hospital_admissions,na.rm=T),na.rm=T)

  ## HOSPITAL ADMISSIONS ####
  plot(data_incidence_sel$sim_date,
       data_incidence_sel$new_hospital_admissions,
       type='l',
       col=alpha(pcolor$H,pcolor$alpha),
       ylab='Hospital admissions',
       xlab='',
       ylim = y_lim,
       yaxt='n',
       xaxt='n')
  add_x_axis(data_incidence_sel$sim_date)
  add_y_axis(y_lim)
  points(hosp_adm_data$date,hosp_adm_data$num_adm,col=pcolor$D,pch=pcolor$pch)
  add_breakpoints()
  add_legend_hosp(pcolor)
  
  # add config tag
  adm_mean_final <- aggregate(new_hospital_admissions ~ sim_date + config_id, data=data_incidence_sel,mean)
  adm_mean_final <- adm_mean_final[adm_mean_final$sim_date == max(adm_mean_final$sim_date),]
  # text(adm_mean_final$sim_date,adm_mean_final$new_hospital_admissions,adm_mean_final$config_id,
  #      pos=4,xpd=TRUE,cex=0.4)
  if(bool_add_axis4){
    if(bool_only_hospital_adm){
      axis(4,adm_mean_final$new_hospital_admissions,adm_mean_final$contact_id,las=2,cex.axis=0.6) # add contact details
    } else{
      axis(4,adm_mean_final$new_hospital_admissions,adm_mean_final$config_id,las=2,cex.axis=0.4) # add full config id
    }  
  }
  
  
  if(bool_only_hospital_adm){ return() } # stop
  
  ## CUMULATIVE: HOSPITAL ####
  y_lim <- range(0,pretty(max(hosp_adm_data$cum_adm,na.rm=T)*2,data_incidence_sel$cumulative_hospital_cases),na.rm=T)
  # plot(data_incidence_sel$sim_date,
  #      data_incidence_sel$cumulative_hospital_cases,
  #      type='l',
  #      col=alpha(pcolor$H,pcolor$alpha),
  #      ylab='Total admissions',
  #      xlab='',
  #      yaxt='n',
  #      xaxt='n',
  #      ylim= range(y_lim))
  colnames_hosp_age <- names(data_incidence_sel)[grepl('cumulative_hospital_cases_age',names(data_incidence_sel))]
  cum_hosp_age     <- aggregate(. ~ sim_date, data = data_incidence_sel[,c('sim_date',colnames_hosp_age)],mean)
  bplot <- barplot(t(as.matrix(cum_hosp_age[,-1])),
                   col=rainbow(9),
                   ylab='Total admissions',
                   xlab='',
                   yaxt='n',
                   xaxt='n',
                   ylim= range(y_lim))
  
  # x axis
  sim_dates <-data_incidence_sel$sim_date 
  date_lim  <- range(sim_dates,na.rm=T)
  bplot_diff <- diff(range(bplot)) / as.numeric(diff(date_lim))
  x_tick_label <- pretty(sim_dates,5)
  x_tick_value <- as.numeric((x_tick_label - min(date_lim)))  * bplot_diff
  axis(1,x_tick_value, format(x_tick_label,'%e %b'),cex.axis=0.9)
  add_y_axis(y_lim)
  add_y_axis_pop(y_lim)
  points(as.numeric((hosp_adm_data$date - min(date_lim)))  * bplot_diff,
         hosp_adm_data$cum_adm,col=pcolor$D,pch=pcolor$pch)
  legend('topleft',
         c('Reported'),
         col=c(pcolor$D),
         pch=c(16),
         lwd=c(NA),
         bg='white',
         cex = 0.6,
         ncol=1)
  # set date format (character vs numeric)
  
  legend('left',
         rev(paste0(seq(0,80,10),'-',
               seq(9,90,10))),
         col=rev(rainbow(9)),
         ncol=1,
         lwd=2,
         title='Age group',
         cex=0.5)
         
  # legend('left',
  #        rev(c('0-18',
  #          '19-59',
  #          '60-79',
  #          '+80')),
  #        col=c(8:5),
  #        lwd=2,
  #        title='Age group',
  #        cex=0.5
  # )
  
  ## INCIDENCE: ALL ####
  y_lim <- range(0,pretty(data_incidence_sel$new_infections),na.rm = T)
  plot(data_incidence_sel$sim_date,
       data_incidence_sel$new_infections,
       type='l',
       col=alpha(pcolor$E,pcolor$alpha),
       ylab='New cases',
       xlab='',
       ylim = y_lim,
       yaxt='n',
       xaxt='n')
  add_x_axis(data_incidence_sel$sim_date)
  add_y_axis(y_lim)
  lines(data_incidence_sel$sim_date,
        data_incidence_sel$new_infectious_cases,
        col=alpha(pcolor$I,pcolor$alpha))
  
  lines(data_incidence_sel$sim_date,
        data_incidence_sel$new_symptomatic_cases,
        col=alpha(pcolor$S,pcolor$alpha))
  
  lines(data_incidence_sel$sim_date,
        data_incidence_sel$new_hospital_admissions,
        col=alpha(pcolor$H,pcolor$alpha))
  # points(hosp_adm_data$date,hosp_adm_data$num_adm,col=pcolor$D,pch=pcolor$pch)
  add_breakpoints()
  add_legend_incidence(pcolor)
  
  if(bool_add_doubling_time){
    flag_dates <- !is.na(data_incidence_sel$doubling_time)
    y_tick    <- quantile(y_lim,0.7)
    x_centre  <- mean(data_incidence_sel$sim_date[flag_dates])
    #doubling_time <- round(mean(data_incidence_sel$doubling_time,na.rm=T),digits = 2)
    doubling_time <- range(round(data_incidence_sel$doubling_time,digits = 2),na.rm=T)
    doubling_time <- unique(doubling_time) # remove duplicates if only one value is present
    
    lines(data_incidence_sel$sim_date[flag_dates],
          rep(y_tick,sum(flag_dates)),
          lty=3)
    text(x_centre,
         y_tick,
         paste('DT:',paste0(doubling_time,collapse='-')),
         pos=3,
         cex=0.8)
  }
  
  
  ## CUMULATIVE: ALL STATES ####
  y_lim <- range(pretty(data_incidence_sel$cumulative_infections))
  plot(data_incidence_sel$sim_date,
       data_incidence_sel$cumulative_infections,
       type='l',
       col=alpha(pcolor$E,pcolor$alpha),
       ylab='Total cases',
       xlab='',
       ylim=y_lim,
       yaxt='n',
       xaxt='n')
  add_x_axis(data_incidence_sel$sim_date)
  add_y_axis(y_lim)
  add_y_axis_pop(y_lim)
  lines(data_incidence_sel$sim_date,
        data_incidence_sel$cumulative_infectious_cases,
        col=alpha(pcolor$I,pcolor$alpha))
  lines(data_incidence_sel$sim_date,
        data_incidence_sel$cumulative_symptomatic_cases,
        col=alpha(pcolor$S,pcolor$alpha))
  lines(data_incidence_sel$sim_date,
        data_incidence_sel$cumulative_hospital_cases,
        col=alpha(pcolor$H,pcolor$alpha))
  #points(hosp_adm_data$date,hosp_adm_data$cum_adm,col=pcolor$D,pch=pcolor$pch)
  
  add_breakpoints()
  if(bool_add_param) {
    add_legend_runinfo(project_summary,input_opt_design,
                       unique(data_incidence_sel$config_id))
  } else {
    add_legend_incidence(pcolor,legend_pos = 'topleft',bool_seroprevalence=TRUE)
  }
  

  ## add reference
  points(prevalence_ref$seroprevalence_date,
         prevalence_ref$point_incidence_mean,
         pch=8)
  add_interval(prevalence_ref$seroprevalence_date,
               y1 = prevalence_ref$point_incidence_low,
               y2 = prevalence_ref$point_incidence_high)
  

} # end function to plot figure


# define the vertical breaks on the plots
add_breakpoints <- function(bool_text=TRUE){
  
  # add start intervention
  add_vertical_line("2020-03-14",bool_text,'Restrictions')
  
  # # add today
  # add_vertical_line(Sys.Date())
  
  # add scenario date (exit wave 1)
  add_vertical_line("2020-05-04",bool_text, "B2B")
  
  # add scenario date (exit wave 2)
  add_vertical_line("2020-05-18",bool_text, 'School')
  
  # add scenario date (exit wave 3)
  add_vertical_line("2020-05-25",bool_text, 'Community')
  
  # add scenario date (summer holiday)
  add_vertical_line("2020-07-01",bool_text,'Holiday')
}

# add vertical line on given date + label on x-axis
add_vertical_line <- function(date_string,bool_text,date_tag = ''){
  
  plot_limits <- par("usr")
  
  v_date <- as.Date(date_string)
  abline(v=v_date,lty=2)
  #axis(1,v_date,format(v_date,'%d/%m'),
       #cex.axis=0.5,padj=-3,tck=-0.005)
  
  if(bool_text)
  {
    v_text <- ifelse(nchar(date_tag)>0,date_tag,format(v_date,'%d/%m'))
    text(x = v_date-1,
         y = mean(plot_limits[3:4]),
         #paste(format(v_date,'%d/%m'),date_tag),
         v_text,
         srt=90, pos=3, offset = +1.5,cex=0.6)
  }
}
  
add_interval <- function(x,y1,y2){
  ## add reference
  arrows(x,y1,
         x,y2,
         angle=90,length=0.05)
  arrows(x,y2,
         x,y1,
         angle=90,length=0.05)
}

add_interval_hor <- function(x1,x2,x_mean,...){
  
  # get mean y-axis
  #y_mean <- mean(par("usr")[3:4])
  y_mean <- par("usr")[4]*3/4
  
  points(x=x_mean,y = y_mean,
         pch = 8,
         ...)
  
  ## add reference
  arrows(x1,y_mean,
         x2,y_mean,
         angle=90,length=0.05,...)
  arrows(x2,y_mean,
         x1,y_mean,
         angle=90,length=0.05,...)
}


# define the legend for hospital(-only) plots
add_legend_hosp <- function(pcolor){
  legend('topleft',
         c('Reported',
           'Simulation'),
         col=c(pcolor$D,pcolor$H),
         pch=c(16,NA),
         lwd=c(NA,2),
         bg='white',
         cex = 0.6,
         ncol=1)
}

# define the legend with all categories
add_legend_incidence <- function(pcolor,legend_pos = 'topright',bool_seroprevalence=FALSE){
  
  legend_text <- c('Exposed (latent)',
                   'Infectious',
                   'Symptomatic',
                   'Hospitalized')
  legend_col <- unlist(pcolor)
  legend_pch <- rep(NA,length(legend_text))
  legend_lwd <- rep(2,length(legend_text))
  
  if(bool_seroprevalence){
    legend_text <- c(legend_text,'Seroprevalence')
    legend_col  <- c(legend_col,1)
    legend_pch  <- c(legend_pch,124)
    legend_lwd  <- c(legend_lwd,NA)
  }
  
  legend(legend_pos,
         legend_text,
         col=legend_col,
         pch=legend_pch,
         lwd=legend_lwd,
         cex=0.6,
         bg='white')
}

add_legend_runinfo <- function(project_summary,input_opt_design,
                               i_config){
  
  # remove NA from i_config
  i_config <- i_config[!is.na(i_config)]
  
  # subset project summary
  project_summary_sel <- project_summary[project_summary$config_id %in% i_config,]
  
  # get run info: population
  run_info <- c(paste0('pop_size = ',unique(project_summary_sel$population_size)/1e3,'k'))

  # add other exp_design parameters (excl. id and non-numeric parameters)
  names_exp_param <- colnames(input_opt_design)
  bool_no_id      <- !names_exp_param %in% c('config_id','contact_id','tracing_id')
  bool_numeric    <- unlist(lapply(input_opt_design[1,],is.numeric))
  names_exp_param <- names_exp_param[bool_no_id & bool_numeric]
  project_summary_sel[,names_exp_param] <- round(project_summary_sel[,names_exp_param],digits = 2)
  
  i_name <- names_exp_param[2]
  if(length(names_exp_param)>0){
    for(i_name in names_exp_param){
    
      if(length(unique(project_summary_sel[,i_name]))<4){
        run_info <- c(run_info,
                      paste0(i_name, ' = ',paste(unique(project_summary_sel[,i_name]),collapse=', '))
        )
      } else{
        run_info <- c(run_info,
                      paste0(i_name, ' = ',paste(range(project_summary_sel[,i_name]),collapse='-'))
        )      
      }

      
      
    }
    
  }
  
  # present info in legend
  legend('topleft',
         run_info,
         cex=0.6,
         bg='white')
}


plot_distancing <- function(project_summary){
  
  opt_distancing <- unique(data.frame(compliance_delay = project_summary$compliance_delay,
                                 cnt_reduction_work = project_summary$cnt_reduction_work,
                                 cnt_reduction_other = project_summary$cnt_reduction_other))


  date_start <- as.Date('2020-03-14')
  
  i_distancing <- 1    
  
  opt_distancing[i_distancing,]
  num_days <- opt_distancing$compliance_delay[i_distancing]
  date_all <- date_start + (1:num_days)-1
  compliance_factor <- 1:num_days / num_days
  
  plot(date_all,compliance_factor*opt_distancing$cnt_reduction_work[i_distancing],
       type='l',
       col=4,
       ylim=0:1)
  lines(date_all,compliance_factor*opt_distancing$cnt_reduction_other[i_distancing])
        
    
  
}

#sim_dates <- data_incidence_scenario$sim_date
add_x_axis <- function(sim_dates,bool_numeric=FALSE,num_ticks = 7,bool_grid = TRUE,las=1){
  
  # set date format (character vs numeric)
  date_format <- ifelse(bool_numeric,'%e/%m','%e %b')
  
  x_ticks <- pretty(sim_dates,num_ticks-2)
  axis(1,x_ticks, format(x_ticks,date_format),cex.axis=0.9,las=las)
  if(bool_grid) { abline(v=pretty(sim_dates,num_ticks),lty=3,col='lightgray') }
  #grid(nx=NA,ny=NULL,lty=3,col='lightgray')
}

add_y_axis <- function(y_lim,bool_grid = TRUE){
  
  y_ticks <- pretty(y_lim,5)
  y_labels <- y_ticks # default
  if(max(y_ticks)>1e4){ y_labels <- paste0(round(y_ticks/1e3),'k')}
  if(max(y_ticks)>1e6){ y_labels <- paste0(round(y_ticks/1e6,digits=2),'M')}
  
  axis(2,y_ticks, y_labels,cex.axis=0.9,las=2)
  if(bool_grid) {abline(h=pretty(y_lim,5),lty=3,col='lightgray')}
}

add_y_axis_pop <- function(y_lim,pop_size = 11e6){
  
  # add axis
  axis(4,pretty(y_lim),paste0(round(pretty(y_lim)/pop_size*100,digits=1),'%'),las=2,cex.axis=0.9)
  mtext('Belgian population (%)',side = 4,line=3,cex=0.7)
}



## HELP FUNCTION TO REFORMAT THE PREVALENCE DATA FROM THE STOCHASTIC MODEL
reformat_prevalence_stochastic_model <- function(){
  
  ref_prevalence <- read.table('./data/prevalence_stochastic_model_20200604.csv',sep=';')
  dim(ref_prevalence)  

  date_steps <- seq(as.Date('2020-03-01'),as.Date('2020-06-01'),1/24)
  num_steps  <- length(date_steps)
  plot(date_steps[-1],colMeans(ref_prevalence[,2:num_steps]),col=2,
       ylim=range(ref_prevalence[,-1]),type='l',xaxt='n',
       xlab='',
       ylab='prevalence stochastic model')
  lines(date_steps[-1],apply(ref_prevalence[,2:num_steps],2,min),col=3)
  lines(date_steps[-1],apply(ref_prevalence[,2:num_steps],2,max),col=3)
  add_breakpoints()
  add_x_axis(date_steps)
  grid(nx=NA,ny=NULL)
  
  date_out <- seq(as.Date('2020-03-01'),as.Date('2020-06-01'),1)
  prevalence_out <- data.frame(mean = approx(date_steps[-num_steps],colMeans(ref_prevalence[,2:num_steps]),date_out)$y,
                               min = approx(date_steps[-num_steps],apply(ref_prevalence[,2:num_steps],2,min),date_out)$y,
                               max = approx(date_steps[-num_steps],apply(ref_prevalence[,2:num_steps],2,max),date_out)$y,
                              date = date_out)
  
  approx(date_steps[-num_steps],colMeans(ref_prevalence[,2:num_steps]),seq(as.Date('2020-04-01'),as.Date('2020-06-01'),15))

  head(prevalence_out)
  saveRDS(prevalence_out,file='./sim_output/prevalence_stochastic_model_20200604.rds')
  
  

}


add_polygon_incidence <- function(data_incidence,colname_burden, scen_color){
  
  if(colname_burden %in% names(data_incidence)){
    hosp_min       <- aggregate(formula(paste(colname_burden,'~ sim_date')), data=  data_incidence, min)
    hosp_max       <- aggregate(formula(paste(colname_burden,'~ sim_date')), data=  data_incidence, max)
    
    hosp_median    <- aggregate(formula(paste(colname_burden,'~ sim_date')), data=  data_incidence, median,na.rm=T)
    hosp_mean      <- aggregate(formula(paste(colname_burden,'~ sim_date')), data=  data_incidence, mean,na.rm=T)
    
    newx           <- c(hosp_min$sim_date,rev(hosp_max$sim_date))
    newy           <- c(hosp_min[,2],rev(hosp_max[,2]))
    
    polygon(newx, newy, col = scen_color, border = scen_color)
  }
    #lines(hosp_median$sim_date,hosp_median$new_hospital_admissions,type='l',lwd=2,col=scen_color)
  #lines(hosp_mean$sim_date,hosp_mean$new_hospital_admissions,type='l',lwd=2,col=scen_color,lty=3)
}

plot_incidence_reproduction <- function(data_incidence,hosp_adm_data,scen_color,plot_main='')
{
  y_lim  <- c(0,700)
  x_lim  <- range(data_incidence$sim_date,na.rm = T)
  
  scen_color <- alpha(scen_color,0.4)
  
  par(fig=c(0,1,0.35,1),mar=c(0,5,2,1))
  plot(hosp_adm_data$date,hosp_adm_data$num_adm,pch=20,
       xaxt='n', yaxt='n',
       xlab = '',ylab='Hospital admissions',
       ylim = y_lim,
       xlim = x_lim,
       main = plot_main)
  sum_scen1 <- add_polygon_incidence(data_incidence,'new_hospital_admissions',scen_color)
  add_y_axis(y_lim)
  add_breakpoints()
  
  par(fig=c(0,1,0,0.34),mar=c(3,5,0,1), new=TRUE)
  plot(0,0,pch=20,
       xaxt='n',yaxt='n',
       xlab = '',ylab='Re',
       ylim = c(0,4.1),
       xlim = x_lim,
       main='')
  add_polygon_incidence(data_incidence,'sec_cases',scen_color)
  add_x_axis(data_incidence$sim_date)
  add_y_axis(0:4)
  add_breakpoints(bool_text=FALSE)
  
}




