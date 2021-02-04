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


############################################################################## #
# EXPLORE TRANSMISSION EVENTS, GENERATION INTERVAL AND REPRODUCTION NUMBER  ####
############################################################################## #
inspect_transmission_dynamics <- function(project_dir,save_pdf = TRUE)
{
  
  # command line message
  smd_print('INSPECT TRANSMISSION DYNAMICS...')
  
  #debug
  if(!exists('save_pdf')) {save_pdf = FALSE}

  # load project summary
  project_summary      <- .rstride$load_project_summary(project_dir)
  
  # retrieve all variable model parameters
  input_opt_design     <- .rstride$get_variable_model_param(project_summary)
  
  # get all transmission output
  data_incidence_all      <- .rstride$load_aggregated_output(project_dir,'data_incidence')
  
  if(length(data_incidence_all) == 1 && is.na(data_incidence_all)){
    smd_print('NO TRANSMISSION DATA AVAILABLE.')
    return(.rstride$no_return_value())
  }
  
  if(!'sec_cases' %in% names(data_incidence_all)){
    smd_print('USE OTHER "event_log_level" TO ANALYSE TRANSMISSION DATA.')
    return(.rstride$no_return_value())
  }
  
  # open pdf stream
  if(save_pdf) .rstride$create_pdf(project_dir,'transmission_inspection',10,7)

  i_config <- 1
  for(i_config in 1:nrow(input_opt_design)){
    
    # reset figure arrangements... and start new plot
    if(save_pdf) par(mfrow=c(3,3))
    
    # subset transmission output corresponding the 'input_opt_design' row
    flag_exp            <- .rstride$get_equal_rows(project_summary,input_opt_design[i_config,])
    data_incidence      <- data_incidence_all[data_incidence_all$exp_id %in% project_summary$exp_id[flag_exp],]
    num_runs_exp        <- sum(flag_exp)
    num_infected_seeds  <- data_incidence$new_infections[1]
  
    # if no incidence data avaiable for this configuration, go to next iteration
    if(nrow(data_incidence)==0){
      next
    }
    
    ## DOUBLING TIME ----
    plot(data_incidence$sim_date,data_incidence$doubling_time,ylim=c(0,14),
         xlab='Date',
         ylab='Doubling time (infections)',
         main='Doubling time (infections)',
         type='l',
         col=alpha(1,0.5))
    add_breakpoints()
    abline(h=3.1,col=4)
    text(max(data_incidence$sim_date),3.1,'3.1',pos=3)
    polygon(x=c(data_incidence$sim_date,rev(data_incidence$sim_date)),
            y = c(rep(2.4,nrow(data_incidence)),rep(4.4,nrow(data_incidence))),
            col=alpha(4,0.1),border = NA)
    
    ## REPRODUCTION NUMBER ----
    data_incidence$sec_cases[data_incidence$sim_date == min(data_incidence$sim_date,na.rm=T)] <- NA
    plot_ymax <- range(c(0,4,data_incidence$sec_cases),na.rm = T)
    plot(data_incidence$sim_date,data_incidence$sec_cases,
         type='l',
         col=alpha(1,0.5),
         xlab='day of infection',ylab='secondary infections',
         ylim=plot_ymax,
         main='reproduction number',
         xaxt='n')
    abline(h=1,lty=3)
    add_x_axis(range(data_incidence$sim_date,na.rm = T))
    add_breakpoints()
    abline(h=3.1,col=4)
    
    text(max(data_incidence$sim_date),3.28,'3.28',pos=3)
    polygon(x=c(data_incidence$sim_date,rev(data_incidence$sim_date)),
            y = c(rep(3.28-1.16,nrow(data_incidence)),rep(3.28+1.16,nrow(data_incidence))),
            col=alpha(4,0.1),border = NA)
    
    if('r0' %in% names(input_opt_design[i_config,])){
      points(min(data_incidence$sim_date),input_opt_design[i_config,'r0'],col=4,pch=9)
      text(min(data_incidence$sim_date),input_opt_design[i_config,'r0'],
           round(input_opt_design[i_config,'r0'],digits=2),pos=3,cex=0.7)
    }
    
    
    
    ## GENERATION INTERVAL   ----
    # note: the generation interval is the time between the infection time of an infected person and the infection time of his or her infector.
    # reference: Kenah et al (2007)
    plot(data_incidence$sim_date,data_incidence$gen_interval,type='l',
         xlab='day of infection',ylab='generation interval [infection]',
         main='generation interval\n[infection]',
         xaxt='n',
         col=alpha(1,0.5))
    add_x_axis(range(data_incidence$sim_date,na.rm=T))
    text(max(data_incidence$sim_date,na.rm=T),5.2,'5.2',pos=3)
    abline(h=5.2,col=4)
    polygon(x=c(data_incidence$sim_date,rev(data_incidence$sim_date)),
            y = c(rep(3.78,nrow(data_incidence)),rep(6.78,nrow(data_incidence))),
            col=alpha(4,0.1),border = NA)
    
    # AGE: INFECTIONS     ----
    age_breaks <- c(0,18,59,79,110)
    age_labels <- levels(cut(0:max(age_breaks),age_breaks,right=F,include.lowest = T))
    
    new_infections_age <- colSums(data_incidence[,paste0('new_infections_age',1:length(age_labels))])
    names(new_infections_age) <- age_labels
    bplot <- barplot(new_infections_age / sum(new_infections_age),
                    xlab='Age category (years)',
                    ylim = c(0,1),
                    ylab= 'Relative incidence',
                    main = 'Total incidence: all')
    grid()
    text(bplot, new_infections_age / sum(new_infections_age),
         round(new_infections_age / sum(new_infections_age),digits=2),
         pos=3)
    
    # AGE: SYMPTOMATIC     ----
    data_incidence$new_symptomatic_cases_age1
    new_symptomatic_cases_age <- colSums(data_incidence[,paste0('new_symptomatic_cases_age',1:length(age_labels))],na.rm = T)
    names(new_symptomatic_cases_age) <- age_labels
    bplot <- barplot(new_symptomatic_cases_age / sum(new_symptomatic_cases_age),
            xlab='Age category (years)',
            ylim = c(0,1),
            ylab= 'Relative incidence symptomatic cases',
            main = 'Total incidence: symptomatic')
    grid()
    text(bplot, new_symptomatic_cases_age / sum(new_symptomatic_cases_age),
         round(new_symptomatic_cases_age / sum(new_symptomatic_cases_age),digits=2),
         pos=3)
    
    # AGE: HOSPITAL ADMISSIONS     ----
    data_incidence$new_hospital_admissions_age1
    new_hospital_admissions_age_all <- data_incidence[,paste0('new_hospital_admissions_age',1:length(age_labels))]
    names(new_hospital_admissions_age_all) <- age_labels
    new_hospital_admissions_age <- colSums(new_hospital_admissions_age_all,na.rm = T)
    bplot <- barplot(new_hospital_admissions_age / sum(new_hospital_admissions_age),
            xlab='Age category (years)',
            ylim = c(0,1),
            ylab= 'Proportion hospital admissions',
            main = 'Total hospital admissions')
    grid()
    text(bplot, new_hospital_admissions_age / sum(new_hospital_admissions_age),
         round(new_hospital_admissions_age / sum(new_hospital_admissions_age),digits=2),
         pos=3)
    
    ## RELATIVE HOSPITAL ADMISSIONS OVER TIME BY AGE
    data_incidence[,paste0('relative_hospital_admissions_age',1:length(age_labels))] <- data_incidence[,paste0('new_hospital_admissions_age',1:length(age_labels))] / data_incidence[,paste0('new_hospital_admissions')]
   if(any(!is.na(data_incidence$relative_hospital_admissions_age1))){
     plot(aggregate(relative_hospital_admissions_age1 ~ sim_date, data = data_incidence,mean),
          ylim=0:1,col=1,lwd=2,type='l',
          xlab='Time',ylab='Proportion hospital admissions')
     lines(aggregate(relative_hospital_admissions_age2 ~ sim_date, data = data_incidence,mean),col=2,lwd=2)
     lines(aggregate(relative_hospital_admissions_age3 ~ sim_date, data = data_incidence,mean),col=3,lwd=2)
     lines(aggregate(relative_hospital_admissions_age4 ~ sim_date, data = data_incidence,mean),col=4,lwd=2)
     legend('top',
            age_labels,
            col = 1:length(age_labels),
            lwd=2,
            ncol=length(age_labels),
            title='Age group (years)',
            cex=0.5)
   }
   
    
    ## LOCATION       ----
    col_location <- names(data_incidence)[grepl('location_',names(data_incidence))]
    loc_names <- col_location
    loc_names <- gsub('location_','',loc_names)
    loc_names <- gsub('Household','HH',loc_names)
    loc_names <- gsub('PrimaryCommunity','Wknd com',loc_names)
    loc_names <- gsub('SecondaryCommunity','Week com',loc_names)
    loc_names <- gsub('K12School','School',loc_names)

    # overall
    summary_location <- colSums(data_incidence[,col_location],na.rm=T) / sum(data_incidence[,col_location],na.rm=T)
    names(summary_location) <- loc_names
    bplot <- barplot(summary_location,las=2,
            ylim=0:1,
            ylab='Relavive incidence')  
    grid()
    text(bplot, summary_location ,
         round(summary_location,digits=2),
         pos=3)
    
    # over time
    y_lim <- range(0,data_incidence[,col_location]*1.1,na.rm=T)
    plot(aggregate(formula(paste(col_location[1],'~ sim_date')), data = data_incidence,mean),
         ylim=y_lim,col=0,
         ylab='Incidence (infections)',
         xlab='Time')
    for(i_col in 1:length(col_location))
    lines(aggregate(formula(paste(col_location[i_col],'~ sim_date')), data = data_incidence,mean),col=i_col,lwd=2)
    
    legend('topright',
           loc_names,
           col = 1:length(loc_names),
           lwd=2,
           ncol=2,
           title='Location',
           cex=0.5)
    
  } # end for-loop to vary the input_opt_design
  
  # close PDF stream
  if(save_pdf) dev.off()
  
  # command line message
  smd_print('INSPECTION OF TRANSMISSION DATA COMPLETE')
  
} # function end


## DEFENSIVE PROGRAMMING
inspect_transmission_data <- function(project_dir,save_pdf = TRUE){
  smd_print('inspect_transmission_data() is depricated... please use inspect_transmission_dynamics',WARNING = T)
  inspect_transmission_dynamics(project_dir,save_pdf)
}

#data_transm <- rstride_out$data_transmission
# get all transmission output
# data_transm_all      <- .rstride$load_aggregated_output(project_dir,'data_transmission')
# data_transm <- data_transm_all[data_transm_all$exp_id == 1,]
get_transmission_statistics <- function(data_transm,
                                        sim_date_range)
{
  
  # 1. Get main statistics
  summary_out <- get_main_transmission_statistics(data_transm,sim_date_range)
  
  
  # 2. Get additional statistics (if possible)
  if(!'infector_id' %in% names(data_transm)){
    return(summary_out)
  }
  
  # setup data.table
  data_transm[,ID := 1:nrow(data_transm),]

  # add recovery date
  data_transm[, date_recovered := infection_date + end_infectiousness]
  
  if(length(unique(data_transm$exp_id))>1){
    smd_print("TRANSMISSION STATISTICS ERROR: MULTIPLE EXPERIMENTS !!", WARNING = T, FORCED = T)
  }
  
  ## REPRODUCTION NUMBER ----
  # day of infection: case
  infection_time <- data.table(local_id       = data_transm$local_id,
                               infector_id    = data_transm$infector_id,
                               infection_date = data_transm$infection_date,
                               recoverd_date  = data_transm$date_recovered)
  
  # day of infection: infector
  infector_time  <- data.table(infector_id            = data_transm$local_id,
                               infector_infection_date = data_transm$infection_date)
  
  
  # set infector_id for infected seeds to -1
  infection_time$infector_id[is.na(infection_time$infector_id)] <- -1

  # merge case and infector timings
  infection_time <- merge(infection_time,infector_time, all.x = TRUE)
  
  # count secondary cases per local_id
  data_infectors <- data_transm[,.(sec_cases = .N),by=infector_id]
  names(data_infectors)[1] <- 'local_id'

  # merge secondary cases with time of infection
  sec_transm    <- merge(infection_time,data_infectors,by='local_id',all=T)
  sec_transm[is.na(sec_cases),sec_cases := 0] # adjust for '0' secondary cases
  sec_transm = sec_transm[!is.na(local_id),]  # remove the "secondary cases as a result of the infected seeding"
  
  ## GENERATION INTERVAL   ----
  # note: the generation interval is the time between the infection time of an infected person and the infection time of his or her infector.
  # reference: Kenah et al (2007)
  # remove the generation intervals counted form the initial infected seed infections
  sec_transm[, gen_interval := as.numeric(infection_date - infector_infection_date)]

  ## RENAME DATE COLUMN
  data_transm[,sim_date := infection_date]
  
  # AGE CATEGORIES     ----
  age_breaks               <- c(seq(0,80,10),110)
  data_transm[,age_cat_num := .(cut(part_age,age_breaks,include.lowest = T,right = T)),]
  data_transm[,age_cat := .(paste0('age',as.numeric(age_cat_num))),]
  
  # set all age categories
  age_cat_all <- paste0('age',1:(length(age_breaks)-1))
  
  ## INCIDENCE ----
  # infections: captured in the "main statistics"   
  
  # infectious
  data_transm[, date_infectiousness := infection_date + start_infectiousness]
  summary_infectious  <- get_summary_table(data_transm,'date_infectiousness','age_cat','new_infectious_cases',age_cat_all)
  
  # symptomatic
  data_transm[, date_symptomatic := infection_date + start_symptoms]
  summary_symptomatic <- get_summary_table(data_transm,'date_symptomatic','age_cat','new_symptomatic_cases',age_cat_all)
  
  # recovered
  summary_recovered   <- get_summary_table(data_transm,'date_recovered','age_cat','new_recovered_cases',age_cat_all)
  
  # hospital admission: captured in the "main statistics"     
  
  ## DOUBLING TIME      ----
  # calculate cumulative cases and double time
  summary_out$doubling_time <- get_longitudinal_doubling_time(summary_out$new_infections)$mean
  # plot(summary_infections$doubling_time)
  
  # add average pre-lockdown doubling time
  ref_dates <- seq(as.Date('2020-02-24'),as.Date('2020-03-08'),1)
  summary_out[,doubling_time_march:=NA_real_]
  summary_out[sim_date %in% ref_dates, doubling_time_march :=  get_doubling_time(new_infections)]
  summary_out[,doubling_time_march]

  ## (A)SYMPTOMATIC TRANSMISSION    ----
  summary_symptomatic_infectors <- data_transm[,.(num_symptomatic_infectors = sum(infector_is_symptomatic,na.rm=T)), by=c('sim_date')]

  
  ## LOCATION ----
  pool_type_opt <- c("College","Household","HouseholdCluster",
                     "K12School","PrimaryCommunity","SecondaryCommunity","Workplace",'Collectivity') # TODO: make this flexible?
  summary_location <- get_summary_table(data_transm,'infection_date','pool_type','location',pool_type_opt) 
  summary_location$location    <- NULL # remove
  if('location_NA' %in% names(summary_location)) {
    summary_location$location_NA <- NULL # remove
  }
  
  # sort
  sorted_names <- (sort(names(summary_location)))
  summary_location <- summary_location[,..sorted_names,]
    
  
  ## AGGREGATE & MERGE ----
  # average generation interval, based on infection date
  sec_transm[,sim_date := infection_date]
  summary_generation     <- sec_transm[,.(gen_interval = mean(gen_interval,na.rm=T)),by=sim_date]
  summary_out    <- merge(summary_out,summary_generation,all.x = TRUE)
  
  # average number of secondary cases, upon recovery
  sec_transm[,sim_date := recoverd_date]
  summary_out    <- merge(summary_out,sec_transm[,.(sec_cases = mean(sec_cases,na.rm=T)),by=sim_date],all.x = TRUE)
  
  # add other details
  summary_out    <- merge(summary_out,summary_infectious,all.x = TRUE,nomatch=0)
  summary_out    <- merge(summary_out,summary_symptomatic,all.x = TRUE)
  summary_out    <- merge(summary_out,summary_symptomatic_infectors,all.x = TRUE)
  summary_out    <- merge(summary_out,summary_location,all.x = TRUE)

  ## CUMULATIVE STATS
  summary_out[,cumulative_infectious_cases := cumsum_na(new_infectious_cases)]
  summary_out[,cumulative_symptomatic_cases := cumsum_na(new_symptomatic_cases)]
  
  # return
  return(summary_out)
}


# separate function with main statistics
get_main_transmission_statistics <- function(data_transm,
                                             sim_date_range)
{
  # setup data.table
  data_transm[,ID := 1:nrow(data_transm),]
  
  ## RENAME DATE COLUMN
  data_transm[,sim_date := infection_date]
  
  # AGE CATEGORIES     ----
  age_breaks               <- c(seq(0,80,10),110)
  data_transm[,age_cat_num := .(cut(part_age,age_breaks,include.lowest = T,right = T)),]
  data_transm[,age_cat := .(paste0('age',as.numeric(age_cat_num))),]
  
  # set all age categories
  age_cat_all <- paste0('age',1:(length(age_breaks)-1))
  
  # set all dates
  summary_date <- data.table(sim_date = seq(min(sim_date_range),max(sim_date_range),1))
  setkey(summary_date,'sim_date')
  
  ## INCIDENCE ----
  # infections
  summary_infections   <- get_summary_table(data_transm,'infection_date','age_cat','new_infections',age_cat_all)
  
  # hospital admission     
  data_transm[, date_hosp_adm := infection_date + hospital_admission_start]
  summary_hospital          <- get_summary_table(data_transm,'date_hosp_adm','age_cat','new_hospital_admissions',age_cat_all)
  
  ## AGGREGATE & MERGE ----
  summary_out    <- merge(summary_date,summary_infections,all.x = TRUE)
  
  # average number of secondary cases, upon recovery
  summary_out    <- merge(summary_out,summary_hospital,all.x = TRUE)
  
  # add exp_id
  summary_out$exp_id <- unique(data_transm$exp_id)
  
  ## CUMULATIVE STATS
  #summary_out[,cumulative_infections := cumsum_na(new_infections)]
  for(i_colname in names(summary_infections)[-1]){
    summary_out[,(gsub('new_infections','cumulative_infections',i_colname)) := cumsum_na(get(i_colname))]
  }
  names(summary_out)
  
  for(i_colname in names(summary_hospital)[-1]){
    summary_out[,(gsub('new_hospital_admissions','cumulative_hospital_cases',i_colname)) := cumsum_na(get(i_colname))]
  }
  names(summary_out)
  
  # return
  return(summary_out)
}

# Function to generate summary tables
# colname_date <- 'infection_date'; colname_value <- 'cnt_location'; prefix <- 'location';
#colname_date <- 'date_hosp_adm'; colname_value <- 'age_cat'; prefix <- 'hosp';
# colname_date <- 'infection_date'; colname_value <-'age_cat';prefix <- 'new_infections'
# data_transm <- ref_case_data_table; colname_date="DATE";colname_value='AGEGROUP';prefix='d'
get_summary_table <- function(data_transm,colname_date,colname_value,prefix,colname_opt=NA){

  # overal summary
  summary_table_general         <- data_transm[,.N,by=colname_date]
  names(summary_table_general)  <- c('sim_date',prefix)

  # specific summary
  summary_table                  <- dcast(data_transm, formula(paste(colname_date, '~' ,colname_value)), value.var='ID', length)
  
  # if not all categories were persent, add column with 0's
  if(all(!is.na(colname_opt)) & any(!colname_opt %in% names(summary_table))){
    summary_table[,c(colname_opt[!colname_opt %in% names(summary_table)]):=0]
    setcolorder(summary_table,c(colname_date,colname_opt))# reorder
  }
  
  # update names
  names(summary_table)           <- c('sim_date',paste(prefix,names(summary_table)[-1],sep='_'))

  # merge
  summary_table <- merge(summary_table,summary_table_general)
  
  # omit NA's (default on LW's MACOS but not on VSC cluster)
  # R versions: 3.5.3 (MACOS) vs. 3.5.1 (VSC)
  summary_table <- na.omit(summary_table,cols='sim_date')

  #check
  head(summary_table)
  
  # return
  return(summary_table)
}

# get 
aggregate_transmission_dynamics <- function(project_dir){
  
  
  # load project summary
  project_summary      <- .rstride$load_project_summary(project_dir)
  
  # disable different start date
  project_summary$start_date <- project_summary$start_date[1]
  
  # retrieve all variable model parameters
  input_opt_design     <- .rstride$get_variable_model_param(project_summary)
  
  # get all transmission output
  data_incidence_all      <- .rstride$load_aggregated_output(project_dir,'data_incidence')
  
  if(length(data_incidence_all) == 1 && is.na(data_incidence_all)){
    smd_print('NO TRANSMISSION DATA AVAILABLE.')
    return(.rstride$no_return_value())
  }
  
  # select one week(7days), after burn in period of 10 days
  sel_period         <- seq(min(data_incidence_all$sim_date,na.rm=T)+9,min(data_incidence_all$sim_date,na.rm=T)+15,1)
  data_incidence_all <- data_incidence_all[data_incidence_all$sim_date %in% sel_period,]
  
  i_config <- 6
  foreach(i_config = 1:nrow(input_opt_design),
          .combine = 'rbind') %do%{
    
    # subset transmission output corresponding the 'input_opt_design' row
    flag_exp            <- .rstride$get_equal_rows(project_summary,input_opt_design[i_config,])
    data_incidence      <- data_incidence_all[data_incidence_all$exp_id %in% project_summary$exp_id[flag_exp],]
    num_runs_exp        <- sum(flag_exp)
    num_infected_seeds  <- data_incidence$new_infections[1]

    return(c(doubling_time = mean(data_incidence$doubling_time,na.rm=TRUE),
             sec_cases     = mean(data_incidence$sec_cases,na.rm=TRUE),
             gen_interval  = mean(data_incidence$gen_interval,na.rm=TRUE),
             r0 = input_opt_design[i_config,]))
    
    
  } -> transmission_results

  return(transmission_results)
}
