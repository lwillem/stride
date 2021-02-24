############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2021, Willem L
############################################################################ #
#
# TO CREATE CALENDAR FILE(S) FOR 2019-2021
#
# CONTAINING:
# 1. Public and school holidays
# 2. Contact reductions (covid19)
#       * pre-, primairy and secondary school
#       * workplace
#       * community
#       * household clusters
# 3. Imported cases (covid19)
# 4. Contact tracing (covid19)
# 5. Universal testing (covid19)
#
############################################################################ #

# debug
if(0==1){
  source('bin/rstride/rStride.R')
  cnt_other_exit_delay <- 21
  show_plots = TRUE
  
  create_calendar_file(show_plots=TRUE)
}


# create calendar files if they do not exist, else re-use them
create_calendar_file <- function(file_name_tag='2020_2021',show_plots = FALSE,file_name=NA)
{
  
  filename_calendar_full <- ifelse(is.na(file_name),
                                   paste('sim_output/calendar_belgium',file_name_tag,'covid19.csv',sep='_'),
                                   file_name)

  ########################################### #
  ## INITIATE DATA                       ####
  ########################################### #
  
  date_start <- as.Date("2019-01-01")
  date_end   <- as.Date("2021-12-31")
  
  # default value in C++ CALENDAR vectors ==>> 0
  
  ########################################### #
  ## 1.a Public holidays                 ####
  ########################################### #
  
  data.table(category = "general",
             date     = as.Date(c(
               '2019-01-01','2019-04-22','2019-05-01','2019-05-30','2019-06-10', # 2019
               '2019-07-21','2019-08-15','2019-11-01','2019-11-11','2019-12-25',
               
               '2020-01-01','2020-04-13','2020-05-01','2020-05-21','2020-06-01', # 2020
               '2020-07-21','2020-08-15','2020-11-01','2020-11-11','2020-12-25',
               
               '2021-01-01','2021-04-05','2021-05-01','2021-05-13','2021-06-24', # 2021
               '2021-07-21','2021-08-15','2021-11-01','2021-11-11','2021-12-25')),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
            ) -> d_calendar_holiday
  
  summary(d_calendar_holiday)
  
  ########################################### #
  ## 1.b School holidays                 ####
  ########################################### #
  
  # start from blanc calendar
  data.table(category = "schools_closed",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0.0,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_school_holidays
  d_college_holidays <- copy(d_school_holidays)
  
  # include K12 school holidays
  d_school_holidays[date %in% c(seq(as.Date('2019-01-01'),as.Date('2019-01-06'),1), # 2019
                                seq(as.Date('2019-03-04'),as.Date('2019-03-10'),1),
                                seq(as.Date('2019-04-08'),as.Date('2019-04-22'),1),
                                seq(as.Date('2019-07-01'),as.Date('2019-08-31'),1),
                                seq(as.Date('2019-10-28'),as.Date('2019-11-03'),1),
                                seq(as.Date('2019-12-23'),as.Date('2019-12-31'),1),
                                
                                seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
                                seq(as.Date('2020-02-24'),as.Date('2020-02-29'),1),
                                seq(as.Date('2020-04-06'),as.Date('2020-04-19'),1),
                                seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),
                                seq(as.Date('2020-11-02'),as.Date('2020-11-08'),1),
                                seq(as.Date('2020-12-21'),as.Date('2020-12-31'),1),
                                
                                seq(as.Date('2021-01-01'),as.Date('2021-01-03'),1), # 2021
                                seq(as.Date('2021-02-15'),as.Date('2021-02-21'),1),
                                seq(as.Date('2021-04-05'),as.Date('2021-04-18'),1),
                                seq(as.Date('2021-07-01'),as.Date('2021-08-31'),1),
                                seq(as.Date('2021-11-01'),as.Date('2021-11-07'),1),
                                seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1)),
             value    := 1.0]
  
  # add college holidays
  d_college_holidays[date %in% c(seq(as.Date('2019-01-01'),as.Date('2019-01-06'),1), # 2019
                          seq(as.Date('2019-03-04'),as.Date('2019-03-10'),1),
                          seq(as.Date('2019-04-08'),as.Date('2019-04-22'),1),
                          seq(as.Date('2019-07-01'),as.Date('2019-09-22'),1), # summer break untill September, 22
                          #seq(as.Date('2019-10-28'),as.Date('2019-11-03'),1), # no fall break
                          seq(as.Date('2019-12-23'),as.Date('2019-12-31'),1),
                          
                          seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
                          seq(as.Date('2020-02-24'),as.Date('2020-02-29'),1),
                          seq(as.Date('2020-04-06'),as.Date('2020-04-19'),1),
                          seq(as.Date('2020-07-01'),as.Date('2020-09-20'),1),# summer break untill September, 20
                          #seq(as.Date('2020-11-02'),as.Date('2020-11-08'),1), # no fall break
                          seq(as.Date('2020-12-21'),as.Date('2020-12-31'),1),
                          
                          seq(as.Date('2021-01-01'),as.Date('2021-01-03'),1), # 2021
                          seq(as.Date('2021-02-15'),as.Date('2021-02-21'),1),
                          seq(as.Date('2021-04-05'),as.Date('2021-04-18'),1),
                          seq(as.Date('2021-07-01'),as.Date('2021-09-19'),1), # summer break untill September, 19
                          #seq(as.Date('2021-11-01'),as.Date('2021-11-07'),1), # no fall break
                          seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1)),
             value    := 1.0]
  
  #K12 school
  tmp_school_holidays <- copy(d_school_holidays)
  tmp_school_holidays[,category:='schools_closed']
  for(i_age in 0:17){
    d_calendar_holiday <- rbind(d_calendar_holiday,copy(tmp_school_holidays[,age:=i_age]))
  }
  
  # College
  tmp_college_holidays <- copy(d_college_holidays)
  tmp_college_holidays[,category:='schools_closed']
  for(i_age in 18:25){
    d_calendar_holiday <- rbind(d_calendar_holiday,copy(tmp_college_holidays[,age:=i_age]))
  }
  
  
  ########################################################### #
  ##  2a. Contact reductions: school closures              ####
  ########################################################### #
  #       * (pre-, primary and secondary school)
  
  # set default school closure
  # school_dates_non_holiday <- date_all[!date_all %in% dcal_school_closure]
  # data.table(category = "schools_closed",
  #            # date     = seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1),
  #            date     = school_dates_non_holiday,
  #            value    = 0.0,
  #            type = 'double',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> d_school_closure
  # 
  # d_school_closure[date %in% seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1),value := 1.0]
  # 
  # tmp_school_closure <- copy(d_school_closure)
  # tmp_school_closure[,category:='schools_closed']
  
  d_calendar_holiday[date %in% seq(as.Date('2020-03-14'),as.Date('2020-06-30'),1), value := 1.0] 

  # set eligible dates for school reopening in May-June 2020
  d_school_reopening <- seq(as.Date('2020-05-18'),as.Date('2020-06-30'),1)
  d_school_reopening_wday <- as.POSIXlt(d_school_reopening)$wday
  
  # preschool (reopens 4d/week)
  d_school_reopening_4d <- d_school_reopening[d_school_reopening_wday %in% 1:4]
  d_calendar_holiday[date %in% d_school_reopening_4d & age %in% c(0,1,2,6,7),value:=0.5]
  
  # primary school (reopens 2d/week)
  d_school_reopening_2d <- d_school_reopening[d_school_reopening_wday %in% 4:5]
  d_school_reopening_2d[1:2] <- d_school_reopening_2d[1:2] - 2 # fix for holidays Thu-Fri in May
  d_calendar_holiday[date %in% d_school_reopening_2d & age %in% c(11),value:=0.5]
  
  
  #secondary school (reopens 1d week)
  d_school_reopening_1d <- d_school_reopening[d_school_reopening_wday %in% 3]
  d_calendar_holiday[date %in% d_school_reopening_1d & age %in% c(17),value:=0.5]
  
  # school reopening September 2020
  # up to primary school
  d_calendar_holiday[date >= as.Date('2020-09-01') &
                       age <= 12 &
                     value == 0, value := 0.5]
  # secondary school
  d_calendar_holiday[date >= as.Date('2020-05-01') &
                       age > 12 & age < 18 &
                       value != 1, value := 0.2]
  # tertiary eduction
  d_calendar_holiday[date >= as.Date('2020-09-01') &
                       age >= 18 &
                       value == 0, value := 0.3]
  
  # tertiary eduction: closed from November 1st
  d_calendar_holiday[date >= as.Date('2020-11-01') &
                       age >= 18 , value := 1]

  ########################################### #
  ##  2b. Contact reductions: other        ####
  ########################################### #
  #       * workplace
  #       * community
  #       * household clusters
  
  # workplace distancing
  data.table(category = "workplace_distancing",
             #date     = seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0.0,
             type     = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_workplace_distancing
  dcal_workplace_distancing[date %in% seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),value := 1.0]
  
  
  # community distancing
  data.table(category = "community_distancing",
             #date     = seq(as.Date('2020-03-14'),as.Date('2020-05-24'),1),
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0.0,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_community_distancing
  dcal_community_distancing[date %in% seq(as.Date('2020-03-14'),as.Date('2020-05-24'),1),value := 1.0]
  
  # collectivity distancing
  data.table(category = "collectivity_distancing",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0.0,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_collectivity_distancing
  
  # household clustering
  data.table(category = "household_clustering",
             date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_household_clustering
  
  ########################################### #
  ## 3. Imported cases                 ####
  ########################################### #
  
  data.table(category = "imported_cases",
             # date     = seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_imported_cases
  dcal_imported_cases[date %in% seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),value := 1]
  
  ########################################### #
  ##  4. Contact tracing                 ####
  ########################################### #
  
  data.table(category = "contact_tracing",
             date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_contact_tracing
  
  ########################################### #
  ## 5. Universal testing                 ####
  ########################################### #

  data.table(category = "universal_testing",
             date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_universal_testing
  
  
  ########################################### #
  ## MERGE HOLIDAYS & OTHER CALENDAR ITEMS ####
  ########################################### #
  
  # get 'dcal_*' variables
  opt_other <- ls(pattern='dcal_')
  
  # combine all 'dcal_*' variable
  d_calendar_all <- foreach(i_other  = opt_other,
                            .init    = d_calendar_holiday,
                            .combine = 'rbind') %do% {
                              get(i_other)
                            } 
  
  # select range
  d_calendar_all <- d_calendar_all[date >= date_start & date <= date_end,]
  range(d_calendar_all$date)
  
  ########################################### #
  ## EXPLORE DATA                         ####
  ########################################### #
  
  plot_calendar(dt_calendar            = d_calendar_all,
                filename_calendar_full = filename_calendar_full,
                show_plots             = show_plots)
  
  ########################################### #
  ## SAVE AS CSV	 	         ####
  ########################################### #
  
  # # format date
  # d_calendar_all[,date:=format(date,'%Y-%m-%d')]
  # format(d_calendar_all$date,'%Y-%m-%d')
  
  # save as csv (all calendar info)
  write.table(d_calendar_all,
              file = filename_calendar_full,sep=',',row.names=F,quote=F)
  
  unique(d_calendar_all$category)
  
  return(filename_calendar_full)
}

# note: variable "b_school_repopening" is not used anymore... but still here for backward compatibility
plot_calendar <- function(dt_calendar, filename_calendar_full, show_plots = TRUE, b_school_repopening=TRUE){
  
  if(show_plots){
    
    # open pdf stream
    pdf(file=gsub('.csv','.pdf',filename_calendar_full),6,6)
    
    category_opt <- unique(dt_calendar$category)
    par(mfrow=c(3,2))
    
    # check if dt_calendar is data.table
    if(!is.data.table(dt_calendar)){
      dt_calendar <- data.table(dt_calendar)
    }
    
    # make sure that "date" is in date format
    dt_calendar$date <- as.Date(dt_calendar$date)
    
    # x_lim      <- range(dt_calendar$date)
    x_lim      <- as.Date(c('2020-02-01','2020-12-31'))
    x_lab_year <- paste(unique(year(dt_calendar$date)),sep='-')
    i_cat <- category_opt[2]
    for(i_cat in category_opt){
      plot(x   = dt_calendar[category == i_cat,date],
           y   = dt_calendar[category == i_cat,value],
           xlim = x_lim,
           ylim = range(0,1,dt_calendar$value[dt_calendar$category == i_cat]),
           col  = 1,
           #type='l',
           pch  = 15,
           #lwd=2,
           main = i_cat,
           bty='n',
           xlab = x_lab_year,
           ylab = unique(dt_calendar[,type]),
           xaxt = 'n'
      )
      add_x_axis(x_lim)
      abline(h=1,lty=3,col='grey')
    }
    
    if("schools_closed" %in% dt_calendar$category){
      
      i_cat <- "schools_closed"

      # convert value into numeric factors (to use as color)
      value_levels            <- c(unique(dt_calendar[category == i_cat & value > 0,'value']))
      dt_calendar$value_level <- factor(dt_calendar$value,levels=unlist(value_levels))
      dt_calendar$value_col   <- as.numeric(dt_calendar$value_level)
      
      plot(x   = dt_calendar[category == i_cat,date],
           y   = dt_calendar[category == i_cat,value],
           xlim = x_lim,
           ylim = range(0,1,dt_calendar$value[dt_calendar$category == i_cat]),
           col  = dt_calendar[category == i_cat,value_col],
           #type='l',
           pch  = 15,
           #lwd=2,
           main = i_cat,
           bty='n',
           xlab = x_lab_year,
           ylab = unique(dt_calendar[,type]),
           xaxt = 'n'
      )
      add_x_axis(x_lim)
      abline(h=1,lty=3,col='grey')
      
     
      
      # plot by age
      plot(x   = dt_calendar[category == i_cat & value == 1,date],
           y   = dt_calendar[category == i_cat & value == 1,age],
           xlim = x_lim,
           ylim = range(0,1,dt_calendar$age,na.rm=T),
           col  = 1,
           pch  = 15,
           main = i_cat,
           bty='n',
           xlab = x_lab_year,
           ylab = 'age',
           xaxt = 'n'
      )
      points(x    = dt_calendar[category == i_cat ,date],
             y    = dt_calendar[category == i_cat ,age],
             col  = dt_calendar[category == i_cat ,value_col],
             pch  = 15
      )
      add_x_axis(x_lim)
    }
    
    # close pdf stream
    dev.off()
  }
}

# db_cat <- "workplace_distancing";date_start<-"2020-03-14";date_end <- "2020-05-02";db_value<-0.75;file_name<-dcal_file
# db_cat <- "workplace_distancing";date_selection <- wp_dist$x; db_value <- wp_dist$y; file_name <- dcal_file
# db_update <- data.frame(c('2020-02-10',0),
#                         c('2020-03-13',0),
#                         c('2020-03-19',0.85),
#                         c('2020-05-02',0.85),
# c('2020-05-03',0.75)); db_cat <- "workplace_distancing";db_age = 'NA'; file_name <- "sim_output/calendar_belgium_wp_fitting_covid19.csv"
# db_category =  "workplace_distancing";db_update = dcal_wp_distancing; file_name = dcal_file
adjust_calendar_file <- function(db_category, db_update, file_name, db_age = 'NA', show_plots=FALSE){
  
  # file_name fix, exclude '../'
  file_name <- gsub('\\.\\.','\\.',file_name)
  
  # read calendar file
  d_calendar_all <- data.table(read.table(file=file_name,sep=',',header=T))
  
  # make sure the value is stored as double
  d_calendar_all$value <- as.double(d_calendar_all$value)
  
  # adjust age == NA
  d_calendar_all[, age_char := as.character(age)]
  d_calendar_all[is.na(age_char), age_char := 'NA']
  
  # check category
  if(!db_category %in% unique(d_calendar_all$category)){
    smd_print("CALENDAR CATEGORY UNKNOWN => STOP CALENDAR ADJUSTMENT")
    smd_print("CALENDAR CATEGORY OPTIONS:", paste0(unique(d_calendar_all$category),collapse=', '))
    
    return(NA)
  }
  
  # extrapolate given dates and values
  db_update  <- data.frame(t(db_update))
  date_out   <- seq(min(as.Date(db_update[,1])),max(as.Date(db_update[,1])),1) 
  date_out   <- date_out[date_out<=max(d_calendar_all$date)]
  db_update  <- approx(x=as.Date(db_update[,1]),
                      y=db_update[,2],
                      xout = date_out)
  names(db_update) <- c('date','value')

  
  # update values for all given date and ages
  for(i_db_age in as.character(db_age)){
    d_calendar_all[as.character(date) %in% as.character(db_update$date) &
                      category == db_category &
                      age_char == i_db_age, 
                   value := db_update$value ]
  }

  d_calendar_all[as.character(date) %in% as.character(db_update$date) &
                   category == db_category &
                   age_char == db_age]
  
  # explore
  plot_calendar(dt_calendar            = d_calendar_all,
                filename_calendar_full = file_name,
                show_plots             = show_plots)
  
  # save as csv 
  write.table(d_calendar_all,
              file = file_name,sep=',',row.names=F,quote=F)
  
}

replace_calendar_value <- function(file_name,db_category,value_orig,value_new,show_plots){
  
  # read calendar file
  d_calendar_all <- data.table(read.table(file=file_name,sep=',',header=T))
  
  # make sure the value is stored as double
  d_calendar_all$value <- as.double(d_calendar_all$value)
  
  d_calendar_all[category == db_category &
                   value == value_orig, value:= value_new]
  
  # explore
  plot_calendar(dt_calendar            = d_calendar_all,
                filename_calendar_full = file_name,
                show_plots             = show_plots)
  
  # save as csv 
  write.table(d_calendar_all,
              file = file_name,sep=',',row.names=F,quote=F)
  
  
  
}

# create calendar file comparable to the original lockdown/exit parameter structure
integrate_lockdown_parameters_into_calendar <- function(config_exp){

  file_name <- smd_file_path(config_exp$output_prefix,'calendar_belgium_covid19_v1_1_param.csv')
  config_exp$holidays_file <- create_calendar_file(file_name = file_name, show_plots = T)
  
  # set dates
  date_start            <- config_exp$start_date
  date_t0               <- as.Date('2020-03-13')
  date_compliance_wp    <- date_t0 + config_exp$compliance_delay_workplace
  date_compliance_other <- date_t0 + config_exp$compliance_delay_other
  date_compliance_collectivity <- date_t0 + config_exp$compliance_delay_collectivity
  date_exit_wp          <- as.Date('2020-05-04')
  date_exit_other       <- as.Date('2020-05-25')
  date_end              <- as.Date('2020-12-31')
  
  # integrate school distancing (general: K6)
  replace_calendar_value(file_name = config_exp$holidays_file,
                         db_category =  "schools_closed",
                         value_orig = 0.5,
                         value_new = config_exp$cnt_reduction_school_exit,
                         show_plots = T)
  
  # integrate school distancing: secondary
  replace_calendar_value(file_name = config_exp$holidays_file,
                         db_category =  "schools_closed",
                         value_orig = 0.2,
                         value_new = ifelse('cnt_reduction_school_exit_secondary' %in% names(config_exp),
                                            config_exp$cnt_reduction_school_exit_secondary,
                                            config_exp$cnt_reduction_school_exit),
                         show_plots = T)

  
  # integrate school distancing: tertiary)
  replace_calendar_value(file_name = config_exp$holidays_file,
                         db_category =  "schools_closed",
                         value_orig = 0.3,
                         value_new = ifelse('cnt_reduction_school_exit_tertiary' %in% names(config_exp),
                                            config_exp$cnt_reduction_school_exit_tertiary,
                                            config_exp$cnt_reduction_school_exit),
                         show_plots = T)

  
  # integreate workplace distancing
  adjust_calendar_file(db_category =  "workplace_distancing",
                       db_update = data.frame(c(as.character(date_t0),0),
                                              c(as.character(date_compliance_wp),config_exp$cnt_reduction_workplace),
                                              c(as.character(date_exit_wp-1),config_exp$cnt_reduction_workplace),
                                              c(as.character(date_exit_wp),config_exp$cnt_reduction_workplace_exit),
                                              c(as.character(date_end),config_exp$cnt_reduction_workplace_exit)),
                       file_name = config_exp$holidays_file )
  # config_exp$cnt_reduction_workplace <- 1
  # config_exp$compliance_delay_workplace <- 0
  # config_exp$cnt_reduction_workplace_exit <- 0
  
  # integreate community distancing
  adjust_calendar_file(db_category =  "community_distancing",
                       db_update = data.frame(c(as.character(date_t0),0),
                                              c(as.character(date_compliance_other),config_exp$cnt_reduction_other),
                                              c(as.character(date_exit_other-1),config_exp$cnt_reduction_other),
                                              c(as.character(date_exit_other),config_exp$cnt_reduction_other_exit),
                                              c(as.character(date_end),config_exp$cnt_reduction_other_exit)),
                       file_name = config_exp$holidays_file,
                       show_plots = T)
  # config_exp$cnt_reduction_other <- 1
  # config_exp$compliance_delay_other <- 0
  # config_exp$cnt_reduction_other_exit <- 0
  
  # integreate collectivity distancing
  if(!any(is.null(c(config_exp$cnt_baseline_collectivity,config_exp$cnt_reduction_collectivity)))){
    adjust_calendar_file(db_category =  "collectivity_distancing",
                         db_update = data.frame(c(as.character(date_start),config_exp$cnt_baseline_collectivity),
                                                c(as.character(date_t0),config_exp$cnt_baseline_collectivity),
                                                c(as.character(date_compliance_collectivity),config_exp$cnt_reduction_collectivity),
                                                c(as.character(date_end),config_exp$cnt_reduction_collectivity)),
                         file_name = config_exp$holidays_file,
                         show_plots = T)
  }
  
  if('temporal_distancing_workplace' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'workplace_distancing',
                                        db_values_char = paste(config_exp$cnt_reduction_workplace,config_exp$temporal_distancing_workplace,sep=','),
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char       = config_exp$dates_distancing_workplace)
  }
  
  if('temporal_distancing_community' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'community_distancing',
                                        db_values_char = paste(config_exp$cnt_reduction_other,config_exp$temporal_distancing_community,sep=','),
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char       = config_exp$dates_distancing_community)
  }
  
  if('temporal_distancing_collectivity' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'collectivity_distancing',
                                        db_values_char = paste(config_exp$cnt_reduction_collectivity,config_exp$temporal_distancing_collectivity,sep=','),
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$dates_distancing_collectivity)
  }
  
  if('temporal_imported_cases' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'imported_cases',
                                        db_values_char = paste(0,config_exp$temporal_imported_cases,sep=','),
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$dates_imported_cases)
  }
  
  

  
  #TODO: add household cluster mixing
  
  # # fix for calendar path
  config_exp$holidays_file <- paste0('../',config_exp$holidays_file)
  
  # return list
  return(config_exp)
}

# db_category <- 'distancing_workplace'
# db_values <- seq(0.8,0.9,length=12)
include_temporal_distancing_factors <- function(db_category,db_values_char,file_name,show_plots=T,db_dates_char=NA){
  
  # split given value string, into numerical values
  db_values <- as.numeric(unlist(strsplit(db_values_char,',')))
  
  if(is.na(db_dates_char)){
    # Start with May 1st using the lockdown measure
    date_start <- as.Date('2020-05-01') # use lockdown values
    
    # Create list with eligible dates for given values (flexible, constant day of the month)
    # note: max number of days per month = 31, max period = 31*number of values
    date_day_select <- 1
    date_all      <- seq(date_start,date_start+(31*length(db_values)),1)
    db_dates      <- date_all[as.numeric(format(date_all,'%d')) == date_day_select]
    
    # safety check!
    if(length(db_values) < length(db_dates)){
      db_dates <- db_dates[1:(length(db_values))]
    }
    
  } else{
    db_dates <- as.Date(unlist(strsplit(db_dates_char,',')))
  }
  
  # add right tail
  db_dates <- c(db_dates,db_dates[length(db_dates)]+365)
  db_values <- c(db_values,db_values[length(db_values)])
  
  adjust_calendar_file(db_category = db_category,
                       db_update   = rbind(as.character(db_dates),db_values),
                       file_name   = file_name,
                       show_plots  = T)
}


if(0==1){ # debug----
  
  dcal_file <- create_calendar_file(file_name_tag = 'wp_fitting',show_plots = T)
  dcal_wp_distancing <- data.frame(c('2020-03-13',0),
                                   c('2020-03-19',0.85),
                                   c('2020-05-02',0.85),
                                   c('2020-05-03',0.5))
  
  adjust_calendar_file(db_category =  "workplace_distancing",db_update = dcal_wp_distancing, file_name = dcal_file)
}

