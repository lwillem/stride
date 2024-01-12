############################################################################ #
#  This file is part of the Stride software. 
#
#  Copyright 2023, Willem L
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
               '2021-07-21','2021-08-15','2021-11-01','2021-11-11','2021-12-25',

               '2022-01-01','2022-04-18','2022-05-01','2022-05-26','2022-06-06', # 2022
               '2022-07-21','2022-08-15','2022-11-01','2022-11-11','2022-12-25')),
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
                                seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1),

                                seq(as.Date('2022-01-01'),as.Date('2022-01-09'),1), # 2022
                                seq(as.Date('2022-02-28'),as.Date('2022-03-06'),1),
                                seq(as.Date('2022-04-04'),as.Date('2022-04-18'),1),
                                seq(as.Date('2022-07-01'),as.Date('2022-08-31'),1),
                                seq(as.Date('2022-10-31'),as.Date('2022-11-06'),1),
                                seq(as.Date('2022-12-26'),as.Date('2022-12-31'),1)),
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
                          seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1),

                          seq(as.Date('2022-01-01'),as.Date('2022-01-09'),1), # 2022
                          seq(as.Date('2022-02-28'),as.Date('2022-03-06'),1),
                          seq(as.Date('2022-04-04'),as.Date('2022-04-18'),1),
                          seq(as.Date('2022-07-01'),as.Date('2022-09-25'),1), ####
                          # seq(as.Date('2022-10-31'),as.Date('2022-11-06'),1), # no fall break
                          seq(as.Date('2022-12-26'),as.Date('2022-12-31'),1)),
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
  
  # # workplace distancing
  # data.table(category = "workplace_distancing",
  #            #date     = seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),
  #            date     = seq(as.Date(date_start),as.Date(date_end),1),
  #            value    = 0.0,
  #            type     = 'double',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_workplace_distancing
  # dcal_workplace_distancing[date %in% seq(as.Date('2020-03-14'),as.Date('2020-05-03'),1),value := 1.0]
  # 
  # 
  # # community distancing
  # data.table(category = "community_distancing",
  #            #date     = seq(as.Date('2020-03-14'),as.Date('2020-05-24'),1),
  #            date     = seq(as.Date(date_start),as.Date(date_end),1),
  #            value    = 0.0,
  #            type = 'double',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_community_distancing
  # dcal_community_distancing[date %in% seq(as.Date('2020-03-14'),as.Date('2020-05-24'),1),value := 1.0]
  
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
  
  # data.table(category = "imported_cases",
  #            # date     = seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),
  #            date     = seq(as.Date(date_start),as.Date(date_end),1),
  #            value    = 0,
  #            type = 'boolean',
  #            age = NA_integer_,
  #            stringsAsFactors = F
  # ) -> dcal_imported_cases
  # dcal_imported_cases[date %in% seq(as.Date('2020-07-01'),as.Date('2020-08-31'),1),value := 1]
  
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
    x_lim      <- c(as.Date('2020-02-01'),max(dt_calendar$date))
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
           xlab = '',
           ylab = unique(dt_calendar[category == i_cat,type]),
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
           xlab = '',
           ylab = unique(dt_calendar[category == i_cat,type]),
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
           xlab = '',
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
  
  d_calendar_categories <- c('general', 
                             'schools_closed', 
                             'workplace_distancing',
                             'collectivity_distancing', 
                             'community_distancing', 
                             'contact_tracing', 
                             'household_clustering', 
                             'imported_cases')
  
  # check category
  if(!db_category %in% d_calendar_categories){
    smd_print("CALENDAR CATEGORY UNKNOWN => STOP CALENDAR ADJUSTMENT")
    smd_print("CALENDAR CATEGORY OPTIONS:", paste0(d_calendar_categories,collapse=', '))
    return(NA)
  }
  
  
  # create data.frame with all information to extrapolate
  df_update  <- data.frame(t(db_update))
  date_out   <- seq(min(as.Date(df_update[,1])),max(as.Date(df_update[,1])),1)
  date_out   <- date_out[date_out<=max(d_calendar_all$date)]
  
  # # if db_category is not present yet, extend first value
  # if(!db_category %in% unique(d_calendar_all$category)){
  #   df_update       <- df_update[c(1,1:nrow(df_update)),]
  #   df_update$V1[1] <- min(d_calendar_all$date)
  #   date_out        <- sort(unique(d_calendar_all$date))
  # }

  # extrapolate given dates and values
  #date_out   <- date_out[date_out<=as.Date("2021-12-31")]
  df_update_full  <- approx(x=as.Date(df_update[,1]),
                            y=df_update[,2],
                            xout = as.Date(date_out),
                            method="linear")
  names(df_update_full) <- c('date','value')

  # integrate (new) values in calendar
  for(i_db_age in as.character(db_age)){
    # remove old values (if any)
    d_calendar_all <- d_calendar_all[!(as.character(date) %in% as.character(df_update_full$date) &
                                       category == db_category &
                                       age_char == i_db_age),]
    # include new values
    dcal_new <- data.table(category = db_category,
                           date     = paste(df_update_full$date),
                           value    = df_update_full$value,
                           type     = 'double',
                           age = ifelse(i_db_age == 'NA', NA_integer_,as.numeric(i_db_age)),
                           age_char = i_db_age,
                           stringsAsFactors = F
    ) 
    d_calendar_all <- rbind(d_calendar_all,dcal_new) 
  }

  # check
  d_calendar_all[as.character(date) %in% as.character(df_update_full$date) &
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

exclude_calendar_category <- function(file_name,db_category,show_plots=FALSE){
  
  # read calendar file
  d_calendar_all <- data.table(read.table(file=file_name,sep=',',header=T))
  
  # remove category (if present)
  d_calendar_all <- d_calendar_all[category != db_category,]
  
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

  # if there are not distancing parameters, return original config_exp
  param_calendar <- config_exp[grepl('cnt_reduction_workplace',names(config_exp)) |   # OR colname contains reduction_workplace
                                   grepl('clustering',names(config_exp)) |            # OR colname contains clustering
                                   grepl('imported',names(config_exp)) |              # OR colname contains imported
                                   grepl('distancing',names(config_exp)) &            # OR colname contains distancing)
                                   !is.na(config_exp)]                                # AND different from NA 
  param_calendar <- unlist(param_calendar)
  param_calendar[is.na(param_calendar)] <- 0

  if(length(param_calendar) == 0 || !any(param_calendar!=0)){
    return(config_exp)
  }
  
  # # else, modify calendar
  file_name_new <- smd_file_path(config_exp$output_prefix,'calendar_belgium_covid19_v1_1_param.csv')
  file_name_exp <- file.path('data',config_exp$holidays_file)
  # config_exp$holidays_file <- create_calendar_file(file_name = file_name, show_plots = T)
  
  if(file.exists(file_name_exp)){
    file.copy(from=file_name_exp,
              to = file_name_new,overwrite = TRUE)
    config_exp$holidays_file <- file_name_new
  } else{
    config_exp$holidays_file <- create_calendar_file(file_name = file_name_new, show_plots = T)
  }
 
  if('distancing_workplace_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'workplace_distancing',
                                        db_values_char = config_exp$distancing_workplace_ratio,
                                        db_delay_char  = config_exp$distancing_workplace_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$distancing_workplace_date)
  }
  
  if('distancing_community_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'community_distancing',
                                        db_values_char = config_exp$distancing_community_ratio,
                                        db_delay_char  = config_exp$distancing_community_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$distancing_community_date)
  }
  
  if('distancing_collectivity_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'collectivity_distancing',
                                        db_values_char = config_exp$distancing_collectivity_ratio,
                                        db_delay_char  = config_exp$distancing_collectivity_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$distancing_collectivity_date)
  }
  
  if('imported_cases_number' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'imported_cases',
                                        db_values_char = config_exp$imported_cases_number,
                                        db_delay_char  = config_exp$imported_cases_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$imported_cases_date)
  }
  
  if('household_clustering_ratio' %in% names(config_exp)){
    include_temporal_distancing_factors(db_category    = 'household_clustering',
                                        db_values_char = config_exp$household_clustering_ratio,
                                        db_delay_char  = config_exp$household_clustering_delay,
                                        file_name      = config_exp$holidays_file,
                                        show_plots     = T,
                                        db_dates_char  = config_exp$household_clustering_date)
  }
  
  # # fix for calendar path
  config_exp$holidays_file <- paste0('../',config_exp$holidays_file)
  
  # return list
  return(config_exp)
}

# db_category <- 'distancing_workplace'
# db_values <- seq(0.8,0.9,length=12)
include_temporal_distancing_factors <- function(db_category,db_values_char,db_delay_char,file_name,show_plots=T,db_dates_char=NA){
  
  # check input parameters
  vector_input_param <- c(db_category,db_values_char,db_delay_char,db_dates_char)
  if(!any(is.na(vector_input_param)) & all(vector_input_param != "NA"))
  {
   
  # make sure input parameters are "character" types
  db_values_char <- c_str(db_values_char)
  db_delay_char  <- c_str(db_delay_char)
  db_dates_char  <- paste(db_dates_char)
  
  # split given value string, into numerical values
  db_values <- as.numeric(unlist(strsplit(db_values_char,',')))
  db_delay  <- as.numeric(unlist(strsplit(db_delay_char,',')))
  db_dates  <- as.Date(unlist(strsplit(db_dates_char,',')))
  
  # account for delay == 0 by using "db_date-1" and "delay 1" 
  bool_delay_zero <- db_delay == 0
  if(any(bool_delay_zero)){
    db_delay[bool_delay_zero] <- 1
    db_dates[bool_delay_zero] <- db_dates[bool_delay_zero] - 1
  }
 
  # account for delay in compliance
  db_dates  <- c(db_dates[1],db_dates + db_delay,db_dates[-1])
  db_values <- c(0,db_values,db_values[-length(db_values)])
  
  # sort
  db_values <- db_values[order(as.Date(db_dates))]
  db_dates  <- db_dates[order(as.Date(db_dates))]
  
  # add right tail
  db_dates  <- c(db_dates,db_dates[length(db_dates)]+365)
  db_values <- c(db_values,db_values[length(db_values)])
  
  adjust_calendar_file(db_category = db_category,
                       db_update   = rbind(as.character(db_dates),db_values),
                       file_name   = file_name,
                       show_plots  = T)
  
  } # end if-clause on is.na
}

# Create a calendar file and fill in contact reduction values
create_new_cnt_calendar_file <- function(file_name, config_exp, end_date="2021-12-31", school_holidays=FALSE) {

  ########################################### #
  ## INITIATE DATA                       ####
  ########################################### #
  date_start <- as.Date(config_exp$start_date)
  date_end   <- as.Date(end_date)

  ####################################### #
  ## Public holidays                 ####
  ####################################### #
  data.table(category = "general",
             date     = as.Date(c(
               '2020-01-01','2020-04-13','2020-05-01','2020-05-21','2020-06-01', # 2020
               '2020-07-21','2020-08-15','2020-11-01','2020-11-11','2020-12-25',

               '2021-01-01','2021-04-05','2021-05-01','2021-05-13','2021-06-24', # 2021
               '2021-07-21','2021-08-15','2021-11-01','2021-11-11','2021-12-25')),
             value    = 1,
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_calendar_holiday
  # summary(d_calendar_holiday)

  ################################################ #
  ## Contact reductions: schools              ####
  ################################################ #
  cnt_reduction_school <- ifelse('cnt_reduction_school' %in% names(config_exp),
                                 config_exp$cnt_reduction_school,
                                 0.0)
  cnt_reduction_school_secondary <- ifelse('cnt_reduction_school_secondary' %in% names(config_exp),
                                           config_exp$cnt_reduction_school_secondary,
                                           cnt_reduction_school)
  cnt_reduction_school_tertiary <- ifelse('cnt_reduction_school_tertiary' %in% names(config_exp),
                                          config_exp$cnt_reduction_school_tertiary,
                                          cnt_reduction_school)

  data.table(category = "schools_closed",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = cnt_reduction_school,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_school_holidays
  data.table(category = "schools_closed",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = cnt_reduction_school_secondary,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_school_holidays_secondary
  data.table(category = "schools_closed",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = cnt_reduction_school_tertiary,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> d_college_holidays

  ####################################### #
  ## School holidays                 ####
  ####################################### #
  # Add school holidays if requested
  if (school_holidays) {
    smd_print("Including school holidays...")
    # include K12 school holidays
    d_school_holidays[date %in% c(seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
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
                                  seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1),

                                  seq(as.Date('2022-01-01'),as.Date('2022-01-09'),1), # 2022
                                  seq(as.Date('2022-02-28'),as.Date('2022-03-06'),1),
                                  seq(as.Date('2022-04-04'),as.Date('2022-04-18'),1),
                                  seq(as.Date('2022-07-01'),as.Date('2022-08-31'),1),
                                  seq(as.Date('2022-10-31'),as.Date('2022-11-06'),1),
                                  seq(as.Date('2022-12-26'),as.Date('2022-12-31'),1)),
                      value    := 1.0]
    d_school_holidays_secondary[date %in% c(seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
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
                                            seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1),

                                            seq(as.Date('2022-01-01'),as.Date('2022-01-09'),1), # 2022
                                            seq(as.Date('2022-02-28'),as.Date('2022-03-06'),1),
                                            seq(as.Date('2022-04-04'),as.Date('2022-04-18'),1),
                                            seq(as.Date('2022-07-01'),as.Date('2022-08-31'),1),
                                            seq(as.Date('2022-10-31'),as.Date('2022-11-06'),1),
                                            seq(as.Date('2022-12-26'),as.Date('2022-12-31'),1)),
                                value    := 1.0]

    # add college holidays
    d_college_holidays[date %in% c(seq(as.Date('2020-01-01'),as.Date('2020-01-05'),1), # 2020
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
                                   seq(as.Date('2021-12-27'),as.Date('2021-12-31'),1),

                                   seq(as.Date('2022-01-01'),as.Date('2022-01-09'),1), # 2022
                                   seq(as.Date('2022-02-28'),as.Date('2022-03-06'),1),
                                   seq(as.Date('2022-04-04'),as.Date('2022-04-18'),1),
                                   seq(as.Date('2022-07-01'),as.Date('2022-09-25'),1), ####
                                   # seq(as.Date('2022-10-31'),as.Date('2022-11-06'),1), # no fall break
                                   seq(as.Date('2022-12-26'),as.Date('2022-12-31'),1)),
                       value    := 1.0]
  }

  # K12 school
  tmp_school_holidays <- copy(d_school_holidays)
  tmp_school_holidays[,category:='schools_closed']
  for (i_age in 0:11) {
    d_calendar_holiday <- rbind(d_calendar_holiday,copy(tmp_school_holidays[,age:=i_age]))
  }
  tmp_school_holidays <- copy(d_school_holidays_secondary)
  tmp_school_holidays[,category:='schools_closed']
  for (i_age in 12:17) {
    d_calendar_holiday <- rbind(d_calendar_holiday,copy(tmp_school_holidays[,age:=i_age]))
  }

  # College
  tmp_college_holidays <- copy(d_college_holidays)
  tmp_college_holidays[,category:='schools_closed']
  for (i_age in 18:25) {
    d_calendar_holiday <- rbind(d_calendar_holiday,copy(tmp_college_holidays[,age:=i_age]))
  }

  ######################################### #
  ##  Contact reductions: other        ####
  ######################################### #
  #       * workplace
  #       * community
  #       * household clusters

  # workplace distancing
  cnt_reduction_workplace <- ifelse('cnt_reduction_workplace' %in% names(config_exp),
                                    config_exp$cnt_reduction_workplace,
                                    0.0)

  data.table(category = "workplace_distancing",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = cnt_reduction_workplace,
             type     = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_workplace_distancing

  # community distancing
  cnt_reduction_other <- ifelse('cnt_reduction_other' %in% names(config_exp),
                                config_exp$cnt_reduction_other,
                                0.0)

  data.table(category = "community_distancing",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = cnt_reduction_other,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_community_distancing

  # collectivity distancing
  cnt_reduction_collectivity <- ifelse('cnt_reduction_collectivity' %in% names(config_exp),
                                       config_exp$cnt_reduction_collectivity,
                                       0.0)

  data.table(category = "collectivity_distancing",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = cnt_reduction_collectivity,
             type = 'double',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_collectivity_distancing

  # household clustering
  data.table(category = "household_clustering",
             # date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),  # TODO: fixed days?
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0,  # TODO: from config
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_household_clustering

  ###################################### #
  ## Imported cases                 ####
  ###################################### #
  data.table(category = "imported_cases",
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 0,  # TODO: from config
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_imported_cases

  ######################################## #
  ##  Contact tracing                 ####
  ######################################## #

  data.table(category = "contact_tracing",
             # date     = seq(as.Date('2020-05-11'),as.Date('2020-08-31'),1),  # TODO: fixed days?
             date     = seq(as.Date(date_start),as.Date(date_end),1),
             value    = 1,  # TODO: from config
             type = 'boolean',
             age = NA_integer_,
             stringsAsFactors = F
  ) -> dcal_contact_tracing

  ############################################# #
  ## MERGE HOLIDAYS & OTHER CALENDAR ITEMS ####
  ############################################# #

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
  ## EXPLORE DATA                        ####
  ########################################### #

  plot_calendar_new(dt_calendar            = d_calendar_all,
                    filename_calendar_full = file_name)

  ########################################### #
  ## SAVE AS CSV	 	                 ####
  ########################################### #

  # # format date
  # d_calendar_all[,date:=format(date,'%Y-%m-%d')]
  # format(d_calendar_all$date,'%Y-%m-%d')

  # save as csv (all calendar info)
  write.table(d_calendar_all,
              file = file_name,sep=',',row.names=F,quote=F)

  unique(d_calendar_all$category)

  return(file_name)
}


plot_calendar_new <- function(dt_calendar, filename_calendar_full){
  smd_print("Plotting calendar...")


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
  x_lim      <- as.Date(c('2021-01-01','2021-12-31'))  # TODO: abstract
  x_lab_year <- paste(unique(year(dt_calendar$date)),sep='-')
  i_cat <- category_opt[2]
  smd_print("cat opt. ", category_opt)

  for(i_cat in category_opt){
    smd_print("cat option ", i_cat)
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

if(0==1){ # debug----
  
  dcal_file <- create_calendar_file(file_name_tag = 'wp_fitting',show_plots = T)
  dcal_wp_distancing <- data.frame(c('2020-03-13',0),
                                   c('2020-03-19',0.85),
                                   c('2020-05-02',0.85),
                                   c('2020-05-03',0.5))
  
  adjust_calendar_file(db_category =  "workplace_distancing",db_update = dcal_wp_distancing, file_name = dcal_file)
}

